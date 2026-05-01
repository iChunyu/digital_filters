#include "static_cheby_filter.h"
#include <math.h>
#include <string.h>

/* ================================================================== */
/*  Chebyshev prototypes (runtime computation — depend on ripple)      */
/* ================================================================== */

static void cheby1_proto(complex_t *poles, uint8_t n, float epsilon)
{
    float mu = asinhf(1.0f / epsilon) / (float)n;
    float sinh_mu = sinhf(mu);
    float cosh_mu = coshf(mu);

    for (uint8_t k = 0; k < n; k++) {
        float theta = (float)(2 * (k + 1) + n - 1)
                    / (2.0f * (float)n) * (float)M_PI;
        poles[k].re = sinh_mu * cosf(theta);
        poles[k].im = cosh_mu * sinf(theta);
    }
}

/*
 * @return Number of finite zeros (n for even order, n−1 for odd).
 */
static uint8_t cheby2_proto(complex_t *poles, complex_t *zeros, uint8_t n,
                             float epsilon)
{
    float mu = asinhf(1.0f / epsilon) / (float)n;
    float sinh_mu = sinhf(mu);
    float cosh_mu = coshf(mu);
    uint8_t nz = 0;

    for (uint8_t k = 0; k < n; k++) {
        float theta = (float)(2 * (k + 1) + n - 1)
                    / (2.0f * (float)n) * (float)M_PI;

        /* Pole: 1 / (Chebyshev I pole) */
        float den = sinh_mu * sinh_mu * cosf(theta) * cosf(theta)
                  + cosh_mu * cosh_mu * sinf(theta) * sinf(theta);
        poles[k].re =  sinh_mu * cosf(theta) / den;
        poles[k].im = -cosh_mu * sinf(theta) / den;

        /* Zero: j / sin(theta).  Skips θ where sin(θ) ≈ 0 (odd N). */
        float s = sinf(theta);
        if (fabsf(s) > 1e-7f) {
            zeros[nz].re = 0.0f;
            zeros[nz].im = 1.0f / s;
            nz++;
        }
    }
    return nz;
}

/* ================================================================== */
/*  Shared design pipeline                                             */
/* ================================================================== */

#define STATIC_CHEBY_MAX_NP 16
#define STATIC_CHEBY_MAX_NS 8

/**
 * @brief Run the Chebyshev design pipeline and deploy biquad coefficients.
 *
 * For Chebyshev I:  pass nz = 0 (prototype has no finite zeros).
 * For Chebyshev II: pass nz = cheby2 finite-zero count.
 *
 * @return Number of sections deployed, or 0 on failure.
 */
static uint8_t static_cheby_design(biquad_filter_t *sections,
                                    uint8_t max_sections,
                                    uint8_t order, uint8_t type,
                                    float wc1, float wc2, float fs,
                                    float w0_norm,
                                    const complex_t *proto_poles,
                                    uint8_t np, uint8_t nz,
                                    const complex_t *proto_zeros)
{
    complex_t poles[STATIC_CHEBY_MAX_NP];
    complex_t zeros[STATIC_CHEBY_MAX_NP];

    memcpy(poles, proto_poles, (size_t)np * sizeof(complex_t));
    if (nz > 0 && proto_zeros != NULL) {
        memcpy(zeros, proto_zeros, (size_t)nz * sizeof(complex_t));
    }

    /* 1. Analog frequency transform */
    switch (type) {
    case FILTER_LOWPASS:
        analog_lp_transform(poles, np, zeros, nz, wc1);
        break;
    case FILTER_HIGHPASS:
        analog_hp_transform(poles, np, zeros, nz, wc1);
        for (uint8_t i = nz; i < np; i++) {
            zeros[i].re = 0.0f;
            zeros[i].im = 0.0f;
        }
        nz = np;
        break;
    case FILTER_BANDPASS: {
        float w0 = sqrtf(wc1 * wc2);
        float xi = wc2 - wc1;
        analog_bp_transform(poles, &np, zeros, &nz, w0, xi);
        break;
    }
    case FILTER_BANDSTOP: {
        float w0 = sqrtf(wc1 * wc2);
        float xi = wc2 - wc1;
        analog_bs_transform(poles, &np, zeros, &nz, w0, xi);
        break;
    }
    default:
        return 0;
    }

    /* 2. Bilinear transform: s → z */
    bilinear_transform(poles, np, fs);
    bilinear_transform(zeros, nz, fs);

    /* 3. Zero-pad: prototype zeros at s=∞ → z = -1 (not for BS). */
    if (type != FILTER_BANDSTOP) {
        for (uint8_t i = nz; i < np; i++) {
            zeros[i].re = -1.0f;
            zeros[i].im =  0.0f;
        }
        nz = np;
    }

    /* 4. Pair poles and zeros → SOS coefficients. */
    uint8_t ns = (np + 1) / 2;
    if (ns > max_sections) return 0;

    float sos[STATIC_CHEBY_MAX_NS][6];
    uint8_t n_sections = zpk2sos(zeros, poles, np, sos, w0_norm);
    if (n_sections > max_sections) return 0;

    /* 5. Deploy to biquad sections. */
    for (uint8_t i = 0; i < n_sections; i++) {
        float num[3] = {sos[i][0], sos[i][1], sos[i][2]};
        float den[3] = {sos[i][3], sos[i][4], sos[i][5]};
        biquad_filter_init(&sections[i], num, den);
    }

    return n_sections;
}

/* ================================================================== */
/*  Per-type init helpers                                              */
/* ================================================================== */

/* ── Chebyshev I ──────────────────────────────────────────────────── */

static uint8_t static_cheby1_lp_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order, float fc, float fs,
                                      float ripple_db)
{
    if (order == 0 || fc <= 0.0f || fc >= fs * 0.5f || ripple_db <= 0.0f)
        return 0;

    float epsilon = sqrtf(powf(10.0f, ripple_db / 10.0f) - 1.0f);
    float wc = 2.0f * (float)M_PI * prewarp(fc, fs);

    complex_t poles[8];
    cheby1_proto(poles, order, epsilon);

    return static_cheby_design(sections, max_sections,
                                order, FILTER_LOWPASS,
                                wc, 0.0f, fs, 0.0f,
                                poles, order, 0, NULL);
}

static uint8_t static_cheby1_hp_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order, float fc, float fs,
                                      float ripple_db)
{
    if (order == 0 || fc <= 0.0f || fc >= fs * 0.5f || ripple_db <= 0.0f)
        return 0;

    float epsilon = sqrtf(powf(10.0f, ripple_db / 10.0f) - 1.0f);
    float wc = 2.0f * (float)M_PI * prewarp(fc, fs);

    complex_t poles[8];
    cheby1_proto(poles, order, epsilon);

    return static_cheby_design(sections, max_sections,
                                order, FILTER_HIGHPASS,
                                wc, 0.0f, fs, (float)M_PI,
                                poles, order, 0, NULL);
}

static uint8_t static_cheby1_bp_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order,
                                      float fc1, float fc2, float fs,
                                      float ripple_db)
{
    if (order == 0 || fc1 <= 0.0f || fc1 >= fs * 0.5f || ripple_db <= 0.0f)
        return 0;
    if (fc2 <= fc1 || fc2 >= fs * 0.5f) return 0;

    float epsilon = sqrtf(powf(10.0f, ripple_db / 10.0f) - 1.0f);
    float wc1 = 2.0f * (float)M_PI * prewarp(fc1, fs);
    float wc2 = 2.0f * (float)M_PI * prewarp(fc2, fs);
    float w0_norm = 2.0f * (float)M_PI * sqrtf(fc1 * fc2) / fs;

    complex_t poles[8];
    cheby1_proto(poles, order, epsilon);

    return static_cheby_design(sections, max_sections,
                                order, FILTER_BANDPASS,
                                wc1, wc2, fs, w0_norm,
                                poles, order, 0, NULL);
}

static uint8_t static_cheby1_bs_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order,
                                      float fc1, float fc2, float fs,
                                      float ripple_db)
{
    if (order == 0 || fc1 <= 0.0f || fc1 >= fs * 0.5f || ripple_db <= 0.0f)
        return 0;
    if (fc2 <= fc1 || fc2 >= fs * 0.5f) return 0;

    float epsilon = sqrtf(powf(10.0f, ripple_db / 10.0f) - 1.0f);
    float wc1 = 2.0f * (float)M_PI * prewarp(fc1, fs);
    float wc2 = 2.0f * (float)M_PI * prewarp(fc2, fs);

    complex_t poles[8];
    cheby1_proto(poles, order, epsilon);

    return static_cheby_design(sections, max_sections,
                                order, FILTER_BANDSTOP,
                                wc1, wc2, fs, 0.0f,
                                poles, order, 0, NULL);
}

/* ── Chebyshev II ─────────────────────────────────────────────────── */

static uint8_t static_cheby2_lp_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order, float fc, float fs,
                                      float ripple_db)
{
    if (order == 0 || fc <= 0.0f || fc >= fs * 0.5f || ripple_db <= 0.0f)
        return 0;

    float epsilon = 1.0f / sqrtf(powf(10.0f, ripple_db / 10.0f) - 1.0f);
    float wc = 2.0f * (float)M_PI * prewarp(fc, fs);

    complex_t poles[8], zeros[8];
    uint8_t nz = cheby2_proto(poles, zeros, order, epsilon);

    return static_cheby_design(sections, max_sections,
                                order, FILTER_LOWPASS,
                                wc, 0.0f, fs, 0.0f,
                                poles, order, nz, zeros);
}

static uint8_t static_cheby2_hp_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order, float fc, float fs,
                                      float ripple_db)
{
    if (order == 0 || fc <= 0.0f || fc >= fs * 0.5f || ripple_db <= 0.0f)
        return 0;

    float epsilon = 1.0f / sqrtf(powf(10.0f, ripple_db / 10.0f) - 1.0f);
    float wc = 2.0f * (float)M_PI * prewarp(fc, fs);

    complex_t poles[8], zeros[8];
    uint8_t nz = cheby2_proto(poles, zeros, order, epsilon);

    return static_cheby_design(sections, max_sections,
                                order, FILTER_HIGHPASS,
                                wc, 0.0f, fs, (float)M_PI,
                                poles, order, nz, zeros);
}

static uint8_t static_cheby2_bp_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order,
                                      float fc1, float fc2, float fs,
                                      float ripple_db)
{
    if (order == 0 || fc1 <= 0.0f || fc1 >= fs * 0.5f || ripple_db <= 0.0f)
        return 0;
    if (fc2 <= fc1 || fc2 >= fs * 0.5f) return 0;

    float epsilon = 1.0f / sqrtf(powf(10.0f, ripple_db / 10.0f) - 1.0f);
    float wc1 = 2.0f * (float)M_PI * prewarp(fc1, fs);
    float wc2 = 2.0f * (float)M_PI * prewarp(fc2, fs);
    float w0_norm = 2.0f * (float)M_PI * sqrtf(fc1 * fc2) / fs;

    complex_t poles[8], zeros[8];
    uint8_t nz = cheby2_proto(poles, zeros, order, epsilon);

    return static_cheby_design(sections, max_sections,
                                order, FILTER_BANDPASS,
                                wc1, wc2, fs, w0_norm,
                                poles, order, nz, zeros);
}

static uint8_t static_cheby2_bs_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order,
                                      float fc1, float fc2, float fs,
                                      float ripple_db)
{
    if (order == 0 || fc1 <= 0.0f || fc1 >= fs * 0.5f || ripple_db <= 0.0f)
        return 0;
    if (fc2 <= fc1 || fc2 >= fs * 0.5f) return 0;

    float epsilon = 1.0f / sqrtf(powf(10.0f, ripple_db / 10.0f) - 1.0f);
    float wc1 = 2.0f * (float)M_PI * prewarp(fc1, fs);
    float wc2 = 2.0f * (float)M_PI * prewarp(fc2, fs);

    complex_t poles[8], zeros[8];
    uint8_t nz = cheby2_proto(poles, zeros, order, epsilon);

    return static_cheby_design(sections, max_sections,
                                order, FILTER_BANDSTOP,
                                wc1, wc2, fs, 0.0f,
                                poles, order, nz, zeros);
}

/* ================================================================== */
/*  Internal helpers                                                    */
/* ================================================================== */

static float static_cheby_update(biquad_filter_t *sections,
                                  uint8_t num_sections, float input)
{
    float x = input;
    for (uint8_t i = 0; i < num_sections; i++) {
        x = biquad_filter_update(&sections[i], x);
    }
    return x;
}

static void static_cheby_reset(biquad_filter_t *sections,
                                uint8_t num_sections, float equilibrium)
{
    float x = equilibrium;
    for (uint8_t i = 0; i < num_sections; i++) {
        biquad_filter_reset(&sections[i], x);
        x = biquad_filter_get_output(&sections[i]);
    }
}

/* ================================================================== */
/*  Per-order init functions (macro-generated)                         */
/* ================================================================== */

/* Chebyshev I — lowpass */
#define X(ord, ns, ol) \
    void cheby1_lp_##ol##_init(cheby1_lp_##ol##_t *f, float fc, float fs, float ripple_db) { \
        f->type = FILTER_LOWPASS; \
        f->order = ord; \
        f->fc1 = fc; \
        f->fc2 = 0.0f; \
        f->fs = fs; \
        f->ripple_db = ripple_db; \
        f->valid = 0; \
        uint8_t n = static_cheby1_lp_init(f->sections, ns, ord, fc, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_STATIC_CHEBY_LP_ORDER
#undef X

/* Chebyshev I — highpass */
#define X(ord, ns, ol) \
    void cheby1_hp_##ol##_init(cheby1_hp_##ol##_t *f, float fc, float fs, float ripple_db) { \
        f->type = FILTER_HIGHPASS; \
        f->order = ord; \
        f->fc1 = fc; \
        f->fc2 = 0.0f; \
        f->fs = fs; \
        f->ripple_db = ripple_db; \
        f->valid = 0; \
        uint8_t n = static_cheby1_hp_init(f->sections, ns, ord, fc, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_STATIC_CHEBY_LP_ORDER
#undef X

/* Chebyshev I — bandpass */
#define X(ord, ns, ol) \
    void cheby1_bp_##ol##_init(cheby1_bp_##ol##_t *f, float fc1, float fc2, float fs, float ripple_db) { \
        f->type = FILTER_BANDPASS; \
        f->order = ord; \
        f->fc1 = fc1; \
        f->fc2 = fc2; \
        f->fs = fs; \
        f->ripple_db = ripple_db; \
        f->valid = 0; \
        uint8_t n = static_cheby1_bp_init(f->sections, ns, ord, fc1, fc2, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_STATIC_CHEBY_BP_ORDER
#undef X

/* Chebyshev I — bandstop */
#define X(ord, ns, ol) \
    void cheby1_bs_##ol##_init(cheby1_bs_##ol##_t *f, float fc1, float fc2, float fs, float ripple_db) { \
        f->type = FILTER_BANDSTOP; \
        f->order = ord; \
        f->fc1 = fc1; \
        f->fc2 = fc2; \
        f->fs = fs; \
        f->ripple_db = ripple_db; \
        f->valid = 0; \
        uint8_t n = static_cheby1_bs_init(f->sections, ns, ord, fc1, fc2, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_STATIC_CHEBY_BP_ORDER
#undef X

/* Chebyshev II — lowpass */
#define X(ord, ns, ol) \
    void cheby2_lp_##ol##_init(cheby2_lp_##ol##_t *f, float fc, float fs, float ripple_db) { \
        f->type = FILTER_LOWPASS; \
        f->order = ord; \
        f->fc1 = fc; \
        f->fc2 = 0.0f; \
        f->fs = fs; \
        f->ripple_db = ripple_db; \
        f->valid = 0; \
        uint8_t n = static_cheby2_lp_init(f->sections, ns, ord, fc, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_STATIC_CHEBY_LP_ORDER
#undef X

/* Chebyshev II — highpass */
#define X(ord, ns, ol) \
    void cheby2_hp_##ol##_init(cheby2_hp_##ol##_t *f, float fc, float fs, float ripple_db) { \
        f->type = FILTER_HIGHPASS; \
        f->order = ord; \
        f->fc1 = fc; \
        f->fc2 = 0.0f; \
        f->fs = fs; \
        f->ripple_db = ripple_db; \
        f->valid = 0; \
        uint8_t n = static_cheby2_hp_init(f->sections, ns, ord, fc, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_STATIC_CHEBY_LP_ORDER
#undef X

/* Chebyshev II — bandpass */
#define X(ord, ns, ol) \
    void cheby2_bp_##ol##_init(cheby2_bp_##ol##_t *f, float fc1, float fc2, float fs, float ripple_db) { \
        f->type = FILTER_BANDPASS; \
        f->order = ord; \
        f->fc1 = fc1; \
        f->fc2 = fc2; \
        f->fs = fs; \
        f->ripple_db = ripple_db; \
        f->valid = 0; \
        uint8_t n = static_cheby2_bp_init(f->sections, ns, ord, fc1, fc2, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_STATIC_CHEBY_BP_ORDER
#undef X

/* Chebyshev II — bandstop */
#define X(ord, ns, ol) \
    void cheby2_bs_##ol##_init(cheby2_bs_##ol##_t *f, float fc1, float fc2, float fs, float ripple_db) { \
        f->type = FILTER_BANDSTOP; \
        f->order = ord; \
        f->fc1 = fc1; \
        f->fc2 = fc2; \
        f->fs = fs; \
        f->ripple_db = ripple_db; \
        f->valid = 0; \
        uint8_t n = static_cheby2_bs_init(f->sections, ns, ord, fc1, fc2, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_STATIC_CHEBY_BP_ORDER
#undef X

/* ================================================================== */
/*  Per-order update / reset (macro-generated)                          */
/* ================================================================== */

/* Chebyshev I — lowpass */
#define X(ord, ns, ol) \
    float cheby1_lp_##ol##_update(cheby1_lp_##ol##_t *f, float input) { \
        if (!f->valid) return input; \
        return static_cheby_update(f->sections, f->num_sections, input); \
    } \
    void cheby1_lp_##ol##_reset(cheby1_lp_##ol##_t *f, float equilibrium) { \
        if (!f->valid) return; \
        static_cheby_reset(f->sections, f->num_sections, equilibrium); \
    }
FOR_EACH_STATIC_CHEBY_LP_ORDER
#undef X

/* Chebyshev I — highpass */
#define X(ord, ns, ol) \
    float cheby1_hp_##ol##_update(cheby1_hp_##ol##_t *f, float input) { \
        if (!f->valid) return input; \
        return static_cheby_update(f->sections, f->num_sections, input); \
    } \
    void cheby1_hp_##ol##_reset(cheby1_hp_##ol##_t *f, float equilibrium) { \
        if (!f->valid) return; \
        static_cheby_reset(f->sections, f->num_sections, equilibrium); \
    }
FOR_EACH_STATIC_CHEBY_LP_ORDER
#undef X

/* Chebyshev I — bandpass */
#define X(ord, ns, ol) \
    float cheby1_bp_##ol##_update(cheby1_bp_##ol##_t *f, float input) { \
        if (!f->valid) return input; \
        return static_cheby_update(f->sections, f->num_sections, input); \
    } \
    void cheby1_bp_##ol##_reset(cheby1_bp_##ol##_t *f, float equilibrium) { \
        if (!f->valid) return; \
        static_cheby_reset(f->sections, f->num_sections, equilibrium); \
    }
FOR_EACH_STATIC_CHEBY_BP_ORDER
#undef X

/* Chebyshev I — bandstop */
#define X(ord, ns, ol) \
    float cheby1_bs_##ol##_update(cheby1_bs_##ol##_t *f, float input) { \
        if (!f->valid) return input; \
        return static_cheby_update(f->sections, f->num_sections, input); \
    } \
    void cheby1_bs_##ol##_reset(cheby1_bs_##ol##_t *f, float equilibrium) { \
        if (!f->valid) return; \
        static_cheby_reset(f->sections, f->num_sections, equilibrium); \
    }
FOR_EACH_STATIC_CHEBY_BP_ORDER
#undef X

/* Chebyshev II — lowpass */
#define X(ord, ns, ol) \
    float cheby2_lp_##ol##_update(cheby2_lp_##ol##_t *f, float input) { \
        if (!f->valid) return input; \
        return static_cheby_update(f->sections, f->num_sections, input); \
    } \
    void cheby2_lp_##ol##_reset(cheby2_lp_##ol##_t *f, float equilibrium) { \
        if (!f->valid) return; \
        static_cheby_reset(f->sections, f->num_sections, equilibrium); \
    }
FOR_EACH_STATIC_CHEBY_LP_ORDER
#undef X

/* Chebyshev II — highpass */
#define X(ord, ns, ol) \
    float cheby2_hp_##ol##_update(cheby2_hp_##ol##_t *f, float input) { \
        if (!f->valid) return input; \
        return static_cheby_update(f->sections, f->num_sections, input); \
    } \
    void cheby2_hp_##ol##_reset(cheby2_hp_##ol##_t *f, float equilibrium) { \
        if (!f->valid) return; \
        static_cheby_reset(f->sections, f->num_sections, equilibrium); \
    }
FOR_EACH_STATIC_CHEBY_LP_ORDER
#undef X

/* Chebyshev II — bandpass */
#define X(ord, ns, ol) \
    float cheby2_bp_##ol##_update(cheby2_bp_##ol##_t *f, float input) { \
        if (!f->valid) return input; \
        return static_cheby_update(f->sections, f->num_sections, input); \
    } \
    void cheby2_bp_##ol##_reset(cheby2_bp_##ol##_t *f, float equilibrium) { \
        if (!f->valid) return; \
        static_cheby_reset(f->sections, f->num_sections, equilibrium); \
    }
FOR_EACH_STATIC_CHEBY_BP_ORDER
#undef X

/* Chebyshev II — bandstop */
#define X(ord, ns, ol) \
    float cheby2_bs_##ol##_update(cheby2_bs_##ol##_t *f, float input) { \
        if (!f->valid) return input; \
        return static_cheby_update(f->sections, f->num_sections, input); \
    } \
    void cheby2_bs_##ol##_reset(cheby2_bs_##ol##_t *f, float equilibrium) { \
        if (!f->valid) return; \
        static_cheby_reset(f->sections, f->num_sections, equilibrium); \
    }
FOR_EACH_STATIC_CHEBY_BP_ORDER
#undef X
