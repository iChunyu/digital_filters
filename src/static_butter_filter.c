#include "static_butter_filter.h"
#include <math.h>
#include <string.h>

/* ================================================================== */
/*  Pre-computed Butterworth prototype poles (orders 1–8)              */
/* ================================================================== */
/*
 * Formula:  θ_k = π · (2k + N + 1) / (2N)
 *           p_k = cos(θ_k) − j·sin(θ_k)     for k = 0..N−1
 *
 * All poles lie in the left half-plane (negative real part).
 * Index: butter_proto[order − 1][pole_index].
 */

static const complex_t butter_proto[8][8] = {
    /* N=1 */ {{-1.0000000000f, -0.0000000000f}},
    /* N=2 */ {{-0.7071067812f, -0.7071067812f},
               {-0.7071067812f,  0.7071067812f}},
    /* N=3 */ {{-0.5000000000f, -0.8660254038f},
               {-1.0000000000f, -0.0000000000f},
               {-0.5000000000f,  0.8660254038f}},
    /* N=4 */ {{-0.3826834324f, -0.9238795325f},
               {-0.9238795325f, -0.3826834324f},
               {-0.9238795325f,  0.3826834324f},
               {-0.3826834324f,  0.9238795325f}},
    /* N=5 */ {{-0.3090169944f, -0.9510565163f},
               {-0.8090169944f, -0.5877852523f},
               {-1.0000000000f, -0.0000000000f},
               {-0.8090169944f,  0.5877852523f},
               {-0.3090169944f,  0.9510565163f}},
    /* N=6 */ {{-0.2588190451f, -0.9659258263f},
               {-0.7071067812f, -0.7071067812f},
               {-0.9659258263f, -0.2588190451f},
               {-0.9659258263f,  0.2588190451f},
               {-0.7071067812f,  0.7071067812f},
               {-0.2588190451f,  0.9659258263f}},
    /* N=7 */ {{-0.2225209340f, -0.9749279122f},
               {-0.6234898019f, -0.7818314825f},
               {-0.9009688679f, -0.4338837391f},
               {-1.0000000000f, -0.0000000000f},
               {-0.9009688679f,  0.4338837391f},
               {-0.6234898019f,  0.7818314825f},
               {-0.2225209340f,  0.9749279122f}},
    /* N=8 */ {{-0.1950903220f, -0.9807852804f},
               {-0.5555702330f, -0.8314696123f},
               {-0.8314696123f, -0.5555702330f},
               {-0.9807852804f, -0.1950903220f},
               {-0.9807852804f,  0.1950903220f},
               {-0.8314696123f,  0.5555702330f},
               {-0.5555702330f,  0.8314696123f},
               {-0.1950903220f,  0.9807852804f}},
};

/* ================================================================== */
/*  Shared design pipeline                                             */
/* ================================================================== */

/*
 * Maximum number of poles after BP/BS transform for 8th order: 2 × 8 = 16.
 * Stack usage during init: poles (128 B) + zeros (128 B) + sos (192 B) ≈ 448 B.
 */
#define STATIC_BUTTER_MAX_NP 16
#define STATIC_BUTTER_MAX_NS 8

/**
 * @brief Run the Butterworth design pipeline and deploy biquad coefficients.
 *
 * @param sections      Output biquad array.
 * @param max_sections  Capacity of @p sections.
 * @param order         Prototype filter order N.
 * @param type          filter_type_e.
 * @param wc1           Pre-warped lower cutoff rad/s (2π·prewarp(fc1, fs)).
 * @param wc2           Pre-warped upper cutoff rad/s for BP/BS (unused for LP/HP).
 * @param fs            Sampling frequency in Hz.
 * @param w0_norm       Digital normalisation frequency (rad/sample).
 * @return              Number of sections deployed, or 0 on failure.
 */
static uint8_t butter_design_static(biquad_filter_t *sections,
                                     uint8_t max_sections,
                                     uint8_t order, uint8_t type,
                                     float wc1, float wc2, float fs,
                                     float w0_norm)
{
    complex_t poles[STATIC_BUTTER_MAX_NP];
    complex_t zeros[STATIC_BUTTER_MAX_NP];
    uint8_t np = order;
    uint8_t nz = 0;

    /* 1. Copy prototype poles from pre-computed table */
    memcpy(poles, butter_proto[order - 1], (size_t)order * sizeof(complex_t));

    /* 2. Analog frequency transform */
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

    /* 3. Bilinear transform: s → z */
    bilinear_transform(poles, np, fs);
    bilinear_transform(zeros, nz, fs);

    /* 4. Zero-pad: prototype zeros at s=∞ → z = -1 (not for BS). */
    if (type != FILTER_BANDSTOP) {
        for (uint8_t i = nz; i < np; i++) {
            zeros[i].re = -1.0f;
            zeros[i].im =  0.0f;
        }
        nz = np;
    }

    /* 5. Pair poles and zeros → SOS coefficients. */
    uint8_t ns = (np + 1) / 2;
    if (ns > max_sections) return 0;

    float sos[STATIC_BUTTER_MAX_NS][6];
    uint8_t n_sections = zpk2sos(zeros, poles, np, sos, w0_norm);
    if (n_sections > max_sections) return 0;

    /* 6. Deploy to biquad sections. */
    for (uint8_t i = 0; i < n_sections; i++) {
        float num[3] = {sos[i][0], sos[i][1], sos[i][2]};
        float den[3] = {sos[i][3], sos[i][4], sos[i][5]};
        biquad_filter_init(&sections[i], num, den);
    }

    return n_sections;
}

/* ================================================================== */
/*  Per-order init helpers (called from macro-generated functions)     */
/* ================================================================== */

static void butter_static_init_lp(static_butter_t *f, uint8_t order,
                                   float fc, float fs)
{
    f->valid = 0;
    f->type  = FILTER_LOWPASS;
    f->order = order;
    f->fc1   = fc;
    f->fc2   = 0.0f;
    f->fs    = fs;

    if (order == 0 || fc <= 0.0f || fc >= fs * 0.5f) return;

    float wc = 2.0f * (float)M_PI * prewarp(fc, fs);
    uint8_t ns = butter_design_static(f->sections, f->num_sections,
                                       order, FILTER_LOWPASS,
                                       wc, 0.0f, fs, 0.0f);
    if (ns == 0) return;
    f->num_sections = ns;
    f->valid = 1;
}

static void butter_static_init_hp(static_butter_t *f, uint8_t order,
                                   float fc, float fs)
{
    f->valid = 0;
    f->type  = FILTER_HIGHPASS;
    f->order = order;
    f->fc1   = fc;
    f->fc2   = 0.0f;
    f->fs    = fs;

    if (order == 0 || fc <= 0.0f || fc >= fs * 0.5f) return;

    float wc = 2.0f * (float)M_PI * prewarp(fc, fs);
    uint8_t ns = butter_design_static(f->sections, f->num_sections,
                                       order, FILTER_HIGHPASS,
                                       wc, 0.0f, fs, (float)M_PI);
    if (ns == 0) return;
    f->num_sections = ns;
    f->valid = 1;
}

static void butter_static_init_bp(static_butter_t *f, uint8_t order,
                                   float fc1, float fc2, float fs)
{
    f->valid = 0;
    f->type  = FILTER_BANDPASS;
    f->order = order;
    f->fc1   = fc1;
    f->fc2   = fc2;
    f->fs    = fs;

    if (order == 0 || fc1 <= 0.0f || fc1 >= fs * 0.5f) return;
    if (fc2 <= fc1 || fc2 >= fs * 0.5f) return;

    float wc1 = 2.0f * (float)M_PI * prewarp(fc1, fs);
    float wc2 = 2.0f * (float)M_PI * prewarp(fc2, fs);
    float w0_norm = 2.0f * (float)M_PI * sqrtf(fc1 * fc2) / fs;

    uint8_t ns = butter_design_static(f->sections, f->num_sections,
                                       order, FILTER_BANDPASS,
                                       wc1, wc2, fs, w0_norm);
    if (ns == 0) return;
    f->num_sections = ns;
    f->valid = 1;
}

static void butter_static_init_bs(static_butter_t *f, uint8_t order,
                                   float fc1, float fc2, float fs)
{
    f->valid = 0;
    f->type  = FILTER_BANDSTOP;
    f->order = order;
    f->fc1   = fc1;
    f->fc2   = fc2;
    f->fs    = fs;

    if (order == 0 || fc1 <= 0.0f || fc1 >= fs * 0.5f) return;
    if (fc2 <= fc1 || fc2 >= fs * 0.5f) return;

    float wc1 = 2.0f * (float)M_PI * prewarp(fc1, fs);
    float wc2 = 2.0f * (float)M_PI * prewarp(fc2, fs);

    uint8_t ns = butter_design_static(f->sections, f->num_sections,
                                       order, FILTER_BANDSTOP,
                                       wc1, wc2, fs, 0.0f);
    if (ns == 0) return;
    f->num_sections = ns;
    f->valid = 1;
}

/* ================================================================== */
/*  Per-order init functions (macro-generated)                         */
/* ================================================================== */

/*
 * LP init:  butter_lp_Nth_init(f, fc, fs)
 * HP init:  butter_hp_Nth_init(f, fc, fs)
 * BP init:  butter_bp_Nth_init(f, fc1, fc2, fs)
 * BS init:  butter_bs_Nth_init(f, fc1, fc2, fs)
 *
 * Each sets num_sections from the struct's array size before calling the
 * shared helper, so the helper can validate capacity.
 */

/* Lowpass */
#define X(order, ns, ol) \
    void butter_lp_##ol##_init(butter_lp_##ol##_t *f, float fc, float fs) { \
        f->num_sections = ns; \
        butter_static_init_lp((static_butter_t *)f, order, fc, fs); \
    }
FOR_EACH_STATIC_BUTTER_LP_ORDER
#undef X

/* Highpass */
#define X(order, ns, ol) \
    void butter_hp_##ol##_init(butter_hp_##ol##_t *f, float fc, float fs) { \
        f->num_sections = ns; \
        butter_static_init_hp((static_butter_t *)f, order, fc, fs); \
    }
FOR_EACH_STATIC_BUTTER_LP_ORDER
#undef X

/* Bandpass */
#define X(order, ns, ol) \
    void butter_bp_##ol##_init(butter_bp_##ol##_t *f, float fc1, float fc2, float fs) { \
        f->num_sections = ns; \
        butter_static_init_bp((static_butter_t *)f, order, fc1, fc2, fs); \
    }
FOR_EACH_STATIC_BUTTER_BP_ORDER
#undef X

/* Bandstop */
#define X(order, ns, ol) \
    void butter_bs_##ol##_init(butter_bs_##ol##_t *f, float fc1, float fc2, float fs) { \
        f->num_sections = ns; \
        butter_static_init_bs((static_butter_t *)f, order, fc1, fc2, fs); \
    }
FOR_EACH_STATIC_BUTTER_BP_ORDER
#undef X

/* ================================================================== */
/*  Unified update / reset                                             */
/* ================================================================== */

float static_butter_update(static_butter_t *f, float input)
{
    if (!f->valid) return input;

    float x = input;
    for (uint8_t i = 0; i < f->num_sections; i++) {
        x = biquad_filter_update(&f->sections[i], x);
    }
    return x;
}

void static_butter_reset(static_butter_t *f, float equilibrium)
{
    if (!f->valid) return;

    float x = equilibrium;
    for (uint8_t i = 0; i < f->num_sections; i++) {
        biquad_filter_reset(&f->sections[i], x);
        x = biquad_filter_get_output(&f->sections[i]);
    }
}
