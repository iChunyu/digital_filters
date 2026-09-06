#include "butter_filter.h"
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
#define BUTTER_MAX_NP 16
#define BUTTER_MAX_NS 8

/**
 * @brief Run the Butterworth design pipeline and deploy biquad coefficients.
 *
 * @param sections      Output biquad array.
 * @param max_sections  Capacity of @p sections.
 * @param order         Prototype filter order N.
 * @param type          filter_type_e.
 * @param wc1           Pre-warped lower cutoff rad/s (2π·prewarp(fc, fs)).
 * @param wc2           Pre-warped upper cutoff rad/s for BP/BS (unused for LP/HP).
 * @param fs            Sampling frequency in Hz.
 * @return              Number of sections deployed, or 0 on failure.
 */
static uint8_t butter_design(biquad_filter_t *sections,
                                     uint8_t max_sections,
                                     uint8_t order, uint8_t type,
                                     float wc1, float wc2, float fs)
{
    complex_t poles[BUTTER_MAX_NP];
    complex_t zeros[BUTTER_MAX_NP];
    uint8_t np = order;
    uint8_t nz = 0;
    float k = 1.0f; /* Butterworth prototype gain */

    /* 1. Copy prototype poles from pre-computed table */
    memcpy(poles, butter_proto[order - 1], (size_t)order * sizeof(complex_t));

    /* LP/BP gain scaling: wc^degree (LP) or xi^degree (BP), folded into the
       bilinear gain below (see bilinear_zpk_gain_scaled). */
    float gs = 0.0f;
    uint8_t gdeg = 0;

    /* 2. Analog frequency transform */
    switch (type) {
    case FILTER_LOWPASS:
        analog_lp_transform(poles, np, zeros, nz, wc1);
        gs = wc1;
        gdeg = order; /* Butterworth: no finite prototype zeros */
        break;
    case FILTER_HIGHPASS:
        k = zpk_hp_bs_gain(k, zeros, nz, poles, np);
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
        gs = xi;
        gdeg = order;
        break;
    }
    case FILTER_BANDSTOP: {
        float w0 = sqrtf(wc1 * wc2);
        float xi = wc2 - wc1;
        k = zpk_hp_bs_gain(k, zeros, nz, poles, np);
        analog_bs_transform(poles, &np, zeros, &nz, w0, xi);
        break;
    }
    default:
        return 0;
    }

    /* 3. Bilinear gain (on s-domain zp, before bilinear transform clobbers
       them).  k == 0 deploys an all-zero-numerator (silence) filter — a
       degenerate outcome, not a valid design. */
    if (gs != 0.0f) {
        k = bilinear_zpk_gain_scaled(k, gs, gdeg, zeros, nz, poles, np,
                                     2.0f * fs);
    } else {
        k = bilinear_zpk_gain(k, zeros, nz, poles, np, 2.0f * fs);
    }
    if (!isfinite(k) || k == 0.0f) return 0;

    /* 4. Bilinear transform: s → z */
    bilinear_transform(poles, np, fs);
    bilinear_transform(zeros, nz, fs);

    /* 5. Zero-pad: prototype zeros at s=∞ → z = -1 (not for BS). */
    if (type != FILTER_BANDSTOP) {
        for (uint8_t i = nz; i < np; i++) {
            zeros[i].re = -1.0f;
            zeros[i].im =  0.0f;
        }
        nz = np;
    }

    /* 6. Pair poles and zeros → SOS coefficients. */
    uint8_t ns = (np + 1) / 2;
    if (ns > max_sections) return 0;

    float sos[BUTTER_MAX_NS][6];
    uint8_t n_sections = zpk2sos_impl(zeros, poles, np, sos, k);
    if (n_sections != ns) return 0;

    /* 7. Deploy to biquad sections. */
    for (uint8_t i = 0; i < n_sections; i++) {
        float num[3] = {sos[i][0], sos[i][1], sos[i][2]};
        float den[3] = {sos[i][3], sos[i][4], sos[i][5]};
        if (!biquad_filter_init(&sections[i], num, den)) return 0;
    }

    return n_sections;
}

/* ================================================================== */
/*  Per-type init helpers                                              */
/* ================================================================== */

/*
 * Analytic cascade-level DC and Nyquist gain check.  Wrong-but-stable
 * pole/zero pairings and gain-scale errors pass every per-section check
 * (each section is finite and Jury-stable) yet wreck the response shape —
 * the confirmed defect measured DC gain 97.9 and a 350x resonance on a
 * bandstop that must sit at ~1.  H(0) and H(π) are exact rational
 * evaluations (no sampling, no trig, O(sections) cost at init time).
 */
static uint8_t butter_check_gains(const biquad_filter_t *sections,
                                  uint8_t num_sections,
                                  float dc_exp, float ny_exp)
{
    float h0 = 1.0f, hn = 1.0f;
    for (uint8_t i = 0; i < num_sections; i++) {
        const biquad_filter_t *b = &sections[i];
        h0 *= (b->num_z[0] + b->num_z[1] + b->num_z[2])
            / (1.0f + b->den_z[1] + b->den_z[2]);
        hn *= (b->num_z[0] - b->num_z[1] + b->num_z[2])
            / (1.0f - b->den_z[1] + b->den_z[2]);
    }
    /* Exact structural gains (0 or 1, enforced by zeros at ±1) tolerate
       ±0.1 — f32 design error is ~1e-3.  Ripple-edge gains (10^(−rp/rs/20)
       for cheby2 even orders) sit on the steepest part of the response
       near the band edges, where f32 bilinear warping shifts the ripple
       pattern; widen those to ±0.25.  The confirmed pairing defects
       deviate by ≥ 0.45, still far outside either window. */
    float tol0 = (dc_exp == 0.0f || dc_exp == 1.0f) ? 0.1f : 0.25f;
    float toln = (ny_exp == 0.0f || ny_exp == 1.0f) ? 0.1f : 0.25f;
    return fabsf(h0 - dc_exp) <= tol0 && fabsf(hn - ny_exp) <= toln;
}

static uint8_t butter_lp_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order, float fc, float fs)
{
    if (order == 0 || order > 8 || fc <= 0.0f || fc >= fs * 0.5f)
        return 0;

    float wc = 2.0f * (float)M_PI * prewarp(fc, fs);
    uint8_t n = butter_design(sections, max_sections,
                              order, FILTER_LOWPASS,
                              wc, 0.0f, fs);
    if (n == 0) return 0;
    /* LP: DC gain 1, Nyquist gain 0 */
    if (!butter_check_gains(sections, n, 1.0f, 0.0f)) return 0;
    return n;
}

static uint8_t butter_hp_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order, float fc, float fs)
{
    if (order == 0 || order > 8 || fc <= 0.0f || fc >= fs * 0.5f)
        return 0;

    float wc = 2.0f * (float)M_PI * prewarp(fc, fs);
    uint8_t n = butter_design(sections, max_sections,
                              order, FILTER_HIGHPASS,
                              wc, 0.0f, fs);
    if (n == 0) return 0;
    /* HP: DC gain 0, Nyquist gain 1 */
    if (!butter_check_gains(sections, n, 0.0f, 1.0f)) return 0;
    return n;
}

static uint8_t butter_bp_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order,
                                      float fc1, float fc2, float fs)
{
    if (order == 0 || order > 8 || fc1 <= 0.0f || fc1 >= fs * 0.5f)
        return 0;
    if (fc2 <= fc1 || fc2 >= fs * 0.5f) return 0;

    float wc1 = 2.0f * (float)M_PI * prewarp(fc1, fs);
    float wc2 = 2.0f * (float)M_PI * prewarp(fc2, fs);

    uint8_t n = butter_design(sections, max_sections,
                              order, FILTER_BANDPASS,
                              wc1, wc2, fs);
    if (n == 0) return 0;
    /* BP: DC and Nyquist gains 0 */
    if (!butter_check_gains(sections, n, 0.0f, 0.0f)) return 0;
    return n;
}

static uint8_t butter_bs_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order,
                                      float fc1, float fc2, float fs)
{
    if (order == 0 || order > 8 || fc1 <= 0.0f || fc1 >= fs * 0.5f)
        return 0;
    if (fc2 <= fc1 || fc2 >= fs * 0.5f) return 0;

    float wc1 = 2.0f * (float)M_PI * prewarp(fc1, fs);
    float wc2 = 2.0f * (float)M_PI * prewarp(fc2, fs);

    uint8_t n = butter_design(sections, max_sections,
                              order, FILTER_BANDSTOP,
                              wc1, wc2, fs);
    if (n == 0) return 0;
    /* BS: DC and Nyquist gains 1 */
    if (!butter_check_gains(sections, n, 1.0f, 1.0f)) return 0;
    return n;
}

/* ================================================================== */
/*  Internal helpers                                                    */
/* ================================================================== */

static float butter_update(biquad_filter_t *sections,
                                   uint8_t num_sections, float input)
{
    float x = input;
    for (uint8_t i = 0; i < num_sections; i++) {
        x = biquad_filter_update(&sections[i], x);
    }
    return x;
}

static void butter_reset(biquad_filter_t *sections,
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

/*
 * LP init:  butter_lp_Nth_init(f, fc, fs)
 * HP init:  butter_hp_Nth_init(f, fc, fs)
 * BP init:  butter_bp_Nth_init(f, fc1, fc2, fs)
 * BS init:  butter_bs_Nth_init(f, fc1, fc2, fs)
 *
 * Each sets metadata fields, validates via the per-type helper, and sets
 * num_sections and valid from the helper's return value.
 */

/* Lowpass */
#define X(ord, ns, ol) \
    void butter_lp_##ol##_init(butter_lp_##ol##_t *f, float fc, float fs) { \
        f->type = FILTER_LOWPASS; \
        f->order = ord; \
        f->fc1 = fc; \
        f->fc2 = 0.0f; \
        f->fs = fs; \
        f->valid = 0; \
        uint8_t n = butter_lp_init(f->sections, ns, ord, fc, fs); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_BUTTER_LP_ORDER
#undef X

/* Highpass */
#define X(ord, ns, ol) \
    void butter_hp_##ol##_init(butter_hp_##ol##_t *f, float fc, float fs) { \
        f->type = FILTER_HIGHPASS; \
        f->order = ord; \
        f->fc1 = fc; \
        f->fc2 = 0.0f; \
        f->fs = fs; \
        f->valid = 0; \
        uint8_t n = butter_hp_init(f->sections, ns, ord, fc, fs); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_BUTTER_LP_ORDER
#undef X

/* Bandpass */
#define X(ord, ns, ol) \
    void butter_bp_##ol##_init(butter_bp_##ol##_t *f, float fc1, float fc2, float fs) { \
        f->type = FILTER_BANDPASS; \
        f->order = ord; \
        f->fc1 = fc1; \
        f->fc2 = fc2; \
        f->fs = fs; \
        f->valid = 0; \
        uint8_t n = butter_bp_init(f->sections, ns, ord, fc1, fc2, fs); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_BUTTER_BP_ORDER
#undef X

/* Bandstop */
#define X(ord, ns, ol) \
    void butter_bs_##ol##_init(butter_bs_##ol##_t *f, float fc1, float fc2, float fs) { \
        f->type = FILTER_BANDSTOP; \
        f->order = ord; \
        f->fc1 = fc1; \
        f->fc2 = fc2; \
        f->fs = fs; \
        f->valid = 0; \
        uint8_t n = butter_bs_init(f->sections, ns, ord, fc1, fc2, fs); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_BUTTER_BP_ORDER
#undef X

/* ================================================================== */
/*  Per-order update / reset (macro-generated)                          */
/* ================================================================== */

/* Lowpass */
#define X(ord, ns, ol) \
    float butter_lp_##ol##_update(butter_lp_##ol##_t *f, float input) { \
        if (!f->valid) return input; \
        return butter_update(f->sections, f->num_sections, input); \
    } \
    void butter_lp_##ol##_reset(butter_lp_##ol##_t *f, float equilibrium) { \
        if (!f->valid) return; \
        butter_reset(f->sections, f->num_sections, equilibrium); \
    }
FOR_EACH_BUTTER_LP_ORDER
#undef X

/* Highpass */
#define X(ord, ns, ol) \
    float butter_hp_##ol##_update(butter_hp_##ol##_t *f, float input) { \
        if (!f->valid) return input; \
        return butter_update(f->sections, f->num_sections, input); \
    } \
    void butter_hp_##ol##_reset(butter_hp_##ol##_t *f, float equilibrium) { \
        if (!f->valid) return; \
        butter_reset(f->sections, f->num_sections, equilibrium); \
    }
FOR_EACH_BUTTER_LP_ORDER
#undef X

/* Bandpass */
#define X(ord, ns, ol) \
    float butter_bp_##ol##_update(butter_bp_##ol##_t *f, float input) { \
        if (!f->valid) return input; \
        return butter_update(f->sections, f->num_sections, input); \
    } \
    void butter_bp_##ol##_reset(butter_bp_##ol##_t *f, float equilibrium) { \
        if (!f->valid) return; \
        butter_reset(f->sections, f->num_sections, equilibrium); \
    }
FOR_EACH_BUTTER_BP_ORDER
#undef X

/* Bandstop */
#define X(ord, ns, ol) \
    float butter_bs_##ol##_update(butter_bs_##ol##_t *f, float input) { \
        if (!f->valid) return input; \
        return butter_update(f->sections, f->num_sections, input); \
    } \
    void butter_bs_##ol##_reset(butter_bs_##ol##_t *f, float equilibrium) { \
        if (!f->valid) return; \
        butter_reset(f->sections, f->num_sections, equilibrium); \
    }
FOR_EACH_BUTTER_BP_ORDER
#undef X
