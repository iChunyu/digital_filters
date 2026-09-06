#include "butter_filter.h"
#include <math.h>
#include <stddef.h>

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
/*  Per-type init helpers                                              */
/* ================================================================== */

static uint8_t butter_lp_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order, float fc, float fs)
{
    if (order == 0 || order > 8 || fc <= 0.0f || fc >= fs * 0.5f)
        return 0;

    float wc = 2.0f * (float)M_PI * prewarp(fc, fs);
    uint8_t n = design_filter(sections, max_sections,
                              FILTER_LOWPASS,
                              wc, 0.0f, fs, 1.0f,
                              butter_proto[order - 1], order, NULL, 0);
    if (n == 0) return 0;
    /* LP: DC gain 1, Nyquist gain 0 */
    if (!check_cascade_gains(sections, n, 1.0f, 0.0f)) return 0;
    return n;
}

static uint8_t butter_hp_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order, float fc, float fs)
{
    if (order == 0 || order > 8 || fc <= 0.0f || fc >= fs * 0.5f)
        return 0;

    float wc = 2.0f * (float)M_PI * prewarp(fc, fs);
    uint8_t n = design_filter(sections, max_sections,
                              FILTER_HIGHPASS,
                              wc, 0.0f, fs, 1.0f,
                              butter_proto[order - 1], order, NULL, 0);
    if (n == 0) return 0;
    /* HP: DC gain 0, Nyquist gain 1 */
    if (!check_cascade_gains(sections, n, 0.0f, 1.0f)) return 0;
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

    uint8_t n = design_filter(sections, max_sections,
                              FILTER_BANDPASS,
                              wc1, wc2, fs, 1.0f,
                              butter_proto[order - 1], order, NULL, 0);
    if (n == 0) return 0;
    /* BP: DC and Nyquist gains 0 */
    if (!check_cascade_gains(sections, n, 0.0f, 0.0f)) return 0;
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

    uint8_t n = design_filter(sections, max_sections,
                              FILTER_BANDSTOP,
                              wc1, wc2, fs, 1.0f,
                              butter_proto[order - 1], order, NULL, 0);
    if (n == 0) return 0;
    /* BS: DC and Nyquist gains 1 */
    if (!check_cascade_gains(sections, n, 1.0f, 1.0f)) return 0;
    return n;
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
