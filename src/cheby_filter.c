#include "cheby_filter.h"
#include <math.h>
#include <string.h>

/* ================================================================== */
/*  Chebyshev prototypes (runtime computation — depend on ripple)      */
/* ================================================================== */

static void cheby1_proto(complex_t *poles, uint8_t n, float epsilon)
{
    float inv_eps = 1.0f / epsilon;
    float mu = logf(inv_eps + sqrtf(inv_eps * inv_eps + 1.0f)) / (float)n;
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
    float inv_eps = 1.0f / epsilon;
    float mu = logf(inv_eps + sqrtf(inv_eps * inv_eps + 1.0f)) / (float)n;
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

#define CHEBY_MAX_NP 16
#define CHEBY_MAX_NS 8

/**
 * @brief Run the Chebyshev design pipeline and deploy biquad coefficients.
 *
 * For Chebyshev I:  pass nz = 0 (prototype has no finite zeros).
 * For Chebyshev II: pass nz = cheby2 finite-zero count.
 *
 * @param k  Prototype gain (prod(-p) for cheby1, prod(-p)/prod(-z) for cheby2).
 * @return Number of sections deployed, or 0 on failure.
 */
static uint8_t cheby_design(biquad_filter_t *sections,
                                    uint8_t max_sections,
                                    uint8_t order, uint8_t type,
                                    float wc1, float wc2, float fs,
                                    float k,
                                    const complex_t *proto_poles,
                                    uint8_t np, uint8_t nz,
                                    const complex_t *proto_zeros)
{
    complex_t poles[CHEBY_MAX_NP];
    complex_t zeros[CHEBY_MAX_NP];
    uint8_t degree = np - nz; /* prototype relative degree */

    memcpy(poles, proto_poles, (size_t)np * sizeof(complex_t));
    if (nz > 0 && proto_zeros != NULL) {
        memcpy(zeros, proto_zeros, (size_t)nz * sizeof(complex_t));
    }

    /* LP/BP gain scaling: wc^degree (LP) or xi^degree (BP), folded into the
       bilinear gain below (see bilinear_zpk_gain_scaled). */
    float gs = 0.0f;
    uint8_t gdeg = degree;

    /* 1. Analog frequency transform */
    switch (type) {
    case FILTER_LOWPASS:
        analog_lp_transform(poles, np, zeros, nz, wc1);
        gs = wc1;
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

    /* 2. Bilinear gain (on s-domain zp, before bilinear transform clobbers
       them).  k == 0 deploys an all-zero-numerator (silence) filter — a
       degenerate outcome, not a valid design. */
    if (gs != 0.0f) {
        k = bilinear_zpk_gain_scaled(k, gs, gdeg, zeros, nz, poles, np,
                                     2.0f * fs);
    } else {
        k = bilinear_zpk_gain(k, zeros, nz, poles, np, 2.0f * fs);
    }
    if (!isfinite(k) || k == 0.0f) return 0;

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

    float sos[CHEBY_MAX_NS][6];
    uint8_t n_sections = zpk2sos_impl(zeros, poles, np, sos, k);
    if (n_sections != ns) return 0;

    /* 6. Deploy to biquad sections. */
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
static uint8_t cheby_check_gains(const biquad_filter_t *sections,
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

/*
 * Expected DC/Nyquist gains per family, type and order parity
 * (verified against scipy):
 *   cheby1: passband edges sit at 10^(−rp/20) for even orders, 1 for odd.
 *   cheby2: stopband edges sit at 10^(−rs/20) for even orders, 0 for odd.
 */
static float cheby1_edge_gain(uint8_t order, float ripple_db)
{
    return (order % 2 == 0) ? powf(10.0f, -ripple_db / 20.0f) : 1.0f;
}

static float cheby2_edge_gain(uint8_t order, float ripple_db)
{
    return (order % 2 == 0) ? powf(10.0f, -ripple_db / 20.0f) : 0.0f;
}

/* ── Chebyshev I ──────────────────────────────────────────────────── */

static uint8_t cheby1_lp_init(biquad_filter_t *sections,
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

    float k = 1.0f / zpk_hp_bs_gain(1.0f, NULL, 0, poles, order);
    if (order % 2 == 0) k /= sqrtf(1.0f + epsilon * epsilon);

    uint8_t n = cheby_design(sections, max_sections,
                             order, FILTER_LOWPASS,
                             wc, 0.0f, fs, k,
                             poles, order, 0, NULL);
    if (n == 0) return 0;
    /* cheby1 LP: DC gain 1 (odd) / 10^(−rp/20) (even), Nyquist 0 */
    if (!cheby_check_gains(sections, n,
                           cheby1_edge_gain(order, ripple_db), 0.0f)) return 0;
    return n;
}

static uint8_t cheby1_hp_init(biquad_filter_t *sections,
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

    float k = 1.0f / zpk_hp_bs_gain(1.0f, NULL, 0, poles, order);
    if (order % 2 == 0) k /= sqrtf(1.0f + epsilon * epsilon);

    uint8_t n = cheby_design(sections, max_sections,
                             order, FILTER_HIGHPASS,
                             wc, 0.0f, fs, k,
                             poles, order, 0, NULL);
    if (n == 0) return 0;
    /* cheby1 HP: DC 0, Nyquist gain 1 (odd) / 10^(−rp/20) (even) */
    if (!cheby_check_gains(sections, n, 0.0f,
                           cheby1_edge_gain(order, ripple_db))) return 0;
    return n;
}

static uint8_t cheby1_bp_init(biquad_filter_t *sections,
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

    float k = 1.0f / zpk_hp_bs_gain(1.0f, NULL, 0, poles, order);
    if (order % 2 == 0) k /= sqrtf(1.0f + epsilon * epsilon);

    uint8_t n = cheby_design(sections, max_sections,
                             order, FILTER_BANDPASS,
                             wc1, wc2, fs, k,
                             poles, order, 0, NULL);
    if (n == 0) return 0;
    /* cheby1 BP: DC and Nyquist gains 0 */
    if (!cheby_check_gains(sections, n, 0.0f, 0.0f)) return 0;
    return n;
}

static uint8_t cheby1_bs_init(biquad_filter_t *sections,
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

    float k = 1.0f / zpk_hp_bs_gain(1.0f, NULL, 0, poles, order);
    if (order % 2 == 0) k /= sqrtf(1.0f + epsilon * epsilon);

    uint8_t n = cheby_design(sections, max_sections,
                             order, FILTER_BANDSTOP,
                             wc1, wc2, fs, k,
                             poles, order, 0, NULL);
    if (n == 0) return 0;
    /* cheby1 BS: DC and Nyquist gains 1 (odd) / 10^(−rp/20) (even) */
    if (!cheby_check_gains(sections, n,
                           cheby1_edge_gain(order, ripple_db),
                           cheby1_edge_gain(order, ripple_db))) return 0;
    return n;
}

/* ── Chebyshev II ─────────────────────────────────────────────────── */

static uint8_t cheby2_lp_init(biquad_filter_t *sections,
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

    float k = 1.0f / zpk_hp_bs_gain(1.0f, zeros, nz, poles, order);

    uint8_t n = cheby_design(sections, max_sections,
                             order, FILTER_LOWPASS,
                             wc, 0.0f, fs, k,
                             poles, order, nz, zeros);
    if (n == 0) return 0;
    /* cheby2 LP: DC gain 1, Nyquist 0 (odd) / 10^(−rs/20) (even) */
    if (!cheby_check_gains(sections, n, 1.0f,
                           cheby2_edge_gain(order, ripple_db))) return 0;
    return n;
}

static uint8_t cheby2_hp_init(biquad_filter_t *sections,
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

    float k = 1.0f / zpk_hp_bs_gain(1.0f, zeros, nz, poles, order);

    uint8_t n = cheby_design(sections, max_sections,
                             order, FILTER_HIGHPASS,
                             wc, 0.0f, fs, k,
                             poles, order, nz, zeros);
    if (n == 0) return 0;
    /* cheby2 HP: DC 0 (odd) / 10^(−rs/20) (even), Nyquist gain 1 */
    if (!cheby_check_gains(sections, n,
                           cheby2_edge_gain(order, ripple_db), 1.0f)) return 0;
    return n;
}

static uint8_t cheby2_bp_init(biquad_filter_t *sections,
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

    float k = 1.0f / zpk_hp_bs_gain(1.0f, zeros, nz, poles, order);

    uint8_t n = cheby_design(sections, max_sections,
                             order, FILTER_BANDPASS,
                             wc1, wc2, fs, k,
                             poles, order, nz, zeros);
    if (n == 0) return 0;
    /* cheby2 BP: DC and Nyquist gains 0 (odd) / 10^(−rs/20) (even) */
    if (!cheby_check_gains(sections, n,
                           cheby2_edge_gain(order, ripple_db),
                           cheby2_edge_gain(order, ripple_db))) return 0;
    return n;
}

static uint8_t cheby2_bs_init(biquad_filter_t *sections,
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

    float k = 1.0f / zpk_hp_bs_gain(1.0f, zeros, nz, poles, order);

    uint8_t n = cheby_design(sections, max_sections,
                             order, FILTER_BANDSTOP,
                             wc1, wc2, fs, k,
                             poles, order, nz, zeros);
    if (n == 0) return 0;
    /* cheby2 BS: DC and Nyquist gains 1 */
    if (!cheby_check_gains(sections, n, 1.0f, 1.0f)) return 0;
    return n;
}

/* ================================================================== */
/*  Internal helpers                                                    */
/* ================================================================== */

static float cheby_update(biquad_filter_t *sections,
                                  uint8_t num_sections, float input)
{
    float x = input;
    for (uint8_t i = 0; i < num_sections; i++) {
        x = biquad_filter_update(&sections[i], x);
    }
    return x;
}

static void cheby_reset(biquad_filter_t *sections,
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
        uint8_t n = cheby1_lp_init(f->sections, ns, ord, fc, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_CHEBY_LP_ORDER
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
        uint8_t n = cheby1_hp_init(f->sections, ns, ord, fc, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_CHEBY_LP_ORDER
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
        uint8_t n = cheby1_bp_init(f->sections, ns, ord, fc1, fc2, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_CHEBY_BP_ORDER
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
        uint8_t n = cheby1_bs_init(f->sections, ns, ord, fc1, fc2, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_CHEBY_BP_ORDER
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
        uint8_t n = cheby2_lp_init(f->sections, ns, ord, fc, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_CHEBY_LP_ORDER
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
        uint8_t n = cheby2_hp_init(f->sections, ns, ord, fc, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_CHEBY_LP_ORDER
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
        uint8_t n = cheby2_bp_init(f->sections, ns, ord, fc1, fc2, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_CHEBY_BP_ORDER
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
        uint8_t n = cheby2_bs_init(f->sections, ns, ord, fc1, fc2, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_CHEBY_BP_ORDER
#undef X

/* ================================================================== */
/*  Per-order update / reset (macro-generated)                          */
/* ================================================================== */

/* Chebyshev I — lowpass */
#define X(ord, ns, ol) \
    float cheby1_lp_##ol##_update(cheby1_lp_##ol##_t *f, float input) { \
        if (!f->valid) return input; \
        return cheby_update(f->sections, f->num_sections, input); \
    } \
    void cheby1_lp_##ol##_reset(cheby1_lp_##ol##_t *f, float equilibrium) { \
        if (!f->valid) return; \
        cheby_reset(f->sections, f->num_sections, equilibrium); \
    }
FOR_EACH_CHEBY_LP_ORDER
#undef X

/* Chebyshev I — highpass */
#define X(ord, ns, ol) \
    float cheby1_hp_##ol##_update(cheby1_hp_##ol##_t *f, float input) { \
        if (!f->valid) return input; \
        return cheby_update(f->sections, f->num_sections, input); \
    } \
    void cheby1_hp_##ol##_reset(cheby1_hp_##ol##_t *f, float equilibrium) { \
        if (!f->valid) return; \
        cheby_reset(f->sections, f->num_sections, equilibrium); \
    }
FOR_EACH_CHEBY_LP_ORDER
#undef X

/* Chebyshev I — bandpass */
#define X(ord, ns, ol) \
    float cheby1_bp_##ol##_update(cheby1_bp_##ol##_t *f, float input) { \
        if (!f->valid) return input; \
        return cheby_update(f->sections, f->num_sections, input); \
    } \
    void cheby1_bp_##ol##_reset(cheby1_bp_##ol##_t *f, float equilibrium) { \
        if (!f->valid) return; \
        cheby_reset(f->sections, f->num_sections, equilibrium); \
    }
FOR_EACH_CHEBY_BP_ORDER
#undef X

/* Chebyshev I — bandstop */
#define X(ord, ns, ol) \
    float cheby1_bs_##ol##_update(cheby1_bs_##ol##_t *f, float input) { \
        if (!f->valid) return input; \
        return cheby_update(f->sections, f->num_sections, input); \
    } \
    void cheby1_bs_##ol##_reset(cheby1_bs_##ol##_t *f, float equilibrium) { \
        if (!f->valid) return; \
        cheby_reset(f->sections, f->num_sections, equilibrium); \
    }
FOR_EACH_CHEBY_BP_ORDER
#undef X

/* Chebyshev II — lowpass */
#define X(ord, ns, ol) \
    float cheby2_lp_##ol##_update(cheby2_lp_##ol##_t *f, float input) { \
        if (!f->valid) return input; \
        return cheby_update(f->sections, f->num_sections, input); \
    } \
    void cheby2_lp_##ol##_reset(cheby2_lp_##ol##_t *f, float equilibrium) { \
        if (!f->valid) return; \
        cheby_reset(f->sections, f->num_sections, equilibrium); \
    }
FOR_EACH_CHEBY_LP_ORDER
#undef X

/* Chebyshev II — highpass */
#define X(ord, ns, ol) \
    float cheby2_hp_##ol##_update(cheby2_hp_##ol##_t *f, float input) { \
        if (!f->valid) return input; \
        return cheby_update(f->sections, f->num_sections, input); \
    } \
    void cheby2_hp_##ol##_reset(cheby2_hp_##ol##_t *f, float equilibrium) { \
        if (!f->valid) return; \
        cheby_reset(f->sections, f->num_sections, equilibrium); \
    }
FOR_EACH_CHEBY_LP_ORDER
#undef X

/* Chebyshev II — bandpass */
#define X(ord, ns, ol) \
    float cheby2_bp_##ol##_update(cheby2_bp_##ol##_t *f, float input) { \
        if (!f->valid) return input; \
        return cheby_update(f->sections, f->num_sections, input); \
    } \
    void cheby2_bp_##ol##_reset(cheby2_bp_##ol##_t *f, float equilibrium) { \
        if (!f->valid) return; \
        cheby_reset(f->sections, f->num_sections, equilibrium); \
    }
FOR_EACH_CHEBY_BP_ORDER
#undef X

/* Chebyshev II — bandstop */
#define X(ord, ns, ol) \
    float cheby2_bs_##ol##_update(cheby2_bs_##ol##_t *f, float input) { \
        if (!f->valid) return input; \
        return cheby_update(f->sections, f->num_sections, input); \
    } \
    void cheby2_bs_##ol##_reset(cheby2_bs_##ol##_t *f, float equilibrium) { \
        if (!f->valid) return; \
        cheby_reset(f->sections, f->num_sections, equilibrium); \
    }
FOR_EACH_CHEBY_BP_ORDER
#undef X
