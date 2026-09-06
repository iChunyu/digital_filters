#include "cheby_filter.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures = 0;

#define CHECK(cond, msg)                                         \
    do {                                                         \
        if (!(cond)) {                                           \
            failures++;                                          \
            fprintf(stderr, "FAIL: %s\n", msg);                 \
        }                                                        \
    } while (0)

#define CLOSE(a, b, eps) (fabsf((a) - (b)) <= (eps))

/* Classify the zeros of all sections: counts of +1, -1, unit-circle
   conjugate pairs, and anything else (stray zeros). */
static void section_zeros_stats(biquad_filter_t *secs, int ns,
                                int *n_plus1, int *n_minus1,
                                int *n_unit, int *n_other)
{
    *n_plus1 = *n_minus1 = *n_unit = *n_other = 0;
    for (int i = 0; i < ns; i++) {
        double b0 = secs[i].num_z[0];
        double b1 = secs[i].num_z[1];
        double b2 = secs[i].num_z[2];
        double s = -b1 / b0;
        double c = b2 / b0;
        double disc = s * s - 4.0 * c;
        if (disc >= 0.0) {
            double r = sqrt(disc);
            double zeros[2] = {(s + r) / 2.0, (s - r) / 2.0};
            for (int j = 0; j < 2; j++) {
                if (fabs(zeros[j] - 1.0) < 1e-3) (*n_plus1)++;
                else if (fabs(zeros[j] + 1.0) < 1e-3) (*n_minus1)++;
                else (*n_other)++;
            }
        } else {
            if (fabs(c - 1.0) < 1e-3) (*n_unit)++;
            else (*n_other)++;
        }
    }
}

/* Steady-state amplitude of a sine at freq through the cascade. */
static float measure_gain_at(biquad_filter_t *secs, uint8_t ns,
                             float freq, float fs, int steps)
{
    float x = 0.0f;
    for (uint8_t i = 0; i < ns; i++) {
        biquad_filter_reset(&secs[i], x);
        x = biquad_filter_get_output(&secs[i]);
    }
    float max_out = 0.0f;
    for (int n = 0; n < steps; n++) {
        float t = (float)n / fs;
        float in_val = sinf(2.0f * (float)M_PI * freq * t);
        x = in_val;
        for (uint8_t i = 0; i < ns; i++) {
            x = biquad_filter_update(&secs[i], x);
        }
        if (n > steps * 3 / 4) {
            float a = fabsf(x);
            if (a > max_out) max_out = a;
        }
    }
    return max_out;
}

/* Steady-state gain for the Nyquist tone (alternating ±1). */
static float measure_nyquist_gain(biquad_filter_t *secs, uint8_t ns, int steps)
{
    float x = 0.0f;
    for (uint8_t i = 0; i < ns; i++) {
        biquad_filter_reset(&secs[i], 0.0f);
    }
    float max_out = 0.0f;
    for (int n = 0; n < steps; n++) {
        float in_val = (n % 2 == 0) ? 1.0f : -1.0f;
        x = in_val;
        for (uint8_t i = 0; i < ns; i++) {
            x = biquad_filter_update(&secs[i], x);
        }
        if (n > steps * 3 / 4) {
            float a = fabsf(x);
            if (a > max_out) max_out = a;
        }
    }
    return max_out;
}

/* Analytic cascade DC/Nyquist gains from section coefficients. */
static float cascade_dc_gain(const biquad_filter_t *secs, uint8_t ns)
{
    float h = 1.0f;
    for (uint8_t i = 0; i < ns; i++)
        h *= (secs[i].num_z[0] + secs[i].num_z[1] + secs[i].num_z[2])
           / (1.0f + secs[i].den_z[1] + secs[i].den_z[2]);
    return h;
}

static float cascade_nyq_gain(const biquad_filter_t *secs, uint8_t ns)
{
    float h = 1.0f;
    for (uint8_t i = 0; i < ns; i++)
        h *= (secs[i].num_z[0] - secs[i].num_z[1] + secs[i].num_z[2])
           / (1.0f - secs[i].den_z[1] + secs[i].den_z[2]);
    return h;
}

/* Max steady-state gain over a fixed frequency grid. */
static float max_gain_over(biquad_filter_t *secs, uint8_t ns,
                           float fs, int steps)
{
    const float grid[] = {20.0f, 24.5f, 50.0f, 100.0f, 200.0f, 300.0f, 480.0f};
    float m = 0.0f;
    for (unsigned i = 0; i < sizeof(grid) / sizeof(grid[0]); i++) {
        float g = measure_gain_at(secs, ns, grid[i], fs, steps);
        if (g > m) m = g;
    }
    return m;
}

int main(void)
{
    float y;

    /* ── Chebyshev I LP 2nd order, fc=2 Hz, fs=20 Hz, 1 dB ripple ────── */

    cheby1_lp_2nd_t c1;
    cheby1_lp_2nd_init(&c1, 2.0f, 20.0f, 1.0f);
    CHECK(c1.valid == 1, "cheby1 LP 2nd init valid");
    CHECK(c1.num_sections == 1, "cheby1 LP 2nd → 1 section");

    cheby1_lp_2nd_reset(&c1, 1.0f);
    y = cheby1_lp_2nd_update(&c1, 1.0f);
    CHECK(CLOSE(y, 0.891251f, 1e-4f), "cheby1 LP 2nd DC gain ~ -1 dB (even order)");

    /* ── Chebyshev I LP 3rd order (odd → 2 sections) ──────────────────── */

    cheby1_lp_3rd_t c1_3;
    cheby1_lp_3rd_init(&c1_3, 3.0f, 20.0f, 0.5f);
    CHECK(c1_3.valid == 1, "cheby1 LP 3rd init valid");
    CHECK(c1_3.num_sections == 2, "cheby1 LP 3rd → 2 sections");

    cheby1_lp_3rd_reset(&c1_3, 1.0f);
    y = cheby1_lp_3rd_update(&c1_3, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "cheby1 LP 3rd DC gain ~ 1");

    /* ── Chebyshev I HP 2nd order ────────────────────────────────────── */

    cheby1_hp_2nd_t c1hp;
    cheby1_hp_2nd_init(&c1hp, 5.0f, 40.0f, 1.0f);
    CHECK(c1hp.valid == 1, "cheby1 HP 2nd init valid");

    cheby1_hp_2nd_reset(&c1hp, 1.0f);
    y = cheby1_hp_2nd_update(&c1hp, 1.0f);
    CHECK(CLOSE(y, 0.0f, 1e-3f), "cheby1 HP 2nd blocks DC");

    /* ── Chebyshev II LP 2nd order, 40 dB stopband ────────────────────── */

    cheby2_lp_2nd_t c2;
    cheby2_lp_2nd_init(&c2, 2.0f, 20.0f, 40.0f);
    CHECK(c2.valid == 1, "cheby2 LP 2nd init valid");
    CHECK(c2.num_sections == 1, "cheby2 LP 2nd → 1 section");

    cheby2_lp_2nd_reset(&c2, 1.0f);
    y = cheby2_lp_2nd_update(&c2, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "cheby2 LP 2nd DC gain ~ 1");

    /* ── Chebyshev II LP 3rd order (odd → 2 sections) ─────────────────── */

    cheby2_lp_3rd_t c2_3;
    cheby2_lp_3rd_init(&c2_3, 3.0f, 20.0f, 40.0f);
    CHECK(c2_3.valid == 1, "cheby2 LP 3rd init valid");
    CHECK(c2_3.num_sections == 2, "cheby2 LP 3rd → 2 sections");

    cheby2_lp_3rd_reset(&c2_3, 1.0f);
    y = cheby2_lp_3rd_update(&c2_3, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "cheby2 LP 3rd DC gain ~ 1");

    /* ── Chebyshev I BP 2nd order, fc1=2, fc2=5, fs=40 Hz ────────────── */

    cheby1_bp_2nd_t c1bp;
    cheby1_bp_2nd_init(&c1bp, 2.0f, 5.0f, 40.0f, 1.0f);
    CHECK(c1bp.valid == 1, "cheby1 BP 2nd init valid");
    CHECK(c1bp.num_sections == 2, "cheby1 BP 2nd → 2 sections");

    cheby1_bp_2nd_reset(&c1bp, 1.0f);
    y = cheby1_bp_2nd_update(&c1bp, 1.0f);
    CHECK(CLOSE(y, 0.0f, 1e-3f), "cheby1 BP 2nd blocks DC");

    /* Centre frequency gain ≈ 1.0 */
    float f0 = sqrtf(2.0f * 5.0f);
    cheby1_bp_2nd_reset(&c1bp, 0.0f);
    float bp_max = 0.0f;
    for (int n = 0; n < 800; n++) {
        float t = (float)n / 40.0f;
        float in_val = sinf(2.0f * (float)M_PI * f0 * t);
        y = cheby1_bp_2nd_update(&c1bp, in_val);
        if (n > 400 && fabsf(y) > bp_max) bp_max = fabsf(y);
    }
    CHECK(CLOSE(bp_max, 1.0f, 0.15f), "cheby1 BP 2nd centre freq gain ~ 1");

    /* ── Chebyshev I BS 2nd order ─────────────────────────────────────── */

    cheby1_bs_2nd_t c1bs;
    cheby1_bs_2nd_init(&c1bs, 2.0f, 5.0f, 40.0f, 1.0f);
    CHECK(c1bs.valid == 1, "cheby1 BS 2nd init valid");

    cheby1_bs_2nd_reset(&c1bs, 1.0f);
    y = cheby1_bs_2nd_update(&c1bs, 1.0f);
    CHECK(CLOSE(y, 0.891251f, 1e-4f), "cheby1 BS 2nd DC gain ~ -1 dB (even order)");

    /* ── Chebyshev II HP 2nd order ────────────────────────────────────── */

    cheby2_hp_2nd_t c2hp;
    cheby2_hp_2nd_init(&c2hp, 5.0f, 40.0f, 40.0f);
    CHECK(c2hp.valid == 1, "cheby2 HP 2nd init valid");

    cheby2_hp_2nd_reset(&c2hp, 0.0f);
    float nyq_gain = 0.0f;
    for (int n = 0; n < 200; n++) {
        float x = (n % 2 == 0) ? 1.0f : -1.0f;
        y = cheby2_hp_2nd_update(&c2hp, x);
        if (n > 100 && fabsf(y) > nyq_gain) nyq_gain = fabsf(y);
    }
    CHECK(CLOSE(nyq_gain, 1.0f, 1e-2f), "cheby2 HP 2nd Nyquist gain ~ 1");

    /* ── Chebyshev II BP 3rd order (odd, tests zero-padding) ──────────── */

    cheby2_bp_3rd_t c2bp3;
    cheby2_bp_3rd_init(&c2bp3, 3.0f, 8.0f, 40.0f, 40.0f);
    CHECK(c2bp3.valid == 1, "cheby2 BP 3rd init valid (odd)");
    CHECK(c2bp3.num_sections == 3, "cheby2 BP 3rd → 3 sections");

    cheby2_bp_3rd_reset(&c2bp3, 1.0f);
    y = cheby2_bp_3rd_update(&c2bp3, 1.0f);
    CHECK(CLOSE(y, 0.0f, 1e-3f), "cheby2 BP 3rd blocks DC");

    /* ── Chebyshev II BS 2nd order ────────────────────────────────────── */

    cheby2_bs_2nd_t c2bs;
    cheby2_bs_2nd_init(&c2bs, 3.0f, 8.0f, 40.0f, 40.0f);
    CHECK(c2bs.valid == 1, "cheby2 BS 2nd init valid");

    cheby2_bs_2nd_reset(&c2bs, 1.0f);
    y = cheby2_bs_2nd_update(&c2bs, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "cheby2 BS 2nd DC gain ~ 1");

    /* ── Invalid params → valid = 0 ────────────────────────────────────── */

    cheby1_lp_2nd_t ci1;
    ci1.valid = 0;
    cheby1_lp_2nd_init(&ci1, 2.0f, 20.0f, 0.0f);
    CHECK(ci1.valid == 0, "cheby1 ripple=0 invalid");

    cheby2_lp_2nd_t ci2;
    ci2.valid = 0;
    cheby2_lp_2nd_init(&ci2, 2.0f, 20.0f, 0.0f);
    CHECK(ci2.valid == 0, "cheby2 ripple=0 invalid");

    /* ── Invalid filter passthrough ────────────────────────────────────── */

    y = cheby1_lp_2nd_update(&ci1, 0.5f);
    CHECK(y == 0.5f, "cheby1 invalid filter passthrough");

    /* ── Regression: cheby2 BP 5th [50,120]@1000, rs=0.5 — mixed real /
          complex zeros must pair by type (zpk2sos pairing bug) ─────────── */

    cheby2_bp_5th_t c2bp_reg;
    cheby2_bp_5th_init(&c2bp_reg, 50.0f, 120.0f, 1000.0f, 0.5f);
    CHECK(c2bp_reg.valid == 1, "cheby2 BP 5th reg init valid");
    CHECK(c2bp_reg.num_sections == 5, "cheby2 BP 5th reg → 5 sections");

    int np1, nm1, nu, no;
    section_zeros_stats(c2bp_reg.sections, c2bp_reg.num_sections,
                        &np1, &nm1, &nu, &no);
    CHECK(np1 == 1 && nm1 == 1,
          "cheby2 BP 5th reg: exactly one zero at +1 and one at -1");
    CHECK(nu == 4, "cheby2 BP 5th reg: 4 conjugate pairs on unit circle");
    CHECK(no == 0, "cheby2 BP 5th reg: no stray zeros");

    y = measure_gain_at(c2bp_reg.sections, c2bp_reg.num_sections,
                        77.46f, 1000.0f, 5000);
    CHECK(CLOSE(y, 1.0f, 0.1f), "cheby2 BP 5th reg: centre gain ~ 1");
    y = measure_nyquist_gain(c2bp_reg.sections, c2bp_reg.num_sections, 2000);
    CHECK(y < 1e-3f, "cheby2 BP 5th reg: Nyquist gain ~ 0");

    /* ── Regression: cheby1 BS 5th [50,120]@1000, rp=3 — real pole must
          pair with real pole (zpk2sos pairing bug) ────────────────────── */

    cheby1_bs_5th_t c1bs_reg;
    cheby1_bs_5th_init(&c1bs_reg, 50.0f, 120.0f, 1000.0f, 3.0f);
    CHECK(c1bs_reg.valid == 1, "cheby1 BS 5th reg init valid");
    CHECK(c1bs_reg.num_sections == 5, "cheby1 BS 5th reg → 5 sections");

    y = measure_gain_at(c1bs_reg.sections, c1bs_reg.num_sections,
                        30.0f, 1000.0f, 8000);
    CHECK(CLOSE(y, 0.768f, 0.1f), "cheby1 BS 5th reg: gain@30Hz ~ 0.768");
    y = measure_gain_at(c1bs_reg.sections, c1bs_reg.num_sections,
                        77.46f, 1000.0f, 8000);
    CHECK(y < 0.01f, "cheby1 BS 5th reg: notch at 77.46Hz");
    y = measure_nyquist_gain(c1bs_reg.sections, c1bs_reg.num_sections, 2000);
    CHECK(CLOSE(y, 1.0f, 0.1f), "cheby1 BS 5th reg: Nyquist gain ~ 1");

    /* ── Analytic cascade gains (no transient effects) ─────────────────── */

    /* cascade_dc_gain / cascade_nyq_gain / max_gain_over are defined
       below via the helpers at the bottom of this file. */

    /* ── Regression: wideband cheby2 configs that used to deploy silently
          wrong filters (mispaired poles → DC 97.9, 350x @ 24.5 Hz) ─────── */

    cheby2_bs_7th_t c2bs7;
    cheby2_bs_7th_init(&c2bs7, 20.0f, 480.0f, 1000.0f, 0.5f);
    CHECK(c2bs7.valid == 1, "cheby2 BS 7th [20,480] rs=0.5 valid");
    if (c2bs7.valid) {
        y = cascade_dc_gain(c2bs7.sections, c2bs7.num_sections);
        CHECK(CLOSE(y, 1.0f, 0.1f), "cheby2 BS 7th [20,480] DC gain ~ 1 (was 97.9)");
        /* 24.5 Hz sits in the stopband: legitimate gain is the rs=0.5
           stopband floor ~0.94; the old defect amplified ~350x here. */
        y = measure_gain_at(c2bs7.sections, c2bs7.num_sections,
                            24.5f, 1000.0f, 4000);
        CHECK(y < 2.0f, "cheby2 BS 7th [20,480] no 350x resonance (was 350x)");
    }

    cheby2_bp_3rd_t c2bp3w;
    cheby2_bp_3rd_init(&c2bp3w, 20.0f, 480.0f, 1000.0f, 40.0f);
    CHECK(c2bp3w.valid == 1, "cheby2 BP 3rd [20,480] rs=40 valid");
    if (c2bp3w.valid) {
        y = cascade_dc_gain(c2bp3w.sections, c2bp3w.num_sections);
        CHECK(fabsf(y) < 0.1f, "cheby2 BP 3rd [20,480] blocks DC (was 0.497)");
    }

    /* ── Near-Nyquist LP rescue (gain-chain f32 overflow fixed) ────────── */

    cheby1_lp_8th_t c1nn;
    c1nn.valid = 0;
    cheby1_lp_8th_init(&c1nn, 480.0f, 1000.0f, 1.0f);
    CHECK(c1nn.valid == 1, "cheby1 LP 8th near-Nyquist valid");
    if (c1nn.valid) {
        y = cascade_dc_gain(c1nn.sections, c1nn.num_sections);
        CHECK(CLOSE(y, 0.891251f, 0.05f), "cheby1 LP 8th near-Nyquist DC = 10^(-rp/20)");
    }

    /* ── Degenerate ripple → fail closed (passthrough) ─────────────────── */

    cheby1_lp_2nd_t c1d;
    c1d.valid = 0;
    cheby1_lp_2nd_init(&c1d, 1000.0f, 48000.0f, 200.0f);
    CHECK(c1d.valid == 0, "cheby1 rp=200 rejected (stability margin)");
    cheby1_lp_2nd_init(&c1d, 1000.0f, 48000.0f, 385.0f);
    CHECK(c1d.valid == 0, "cheby1 rp=385 rejected (k=0)");
    y = cheby1_lp_2nd_update(&c1d, 0.5f);
    CHECK(y == 0.5f, "degenerate cheby1 passthrough");

    /* ── Sweep: cheby2 BP/BS matrix must deploy with sane responses ────── */

    static const float sw_bands[3][2] = {{50.0f, 120.0f}, {20.0f, 480.0f},
                                         {10.0f, 499.0f}};
    static const float sw_rs[2] = {0.5f, 40.0f};

    #define X(ord, ns, ol) \
        for (int sw_bi = 0; sw_bi < 3; sw_bi++) \
        for (int sw_ri = 0; sw_ri < 2; sw_ri++) { \
            cheby2_bp_##ol##_t swf; \
            swf.valid = 0; \
            cheby2_bp_##ol##_init(&swf, sw_bands[sw_bi][0], sw_bands[sw_bi][1], \
                                  1000.0f, sw_rs[sw_ri]); \
            CHECK(swf.valid == 1, "sweep cheby2 BP " #ol " valid"); \
            if (swf.valid) { \
                float swg = ((ord) % 2 == 0) \
                          ? powf(10.0f, -sw_rs[sw_ri] / 20.0f) : 0.0f; \
                float swt = (swg == 0.0f) ? 0.1f : 0.25f; \
                y = cascade_dc_gain(swf.sections, swf.num_sections); \
                CHECK(fabsf(y - swg) <= swt, "sweep cheby2 BP " #ol " DC window"); \
                y = cascade_nyq_gain(swf.sections, swf.num_sections); \
                CHECK(fabsf(y - swg) <= swt, "sweep cheby2 BP " #ol " Nyq window"); \
                y = max_gain_over(swf.sections, swf.num_sections, 1000.0f, 2000); \
                CHECK(y < 2.0f, "sweep cheby2 BP " #ol " max|H| sane"); \
            } \
        }
    FOR_EACH_CHEBY_BP_ORDER
    #undef X

    #define X(ord, ns, ol) \
        for (int sw_bi = 0; sw_bi < 3; sw_bi++) \
        for (int sw_ri = 0; sw_ri < 2; sw_ri++) { \
            cheby2_bs_##ol##_t swf; \
            swf.valid = 0; \
            cheby2_bs_##ol##_init(&swf, sw_bands[sw_bi][0], sw_bands[sw_bi][1], \
                                  1000.0f, sw_rs[sw_ri]); \
            CHECK(swf.valid == 1, "sweep cheby2 BS " #ol " valid"); \
            if (swf.valid) { \
                y = cascade_dc_gain(swf.sections, swf.num_sections); \
                CHECK(CLOSE(y, 1.0f, 0.1f), "sweep cheby2 BS " #ol " DC window"); \
                y = cascade_nyq_gain(swf.sections, swf.num_sections); \
                CHECK(CLOSE(y, 1.0f, 0.1f), "sweep cheby2 BS " #ol " Nyq window"); \
                y = max_gain_over(swf.sections, swf.num_sections, 1000.0f, 2000); \
                CHECK(y < 2.0f, "sweep cheby2 BS " #ol " max|H| sane"); \
            } \
        }
    FOR_EACH_CHEBY_BP_ORDER
    #undef X

    /* ── Sweep: cheby1 LP matrix ───────────────────────────────────────── */

    static const float sw_lp_fc[2] = {100.0f, 480.0f};
    static const float sw_lp_rp[2] = {0.5f, 3.0f};

    #define X(ord, ns, ol) \
        for (int sw_fi = 0; sw_fi < 2; sw_fi++) \
        for (int sw_ri = 0; sw_ri < 2; sw_ri++) { \
            cheby1_lp_##ol##_t swf; \
            swf.valid = 0; \
            cheby1_lp_##ol##_init(&swf, sw_lp_fc[sw_fi], 1000.0f, sw_lp_rp[sw_ri]); \
            CHECK(swf.valid == 1, "sweep cheby1 LP " #ol " valid"); \
            if (swf.valid) { \
                float swg = ((ord) % 2 == 0) \
                          ? powf(10.0f, -sw_lp_rp[sw_ri] / 20.0f) : 1.0f; \
                y = cascade_dc_gain(swf.sections, swf.num_sections); \
                CHECK(fabsf(y - swg) <= 0.1f, "sweep cheby1 LP " #ol " DC window"); \
                y = cascade_nyq_gain(swf.sections, swf.num_sections); \
                CHECK(fabsf(y) <= 0.1f, "sweep cheby1 LP " #ol " Nyq window"); \
                float swb = 2.0f * powf(10.0f, sw_lp_rp[sw_ri] / 20.0f); \
                y = max_gain_over(swf.sections, swf.num_sections, 1000.0f, 2000); \
                CHECK(y < swb, "sweep cheby1 LP " #ol " max|H| sane"); \
            } \
        }
    FOR_EACH_CHEBY_LP_ORDER
    #undef X

    /* ── Report ───────────────────────────────────────────────────────── */

    if (failures) {
        fprintf(stderr, "%d test(s) FAILED.\n", failures);
        return EXIT_FAILURE;
    }

    printf("All tests passed.\n");
    return EXIT_SUCCESS;
}
