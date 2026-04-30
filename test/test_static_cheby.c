#include "static_cheby_filter.h"
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

int main(void)
{
    float y;

    /* ── Chebyshev I LP 2nd order, fc=2 Hz, fs=20 Hz, 1 dB ripple ────── */

    cheby1_lp_2nd_t c1;
    cheby1_lp_2nd_init(&c1, 2.0f, 20.0f, 1.0f);
    CHECK(c1.valid == 1, "cheby1 LP 2nd init valid");
    CHECK(c1.num_sections == 1, "cheby1 LP 2nd → 1 section");

    static_cheby_reset((static_cheby_t *)&c1, 1.0f);
    y = static_cheby_update((static_cheby_t *)&c1, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "cheby1 LP 2nd DC gain ~ 1");

    /* ── Chebyshev I LP 3rd order (odd → 2 sections) ──────────────────── */

    cheby1_lp_3rd_t c1_3;
    cheby1_lp_3rd_init(&c1_3, 3.0f, 20.0f, 0.5f);
    CHECK(c1_3.valid == 1, "cheby1 LP 3rd init valid");
    CHECK(c1_3.num_sections == 2, "cheby1 LP 3rd → 2 sections");

    static_cheby_reset((static_cheby_t *)&c1_3, 1.0f);
    y = static_cheby_update((static_cheby_t *)&c1_3, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "cheby1 LP 3rd DC gain ~ 1");

    /* ── Chebyshev I HP 2nd order ────────────────────────────────────── */

    cheby1_hp_2nd_t c1hp;
    cheby1_hp_2nd_init(&c1hp, 5.0f, 40.0f, 1.0f);
    CHECK(c1hp.valid == 1, "cheby1 HP 2nd init valid");

    static_cheby_reset((static_cheby_t *)&c1hp, 1.0f);
    y = static_cheby_update((static_cheby_t *)&c1hp, 1.0f);
    CHECK(CLOSE(y, 0.0f, 1e-3f), "cheby1 HP 2nd blocks DC");

    /* ── Chebyshev II LP 2nd order, 40 dB stopband ────────────────────── */

    cheby2_lp_2nd_t c2;
    cheby2_lp_2nd_init(&c2, 2.0f, 20.0f, 40.0f);
    CHECK(c2.valid == 1, "cheby2 LP 2nd init valid");
    CHECK(c2.num_sections == 1, "cheby2 LP 2nd → 1 section");

    static_cheby_reset((static_cheby_t *)&c2, 1.0f);
    y = static_cheby_update((static_cheby_t *)&c2, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "cheby2 LP 2nd DC gain ~ 1");

    /* ── Chebyshev II LP 3rd order (odd → 2 sections) ─────────────────── */

    cheby2_lp_3rd_t c2_3;
    cheby2_lp_3rd_init(&c2_3, 3.0f, 20.0f, 40.0f);
    CHECK(c2_3.valid == 1, "cheby2 LP 3rd init valid");
    CHECK(c2_3.num_sections == 2, "cheby2 LP 3rd → 2 sections");

    static_cheby_reset((static_cheby_t *)&c2_3, 1.0f);
    y = static_cheby_update((static_cheby_t *)&c2_3, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "cheby2 LP 3rd DC gain ~ 1");

    /* ── Chebyshev I BP 2nd order, fc1=2, fc2=5, fs=40 Hz ────────────── */

    cheby1_bp_2nd_t c1bp;
    cheby1_bp_2nd_init(&c1bp, 2.0f, 5.0f, 40.0f, 1.0f);
    CHECK(c1bp.valid == 1, "cheby1 BP 2nd init valid");
    CHECK(c1bp.num_sections == 2, "cheby1 BP 2nd → 2 sections");

    static_cheby_reset((static_cheby_t *)&c1bp, 1.0f);
    y = static_cheby_update((static_cheby_t *)&c1bp, 1.0f);
    CHECK(CLOSE(y, 0.0f, 1e-3f), "cheby1 BP 2nd blocks DC");

    /* Centre frequency gain ≈ 1.0 */
    float f0 = sqrtf(2.0f * 5.0f);
    static_cheby_reset((static_cheby_t *)&c1bp, 0.0f);
    float bp_max = 0.0f;
    for (int n = 0; n < 800; n++) {
        float t = (float)n / 40.0f;
        float in_val = sinf(2.0f * (float)M_PI * f0 * t);
        y = static_cheby_update((static_cheby_t *)&c1bp, in_val);
        if (n > 400 && fabsf(y) > bp_max) bp_max = fabsf(y);
    }
    CHECK(CLOSE(bp_max, 1.0f, 0.15f), "cheby1 BP 2nd centre freq gain ~ 1");

    /* ── Chebyshev I BS 2nd order ─────────────────────────────────────── */

    cheby1_bs_2nd_t c1bs;
    cheby1_bs_2nd_init(&c1bs, 2.0f, 5.0f, 40.0f, 1.0f);
    CHECK(c1bs.valid == 1, "cheby1 BS 2nd init valid");

    static_cheby_reset((static_cheby_t *)&c1bs, 1.0f);
    y = static_cheby_update((static_cheby_t *)&c1bs, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "cheby1 BS 2nd DC gain ~ 1");

    /* ── Chebyshev II HP 2nd order ────────────────────────────────────── */

    cheby2_hp_2nd_t c2hp;
    cheby2_hp_2nd_init(&c2hp, 5.0f, 40.0f, 40.0f);
    CHECK(c2hp.valid == 1, "cheby2 HP 2nd init valid");

    static_cheby_reset((static_cheby_t *)&c2hp, 0.0f);
    float nyq_gain = 0.0f;
    for (int n = 0; n < 200; n++) {
        float x = (n % 2 == 0) ? 1.0f : -1.0f;
        y = static_cheby_update((static_cheby_t *)&c2hp, x);
        if (n > 100 && fabsf(y) > nyq_gain) nyq_gain = fabsf(y);
    }
    CHECK(CLOSE(nyq_gain, 1.0f, 1e-2f), "cheby2 HP 2nd Nyquist gain ~ 1");

    /* ── Chebyshev II BP 3rd order (odd, tests zero-padding) ──────────── */

    cheby2_bp_3rd_t c2bp3;
    cheby2_bp_3rd_init(&c2bp3, 3.0f, 8.0f, 40.0f, 40.0f);
    CHECK(c2bp3.valid == 1, "cheby2 BP 3rd init valid (odd)");
    CHECK(c2bp3.num_sections == 3, "cheby2 BP 3rd → 3 sections");

    static_cheby_reset((static_cheby_t *)&c2bp3, 1.0f);
    y = static_cheby_update((static_cheby_t *)&c2bp3, 1.0f);
    CHECK(CLOSE(y, 0.0f, 1e-3f), "cheby2 BP 3rd blocks DC");

    /* ── Chebyshev II BS 2nd order ────────────────────────────────────── */

    cheby2_bs_2nd_t c2bs;
    cheby2_bs_2nd_init(&c2bs, 3.0f, 8.0f, 40.0f, 40.0f);
    CHECK(c2bs.valid == 1, "cheby2 BS 2nd init valid");

    static_cheby_reset((static_cheby_t *)&c2bs, 1.0f);
    y = static_cheby_update((static_cheby_t *)&c2bs, 1.0f);
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

    y = static_cheby_update((static_cheby_t *)&ci1, 0.5f);
    CHECK(y == 0.5f, "cheby1 invalid filter passthrough");

    /* ── Report ───────────────────────────────────────────────────────── */

    if (failures) {
        fprintf(stderr, "%d test(s) FAILED.\n", failures);
        return EXIT_FAILURE;
    }

    printf("All tests passed.\n");
    return EXIT_SUCCESS;
}
