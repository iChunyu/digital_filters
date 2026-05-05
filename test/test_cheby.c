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

int main(void)
{
    float y;

    /* ── Chebyshev I LP 2nd order, fc=2 Hz, fs=20 Hz, 1 dB ripple ──── */

    cheby1_t c1;
    cheby1_init(&c1, FILTER_LOWPASS, 2, 2.0f, 0.0f, 20.0f, 1.0f);
    CHECK(c1.valid == 1, "cheby1 LP init");
    CHECK(c1.sos->num_sections == 1, "cheby1 2nd-order → 1 section");

    cheby1_reset(&c1, 1.0f);
    y = cheby1_update(&c1, 1.0f);
    CHECK(CLOSE(y, 0.891251f, 1e-4f), "cheby1 LP DC gain ~ -1 dB (even order)");

    cheby1_destroy(&c1);
    CHECK(c1.valid == 0, "cheby1 valid cleared");

    /* ── Chebyshev I LP 3rd order (odd → 1st+2nd sections) ─────────── */

    cheby1_t c1_3;
    cheby1_init(&c1_3, FILTER_LOWPASS, 3, 3.0f, 0.0f, 20.0f, 0.5f);
    CHECK(c1_3.valid == 1, "cheby1 3rd-order LP init");
    CHECK(c1_3.sos->num_sections == 2, "cheby1 3rd-order → 2 sections");

    cheby1_reset(&c1_3, 1.0f);
    y = cheby1_update(&c1_3, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "cheby1 3rd-order DC gain ~ 1");

    cheby1_destroy(&c1_3);

    /* ── Chebyshev I HP 2nd order ──────────────────────────────────── */

    cheby1_t c1hp;
    cheby1_init(&c1hp, FILTER_HIGHPASS, 2, 5.0f, 0.0f, 40.0f, 1.0f);
    CHECK(c1hp.valid == 1, "cheby1 HP init");

    cheby1_reset(&c1hp, 1.0f);
    y = cheby1_update(&c1hp, 1.0f);
    CHECK(CLOSE(y, 0.0f, 1e-3f), "cheby1 HP blocks DC");

    cheby1_destroy(&c1hp);

    /* ── Chebyshev II LP 2nd order, 40 dB stopband ─────────────────── */

    cheby2_t c2;
    cheby2_init(&c2, FILTER_LOWPASS, 2, 2.0f, 0.0f, 20.0f, 40.0f);
    CHECK(c2.valid == 1, "cheby2 LP init");
    CHECK(c2.sos->num_sections == 1, "cheby2 2nd-order → 1 section");

    cheby2_reset(&c2, 1.0f);
    y = cheby2_update(&c2, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "cheby2 LP DC gain ~ 1");

    cheby2_destroy(&c2);

    /* ── Chebyshev II LP 3rd order (odd → tests infinite-zero path) ── */

    cheby2_t c2_3;
    cheby2_init(&c2_3, FILTER_LOWPASS, 3, 3.0f, 0.0f, 20.0f, 40.0f);
    CHECK(c2_3.valid == 1, "cheby2 3rd-order LP init");
    CHECK(c2_3.sos->num_sections == 2, "cheby2 3rd-order → 2 sections");

    cheby2_reset(&c2_3, 1.0f);
    y = cheby2_update(&c2_3, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "cheby2 3rd-order DC gain ~ 1");

    cheby2_destroy(&c2_3);

    /* ── Chebyshev I BP 2nd order, fc1=2 Hz, fc2=5 Hz, fs=40 Hz ─────── */

    cheby1_t c1bp;
    cheby1_init(&c1bp, FILTER_BANDPASS, 2, 2.0f, 5.0f, 40.0f, 1.0f);
    CHECK(c1bp.valid == 1, "cheby1 BP init");
    CHECK(c1bp.sos->num_sections == 2, "cheby1 BP → 2 sections");

    cheby1_reset(&c1bp, 1.0f);
    y = cheby1_update(&c1bp, 1.0f);
    CHECK(CLOSE(y, 0.0f, 1e-3f), "cheby1 BP blocks DC");

    /* Centre frequency gain ≈ 1.0 */
    float f0 = sqrtf(2.0f * 5.0f); /* ~3.16 Hz */
    cheby1_reset(&c1bp, 0.0f);
    float bp_max = 0.0f;
    for (int n = 0; n < 800; n++) {
        float t = (float)n / 40.0f;
        float in_val = sinf(2.0f * (float)M_PI * f0 * t);
        y = cheby1_update(&c1bp, in_val);
        if (n > 400 && fabsf(y) > bp_max) bp_max = fabsf(y);
    }
    CHECK(CLOSE(bp_max, 1.0f, 0.15f), "cheby1 BP centre freq gain ~ 1");

    cheby1_destroy(&c1bp);

    /* ── Chebyshev I BS 2nd order ───────────────────────────────────── */

    cheby1_t c1bs;
    cheby1_init(&c1bs, FILTER_BANDSTOP, 2, 2.0f, 5.0f, 40.0f, 1.0f);
    CHECK(c1bs.valid == 1, "cheby1 BS init");

    cheby1_reset(&c1bs, 1.0f);
    y = cheby1_update(&c1bs, 1.0f);
    CHECK(CLOSE(y, 0.891251f, 1e-4f), "cheby1 BS DC gain ~ -1 dB (even order)");

    cheby1_destroy(&c1bs);

    /* ── Chebyshev II HP 2nd order ─────────────────────────────────── */

    cheby2_t c2hp;
    cheby2_init(&c2hp, FILTER_HIGHPASS, 2, 5.0f, 0.0f, 40.0f, 40.0f);
    CHECK(c2hp.valid == 1, "cheby2 HP init");

    /* Chebyshev II HP zeros are not all at DC — test passband (Nyquist) */
    cheby2_reset(&c2hp, 0.0f);
    float nyq_gain = 0.0f;
    for (int n = 0; n < 200; n++) {
        float x = (n % 2 == 0) ? 1.0f : -1.0f;
        y = cheby2_update(&c2hp, x);
        if (n > 100 && fabsf(y) > nyq_gain) nyq_gain = fabsf(y);
    }
    CHECK(CLOSE(nyq_gain, 1.0f, 1e-2f), "cheby2 HP Nyquist gain ~ 1");

    cheby2_destroy(&c2hp);

    /* ── Chebyshev II BP 3rd order (odd, nz=np-1, exercises padding) ── */
    /* NB: even-order Chebyshev II has np=nz → no infinite prototype zeros,
       so BP has non-zero DC gain.  Odd order is needed for DC rejection. */

    cheby2_t c2bp3;
    cheby2_init(&c2bp3, FILTER_BANDPASS, 3, 3.0f, 8.0f, 40.0f, 40.0f);
    CHECK(c2bp3.valid == 1, "cheby2 BP 3rd-order init (odd)");
    CHECK(c2bp3.sos->num_sections == 3, "cheby2 BP 3rd → 3 sections");

    cheby2_reset(&c2bp3, 1.0f);
    y = cheby2_update(&c2bp3, 1.0f);
    CHECK(CLOSE(y, 0.0f, 1e-3f), "cheby2 BP 3rd-order blocks DC");

    cheby2_destroy(&c2bp3);

    /* ── Chebyshev II BS 2nd order ──────────────────────────────────── */

    cheby2_t c2bs;
    cheby2_init(&c2bs, FILTER_BANDSTOP, 2, 3.0f, 8.0f, 40.0f, 40.0f);
    CHECK(c2bs.valid == 1, "cheby2 BS init");

    cheby2_reset(&c2bs, 1.0f);
    y = cheby2_update(&c2bs, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "cheby2 BS DC gain ~ 1");

    cheby2_destroy(&c2bs);

    /* ── Invalid params → valid = 0 ────────────────────────────────── */
    cheby1_t ci1;
    cheby1_init(&ci1, FILTER_LOWPASS, 0, 2.0f, 0.0f, 20.0f, 1.0f);
    CHECK(ci1.valid == 0, "cheby1 order=0 invalid");

    cheby1_init(&ci1, FILTER_LOWPASS, 2, 2.0f, 0.0f, 20.0f, 0.0f);
    CHECK(ci1.valid == 0, "cheby1 ripple=0 invalid");

    cheby2_t ci2;
    cheby2_init(&ci2, FILTER_LOWPASS, 0, 2.0f, 0.0f, 20.0f, 40.0f);
    CHECK(ci2.valid == 0, "cheby2 order=0 invalid");

    /* ── Report ────────────────────────────────────────────────────── */

    if (failures) {
        fprintf(stderr, "%d test(s) FAILED.\n", failures);
        return EXIT_FAILURE;
    }
    printf("All tests passed.\n");
    return EXIT_SUCCESS;
}
