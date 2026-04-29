#include "butter_filter.h"
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

/* Run a tone through the filter and measure steady-state amplitude. */
static float measure_gain(butter_t *b, float freq, float fs, int steps)
{
    butter_reset(b, 0.0f);
    float max_out = 0.0f;

    for (int n = 0; n < steps; n++) {
        float t = (float)n / fs;
        float in_val = sinf(2.0f * (float)M_PI * freq * t);
        float out_val = butter_update(b, in_val);
        if (n > steps / 2) { /* skip transient */
            if (fabsf(out_val) > max_out) {
                max_out = fabsf(out_val);
            }
        }
    }
    return max_out; /* approximate amplitude gain (input amplitude = 1) */
}

int main(void)
{
    /* ── LP Butterworth 2nd order, fc=2 Hz, fs=20 Hz ──────────────── */

    butter_t blp;
    butter_init(&blp, FILTER_LOWPASS, 2, 2.0f, 0.0f, 20.0f);
    CHECK(blp.valid == 1, "LP init valid");
    CHECK(blp.sos != NULL, "LP sos allocated");
    CHECK(blp.sos->num_sections == 1, "LP 2nd-order → 1 section");

    /* DC gain ≈ 1.0 */
    butter_reset(&blp, 1.0f);
    float y = butter_update(&blp, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "LP DC gain ~ 1");

    /* Attenuation at Nyquist (10 Hz): gain should be ≪ 1 */
    float gn = measure_gain(&blp, 10.0f, 20.0f, 400);
    CHECK(gn < 0.15f, "LP Nyquist attenuation");

    butter_destroy(&blp);
    CHECK(blp.valid == 0, "LP valid cleared after destroy");
    CHECK(blp.sos == NULL, "LP sos NULL after destroy");

    /* ── HP Butterworth 2nd order, fc=5 Hz, fs=40 Hz ──────────────── */

    butter_t bhp;
    butter_init(&bhp, FILTER_HIGHPASS, 2, 5.0f, 0.0f, 40.0f);
    CHECK(bhp.valid == 1, "HP init valid");

    /* DC gain ≈ 0 */
    butter_reset(&bhp, 1.0f);
    y = butter_update(&bhp, 1.0f);
    CHECK(CLOSE(y, 0.0f, 1e-3f), "HP blocks DC");

    /* Nyquist gain ≈ 1.0 */
    butter_reset(&bhp, 0.0f);
    float nyq_gain = 0.0f;
    for (int n = 0; n < 200; n++) {
        float x = (n % 2 == 0) ? 1.0f : -1.0f; /* Nyquist tone */
        y = butter_update(&bhp, x);
        if (n > 100 && fabsf(y) > nyq_gain) nyq_gain = fabsf(y);
    }
    CHECK(CLOSE(nyq_gain, 1.0f, 1e-3f), "HP Nyquist gain ~ 1");

    butter_destroy(&bhp);

    /* ── Invalid parameters → valid = 0, passthrough ──────────────── */

    butter_t binv;

    butter_init(&binv, FILTER_LOWPASS, 0, 2.0f, 0.0f, 20.0f);
    CHECK(binv.valid == 0, "order=0 invalid");
    y = butter_update(&binv, 0.5f);
    CHECK(y == 0.5f, "invalid filter passthrough");

    butter_init(&binv, FILTER_LOWPASS, 2, 0.0f, 0.0f, 20.0f);
    CHECK(binv.valid == 0, "fc=0 invalid");

    butter_init(&binv, FILTER_LOWPASS, 2, 10.0f, 0.0f, 20.0f);
    CHECK(binv.valid == 0, "fc=fs/2 invalid");

    /* ── LP 4th-order → 2 sections ────────────────────────────────── */

    butter_t b4;
    butter_init(&b4, FILTER_LOWPASS, 4, 3.0f, 0.0f, 20.0f);
    CHECK(b4.valid == 1, "4th-order LP init");
    CHECK(b4.sos->num_sections == 2, "4th-order → 2 sections");

    butter_reset(&b4, 1.0f);
    y = butter_update(&b4, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "4th-order LP DC gain ~ 1");

    butter_destroy(&b4);

    /* ── LP 1st-order ─────────────────────────────────────────────── */

    butter_t b1;
    butter_init(&b1, FILTER_LOWPASS, 1, 2.0f, 0.0f, 20.0f);
    CHECK(b1.valid == 1, "1st-order LP init");
    CHECK(b1.sos->num_sections == 1, "1st-order → 1 section");

    butter_reset(&b1, 1.0f);
    y = butter_update(&b1, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "1st-order LP DC gain ~ 1");

    butter_destroy(&b1);

    /* ── BP 2nd-order, fc1=2 Hz, fc2=5 Hz, fs=40 Hz ─────────────────── */

    butter_t bbp;
    butter_init(&bbp, FILTER_BANDPASS, 2, 2.0f, 5.0f, 40.0f);
    CHECK(bbp.valid == 1, "BP init valid");
    CHECK(bbp.sos->num_sections == 2, "BP 2nd-order (4 poles) → 2 sections");

    /* DC gain ≈ 0 */
    butter_reset(&bbp, 1.0f);
    y = butter_update(&bbp, 1.0f);
    CHECK(CLOSE(y, 0.0f, 1e-3f), "BP blocks DC");

    /* Centre frequency gain ≈ 1.0 */
    float f0 = sqrtf(2.0f * 5.0f); /* ~3.16 Hz */
    gn = measure_gain(&bbp, f0, 40.0f, 800);
    CHECK(CLOSE(gn, 1.0f, 0.1f), "BP centre freq gain ~ 1");

    /* Out-of-band attenuation at 20 Hz (Nyquist). */
    gn = measure_gain(&bbp, 20.0f, 40.0f, 800);
    CHECK(gn < 0.15f, "BP Nyquist attenuation");

    butter_destroy(&bbp);

    /* ── BP 1st-order (2 poles after doubling) ──────────────────────── */

    butter_t bbp1;
    butter_init(&bbp1, FILTER_BANDPASS, 1, 3.0f, 6.0f, 40.0f);
    CHECK(bbp1.valid == 1, "BP 1st-order init");
    CHECK(bbp1.sos->num_sections == 1, "BP 1st-order → 1 section");
    butter_destroy(&bbp1);

    /* ── BS 2nd-order, fc1=2 Hz, fc2=5 Hz, fs=40 Hz ─────────────────── */

    butter_t bbs;
    butter_init(&bbs, FILTER_BANDSTOP, 2, 2.0f, 5.0f, 40.0f);
    CHECK(bbs.valid == 1, "BS init valid");
    CHECK(bbs.sos->num_sections == 2, "BS 2nd-order (4 poles) → 2 sections");

    /* DC gain ≈ 1 */
    butter_reset(&bbs, 1.0f);
    y = butter_update(&bbs, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "BS DC gain ~ 1");

    /* Centre frequency attenuation (notch). */
    gn = measure_gain(&bbs, f0, 40.0f, 800);
    CHECK(gn < 0.15f, "BS notch at centre freq");

    /* Nyquist gain ≈ 1.0 */
    butter_reset(&bbs, 0.0f);
    nyq_gain = 0.0f;
    for (int n = 0; n < 200; n++) {
        float x = (n % 2 == 0) ? 1.0f : -1.0f;
        y = butter_update(&bbs, x);
        if (n > 100 && fabsf(y) > nyq_gain) nyq_gain = fabsf(y);
    }
    CHECK(CLOSE(nyq_gain, 1.0f, 1e-3f), "BS Nyquist gain ~ 1");

    butter_destroy(&bbs);

    /* ── BP/BS invalid: fc1 >= fc2 ──────────────────────────────────── */

    butter_init(&bbp, FILTER_BANDPASS, 2, 5.0f, 2.0f, 40.0f);
    CHECK(bbp.valid == 0, "BP fc1>=fc2 invalid");

    butter_init(&bbs, FILTER_BANDSTOP, 2, 5.0f, 5.0f, 40.0f);
    CHECK(bbs.valid == 0, "BS fc1==fc2 invalid");

    /* ── BP/BS invalid: fc2 >= fs/2 ─────────────────────────────────── */

    butter_init(&bbp, FILTER_BANDPASS, 2, 2.0f, 20.0f, 40.0f);
    CHECK(bbp.valid == 0, "BP fc2=fs/2 invalid");

    butter_init(&bbs, FILTER_BANDSTOP, 2, 2.0f, 25.0f, 40.0f);
    CHECK(bbs.valid == 0, "BS fc2>=fs/2 invalid");

    /* ── Report ────────────────────────────────────────────────────── */

    if (failures) {
        fprintf(stderr, "%d test(s) FAILED.\n", failures);
        return EXIT_FAILURE;
    }

    printf("All tests passed.\n");
    return EXIT_SUCCESS;
}
