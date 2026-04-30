#include "static_butter_filter.h"
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
static float measure_gain(biquad_filter_t *sections, uint8_t ns, uint8_t valid,
                          float freq, float fs, int steps)
{
    if (!valid) return 0.0f;

    float x = 0.0f;
    for (uint8_t i = 0; i < ns; i++) {
        biquad_filter_reset(&sections[i], x);
        x = biquad_filter_get_output(&sections[i]);
    }

    float max_out = 0.0f;
    for (int n = 0; n < steps; n++) {
        float t = (float)n / fs;
        float in_val = sinf(2.0f * (float)M_PI * freq * t);
        x = in_val;
        for (uint8_t i = 0; i < ns; i++) {
            x = biquad_filter_update(&sections[i], x);
        }
        if (n > steps / 2) {
            if (fabsf(x) > max_out) {
                max_out = fabsf(x);
            }
        }
    }
    return max_out;
}

int main(void)
{
    float y;

    /* ── LP 2nd order, fc=2 Hz, fs=20 Hz ──────────────────────────────── */

    butter_lp_2nd_t blp;
    butter_lp_2nd_init(&blp, 2.0f, 20.0f);
    CHECK(blp.valid == 1, "LP 2nd init valid");
    CHECK(blp.type == FILTER_LOWPASS, "LP 2nd type");
    CHECK(blp.order == 2, "LP 2nd order");
    CHECK(blp.num_sections == 1, "LP 2nd → 1 section");

    /* DC gain ≈ 1.0 */
    butter_lp_2nd_reset(&blp, 1.0f);
    y = butter_lp_2nd_update(&blp, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "LP 2nd DC gain ~ 1");

    /* Attenuation at Nyquist (10 Hz) */
    float gn = measure_gain(blp.sections, blp.num_sections, blp.valid, 10.0f, 20.0f, 400);
    CHECK(gn < 0.15f, "LP 2nd Nyquist attenuation");

    /* ── HP 2nd order, fc=5 Hz, fs=40 Hz ──────────────────────────────── */

    butter_hp_2nd_t bhp;
    butter_hp_2nd_init(&bhp, 5.0f, 40.0f);
    CHECK(bhp.valid == 1, "HP 2nd init valid");
    CHECK(bhp.num_sections == 1, "HP 2nd → 1 section");

    /* DC gain ≈ 0 */
    butter_hp_2nd_reset(&bhp, 1.0f);
    y = butter_hp_2nd_update(&bhp, 1.0f);
    CHECK(CLOSE(y, 0.0f, 1e-3f), "HP 2nd blocks DC");

    /* Nyquist gain ≈ 1.0 */
    butter_hp_2nd_reset(&bhp, 0.0f);
    float nyq_gain = 0.0f;
    for (int n = 0; n < 200; n++) {
        float x = (n % 2 == 0) ? 1.0f : -1.0f;
        y = butter_hp_2nd_update(&bhp, x);
        if (n > 100 && fabsf(y) > nyq_gain) nyq_gain = fabsf(y);
    }
    CHECK(CLOSE(nyq_gain, 1.0f, 1e-3f), "HP 2nd Nyquist gain ~ 1");

    /* ── Invalid parameters → valid = 0, passthrough ──────────────────── */

    butter_lp_2nd_t binv;

    /* fc=0 invalid */
    binv.valid = 0;  /* suppress uninit warning in macro */
    butter_lp_2nd_init(&binv, 0.0f, 20.0f);
    CHECK(binv.valid == 0, "LP 2nd fc=0 invalid");
    y = butter_lp_2nd_update(&binv, 0.5f);
    CHECK(y == 0.5f, "invalid filter passthrough");

    /* fc=fs/2 invalid */
    butter_lp_2nd_init(&binv, 10.0f, 20.0f);
    CHECK(binv.valid == 0, "LP 2nd fc=fs/2 invalid");

    /* ── LP 4th-order → 2 sections ────────────────────────────────────── */

    butter_lp_4th_t b4;
    butter_lp_4th_init(&b4, 3.0f, 20.0f);
    CHECK(b4.valid == 1, "LP 4th init valid");
    CHECK(b4.num_sections == 2, "LP 4th → 2 sections");

    butter_lp_4th_reset(&b4, 1.0f);
    y = butter_lp_4th_update(&b4, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "LP 4th DC gain ~ 1");

    /* ── LP 1st-order ─────────────────────────────────────────────────── */

    butter_lp_1st_t b1;
    butter_lp_1st_init(&b1, 2.0f, 20.0f);
    CHECK(b1.valid == 1, "LP 1st init valid");
    CHECK(b1.num_sections == 1, "LP 1st → 1 section");

    butter_lp_1st_reset(&b1, 1.0f);
    y = butter_lp_1st_update(&b1, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "LP 1st DC gain ~ 1");

    /* ── HP 3rd-order → 2 sections ────────────────────────────────────── */

    butter_hp_3rd_t bhp3;
    butter_hp_3rd_init(&bhp3, 5.0f, 40.0f);
    CHECK(bhp3.valid == 1, "HP 3rd init valid");
    CHECK(bhp3.num_sections == 2, "HP 3rd → 2 sections");

    butter_hp_3rd_reset(&bhp3, 1.0f);
    y = butter_hp_3rd_update(&bhp3, 1.0f);
    CHECK(CLOSE(y, 0.0f, 1e-3f), "HP 3rd blocks DC");

    /* ── BP 2nd order, fc1=2 Hz, fc2=5 Hz, fs=40 Hz ──────────────────── */

    butter_bp_2nd_t bbp;
    butter_bp_2nd_init(&bbp, 2.0f, 5.0f, 40.0f);
    CHECK(bbp.valid == 1, "BP 2nd init valid");
    CHECK(bbp.num_sections == 2, "BP 2nd → 2 sections");

    /* DC gain ≈ 0 */
    butter_bp_2nd_reset(&bbp, 1.0f);
    y = butter_bp_2nd_update(&bbp, 1.0f);
    CHECK(CLOSE(y, 0.0f, 1e-3f), "BP 2nd blocks DC");

    /* Centre frequency gain ≈ 1.0 */
    float f0 = sqrtf(2.0f * 5.0f);
    gn = measure_gain(bbp.sections, bbp.num_sections, bbp.valid, f0, 40.0f, 800);
    CHECK(CLOSE(gn, 1.0f, 0.1f), "BP 2nd centre freq gain ~ 1");

    /* Out-of-band attenuation at Nyquist */
    gn = measure_gain(bbp.sections, bbp.num_sections, bbp.valid, 20.0f, 40.0f, 800);
    CHECK(gn < 0.15f, "BP 2nd Nyquist attenuation");

    /* ── BP 1st-order → 1 section ─────────────────────────────────────── */

    butter_bp_1st_t bbp1;
    butter_bp_1st_init(&bbp1, 3.0f, 6.0f, 40.0f);
    CHECK(bbp1.valid == 1, "BP 1st init valid");
    CHECK(bbp1.num_sections == 1, "BP 1st → 1 section");

    /* ── BS 2nd order, fc1=2 Hz, fc2=5 Hz, fs=40 Hz ──────────────────── */

    butter_bs_2nd_t bbs;
    butter_bs_2nd_init(&bbs, 2.0f, 5.0f, 40.0f);
    CHECK(bbs.valid == 1, "BS 2nd init valid");
    CHECK(bbs.num_sections == 2, "BS 2nd → 2 sections");

    /* DC gain ≈ 1 */
    butter_bs_2nd_reset(&bbs, 1.0f);
    y = butter_bs_2nd_update(&bbs, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "BS 2nd DC gain ~ 1");

    /* Centre frequency notch */
    gn = measure_gain(bbs.sections, bbs.num_sections, bbs.valid, f0, 40.0f, 800);
    CHECK(gn < 0.15f, "BS 2nd notch at centre freq");

    /* Nyquist gain ≈ 1.0 */
    butter_bs_2nd_reset(&bbs, 0.0f);
    nyq_gain = 0.0f;
    for (int n = 0; n < 200; n++) {
        float x = (n % 2 == 0) ? 1.0f : -1.0f;
        y = butter_bs_2nd_update(&bbs, x);
        if (n > 100 && fabsf(y) > nyq_gain) nyq_gain = fabsf(y);
    }
    CHECK(CLOSE(nyq_gain, 1.0f, 1e-3f), "BS 2nd Nyquist gain ~ 1");

    /* ── BP/BS invalid: fc1 >= fc2 ────────────────────────────────────── */

    butter_bp_2nd_init(&bbp, 5.0f, 2.0f, 40.0f);
    CHECK(bbp.valid == 0, "BP fc1>=fc2 invalid");

    butter_bs_2nd_t bbs_inv;
    butter_bs_2nd_init(&bbs_inv, 5.0f, 5.0f, 40.0f);
    CHECK(bbs_inv.valid == 0, "BS fc1==fc2 invalid");

    /* ── BP/BS invalid: fc2 >= fs/2 ───────────────────────────────────── */

    butter_bp_2nd_init(&bbp, 2.0f, 20.0f, 40.0f);
    CHECK(bbp.valid == 0, "BP fc2=fs/2 invalid");

    butter_bs_2nd_init(&bbs, 2.0f, 25.0f, 40.0f);
    CHECK(bbs.valid == 0, "BS fc2>=fs/2 invalid");

    /* ── LP 8th-order → 4 sections ────────────────────────────────────── */

    butter_lp_8th_t b8;
    butter_lp_8th_init(&b8, 3.0f, 20.0f);
    CHECK(b8.valid == 1, "LP 8th init valid");
    CHECK(b8.num_sections == 4, "LP 8th → 4 sections");

    butter_lp_8th_reset(&b8, 1.0f);
    y = butter_lp_8th_update(&b8, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "LP 8th DC gain ~ 1");

    /* ── Struct sizes ─────────────────────────────────────────────────── */

    CHECK(sizeof(butter_lp_3rd_t) > sizeof(butter_lp_1st_t),
          "LP 3rd bigger than LP 1st");
    CHECK(sizeof(butter_bp_8th_t) > sizeof(butter_bp_1st_t),
          "BP 8th bigger than BP 1st");

    /* ── Report ───────────────────────────────────────────────────────── */

    if (failures) {
        fprintf(stderr, "%d test(s) FAILED.\n", failures);
        return EXIT_FAILURE;
    }

    printf("All tests passed.\n");
    return EXIT_SUCCESS;
}
