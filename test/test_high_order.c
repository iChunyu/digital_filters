#include "butter_filter.h"
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
    /* ── 8th-order Butterworth LP ───────────────────────────────────── */
    butter_t b8;
    butter_init(&b8, FILTER_LOWPASS, 8, 50.0f, 0.0f, 400.0f);
    CHECK(b8.valid == 1, "8th-order LP valid");
    CHECK(b8.sos->num_sections == 4, "8th-order → 4 sections");
    butter_reset(&b8, 1.0f);
    float y = butter_update(&b8, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "8th-order LP DC gain ~ 1");
    butter_destroy(&b8);

    /* ── 10th-order Butterworth LP ──────────────────────────────────── */
    butter_t b10;
    butter_init(&b10, FILTER_LOWPASS, 10, 50.0f, 0.0f, 400.0f);
    CHECK(b10.valid == 1, "10th-order LP valid");
    CHECK(b10.sos->num_sections == 5, "10th-order → 5 sections");
    butter_reset(&b10, 1.0f);
    y = butter_update(&b10, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "10th-order LP DC gain ~ 1");
    butter_destroy(&b10);

    /* ── 12th-order Butterworth LP ──────────────────────────────────── */
    butter_t b12;
    butter_init(&b12, FILTER_LOWPASS, 12, 50.0f, 0.0f, 400.0f);
    CHECK(b12.valid == 1, "12th-order LP valid");
    CHECK(b12.sos->num_sections == 6, "12th-order → 6 sections");
    butter_reset(&b12, 1.0f);
    y = butter_update(&b12, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "12th-order LP DC gain ~ 1");
    butter_destroy(&b12);

    /* ── 12th-order Butterworth BP (24 effective order) ──────────────── */
    butter_t b12bp;
    butter_init(&b12bp, FILTER_BANDPASS, 12, 20.0f, 50.0f, 400.0f);
    CHECK(b12bp.valid == 1, "12th-order BP valid");
    CHECK(b12bp.sos->num_sections == 12, "12th-order BP → 12 sections");
    butter_reset(&b12bp, 1.0f);
    y = butter_update(&b12bp, 1.0f);
    CHECK(CLOSE(y, 0.0f, 1e-2f), "12th-order BP blocks DC");
    butter_destroy(&b12bp);

    /* ── 12th-order Chebyshev I LP ───────────────────────────────────── */
    cheby1_t c1;
    cheby1_init(&c1, FILTER_LOWPASS, 12, 50.0f, 0.0f, 400.0f, 1.0f);
    CHECK(c1.valid == 1, "12th-order Cheby1 LP valid");
    CHECK(c1.sos->num_sections == 6, "12th-order Cheby1 → 6 sections");
    /* Even-order: DC gain = 1/√(1+ε²) */
    cheby1_reset(&c1, 1.0f);
    y = cheby1_update(&c1, 1.0f);
    {
        float eps = sqrtf(powf(10.0f, 0.1f) - 1.0f);
        float expected = 1.0f / sqrtf(1.0f + eps * eps);
        CHECK(CLOSE(y, expected, 1e-3f), "12th-order Cheby1 LP DC gain correct");
    }
    cheby1_destroy(&c1);

    /* ── 12th-order Chebyshev II LP ──────────────────────────────────── */
    cheby2_t c2;
    cheby2_init(&c2, FILTER_LOWPASS, 12, 50.0f, 0.0f, 400.0f, 40.0f);
    CHECK(c2.valid == 1, "12th-order Cheby2 LP valid");
    CHECK(c2.sos->num_sections == 6, "12th-order Cheby2 → 6 sections");
    cheby2_reset(&c2, 1.0f);
    y = cheby2_update(&c2, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-2f), "12th-order Cheby2 LP DC gain ~ 1");
    cheby2_destroy(&c2);

    /* ── Order 13 → invalid (>12) ────────────────────────────────────── */
    butter_t b13;
    butter_init(&b13, FILTER_LOWPASS, 13, 50.0f, 0.0f, 400.0f);
    CHECK(b13.valid == 0, "13th-order invalid (max 12)");

    /* ── Report ──────────────────────────────────────────────────────── */
    if (failures) {
        fprintf(stderr, "%d test(s) FAILED.\n", failures);
        return EXIT_FAILURE;
    }
    printf("All high-order tests passed.\n");
    return EXIT_SUCCESS;
}
