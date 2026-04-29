#include "sos_filter.h"
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
    /* ── single-section SOS = single biquad ──────────────────────── */

    sos_filter_t s1;
    sos_filter_init(&s1, 1, 10.0f);

    float num[3] = {0.2f, 0.4f, 0.2f};
    float den[3] = {1.0f, 0.0f, 0.0f};
    sos_filter_set_section(&s1, 0, num, den);

    float y = sos_filter_update(&s1, 1.0f);
    y = sos_filter_update(&s1, 1.0f);
    y = sos_filter_update(&s1, 1.0f);
    CHECK(CLOSE(y, 0.8f, 1e-5f), "single-section SOS DC gain");

    sos_filter_destroy(&s1);

    /* ── identity cascade   ─────────────────────────────────────── */

    sos_filter_t id2;
    sos_filter_init(&id2, 2, 10.0f);

    y = sos_filter_update(&id2, 0.5f);
    CHECK(y == 0.5f, "two-section identity DC");
    y = sos_filter_update(&id2, -0.25f);
    CHECK(y == -0.25f, "two-section identity negative");

    sos_filter_destroy(&id2);

    /* ── set_section with out-of-bounds index (must not crash) ──── */

    sos_filter_t s3;
    sos_filter_init(&s3, 1, 10.0f);
    sos_filter_set_section(&s3, 5, num, den);  /* OOB, ignored */
    y = sos_filter_update(&s3, 0.5f);
    CHECK(y == 0.5f, "OOB set_section leaves identity");

    sos_filter_destroy(&s3);

    /* ── reset propagates equilibrium through cascade ───────────── */
    /* Section 0: DC gain = 2.0, Section 1: DC gain = 0.5 */
    /* Overall DC gain = 1.0 */

    float num0[3] = {2.0f, 0.0f, 0.0f};
    float den0[3] = {1.0f, 0.0f, 0.0f};
    float num1[3] = {0.5f, 0.0f, 0.0f};
    float den1[3] = {1.0f, 0.0f, 0.0f};

    sos_filter_t s4;
    sos_filter_init(&s4, 2, 10.0f);
    sos_filter_set_section(&s4, 0, num0, den0);
    sos_filter_set_section(&s4, 1, num1, den1);
    sos_filter_reset(&s4, 1.0f);

    y = sos_filter_update(&s4, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-5f), "reset propagates equilibrium through cascade");

    sos_filter_destroy(&s4);

    /* ── get_input / get_output consistency ──────────────────────── */

    sos_filter_t s5;
    sos_filter_init(&s5, 1, 10.0f);
    sos_filter_set_section(&s5, 0, num, den);

    float x = 0.7f;
    sos_filter_update(&s5, x);
    float out = sos_filter_get_output(&s5);
    float in  = sos_filter_get_input(&s5);
    CHECK(CLOSE(in, x, 1e-6f), "sos get_input reconstructs x[n]");

    sos_filter_destroy(&s5);

    /* ── report ───────────────────────────────────────────────────── */

    if (failures) {
        fprintf(stderr, "%d test(s) FAILED.\n", failures);
        return EXIT_FAILURE;
    }

    printf("All tests passed.\n");
    return EXIT_SUCCESS;
}
