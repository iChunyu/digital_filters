#include "biquad_filter.h"
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
    /* ── identity / pass-through ─────────────────────────────────── */

    biquad_filter_t id;
    biquad_filter_set_empty(&id);

    float y = biquad_filter_update(&id, 0.5f);
    CHECK(y == 0.5f, "identity DC");

    y = biquad_filter_update(&id, -0.25f);
    CHECK(y == -0.25f, "identity negative");

    /* ── known filter: 2nd-order LPF  ────────────────────────────── */
    /* s-domain: H(s) = 1 / (s^2 + sqrt(2)*s + 1)   (Butterworth, fc=1 Hz) */
    float num_s[3] = {1.0f, 0.0f, 0.0f};
    float den_s[3] = {1.0f, 1.41421356f, 1.0f};

    float num_z[3], den_z[3];
    biquad_c2d_bilnear(num_z, den_z, num_s, den_s, 10.0f);

    biquad_filter_t lpf;
    biquad_filter_init(&lpf, num_z, den_z);

    /* DC gain should be 1.0 */
    biquad_filter_reset(&lpf, 1.0f);
    y = biquad_filter_update(&lpf, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "LPF DC gain ~ 1");

    /* ── zero-coefficient reset  ──────────────────────────────────── */

    float zero_den[3] = {0.0f, 1.0f, 1.0f};
    biquad_filter_init(&lpf, num_z, zero_den);
    y = biquad_filter_update(&lpf, 0.5f);
    CHECK(y == 0.5f, "zero den_z[0] falls back to identity");

    /* ── get_output / get_input consistency  ──────────────────────── */

    float num[3] = {0.2f, 0.4f, 0.2f};
    float den[3] = {1.0f, 0.0f, 0.0f};
    biquad_filter_t f;
    biquad_filter_init(&f, num, den);

    float x = 0.3f;
    biquad_filter_update(&f, x);
    float out = biquad_filter_get_output(&f);
    float in  = biquad_filter_get_input(&f);
    CHECK(CLOSE(in, x, 1e-6f), "get_input reconstructs x[n]");

    /* ── report ───────────────────────────────────────────────────── */

    if (failures) {
        fprintf(stderr, "%d test(s) FAILED.\n", failures);
        return EXIT_FAILURE;
    }

    printf("All tests passed.\n");
    return EXIT_SUCCESS;
}
