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
    biquad_c2d_bilinear(num_z, den_z, num_s, den_s, 10.0f);

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

    /* ── non-finite numerator → identity fallback ─────────────────── */

    float nan_num[3] = {NAN, 1.0f, 1.0f};
    biquad_filter_t nf;
    uint8_t rc = biquad_filter_init(&nf, nan_num, num_z);
    CHECK(rc == 0, "NaN numerator rejected");
    y = biquad_filter_update(&nf, 0.75f);
    CHECK(y == 0.75f, "NaN numerator → identity passthrough");

    /* ── infinite leading denominator → identity, not silence ─────────── */

    float inf_den[3] = {INFINITY, 0.0f, 0.0f};
    rc = biquad_filter_init(&nf, num_z, inf_den);
    CHECK(rc == 0, "den_z[0] = Inf rejected");
    y = biquad_filter_update(&nf, 0.5f);
    CHECK(y == 0.5f, "den_z[0] = Inf → identity passthrough");

    /* ── stability margin: near-unit-circle poles rejected ────────────── */

    float marg_num[3] = {1.0f, 0.0f, 0.0f};

    float marg_den_bad[3] = {1.0f, 0.0f, 0.99999f};   /* a2 just below 1 */
    rc = biquad_filter_init(&nf, marg_num, marg_den_bad);
    CHECK(rc == 0, "a2 = 0.99999 rejected (margin)");

    float marg_den_ok[3] = {1.0f, 0.0f, 0.9998f};     /* inside margin */
    rc = biquad_filter_init(&nf, marg_num, marg_den_ok);
    CHECK(rc == 1, "a2 = 0.9998 accepted");

    float marg1_bad[3] = {1.0f, 0.99999f, 0.0f};      /* 1st-order |a1| ~ 1 */
    rc = biquad_filter_init(&nf, marg_num, marg1_bad);
    CHECK(rc == 0, "1st-order |a1| = 0.99999 rejected (margin)");

    /* ── reset with non-finite equilibrium → zero state ───────────────── */

    biquad_filter_t fr;
    biquad_filter_init(&fr, num_z, den_z);
    biquad_filter_reset(&fr, NAN);
    CHECK(fr.w[0] == 0.0f && fr.w[1] == 0.0f && fr.w[2] == 0.0f,
          "reset(NaN) zeroes the state");
    y = biquad_filter_get_output(&fr);
    CHECK(y == 0.0f, "reset(NaN) → zero output");

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
    CHECK(CLOSE(out, 0.2f * x, 1e-6f), "get_output matches y[n] = b0·x");

    /* ── margin symmetry: real pairs of opposite sign rejected too ────────
       A one-sided a2 > 0.9999 check missed a2 ≈ −0.99995 (real poles
       ±0.99997) and product-based checks are blind to a dominant pole of
       an unequal real pair. ──────────────────────────────────────────────── */

    float marg_den_neg[3] = {1.0f, 0.0f, -0.99999f};   /* real pair ±0.999995 */
    rc = biquad_filter_init(&nf, marg_num, marg_den_neg);
    CHECK(rc == 0, "a2 = -0.99999 rejected (real pair, symmetric margin)");

    float dom_num[3] = {1.0f, 0.0f, 0.0f};
    /* poles 0.999999 and 0.85: a2 = 0.85 passes a product-based check,
       and the Jury sum 1 + a1 + a2 = 1e-6 is positive — only the radius
       check on the dominant pole can reject this. */
    float dom_den[3] = {1.0f, -1.849999f, 0.85f};
    rc = biquad_filter_init(&nf, dom_num, dom_den);
    CHECK(rc == 0, "dominant pole 0.999999 rejected (product a2=0.85 blind spot)");

    /* ── reset on a pure integrator (1 + a1 + a2 == 0) ───────────────────
       The @note contract: no steady state exists, the state must be forced
       to zero — never inf/NaN.  biquad_filter_reset is public API on a
       public struct; coefficients need not have passed init. ──────────── */

    biquad_filter_t integ;
    biquad_filter_set_empty(&integ);
    integ.den_z[1] = -1.0f;    /* H(z) = 1 / (1 - z^-1): pole at z = 1 */
    biquad_filter_reset(&integ, 1.0f);
    CHECK(integ.w[0] == 0.0f && integ.w[1] == 0.0f && integ.w[2] == 0.0f,
          "reset(integrator) forces zero state (not inf)");
    y = biquad_filter_update(&integ, 0.5f);
    CHECK(isfinite(y), "integrator update stays finite after guarded reset");

    /* ── non-finite input poisons state (documented hot-path semantics) ──
       The per-sample update deliberately has no NaN guard (branch cost on
       the MCU hot path); a NaN input propagates until reset recovers. ──── */

    biquad_filter_t ns;
    biquad_filter_init(&ns, num, den);
    biquad_filter_reset(&ns, 0.0f);
    y = biquad_filter_update(&ns, 1.0f);
    CHECK(CLOSE(y, 0.2f, 1e-6f), "NaN-poison test: sane state before");
    y = biquad_filter_update(&ns, NAN);
    CHECK(isnan(y), "NaN input → NaN output (no hot-path guard)");
    y = biquad_filter_update(&ns, 1.0f);
    CHECK(isnan(y), "NaN input poisons subsequent samples");
    biquad_filter_reset(&ns, 0.0f);
    y = biquad_filter_update(&ns, 1.0f);
    CHECK(CLOSE(y, 0.2f, 1e-6f), "reset recovers from NaN-poisoned state");

    /* ── report ───────────────────────────────────────────────────── */

    if (failures) {
        fprintf(stderr, "%d test(s) FAILED.\n", failures);
        return EXIT_FAILURE;
    }

    printf("All tests passed.\n");
    return EXIT_SUCCESS;
}
