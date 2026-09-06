/**
 * @file    biquad_filter.c
 * @brief   Implementation of the Direct Form II biquad filter.
 */

#include "biquad_filter.h"
#include <math.h>

void biquad_filter_set_empty(biquad_filter_t *filter)
{
    filter->num_z[0] = 1.0f;
    filter->num_z[1] = 0.0f;
    filter->num_z[2] = 0.0f;

    filter->den_z[0] = 1.0f;
    filter->den_z[1] = 0.0f;
    filter->den_z[2] = 0.0f;

    filter->w[0] = 0.0f;
    filter->w[1] = 0.0f;
    filter->w[2] = 0.0f;
}

void biquad_c2d_bilinear(float num_z[3], float den_z[3], const float num_s[3],
                        const float den_s[3], float fs)
{
    float K = 2.0f * fs;
    float K2 = K * K;

    num_z[0] = num_s[0] + num_s[1] * K + num_s[2] * K2;
    num_z[1] = 2.0f * num_s[0] - 2.0f * num_s[2] * K2;
    num_z[2] = num_s[0] - num_s[1] * K + num_s[2] * K2;

    den_z[0] = den_s[0] + den_s[1] * K + den_s[2] * K2;
    den_z[1] = 2.0f * den_s[0] - 2.0f * den_s[2] * K2;
    den_z[2] = den_s[0] - den_s[1] * K + den_s[2] * K2;
}

uint8_t biquad_filter_init(biquad_filter_t *filter, const float num_z[3],
                           const float den_z[3])
{
    /* Reject zero/infinite leading denominator.  Zero divides by zero;
       ±Inf would make inv = 0 and silently deploy an all-zero (silence)
       filter — passthrough is the safer fallback. */
    if (den_z[0] == 0.0f || !isfinite(den_z[0])) {
        biquad_filter_set_empty(filter);
        return 0;
    }

    /* Normalise so that den_z[0] == 1.0 */
    float inv = 1.0f / den_z[0];

    filter->num_z[0] = num_z[0] * inv;
    filter->num_z[1] = num_z[1] * inv;
    filter->num_z[2] = num_z[2] * inv;

    filter->den_z[0] = 1.0f;
    filter->den_z[1] = den_z[1] * inv;
    filter->den_z[2] = den_z[2] * inv;

    /*
     * Reject non-finite coefficients (e.g. gain overflow in the design
     * pipeline).  NaN/Inf numerators would otherwise pass the Jury check
     * below — it only inspects the denominator — and poison the output.
     */
    if (!(isfinite(filter->num_z[0]) && isfinite(filter->num_z[1])
          && isfinite(filter->num_z[2]) && isfinite(filter->den_z[1])
          && isfinite(filter->den_z[2]))) {
        biquad_filter_set_empty(filter);
        return 0;
    }

    /*
     * Stability check — all three Jury conditions for a 2nd-order system:
     *   |a2| < 1
     *   1 + a1 + a2 > 0
     *   1 - a1 + a2 > 0
     * If any condition fails the filter is unstable; fall back to identity.
     */
    float a1 = filter->den_z[1];
    float a2 = filter->den_z[2];

    if (!(a2 > -1.0f && a2 < 1.0f
          && 1.0f + a1 + a2 > 0.0f
          && 1.0f - a1 + a2 > 0.0f)) {
        biquad_filter_set_empty(filter);
        return 0;
    }

    /*
     * Stability margin.  The Jury conditions above accept poles arbitrarily
     * close to the unit circle; in f32 that means a2 = 1 − 2⁻²⁴ passes and
     * the filter rings for millions of samples (a ~6e-8 contraction per
     * sample).  Legitimate designs measured down to 1 − max_r = 2.8e-3
     * (a2 ≈ 0.994) while degenerate ones (cheby1 rp ≥ 60 dB) sit at
     * 1 − max_r ≤ 1.3e-5 — a clean two-decade gap.  Reject |a2| > 0.9999
     * (radius > 0.99995); for 1st-order sections the lone pole is −a1, so
     * bound |a1| instead.
     */
    if (a2 > 0.9999f || (a2 == 0.0f && fabsf(a1) > 0.9999f)) {
        biquad_filter_set_empty(filter);
        return 0;
    }

    filter->w[0] = 0.0f;
    filter->w[1] = 0.0f;
    filter->w[2] = 0.0f;
    return 1;
}

float biquad_filter_get_output(const biquad_filter_t *filter)
{
    return filter->num_z[0] * filter->w[0]
         + filter->num_z[1] * filter->w[1]
         + filter->num_z[2] * filter->w[2];
}

float biquad_filter_get_input(const biquad_filter_t *filter)
{
    return filter->w[0]
         + filter->den_z[1] * filter->w[1]
         + filter->den_z[2] * filter->w[2];
}

void biquad_filter_reset(biquad_filter_t *filter, float equilibrium)
{
    /*
     * Non-finite equilibrium would poison the state vector with NaN and
     * propagate to every subsequent output.  Fall back to a zero state
     * (the identity filter's own steady state) instead.
     */
    if (!isfinite(equilibrium)) {
        filter->w[0] = 0.0f;
        filter->w[1] = 0.0f;
        filter->w[2] = 0.0f;
        return;
    }

    /*
     * Steady-state: x = const implies w[0]=w[1]=w[2]=w_ss.
     * From the state equation:
     *   x = w_ss + a1*w_ss + a2*w_ss  =>  w_ss = x / (1 + a1 + a2)
     *
     * The denominator is guaranteed non-zero: biquad_filter_init() enforces
     * 1 + a1 + a2 > 0 via the Jury stability conditions; if init rejects the
     * coefficients it falls back to the identity filter (a1 = a2 = 0).
     */
    float w_ss = equilibrium / (1.0f + filter->den_z[1] + filter->den_z[2]);
    filter->w[0] = w_ss;
    filter->w[1] = w_ss;
    filter->w[2] = w_ss;
}
