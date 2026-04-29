/**
 * @file    biquad_filter.c
 * @brief   Implementation of the Direct Form II biquad filter.
 */

#include "biquad_filter.h"

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

void biquad_c2d_bilnear(float num_z[3], float den_z[3], const float num_s[3],
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

void biquad_filter_init(biquad_filter_t *filter, const float num_z[3],
                        const float den_z[3])
{
    /* Reject zero leading denominator — would cause division by zero. */
    if (den_z[0] == 0.0f) {
        biquad_filter_set_empty(filter);
        return;
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
        return;
    }

    filter->w[0] = 0.0f;
    filter->w[1] = 0.0f;
    filter->w[2] = 0.0f;
}

float biquad_filter_update(biquad_filter_t *filter, float input)
{
    /* Shift state */
    filter->w[2] = filter->w[1];
    filter->w[1] = filter->w[0];

    /* Compute new state: w[n] = x[n] - a1*w[n-1] - a2*w[n-2] */
    filter->w[0] = input
                   - filter->den_z[1] * filter->w[1]
                   - filter->den_z[2] * filter->w[2];

    /* Output: y[n] = b0*w[n] + b1*w[n-1] + b2*w[n-2] */
    return biquad_filter_get_output(filter);
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
