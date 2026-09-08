/**
 * @file    biquad_filter.c
 * @brief   直接 II 型 biquad 滤波器的实现。
 */

#include "biquad_filter.h"
#include <math.h>

/**
 * @brief 三项补偿求和（TwoSum 式）。
 *
 * 裸 f32 计算 1 + a1 + a2 在真余量为 1.0 的若干 ulp 时会恰好消成
 * 0.0f（窄带设计：a1 ≈ −2，a2 ≈ 1 − ε）——把极点远离单位圆的
 * 稳定滤波器误拒，同样也会误算 reset 的分母。把每次加法的
 * 舍入误差折回即可在 f32 代价下恢复余量。
 *
 * @param x  第一加数。
 * @param y  第二加数。
 * @param z  第三加数。
 * @return   补偿后的和 x + y + z。
 */
static float sum3f(float x, float y, float z)
{
    float s = x;
    float c = 0.0f;
    float t = s + y;
    c += (fabsf(s) >= fabsf(y)) ? (s - t) + y : (y - t) + s;
    s = t;
    t = s + z;
    c += (fabsf(s) >= fabsf(z)) ? (s - t) + z : (z - t) + s;
    s = t;
    return s + c;
}

/**
 * @brief 状态向量清零。
 *
 * @param[out] filter  滤波器对象指针。
 */
static void biquad_zero_state(biquad_filter_t *filter)
{
    filter->w[0] = 0.0f;
    filter->w[1] = 0.0f;
    filter->w[2] = 0.0f;
}

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
     * The two Jury sums use compensated summation (see sum3f): naive f32
     * evaluation cancels to exactly 0.0f for narrowband designs whose true
     * residual is a few ulps of 1.0, spuriously rejecting stable filters.
     * A compensated sum of EXACTLY 0.0f, by contrast, means the f32
     * coefficients place a pole exactly on the unit circle (quantization
     * has collapsed it onto z = ±1) — that one is genuinely rejectable,
     * and the margin check below cannot resolve its radius.
     */
    float a1 = filter->den_z[1];
    float a2 = filter->den_z[2];

    if (!(a2 > -1.0f && a2 < 1.0f
          && sum3f(1.0f, a1, a2) > 0.0f
          && sum3f(1.0f, -a1, a2) > 0.0f)) {
        biquad_filter_set_empty(filter);
        return 0;
    }

    /*
     * Stability margin.  The Jury conditions above accept poles arbitrarily
     * close to the unit circle; in f32 a2 = 1 − 2⁻²⁴ passes and the filter
     * rings for millions of samples (~6e-8 contraction per sample).  Reject
     * any pole with radius r > 0.99995 (1 − r < 5e-5).
     *
     * The check is evaluated on the pole radii themselves, not on a2:
     * a2 = r1·r2 is blind to a dominant pole of an unequal real pair
     * (r1 ≈ 1, r2 ≈ 0.85 → a2 ≈ 0.85 sails through), and a one-sided
     * a2 > 0.9999 test misses real pairs of opposite sign (a2 ≈ −0.99995).
     * For a conjugate pair r² = a2 (so r > 0.99995 ⟺ a2 > 0.9999); for
     * real roots r_max = (|a1| + √(a1² − 4·a2)) / 2, rearranged so no
     * sqrtf is needed — the biquad init path must stay callable from
     * bare-metal firmware that links no libm.  The formula also covers
     * first-order sections (a2 = 0 → r_max = |a1|).
     *
     * a1² − 4·a2 is computed with a Dekker-split compensation: wide-band
     * BP/BS designs land near-real pole pairs within ~1e-4 of the unit
     * circle, where the raw f32 discriminant carries ~5e-7 noise while the
     * true value (e.g. −4·(imag part)² ≈ −1e-7) is no larger — the sign
     * flip between the conjugate/real branches would turn the margin
     * decision into a coin toss on legitimate designs.
     */
    float p = a1 * 4097.0f;          /* Dekker split: a1 = hi + lo */
    float hi = p - (p - a1);
    float lo = a1 - hi;
    float sq = a1 * a1;
    float err = ((hi * hi - sq) + 2.0f * hi * lo) + lo * lo;
    float disc = (sq - 4.0f * a2) + err;
    int reject;
    if (disc < 0.0f) {
        reject = a2 > 0.9999f;                 /* conjugate pair: a2 = r² */
    } else {
        /* r_max > 0.99995  ⟺  √disc > 1.9999 − |a1|  (squared). */
        float rhs = 1.9999f - fabsf(a1);
        reject = (rhs < 0.0f) || (disc > rhs * rhs);
    }
    if (reject) {
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
        biquad_zero_state(filter);
        return;
    }

    /*
     * Steady-state: x = const implies w[0]=w[1]=w[2]=w_ss.
     * From the state equation:
     *   x = w_ss + a1*w_ss + a2*w_ss  =>  w_ss = x / (1 + a1 + a2)
     *
     * This is a public API on a fully public struct — the coefficients need
     * not have passed biquad_filter_init().  Compensated summation keeps the
     * f32 evaluation of the denominator from cancelling to exactly 0.0f on
     * init-valid narrowband filters, and the guards below implement the
     * @note contract: 1 + a1 + a2 == 0 (e.g. a pure integrator) has no
     * steady state → force zero; a denormal-tiny denominator would make
     * w_ss overflow to inf and poison every subsequent update with NaN,
     * so force zero there too.
     */
    float denom = sum3f(1.0f, filter->den_z[1], filter->den_z[2]);
    if (denom == 0.0f || !isfinite(denom)) {
        biquad_zero_state(filter);
        return;
    }
    float w_ss = equilibrium / denom;
    if (!isfinite(w_ss)) {
        biquad_zero_state(filter);
        return;
    }
    filter->w[0] = w_ss;
    filter->w[1] = w_ss;
    filter->w[2] = w_ss;
}
