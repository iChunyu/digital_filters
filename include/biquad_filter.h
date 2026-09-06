/**
 * @file    biquad_filter.h
 * @brief   Second-order IIR filter (Direct Form II / canonical form).
 *
 * Implements the discrete-time transfer function:
 * @f[
 *   H(z) = \frac{b_0 + b_1 z^{-1} + b_2 z^{-2}}
 *                {1   + a_1 z^{-1} + a_2 z^{-2}}
 * @f]
 *
 * The canonical Direct Form II uses only 3 state variables (@p w[0..2]),
 * making it the most memory-efficient realisation for a second-order section.
 */

#ifndef BIQUAD_FILTER_H_
#define BIQUAD_FILTER_H_

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Biquad filter object (Direct Form II).
 *
 * Coefficients are stored normalised so that @p den_z[0] is always 1.0.
 * The state vector @p w[] holds the intermediate values in the canonical
 * realisation:
 *
 * @verbatim
 *   x[n] --->(+) --------> w[n] -----[b0]--->(+)---> y[n]
 *             ^-           |                  ^
 *             |           z^-1                |
 *             |            |                  |
 *             +----[a1]----w[n-1]----[b1]-----+
 *             |            |                  |
 *             |           z^-1                |
 *             |            |                  |
 *             +----[a2]----w[n-2]----[b2]-----+
 * @endverbatim
 */
typedef struct {
    float num_z[3]; /**< Numerator coefficients   b0, b1, b2               */
    float den_z[3]; /**< Denominator coefficients  1.0, a1, a2 (normalised) */
    float w[3];     /**< State variables           w[n], w[n-1], w[n-2]    */
} biquad_filter_t;

/**
 * @brief Convert s-domain biquad coefficients to z-domain via the bilinear
 *        transform.
 *
 * Maps a continuous-time transfer function
 * @f[
 *   H(s) = \frac{B_0 + B_1 s + B_2 s^2}
 *                {A_0 + A_1 s + A_2 s^2}
 * @f]
 * to the discrete-time equivalent
 * @f[
 *   H(z) = \frac{b_0 + b_1 z^{-1} + b_2 z^{-2}}
 *                {a_0 + a_1 z^{-1} + a_2 z^{-2}}
 * @f]
 * using the substitution @f$ s = 2 f_s \frac{1 - z^{-1}}{1 + z^{-1}} @f$.
 *
 * @param[out] num_z  Resulting z-domain numerator   [b0, b1, b2].
 * @param[out] den_z  Resulting z-domain denominator [a0, a1, a2].
 * @param[in]  num_s  s-domain numerator   [B0, B1, B2].
 * @param[in]  den_s  s-domain denominator [A0, A1, A2].
 * @param[in]  fs     Sampling frequency in Hz.
 */
void biquad_c2d_bilinear(float num_z[3], float den_z[3], const float num_s[3],
                        const float den_s[3], float fs);

/**
 * @brief Set a biquad filter to unity pass-through (identity).
 *
 * Initialises the filter so that @f$ H(z) = 1 @f$ — the output equals the
 * input sample-for-sample with no filtering and zero internal state.
 *
 * @param[out] filter  Pointer to the filter object.
 */
void biquad_filter_set_empty(biquad_filter_t *filter);

/**
 * @brief Initialise a biquad filter with given z-domain coefficients.
 *
 * Coefficients are normalised so that @p den_z[0] becomes 1.0.
 * If @p den_z[0] is zero or non-finite, any coefficient is non-finite,
 * the resulting poles lie outside the unit circle (unstable), or any
 * pole radius exceeds 0.99995 (1 − r < 5e-5 — rings for ≥ 10⁴ samples:
 * |a2| > 0.9999 for conjugate pairs, |a1| > 0.99995 for first-order
 * sections, the general root formula for unequal real pairs), the filter
 * is silently replaced by a unity pass-through (identity) and 0 is
 * returned — passthrough is always preferred over divergence.
 *
 * @param[out] filter  Pointer to the filter object.
 * @param[in]  num_z   Numerator coefficients   [b0, b1, b2] in z-domain.
 * @param[in]  den_z   Denominator coefficients [a0, a1, a2] in z-domain.
 * @return             1 if the filter was deployed, 0 if it fell back
 *                     to identity.
 */
uint8_t biquad_filter_init(biquad_filter_t *filter, const float num_z[3],
                           const float den_z[3]);

/**
 * @brief Process one input sample and return the filtered output.
 *
 * This advances the internal state by one time step.
 *
 * Defined static inline (header-only) to eliminate per-section function-call
 * overhead on resource-constrained MCUs.  The output computation is fused
 * with the state update so @p w[] values are loaded only once.
 *
 * @note A non-finite input (NaN/Inf) poisons the state vector and every
 *       subsequent output until the filter is reset.  No guard is performed
 *       here on purpose — it would cost branches on every sample of the MCU
 *       hot path.  Callers feeding untrusted or sensor data should sanitise
 *       at the source.
 *
 * @param[in,out] filter  Pointer to the filter object.
 * @param[in]     input   Current input sample.
 *
 * @return Filtered output sample @f$ y[n] @f$.
 */
static inline float biquad_filter_update(biquad_filter_t *filter, float input)
{
    /* Snapshot current state before shifting (values become w[n-1], w[n-2]). */
    const float w1 = filter->w[0];
    const float w2 = filter->w[1];

    /* Cache coefficients — avoids reloading through pointer on each access. */
    const float a1 = filter->den_z[1];
    const float a2 = filter->den_z[2];
    const float b0 = filter->num_z[0];
    const float b1 = filter->num_z[1];
    const float b2 = filter->num_z[2];

    /* Compute new state: w[n] = x[n] - a1·w[n-1] - a2·w[n-2] */
    const float w0 = input - a1 * w1 - a2 * w2;

    /* Commit state to memory. */
    filter->w[2] = w2;
    filter->w[1] = w1;
    filter->w[0] = w0;

    /* Output: y[n] = b0·w[n] + b1·w[n-1] + b2·w[n-2] */
    return b0 * w0 + b1 * w1 + b2 * w2;
}

/**
 * @brief Return the current output without advancing the state.
 *
 * @param[in] filter  Pointer to the filter object.
 *
 * @return Current output @f$ y[n] @f$.
 */
float biquad_filter_get_output(const biquad_filter_t *filter);

/**
 * @brief Return the current input reconstructed from state.
 *
 * Useful for debugging or cascading — recovers @f$ x[n] @f$ from the
 * internal @p w[] state and denominator coefficients.
 *
 * @param[in] filter  Pointer to the filter object.
 *
 * @return Reconstructed input @f$ x[n] @f$.
 */
float biquad_filter_get_input(const biquad_filter_t *filter);

/**
 * @brief Reset the filter to a steady-state equilibrium.
 *
 * Computes the state vector @p w[] such that a constant input of value
 * @p equilibrium produces the same constant output (DC gain matching).
 *
 * @note If @f$ 1 + a_1 + a_2 = 0 @f$ (e.g. a pure integrator), the
 *       steady-state is undefined; the state is forced to zero in that case.
 *       A non-finite @p equilibrium also forces the state to zero instead
 *       of poisoning it with NaN.
 *
 * @param[in,out] filter      Pointer to the filter object.
 * @param[in]     equilibrium  Constant input value at steady-state.
 */
void biquad_filter_reset(biquad_filter_t *filter, float equilibrium);

/**
 * @brief Process one sample through a cascade of biquad sections.
 *
 * Shared core for the per-order @p _update functions of the higher-order
 * filter families.  Defined static inline so the whole per-sample path is
 * call-free; @p num_sections must be a compile-time literal (the families
 * pass their X-macro section count) so the loop bound folds away.
 *
 * @param[in,out] sections      Cascade of biquad sections.
 * @param[in]     num_sections  Number of sections (compile-time literal).
 * @param[in]     input         Current input sample.
 * @return                      Filtered output sample.
 */
static inline float biquad_cascade_update(biquad_filter_t *sections,
                                          uint8_t num_sections, float input)
{
    float x = input;
    for (uint8_t i = 0; i < num_sections; i++) {
        x = biquad_filter_update(&sections[i], x);
    }
    return x;
}

/**
 * @brief Reset a cascade of biquad sections to steady-state.
 *
 * Shared core for the per-order @p _reset functions of the higher-order
 * filter families.  Each section's steady state is computed for the
 * steady-state output of the preceding section.
 *
 * @param[in,out] sections      Cascade of biquad sections.
 * @param[in]     num_sections  Number of sections (compile-time literal).
 * @param[in]     equilibrium   Constant input value at steady-state.
 */
static inline void biquad_cascade_reset(biquad_filter_t *sections,
                                        uint8_t num_sections, float equilibrium)
{
    float x = equilibrium;
    for (uint8_t i = 0; i < num_sections; i++) {
        biquad_filter_reset(&sections[i], x);
        x = biquad_filter_get_output(&sections[i]);
    }
}


#ifdef __cplusplus
}
#endif

#endif /* BIQUAD_FILTER_H_ */
