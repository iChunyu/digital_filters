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
 * If @p den_z[0] is zero, or if the resulting poles lie outside the unit
 * circle (unstable), the filter is silently replaced by a unity pass-through
 * (identity).
 *
 * @param[out] filter  Pointer to the filter object.
 * @param[in]  num_z   Numerator coefficients   [b0, b1, b2] in z-domain.
 * @param[in]  den_z   Denominator coefficients [a0, a1, a2] in z-domain.
 */
void biquad_filter_init(biquad_filter_t *filter, const float num_z[3],
                        const float den_z[3]);

/**
 * @brief Process one input sample and return the filtered output.
 *
 * This advances the internal state by one time step.
 *
 * @param[in,out] filter  Pointer to the filter object.
 * @param[in]     input   Current input sample.
 *
 * @return Filtered output sample @f$ y[n] @f$.
 */
float biquad_filter_update(biquad_filter_t *filter, float input);

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
 *
 * @param[in,out] filter      Pointer to the filter object.
 * @param[in]     equilibrium  Constant input value at steady-state.
 */
void biquad_filter_reset(biquad_filter_t *filter, float equilibrium);


#ifdef __cplusplus
}
#endif

#endif /* BIQUAD_FILTER_H_ */
