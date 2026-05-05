#ifndef FILTER_UTILS_H_
#define FILTER_UTILS_H_

#include <math.h>
#include <stdint.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    FILTER_LOWPASS,
    FILTER_HIGHPASS,
    FILTER_BANDPASS,
    FILTER_BANDSTOP
} filter_type_e;

typedef struct {
    float re;
    float im;
} complex_t;

/**
 * @brief Pre-warp a digital cutoff frequency to its analog equivalent.
 *
 * Uses the bilinear transform frequency mapping:
 *   f_analog = fs/pi * tan(pi * f_digital / fs)
 *
 * @param fd  Desired digital cutoff frequency in Hz (0 < fd < fs/2).
 * @param fs  Sampling frequency in Hz.
 * @return    Equivalent analog cutoff frequency in Hz.
 */
float prewarp(float fd, float fs);

/**
 * @brief Lowpass analog frequency transform: scale poles and zeros by wc.
 *
 * @param[in,out] poles  Array of np complex poles.
 * @param[in]     np     Number of poles.
 * @param[in,out] zeros  Array of nz complex zeros.
 * @param[in]     nz     Number of zeros.
 * @param[in]     wc     Analog cutoff frequency (rad/s, i.e. 2*pi*prewarp).
 */
void analog_lp_transform(complex_t *poles, uint8_t np,
                         complex_t *zeros, uint8_t nz, float wc);

/**
 * @brief Highpass analog frequency transform: invert poles/zeros and scale.
 *
 * Each pole p becomes wc/p, each zero z becomes wc/z.
 * The caller must append np - nz zeros at the origin (0, 0) after this call
 * to account for the infinite zeros in the prototype.
 *
 * @param[in,out] poles  Array of np complex poles.
 * @param[in]     np     Number of poles.
 * @param[in,out] zeros  Array of nz complex zeros.
 * @param[in]     nz     Number of zeros.
 * @param[in]     wc     Analog cutoff frequency (rad/s).
 */
void analog_hp_transform(complex_t *poles, uint8_t np,
                         complex_t *zeros, uint8_t nz, float wc);

/**
 * @brief Bilinear transform: map s-domain poles/zeros to z-domain.
 *
 * Applies z = (2*fs + s) / (2*fs - s) to each element in-place.
 * The caller must append (np - nz) zeros at z = -1 afterward to account
 * for the excess poles.
 *
 * @param[in,out] zp  Array of n complex poles or zeros (modified in-place).
 * @param[in]     n   Number of elements.
 * @param[in]     fs  Sampling frequency in Hz.
 */
void bilinear_transform(complex_t *zp, uint8_t n, float fs);

/**
 * @brief Analog bandpass frequency transform.
 *
 * Substitutes s → (s² + ω₀²) / (ξ·s), doubling the order.
 * Each pole and zero splits into two via the quadratic formula.
 * The function appends (np_old − nz_old) zeros at the origin (0, 0) to
 * account for infinite prototype zeros.
 *
 * @param[in,out] poles  Array with capacity ≥ 2·*np. Transformed in-place.
 * @param[in,out] np     Input: number of prototype poles. Output: 2 × input.
 * @param[in,out] zeros  Array with capacity ≥ *np + *nz. Transformed in-place.
 * @param[in,out] nz     Input: number of prototype zeros. Output: *np_old + *nz_old.
 * @param[in]     w0     Analog center frequency ω₀ = √(ω₁·ω₂) (rad/s).
 * @param[in]     xi     Bandwidth ξ = ω₂ − ω₁ (rad/s).
 */
void analog_bp_transform(complex_t *poles, uint8_t *np,
                         complex_t *zeros, uint8_t *nz,
                         float w0, float xi);

/**
 * @brief Analog bandstop frequency transform.
 *
 * Substitutes s → ξ·s / (s² + ω₀²), doubling the order.
 * The function appends 2·(np_old − nz_old) zeros at ±jω₀ to account
 * for infinite prototype zeros.
 *
 * @param[in,out] poles  Array with capacity ≥ 2·*np. Transformed in-place.
 * @param[in,out] np     Input: number of prototype poles. Output: 2 × input.
 * @param[in,out] zeros  Array with capacity ≥ 2·*np. Transformed in-place.
 * @param[in,out] nz     Input: number of prototype zeros. Output: 2 × *np_old.
 * @param[in]     w0     Analog center frequency ω₀ = √(ω₁·ω₂) (rad/s).
 * @param[in]     xi     Bandwidth ξ = ω₂ − ω₁ (rad/s).
 */
void analog_bs_transform(complex_t *poles, uint8_t *np,
                         complex_t *zeros, uint8_t *nz,
                         float w0, float xi);

/**
 * @brief Compute the gain adjustment for an analog LP→LP frequency transform.
 *
 * k_lp = k * wo^degree  (scipy lp2lp_zpk convention).
 *
 * @param k       Current system gain.
 * @param wo      Target cutoff angular frequency (rad/s).
 * @param degree  Relative degree = np − nz (number of excess poles).
 * @return        Adjusted gain.
 */
float zpk_lp_gain(float k, float wo, uint8_t degree);

/**
 * @brief Compute the gain adjustment for an analog LP→HP or LP→BS transform.
 *
 * k = k · Re(∏(−z) / ∏(−p))  evaluated on the prototype (pre-transform)
 * poles and zeros.  When nz = 0 the empty product ∏(−z) is 1.
 *
 * @param k   Current system gain.
 * @param z   Prototype zeros (pre-transform), nz elements.
 * @param nz  Number of prototype zeros.
 * @param p   Prototype poles (pre-transform), np elements.
 * @param np  Number of prototype poles.
 * @return    Adjusted gain.
 */
float zpk_hp_bs_gain(float k, const complex_t *z, uint8_t nz,
                     const complex_t *p, uint8_t np);

/**
 * @brief Compute the gain adjustment for an analog LP→BP frequency transform.
 *
 * k_bp = k * bw^degree  (scipy lp2bp_zpk convention).
 *
 * @param k       Current system gain.
 * @param bw      Bandwidth ξ (rad/s).
 * @param degree  Relative degree = np − nz.
 * @return        Adjusted gain.
 */
float zpk_bp_gain(float k, float bw, uint8_t degree);

/**
 * @brief Compute the gain adjustment for the bilinear transform.
 *
 * k_z = k · Re(∏(K − z) / ∏(K − p))  where K = 2·fs and z,p are the
 * analog-domain (s-plane) zeros and poles BEFORE the bilinear substitution.
 *
 * @param k   Current system gain.
 * @param z   Analog-domain zeros, nz elements.
 * @param nz  Number of analog zeros.
 * @param p   Analog-domain poles, np elements.
 * @param np  Number of analog poles.
 * @param K   Bilinear constant = 2·fs.
 * @return    Adjusted gain.
 */
float bilinear_zpk_gain(float k, const complex_t *z, uint8_t nz,
                         const complex_t *p, uint8_t np, float K);

/**
 * @brief Convert z-domain pole/zero arrays to second-order section coefficients.
 *
 * Pairs poles with nearest zeros using the "most unfavorable pole first"
 * algorithm and produces biquad coefficients.  Each section is built with
 * unity numerator gain; the global system gain @p k is applied to the
 * numerator of the first section only (matching scipy zpk2sos convention).
 *
 * @param[in]  zeros    Array of n z-domain zeros.
 * @param[in]  poles    Array of n z-domain poles.
 * @param[in]  n        Number of poles (must equal number of zeros).
 * @param[out] sos      SOS matrix with ceil(n/2) rows, each [b0,b1,b2, 1,a1,a2].
 *                      Caller must allocate ceil(n/2) rows.
 * @param[in]  k        Overall system gain applied to sos[0] numerator.
 * @return              Number of SOS sections = ceil(n/2).
 */
uint8_t zpk2sos(const complex_t *zeros, const complex_t *poles, uint8_t n,
                float (*sos)[6], float k);

#ifdef __cplusplus
}
#endif

#endif /* FILTER_UTILS_H_ */
