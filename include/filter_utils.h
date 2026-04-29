#ifndef FILTER_UTILS_H_
#define FILTER_UTILS_H_

#include <stdint.h>

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
 * @brief Convert z-domain pole/zero arrays to second-order section coefficients.
 *
 * Pairs poles with nearest zeros using the "most unfavorable pole first"
 * algorithm, produces biquad coefficients, and normalises each section for
 * unity gain at the reference frequency (DC for lowpass, Nyquist for highpass).
 *
 * @param[in]  zeros  Array of n z-domain zeros.
 * @param[in]  poles  Array of n z-domain poles.
 * @param[in]  n      Number of poles (must equal number of zeros).
 * @param[out] sos    SOS matrix with ceil(n/2) rows, each [b0,b1,b2, 1,a1,a2].
 *                    Caller must allocate ceil(n/2) rows.
 * @param[in]  type   Filter type (for gain normalisation reference).
 * @return            Number of SOS sections = ceil(n/2).
 */
uint8_t zpk2sos(const complex_t *zeros, const complex_t *poles, uint8_t n,
                float (*sos)[6], filter_type_e type);

#endif /* FILTER_UTILS_H_ */
