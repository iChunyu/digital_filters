#ifndef FILTER_UTILS_H_
#define FILTER_UTILS_H_

#include <math.h>
#include <stdint.h>
#include "biquad_filter.h"

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
 * @brief Compute the combined gain adjustment for an LP/BP frequency
 *        transform followed by the bilinear transform.
 *
 * k = k · s^degree · Re(∏(K − z) / ∏(K − p))  with K = 2·fs.
 *
 * The s factors are interleaved with the (K − p) divisions so no
 * intermediate overflows f32 — computing s^degree standalone overflows
 * for near-Nyquist designs (wc^8 ≈ FLT_MAX) before the ∏(K − p) factors of
 * the same magnitude can cancel it.
 *
 * @param k       Current system gain.
 * @param s       Scale factor per unit of degree: wc (rad/s) for LP,
 *                bandwidth ξ (rad/s) for BP.
 * @param degree  Relative degree = np − nz of the prototype.
 * @param z       Analog-domain zeros (post frequency transform), nz elements.
 * @param nz      Number of analog zeros.
 * @param p       Analog-domain poles (post frequency transform), np elements.
 * @param np      Number of analog poles.
 * @param K       Bilinear constant = 2·fs.
 * @return        Adjusted gain.
 */
float bilinear_zpk_gain_scaled(float k, float s, uint8_t degree,
                               const complex_t *z, uint8_t nz,
                               const complex_t *p, uint8_t np, float K);

/**
 * @brief Convert z-domain pole/zero arrays to second-order section coefficients.
 *
 * Pairs poles with nearest zeros using the "most unfavorable pole first"
 * algorithm and produces biquad coefficients.  Each section is built with
 * unity numerator gain; the global system gain @p k is applied to the
 * numerator of the first section only (matching scipy zpk2sos convention).
 *
 * Works directly on the caller's mutable pole/zero arrays without making
 * internal copies (saves ~256 bytes of stack — significant on MCUs where
 * zpk2sos is called deep in the init call chain); the arrays are read but
 * not modified, ownership tracking uses an internal used[] bitmap.
 *
 * Fail-closed: returns 0 if the pole/zero sets are unbalanced, any element
 * is left unpaired (e.g. a synthesised conjugate), or the sections' roots
 * do not reproduce the input pole/zero multisets (cross-pair misclaim).
 * Callers must treat 0 as "design failed → deploy passthrough".
 *
 * @param[in]  zeros    Array of n z-domain zeros.
 * @param[in]  poles    Array of n z-domain poles.
 * @param[in]  n        Number of poles (must equal number of zeros).
 * @param[out] sos      SOS matrix with ceil(n/2) rows, each [b0,b1,b2, 1,a1,a2].
 *                      Caller must allocate ceil(n/2) rows.
 * @param[in]  k        Overall system gain applied to sos[0] numerator.
 * @return              Number of SOS sections = ceil(n/2), or 0 on failure.
 */
uint8_t zpk2sos(complex_t *zeros, complex_t *poles, uint8_t n,
                float (*sos)[6], float k);

/**
 * @brief Shared IIR design pipeline: analog prototype → frequency
 *        transform → gain folding → bilinear → zpk2sos → deployment.
 *
 * Used by the Butterworth and Chebyshev families (identical pipeline,
 * different prototypes): Butterworth passes its ROM pole table with
 * k = 1, nz = 0; Chebyshev computes its prototype at runtime and passes
 * its own k and finite zeros.  Keeping the pipeline in ONE place means
 * every fail-closed chokepoint (k finite/nonzero, section-count cap,
 * per-section init, gain gate) exists exactly once.
 *
 * Fail-closed: returns 0 (caller deploys passthrough) on any chokepoint.
 *
 * @param[out] sections      Output biquad array.
 * @param[in]  max_sections  Capacity of @p sections.
 * @param[in]  type          filter_type_e.
 * @param[in]  wc1           Pre-warped lower cutoff rad/s (2π·prewarp(fc, fs)).
 * @param[in]  wc2           Pre-warped upper cutoff rad/s for BP/BS (unused for LP/HP).
 * @param[in]  fs            Sampling frequency in Hz.
 * @param[in]  k             Prototype gain.
 * @param[in]  proto_poles   Prototype poles, np elements (copied).
 * @param[in]  np            Number of prototype poles.
 * @param[in]  proto_zeros   Prototype finite zeros, nz elements (may be NULL).
 * @param[in]  nz            Number of prototype finite zeros.
 * @return                   Number of sections deployed, or 0 on failure.
 */
uint8_t design_filter(biquad_filter_t *sections, uint8_t max_sections,
                      uint8_t type,
                      float wc1, float wc2, float fs,
                      float k,
                      const complex_t *proto_poles, uint8_t np,
                      const complex_t *proto_zeros, uint8_t nz);

/**
 * @brief Analytic cascade-level DC and Nyquist gain check.
 *
 * Wrong-but-stable pole/zero pairings and gain-scale errors pass every
 * per-section check (each section is finite and Jury-stable) yet wreck the
 * response shape — the confirmed defect measured DC gain 97.9 and a 350×
 * resonance on a bandstop that must sit at ~1.  H(0) and H(π) are exact
 * rational evaluations (no sampling, no trig, O(sections) cost at init
 * time).
 *
 * Window calibration (measured on deployed f32 coefficients):
 * accepted designs show realized error ≤ ~3.4e-2 (worst case: narrowband
 * fc ≈ 6 Hz @ 48 kHz, where f32 coefficient quantization is genuinely
 * felt); the confirmed pairing defects deviate by ≥ 0.45.  The ±0.1 /
 * ±0.25 windows keep ≥ 3× separation on both sides.
 *
 * @param[in] sections      Deployed biquad cascade.
 * @param[in] num_sections  Number of deployed sections.
 * @param[in] dc_exp        Expected cascade DC gain (0, 1 or ripple-edge).
 * @param[in] ny_exp        Expected cascade Nyquist gain (0, 1 or ripple-edge).
 * @return                  1 if both gains match within tolerance.
 */
uint8_t check_cascade_gains(const biquad_filter_t *sections, uint8_t num_sections,
                            float dc_exp, float ny_exp);

#ifdef __cplusplus
}
#endif

#endif /* FILTER_UTILS_H_ */
