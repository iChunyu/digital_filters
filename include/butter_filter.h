#ifndef BUTTER_FILTER_H_
#define BUTTER_FILTER_H_

#include <stdint.h>
#include "sos_filter.h"
#include "filter_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Butterworth filter object.
 *
 * Wraps a cascade of second-order sections (sos_filter_t) designed from the
 * Butterworth analog prototype via pre-warping, frequency transform, bilinear
 * discretisation, and pole-zero pairing.
 */
typedef struct {
    uint8_t       valid;  /**< 1 if initialisation succeeded, 0 otherwise */
    uint8_t       type;   /**< filter_type_e: LOWPASS, HIGHPASS, etc.       */
    uint8_t       order;  /**< Filter order N                               */
    float         fc1;    /**< Cutoff / lower band-edge in Hz               */
    float         fc2;    /**< Upper band-edge in Hz (unused for LP / HP)   */
    float         fs;     /**< Sampling frequency in Hz                     */
    sos_filter_t *sos;    /**< Cascade of biquad sections (allocated on init) */
} butter_t;

/**
 * @brief Design and initialise a Butterworth filter.
 *
 * Generates the analog prototype, applies the frequency transform and bilinear
 * discretisation, pairs poles and zeros into second-order sections, and
 * deploys them into the internal sos_filter_t cascade.
 *
 * On failure (invalid order, cutoff out of range, etc.) @p valid is set to 0
 * and @p sos is left as NULL.
 *
 * @param[out] b      Pointer to the butter_t object.
 * @param[in]  type   Filter type (FILTER_LOWPASS or FILTER_HIGHPASS).
 * @param[in]  order  Filter order N (≥ 1).
 * @param[in]  fc1    Cutoff frequency in Hz (0 < fc1 < fs/2).
 * @param[in]  fc2    Upper band-edge in Hz (unused for LP/HP, pass 0).
 * @param[in]  fs     Sampling frequency in Hz.
 */
void butter_init(butter_t *b, uint8_t type, uint8_t order,
                 float fc1, float fc2, float fs);

/**
 * @brief Free the internal SOS cascade and mark the filter invalid.
 *
 * Safe to call on an already-destroyed or never-initialised filter
 * (sos == NULL is a no-op).
 *
 * @param[in,out] b  Pointer to the butter_t object.
 */
void butter_destroy(butter_t *b);

/**
 * @brief Process one input sample and return the filtered output.
 *
 * If the filter was not successfully initialised (@p valid == 0) the input is
 * passed through unchanged.
 *
 * @param[in,out] b      Pointer to the butter_t object.
 * @param[in]     input  Current input sample.
 * @return               Filtered output sample.
 */
float butter_update(butter_t *b, float input);

/**
 * @brief Reset the filter to a steady-state equilibrium.
 *
 * Propagates the given constant input through the cascade so that subsequent
 * calls to butter_update() produce the same constant output immediately.
 * No-op if @p valid == 0.
 *
 * @param[in,out] b           Pointer to the butter_t object.
 * @param[in]     equilibrium  Constant input value at steady-state.
 */
void butter_reset(butter_t *b, float equilibrium);

#ifdef __cplusplus
}
#endif

#endif /* BUTTER_FILTER_H_ */
