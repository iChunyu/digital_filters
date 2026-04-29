#ifndef CHEBY_FILTER_H_
#define CHEBY_FILTER_H_

#include <stdint.h>
#include "sos_filter.h"
#include "filter_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================== */

/**
 * @brief Chebyshev Type I filter object (passband ripple).
 *
 * Wraps a cascade of second-order sections designed from the Chebyshev Type I
 * analog prototype.  The passband has equiripple behaviour between 0 and
 * @p ripple_db dB; the stopband rolls off monotonically.
 */
typedef struct {
    uint8_t       valid;      /**< 1 if initialisation succeeded           */
    uint8_t       type;       /**< filter_type_e: LOWPASS, HIGHPASS, etc.  */
    uint8_t       order;      /**< Filter order N                          */
    float         fc1;        /**< Cutoff / lower band-edge in Hz          */
    float         fc2;        /**< Upper band-edge in Hz (unused for LP/HP)*/
    float         fs;         /**< Sampling frequency in Hz                */
    float         ripple_db;  /**< Passband ripple (dB), > 0               */
    sos_filter_t *sos;        /**< Cascade of biquad sections              */
} cheby1_t;

void cheby1_init(cheby1_t *c, uint8_t type, uint8_t order,
                 float fc1, float fc2, float fs, float ripple_db);
void cheby1_destroy(cheby1_t *c);
float cheby1_update(cheby1_t *c, float input);
float cheby1_reset(cheby1_t *c, float equilibrium);

/* ================================================================== */

/**
 * @brief Chebyshev Type II filter object (stopband ripple).
 *
 * Wraps a cascade of second-order sections designed from the Chebyshev Type II
 * analog prototype.  The passband is monotonic; the stopband has equiripple
 * behaviour with minimum attenuation @p ripple_db dB.
 */
typedef struct {
    uint8_t       valid;      /**< 1 if initialisation succeeded           */
    uint8_t       type;       /**< filter_type_e: LOWPASS, HIGHPASS, etc.  */
    uint8_t       order;      /**< Filter order N                          */
    float         fc1;        /**< Cutoff / lower band-edge in Hz          */
    float         fc2;        /**< Upper band-edge in Hz (unused for LP/HP)*/
    float         fs;         /**< Sampling frequency in Hz                */
    float         ripple_db;  /**< Minimum stopband attenuation (dB), > 0  */
    sos_filter_t *sos;        /**< Cascade of biquad sections              */
} cheby2_t;

void cheby2_init(cheby2_t *c, uint8_t type, uint8_t order,
                 float fc1, float fc2, float fs, float ripple_db);
void cheby2_destroy(cheby2_t *c);
float cheby2_update(cheby2_t *c, float input);
float cheby2_reset(cheby2_t *c, float equilibrium);

#ifdef __cplusplus
}
#endif

#endif /* CHEBY_FILTER_H_ */
