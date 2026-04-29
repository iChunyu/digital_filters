#ifndef SOS_FILTER_H_
#define SOS_FILTER_H_

#include <stdint.h>
#include "biquad_filter.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Cascade of second-order sections (biquads) forming a higher-order
 *        IIR filter.
 *
 * Each section is a biquad_filter_t.  Memory for the sections is allocated by
 * sos_filter_init() and freed by sos_filter_destroy().
 */
typedef struct {
    uint8_t num_sections;      /**< Number of biquad sections in the cascade */
    biquad_filter_t *sections; /**< Array of biquad sections (length num_sections) */
} sos_filter_t;

/**
 * @brief Initialise a cascaded SOS filter.
 *
 * Allocates memory for @p num_sections biquad sections and initialises each
 * to unity pass-through (identity).
 *
 * @param[out] filter       Pointer to the SOS filter object.
 * @param[in]  num_sections Number of second-order sections.
 */
void sos_filter_init(sos_filter_t *filter, uint8_t num_sections);

/**
 * @brief Free the memory allocated by sos_filter_init().
 *
 * @param[in,out] filter  Pointer to the SOS filter object.
 */
void sos_filter_destroy(sos_filter_t *filter);

/**
 * @brief Set the coefficients of one biquad section.
 *
 * If @p section_index is out of bounds the call is silently ignored.
 *
 * @param[in,out] filter         Pointer to the SOS filter object.
 * @param[in]     section_index  Zero-based index of the section to configure.
 * @param[in]     num_z          Numerator coefficients   [b0, b1, b2] in z-domain.
 * @param[in]     den_z          Denominator coefficients [a0, a1, a2] in z-domain.
 */
void sos_filter_set_section(sos_filter_t *filter, uint8_t section_index,
                             const float num_z[3], const float den_z[3]);

/**
 * @brief Process one input sample through the cascade and return the output.
 *
 * Samples flow: input -> section[0] -> section[1] -> ... -> section[N-1] -> output.
 *
 * @param[in,out] filter  Pointer to the SOS filter object.
 * @param[in]     input   Current input sample.
 *
 * @return Filtered output sample.
 */
float sos_filter_update(sos_filter_t *filter, float input);

/**
 * @brief Return the current output of the cascade (last section) without
 *        advancing state.
 *
 * @param[in] filter  Pointer to the SOS filter object.
 *
 * @return Current output.
 */
float sos_filter_get_output(const sos_filter_t *filter);

/**
 * @brief Reconstruct the input to the first section from its internal state.
 *
 * This recovers the overall filter input @f$ x[n] @f$, useful for debugging
 * or cascading multiple SOS filters.
 *
 * @param[in] filter  Pointer to the SOS filter object.
 *
 * @return Reconstructed input.
 */
float sos_filter_get_input(const sos_filter_t *filter);

/**
 * @brief Reset all sections to steady-state equilibrium.
 *
 * The given @p equilibrium is the steady-state input to the cascade.  Each
 * section's intermediate steady-state value is computed by propagating the
 * DC gain through the preceding sections, so the overall DC response is
 * preserved.
 *
 * @param[in,out] filter      Pointer to the SOS filter object.
 * @param[in]     equilibrium  Constant input value at steady-state.
 */
void sos_filter_reset(sos_filter_t *filter, float equilibrium);

#ifdef __cplusplus
}
#endif

#endif  /* SOS_FILTER_H_ */
