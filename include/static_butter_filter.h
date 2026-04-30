#ifndef STATIC_BUTTER_FILTER_H_
#define STATIC_BUTTER_FILTER_H_

#include <stdint.h>
#include "biquad_filter.h"
#include "filter_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ── Common prefix for all static Butterworth structs ─────────────────── */

#define STATIC_BUTTER_FIELDS \
    uint8_t  valid;        /* 1 = init succeeded                          */ \
    uint8_t  type;         /* filter_type_e: LOWPASS, HIGHPASS, etc.     */ \
    uint8_t  order;        /* filter order N                              */ \
    uint8_t  num_sections; /* number of active biquad sections            */ \
    float    fc1;          /* cutoff / lower band-edge in Hz              */ \
    float    fc2;          /* upper band-edge in Hz (0 for LP/HP)         */ \
    float    fs;           /* sampling frequency in Hz                    */

/* ── Order tables (X-macro) ───────────────────────────────────────────── */
/* order, sections_for_lp_hp, ordinal_label */

#define FOR_EACH_STATIC_BUTTER_LP_ORDER \
    X(1, 1, 1st) \
    X(2, 1, 2nd) \
    X(3, 2, 3rd) \
    X(4, 2, 4th) \
    X(5, 3, 5th) \
    X(6, 3, 6th) \
    X(7, 4, 7th) \
    X(8, 4, 8th)

/* order, sections_for_bp_bs, ordinal_label */

#define FOR_EACH_STATIC_BUTTER_BP_ORDER \
    X(1, 1, 1st) \
    X(2, 2, 2nd) \
    X(3, 3, 3rd) \
    X(4, 4, 4th) \
    X(5, 5, 5th) \
    X(6, 6, 6th) \
    X(7, 7, 7th) \
    X(8, 8, 8th)

/* ── Per-order struct typedefs ────────────────────────────────────────── */

/* Lowpass */
#define X(order, ns, ol) \
    typedef struct { STATIC_BUTTER_FIELDS biquad_filter_t sections[ns]; } butter_lp_##ol##_t;
FOR_EACH_STATIC_BUTTER_LP_ORDER
#undef X

/* Highpass */
#define X(order, ns, ol) \
    typedef struct { STATIC_BUTTER_FIELDS biquad_filter_t sections[ns]; } butter_hp_##ol##_t;
FOR_EACH_STATIC_BUTTER_LP_ORDER
#undef X

/* Bandpass */
#define X(order, ns, ol) \
    typedef struct { STATIC_BUTTER_FIELDS biquad_filter_t sections[ns]; } butter_bp_##ol##_t;
FOR_EACH_STATIC_BUTTER_BP_ORDER
#undef X

/* Bandstop */
#define X(order, ns, ol) \
    typedef struct { STATIC_BUTTER_FIELDS biquad_filter_t sections[ns]; } butter_bs_##ol##_t;
FOR_EACH_STATIC_BUTTER_BP_ORDER
#undef X

/* ── Per-order init declarations ──────────────────────────────────────── */

/* Lowpass: init(f, fc, fs) */
#define X(order, ns, ol) \
    void butter_lp_##ol##_init(butter_lp_##ol##_t *f, float fc, float fs);
FOR_EACH_STATIC_BUTTER_LP_ORDER
#undef X

/* Highpass: init(f, fc, fs) */
#define X(order, ns, ol) \
    void butter_hp_##ol##_init(butter_hp_##ol##_t *f, float fc, float fs);
FOR_EACH_STATIC_BUTTER_LP_ORDER
#undef X

/* Bandpass: init(f, fc1, fc2, fs) */
#define X(order, ns, ol) \
    void butter_bp_##ol##_init(butter_bp_##ol##_t *f, float fc1, float fc2, float fs);
FOR_EACH_STATIC_BUTTER_BP_ORDER
#undef X

/* Bandstop: init(f, fc1, fc2, fs) */
#define X(order, ns, ol) \
    void butter_bs_##ol##_init(butter_bs_##ol##_t *f, float fc1, float fc2, float fs);
FOR_EACH_STATIC_BUTTER_BP_ORDER
#undef X

/* ── Per-order update / reset declarations ────────────────────────────── */

/**
 * @brief Process one sample through a statically-allocated Butterworth filter.
 *
 * Functions follow the naming convention butter_{lp,hp,bp,bs}_{1st..8th}_update.
 * If the filter is invalid (!valid), returns @p input unchanged (passthrough).
 *
 * @param[in,out] f      Pointer to the filter struct.
 * @param[in]     input  Current input sample.
 * @return               Filtered output.
 */

/**
 * @brief Reset a statically-allocated Butterworth filter to steady-state.
 *
 * Functions follow the naming convention butter_{lp,hp,bp,bs}_{1st..8th}_reset.
 * No-op if the filter is invalid (!valid).
 *
 * @param[in,out] f           Pointer to the filter struct.
 * @param[in]     equilibrium  Constant input value at steady-state.
 */

/* Lowpass */
#define X(order, ns, ol) \
    float butter_lp_##ol##_update(butter_lp_##ol##_t *f, float input); \
    void  butter_lp_##ol##_reset(butter_lp_##ol##_t *f, float equilibrium);
FOR_EACH_STATIC_BUTTER_LP_ORDER
#undef X

/* Highpass */
#define X(order, ns, ol) \
    float butter_hp_##ol##_update(butter_hp_##ol##_t *f, float input); \
    void  butter_hp_##ol##_reset(butter_hp_##ol##_t *f, float equilibrium);
FOR_EACH_STATIC_BUTTER_LP_ORDER
#undef X

/* Bandpass */
#define X(order, ns, ol) \
    float butter_bp_##ol##_update(butter_bp_##ol##_t *f, float input); \
    void  butter_bp_##ol##_reset(butter_bp_##ol##_t *f, float equilibrium);
FOR_EACH_STATIC_BUTTER_BP_ORDER
#undef X

/* Bandstop */
#define X(order, ns, ol) \
    float butter_bs_##ol##_update(butter_bs_##ol##_t *f, float input); \
    void  butter_bs_##ol##_reset(butter_bs_##ol##_t *f, float equilibrium);
FOR_EACH_STATIC_BUTTER_BP_ORDER
#undef X

#ifdef __cplusplus
}
#endif

#endif /* STATIC_BUTTER_FILTER_H_ */
