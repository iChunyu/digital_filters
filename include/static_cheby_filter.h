#ifndef STATIC_CHEBY_FILTER_H_
#define STATIC_CHEBY_FILTER_H_

#include <stdint.h>
#include "biquad_filter.h"
#include "filter_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ── Common prefix for all static Chebyshev structs ───────────────────── */

#define STATIC_CHEBY_FIELDS \
    uint8_t  valid;        /* 1 = init succeeded                          */ \
    uint8_t  type;         /* filter_type_e: LOWPASS, HIGHPASS, etc.     */ \
    uint8_t  order;        /* filter order N                              */ \
    uint8_t  num_sections; /* number of active biquad sections            */ \
    float    fc1;          /* cutoff / lower band-edge in Hz              */ \
    float    fc2;          /* upper band-edge in Hz (0 for LP/HP)         */ \
    float    fs;           /* sampling frequency in Hz                    */ \
    float    ripple_db;    /* passband ripple (Type I) or stopband attn (Type II) */

/* ── Order tables (X-macro) ───────────────────────────────────────────── */
/* order, sections_for_lp_hp, ordinal_label */

#define FOR_EACH_STATIC_CHEBY_LP_ORDER \
    X(1, 1, 1st) \
    X(2, 1, 2nd) \
    X(3, 2, 3rd) \
    X(4, 2, 4th) \
    X(5, 3, 5th) \
    X(6, 3, 6th) \
    X(7, 4, 7th) \
    X(8, 4, 8th)

#define FOR_EACH_STATIC_CHEBY_BP_ORDER \
    X(1, 1, 1st) \
    X(2, 2, 2nd) \
    X(3, 3, 3rd) \
    X(4, 4, 4th) \
    X(5, 5, 5th) \
    X(6, 6, 6th) \
    X(7, 7, 7th) \
    X(8, 8, 8th)

/* ── Per-order struct typedefs ────────────────────────────────────────── */

/* Chebyshev I — lowpass */
#define X(order, ns, ol) \
    typedef struct { STATIC_CHEBY_FIELDS biquad_filter_t sections[ns]; } cheby1_lp_##ol##_t;
FOR_EACH_STATIC_CHEBY_LP_ORDER
#undef X

/* Chebyshev I — highpass */
#define X(order, ns, ol) \
    typedef struct { STATIC_CHEBY_FIELDS biquad_filter_t sections[ns]; } cheby1_hp_##ol##_t;
FOR_EACH_STATIC_CHEBY_LP_ORDER
#undef X

/* Chebyshev I — bandpass */
#define X(order, ns, ol) \
    typedef struct { STATIC_CHEBY_FIELDS biquad_filter_t sections[ns]; } cheby1_bp_##ol##_t;
FOR_EACH_STATIC_CHEBY_BP_ORDER
#undef X

/* Chebyshev I — bandstop */
#define X(order, ns, ol) \
    typedef struct { STATIC_CHEBY_FIELDS biquad_filter_t sections[ns]; } cheby1_bs_##ol##_t;
FOR_EACH_STATIC_CHEBY_BP_ORDER
#undef X

/* Chebyshev II — lowpass */
#define X(order, ns, ol) \
    typedef struct { STATIC_CHEBY_FIELDS biquad_filter_t sections[ns]; } cheby2_lp_##ol##_t;
FOR_EACH_STATIC_CHEBY_LP_ORDER
#undef X

/* Chebyshev II — highpass */
#define X(order, ns, ol) \
    typedef struct { STATIC_CHEBY_FIELDS biquad_filter_t sections[ns]; } cheby2_hp_##ol##_t;
FOR_EACH_STATIC_CHEBY_LP_ORDER
#undef X

/* Chebyshev II — bandpass */
#define X(order, ns, ol) \
    typedef struct { STATIC_CHEBY_FIELDS biquad_filter_t sections[ns]; } cheby2_bp_##ol##_t;
FOR_EACH_STATIC_CHEBY_BP_ORDER
#undef X

/* Chebyshev II — bandstop */
#define X(order, ns, ol) \
    typedef struct { STATIC_CHEBY_FIELDS biquad_filter_t sections[ns]; } cheby2_bs_##ol##_t;
FOR_EACH_STATIC_CHEBY_BP_ORDER
#undef X

/* ── Per-order init declarations ──────────────────────────────────────── */

/* Chebyshev I — lowpass */
#define X(order, ns, ol) \
    void cheby1_lp_##ol##_init(cheby1_lp_##ol##_t *f, float fc, float fs, float ripple_db);
FOR_EACH_STATIC_CHEBY_LP_ORDER
#undef X

/* Chebyshev I — highpass */
#define X(order, ns, ol) \
    void cheby1_hp_##ol##_init(cheby1_hp_##ol##_t *f, float fc, float fs, float ripple_db);
FOR_EACH_STATIC_CHEBY_LP_ORDER
#undef X

/* Chebyshev I — bandpass */
#define X(order, ns, ol) \
    void cheby1_bp_##ol##_init(cheby1_bp_##ol##_t *f, float fc1, float fc2, float fs, float ripple_db);
FOR_EACH_STATIC_CHEBY_BP_ORDER
#undef X

/* Chebyshev I — bandstop */
#define X(order, ns, ol) \
    void cheby1_bs_##ol##_init(cheby1_bs_##ol##_t *f, float fc1, float fc2, float fs, float ripple_db);
FOR_EACH_STATIC_CHEBY_BP_ORDER
#undef X

/* Chebyshev II — lowpass */
#define X(order, ns, ol) \
    void cheby2_lp_##ol##_init(cheby2_lp_##ol##_t *f, float fc, float fs, float ripple_db);
FOR_EACH_STATIC_CHEBY_LP_ORDER
#undef X

/* Chebyshev II — highpass */
#define X(order, ns, ol) \
    void cheby2_hp_##ol##_init(cheby2_hp_##ol##_t *f, float fc, float fs, float ripple_db);
FOR_EACH_STATIC_CHEBY_LP_ORDER
#undef X

/* Chebyshev II — bandpass */
#define X(order, ns, ol) \
    void cheby2_bp_##ol##_init(cheby2_bp_##ol##_t *f, float fc1, float fc2, float fs, float ripple_db);
FOR_EACH_STATIC_CHEBY_BP_ORDER
#undef X

/* Chebyshev II — bandstop */
#define X(order, ns, ol) \
    void cheby2_bs_##ol##_init(cheby2_bs_##ol##_t *f, float fc1, float fc2, float fs, float ripple_db);
FOR_EACH_STATIC_CHEBY_BP_ORDER
#undef X

/* ── Per-order update / reset declarations ────────────────────────────── */

/**
 * @brief Process one sample through a statically-allocated Chebyshev filter.
 *
 * Functions follow the naming convention
 * cheby{1,2}_{lp,hp,bp,bs}_{1st..8th}_update.
 * If the filter is invalid (!valid), returns @p input unchanged (passthrough).
 *
 * @param[in,out] f      Pointer to the filter struct.
 * @param[in]     input  Current input sample.
 * @return               Filtered output.
 */

/**
 * @brief Reset a statically-allocated Chebyshev filter to steady-state.
 *
 * Functions follow the naming convention
 * cheby{1,2}_{lp,hp,bp,bs}_{1st..8th}_reset.
 * No-op if the filter is invalid (!valid).
 *
 * @param[in,out] f           Pointer to the filter struct.
 * @param[in]     equilibrium  Constant input value at steady-state.
 */

/* Chebyshev I — lowpass */
#define X(order, ns, ol) \
    float cheby1_lp_##ol##_update(cheby1_lp_##ol##_t *f, float input); \
    void  cheby1_lp_##ol##_reset(cheby1_lp_##ol##_t *f, float equilibrium);
FOR_EACH_STATIC_CHEBY_LP_ORDER
#undef X

/* Chebyshev I — highpass */
#define X(order, ns, ol) \
    float cheby1_hp_##ol##_update(cheby1_hp_##ol##_t *f, float input); \
    void  cheby1_hp_##ol##_reset(cheby1_hp_##ol##_t *f, float equilibrium);
FOR_EACH_STATIC_CHEBY_LP_ORDER
#undef X

/* Chebyshev I — bandpass */
#define X(order, ns, ol) \
    float cheby1_bp_##ol##_update(cheby1_bp_##ol##_t *f, float input); \
    void  cheby1_bp_##ol##_reset(cheby1_bp_##ol##_t *f, float equilibrium);
FOR_EACH_STATIC_CHEBY_BP_ORDER
#undef X

/* Chebyshev I — bandstop */
#define X(order, ns, ol) \
    float cheby1_bs_##ol##_update(cheby1_bs_##ol##_t *f, float input); \
    void  cheby1_bs_##ol##_reset(cheby1_bs_##ol##_t *f, float equilibrium);
FOR_EACH_STATIC_CHEBY_BP_ORDER
#undef X

/* Chebyshev II — lowpass */
#define X(order, ns, ol) \
    float cheby2_lp_##ol##_update(cheby2_lp_##ol##_t *f, float input); \
    void  cheby2_lp_##ol##_reset(cheby2_lp_##ol##_t *f, float equilibrium);
FOR_EACH_STATIC_CHEBY_LP_ORDER
#undef X

/* Chebyshev II — highpass */
#define X(order, ns, ol) \
    float cheby2_hp_##ol##_update(cheby2_hp_##ol##_t *f, float input); \
    void  cheby2_hp_##ol##_reset(cheby2_hp_##ol##_t *f, float equilibrium);
FOR_EACH_STATIC_CHEBY_LP_ORDER
#undef X

/* Chebyshev II — bandpass */
#define X(order, ns, ol) \
    float cheby2_bp_##ol##_update(cheby2_bp_##ol##_t *f, float input); \
    void  cheby2_bp_##ol##_reset(cheby2_bp_##ol##_t *f, float equilibrium);
FOR_EACH_STATIC_CHEBY_BP_ORDER
#undef X

/* Chebyshev II — bandstop */
#define X(order, ns, ol) \
    float cheby2_bs_##ol##_update(cheby2_bs_##ol##_t *f, float input); \
    void  cheby2_bs_##ol##_reset(cheby2_bs_##ol##_t *f, float equilibrium);
FOR_EACH_STATIC_CHEBY_BP_ORDER
#undef X

#ifdef __cplusplus
}
#endif

#endif /* STATIC_CHEBY_FILTER_H_ */
