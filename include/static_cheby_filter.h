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

/**
 * @brief Base type for unified static_cheby_update() / static_cheby_reset().
 *
 * Cast any cheby1_lp_Nth_t*, cheby2_hp_Nth_t*, etc. pointer to
 * static_cheby_t* — the common prefix guarantees identical memory layout.
 */
typedef struct {
    STATIC_CHEBY_FIELDS
    biquad_filter_t sections[1]; /* placeholder; real count = num_sections */
} static_cheby_t;

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

/* ── Unified update / reset ───────────────────────────────────────────── */

float static_cheby_update(static_cheby_t *f, float input);
void  static_cheby_reset(static_cheby_t *f, float equilibrium);

#ifdef __cplusplus
}
#endif

#endif /* STATIC_CHEBY_FILTER_H_ */
