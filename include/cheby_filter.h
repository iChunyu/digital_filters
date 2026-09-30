/* Chebyshev I / II 族：原型运行时计算（依赖纹波）；依据见 docs/cheby_filter.md。 */

#ifndef CHEBY_FILTER_H_
#define CHEBY_FILTER_H_

#include <stdint.h>
#include "biquad_filter.h"
#include "filter_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

#define CHEBY_FIELDS \
    uint8_t  valid; \
    uint8_t  type; \
    uint8_t  order; \
    uint8_t  num_sections; \
    float    fc1; \
    float    fc2;          /* 上带边（Hz，LP/HP 为 0） */ \
    float    fs; \
    float    ripple_db;    /* 通带纹波（Type I）或阻带衰减（Type II） */


#define FOR_EACH_CHEBY_LP_ORDER \
    X(1, 1, 1st) \
    X(2, 1, 2nd) \
    X(3, 2, 3rd) \
    X(4, 2, 4th) \
    X(5, 3, 5th) \
    X(6, 3, 6th) \
    X(7, 4, 7th) \
    X(8, 4, 8th)

#define FOR_EACH_CHEBY_BP_ORDER \
    X(1, 1, 1st) \
    X(2, 2, 2nd) \
    X(3, 3, 3rd) \
    X(4, 4, 4th) \
    X(5, 5, 5th) \
    X(6, 6, 6th) \
    X(7, 7, 7th) \
    X(8, 8, 8th)

/* 各阶结构体 typedef */

#define X(order, ns, ol) \
    typedef struct { CHEBY_FIELDS biquad_filter_t sections[ns]; } cheby1_lp_##ol##_t;
FOR_EACH_CHEBY_LP_ORDER
#undef X

#define X(order, ns, ol) \
    typedef struct { CHEBY_FIELDS biquad_filter_t sections[ns]; } cheby1_hp_##ol##_t;
FOR_EACH_CHEBY_LP_ORDER
#undef X

#define X(order, ns, ol) \
    typedef struct { CHEBY_FIELDS biquad_filter_t sections[ns]; } cheby1_bp_##ol##_t;
FOR_EACH_CHEBY_BP_ORDER
#undef X

#define X(order, ns, ol) \
    typedef struct { CHEBY_FIELDS biquad_filter_t sections[ns]; } cheby1_bs_##ol##_t;
FOR_EACH_CHEBY_BP_ORDER
#undef X

#define X(order, ns, ol) \
    typedef struct { CHEBY_FIELDS biquad_filter_t sections[ns]; } cheby2_lp_##ol##_t;
FOR_EACH_CHEBY_LP_ORDER
#undef X

#define X(order, ns, ol) \
    typedef struct { CHEBY_FIELDS biquad_filter_t sections[ns]; } cheby2_hp_##ol##_t;
FOR_EACH_CHEBY_LP_ORDER
#undef X

#define X(order, ns, ol) \
    typedef struct { CHEBY_FIELDS biquad_filter_t sections[ns]; } cheby2_bp_##ol##_t;
FOR_EACH_CHEBY_BP_ORDER
#undef X

#define X(order, ns, ol) \
    typedef struct { CHEBY_FIELDS biquad_filter_t sections[ns]; } cheby2_bs_##ol##_t;
FOR_EACH_CHEBY_BP_ORDER
#undef X

/* 各阶 init 声明（末尾多一个 ripple_db） */

#define X(order, ns, ol) \
    void cheby1_lp_##ol##_init(cheby1_lp_##ol##_t *f, float fc, float fs, float ripple_db);
FOR_EACH_CHEBY_LP_ORDER
#undef X

#define X(order, ns, ol) \
    void cheby1_hp_##ol##_init(cheby1_hp_##ol##_t *f, float fc, float fs, float ripple_db);
FOR_EACH_CHEBY_LP_ORDER
#undef X

#define X(order, ns, ol) \
    void cheby1_bp_##ol##_init(cheby1_bp_##ol##_t *f, float fc1, float fc2, float fs, float ripple_db);
FOR_EACH_CHEBY_BP_ORDER
#undef X

#define X(order, ns, ol) \
    void cheby1_bs_##ol##_init(cheby1_bs_##ol##_t *f, float fc1, float fc2, float fs, float ripple_db);
FOR_EACH_CHEBY_BP_ORDER
#undef X

#define X(order, ns, ol) \
    void cheby2_lp_##ol##_init(cheby2_lp_##ol##_t *f, float fc, float fs, float ripple_db);
FOR_EACH_CHEBY_LP_ORDER
#undef X

#define X(order, ns, ol) \
    void cheby2_hp_##ol##_init(cheby2_hp_##ol##_t *f, float fc, float fs, float ripple_db);
FOR_EACH_CHEBY_LP_ORDER
#undef X

#define X(order, ns, ol) \
    void cheby2_bp_##ol##_init(cheby2_bp_##ol##_t *f, float fc1, float fc2, float fs, float ripple_db);
FOR_EACH_CHEBY_BP_ORDER
#undef X

#define X(order, ns, ol) \
    void cheby2_bs_##ol##_init(cheby2_bs_##ol##_t *f, float fc1, float fc2, float fs, float ripple_db);
FOR_EACH_CHEBY_BP_ORDER
#undef X

/* 各阶 update / reset（static inline） */

#define X(order, ns, ol) \
    static inline float cheby1_lp_##ol##_update(cheby1_lp_##ol##_t *f, float input) \
    { \
        if (!f->valid) return input; \
        return biquad_cascade_update(f->sections, ns, input); \
    } \
    static inline void cheby1_lp_##ol##_reset(cheby1_lp_##ol##_t *f, float equilibrium) \
    { \
        if (!f->valid) return; \
        biquad_cascade_reset(f->sections, ns, equilibrium); \
    }
FOR_EACH_CHEBY_LP_ORDER
#undef X

#define X(order, ns, ol) \
    static inline float cheby1_hp_##ol##_update(cheby1_hp_##ol##_t *f, float input) \
    { \
        if (!f->valid) return input; \
        return biquad_cascade_update(f->sections, ns, input); \
    } \
    static inline void cheby1_hp_##ol##_reset(cheby1_hp_##ol##_t *f, float equilibrium) \
    { \
        if (!f->valid) return; \
        biquad_cascade_reset(f->sections, ns, equilibrium); \
    }
FOR_EACH_CHEBY_LP_ORDER
#undef X

#define X(order, ns, ol) \
    static inline float cheby1_bp_##ol##_update(cheby1_bp_##ol##_t *f, float input) \
    { \
        if (!f->valid) return input; \
        return biquad_cascade_update(f->sections, ns, input); \
    } \
    static inline void cheby1_bp_##ol##_reset(cheby1_bp_##ol##_t *f, float equilibrium) \
    { \
        if (!f->valid) return; \
        biquad_cascade_reset(f->sections, ns, equilibrium); \
    }
FOR_EACH_CHEBY_BP_ORDER
#undef X

#define X(order, ns, ol) \
    static inline float cheby1_bs_##ol##_update(cheby1_bs_##ol##_t *f, float input) \
    { \
        if (!f->valid) return input; \
        return biquad_cascade_update(f->sections, ns, input); \
    } \
    static inline void cheby1_bs_##ol##_reset(cheby1_bs_##ol##_t *f, float equilibrium) \
    { \
        if (!f->valid) return; \
        biquad_cascade_reset(f->sections, ns, equilibrium); \
    }
FOR_EACH_CHEBY_BP_ORDER
#undef X

#define X(order, ns, ol) \
    static inline float cheby2_lp_##ol##_update(cheby2_lp_##ol##_t *f, float input) \
    { \
        if (!f->valid) return input; \
        return biquad_cascade_update(f->sections, ns, input); \
    } \
    static inline void cheby2_lp_##ol##_reset(cheby2_lp_##ol##_t *f, float equilibrium) \
    { \
        if (!f->valid) return; \
        biquad_cascade_reset(f->sections, ns, equilibrium); \
    }
FOR_EACH_CHEBY_LP_ORDER
#undef X

#define X(order, ns, ol) \
    static inline float cheby2_hp_##ol##_update(cheby2_hp_##ol##_t *f, float input) \
    { \
        if (!f->valid) return input; \
        return biquad_cascade_update(f->sections, ns, input); \
    } \
    static inline void cheby2_hp_##ol##_reset(cheby2_hp_##ol##_t *f, float equilibrium) \
    { \
        if (!f->valid) return; \
        biquad_cascade_reset(f->sections, ns, equilibrium); \
    }
FOR_EACH_CHEBY_LP_ORDER
#undef X

#define X(order, ns, ol) \
    static inline float cheby2_bp_##ol##_update(cheby2_bp_##ol##_t *f, float input) \
    { \
        if (!f->valid) return input; \
        return biquad_cascade_update(f->sections, ns, input); \
    } \
    static inline void cheby2_bp_##ol##_reset(cheby2_bp_##ol##_t *f, float equilibrium) \
    { \
        if (!f->valid) return; \
        biquad_cascade_reset(f->sections, ns, equilibrium); \
    }
FOR_EACH_CHEBY_BP_ORDER
#undef X

#define X(order, ns, ol) \
    static inline float cheby2_bs_##ol##_update(cheby2_bs_##ol##_t *f, float input) \
    { \
        if (!f->valid) return input; \
        return biquad_cascade_update(f->sections, ns, input); \
    } \
    static inline void cheby2_bs_##ol##_reset(cheby2_bs_##ol##_t *f, float equilibrium) \
    { \
        if (!f->valid) return; \
        biquad_cascade_reset(f->sections, ns, equilibrium); \
    }
FOR_EACH_CHEBY_BP_ORDER
#undef X

#ifdef __cplusplus
}
#endif

#endif
