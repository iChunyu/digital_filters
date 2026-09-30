/* Butterworth 族：ROM 原型表 + 共享管线；结构与闸门依据见 docs/butter_filter.md。 */

#ifndef BUTTER_FILTER_H_
#define BUTTER_FILTER_H_

#include <stdint.h>
#include "biquad_filter.h"
#include "filter_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

#define BUTTER_FIELDS \
    uint8_t  valid; \
    uint8_t  type;         /* filter_type_e */ \
    uint8_t  order; \
    uint8_t  num_sections; \
    float    fc1; \
    float    fc2;          /* 上带边（Hz，LP/HP 为 0） */ \
    float    fs;


#define FOR_EACH_BUTTER_LP_ORDER \
    X(1, 1, 1st) \
    X(2, 1, 2nd) \
    X(3, 2, 3rd) \
    X(4, 2, 4th) \
    X(5, 3, 5th) \
    X(6, 3, 6th) \
    X(7, 4, 7th) \
    X(8, 4, 8th)


#define FOR_EACH_BUTTER_BP_ORDER \
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
    typedef struct { BUTTER_FIELDS biquad_filter_t sections[ns]; } butter_lp_##ol##_t;
FOR_EACH_BUTTER_LP_ORDER
#undef X

#define X(order, ns, ol) \
    typedef struct { BUTTER_FIELDS biquad_filter_t sections[ns]; } butter_hp_##ol##_t;
FOR_EACH_BUTTER_LP_ORDER
#undef X

#define X(order, ns, ol) \
    typedef struct { BUTTER_FIELDS biquad_filter_t sections[ns]; } butter_bp_##ol##_t;
FOR_EACH_BUTTER_BP_ORDER
#undef X

#define X(order, ns, ol) \
    typedef struct { BUTTER_FIELDS biquad_filter_t sections[ns]; } butter_bs_##ol##_t;
FOR_EACH_BUTTER_BP_ORDER
#undef X

/* 各阶 init 声明 */

#define X(order, ns, ol) \
    void butter_lp_##ol##_init(butter_lp_##ol##_t *f, float fc, float fs);
FOR_EACH_BUTTER_LP_ORDER
#undef X

#define X(order, ns, ol) \
    void butter_hp_##ol##_init(butter_hp_##ol##_t *f, float fc, float fs);
FOR_EACH_BUTTER_LP_ORDER
#undef X

#define X(order, ns, ol) \
    void butter_bp_##ol##_init(butter_bp_##ol##_t *f, float fc1, float fc2, float fs);
FOR_EACH_BUTTER_BP_ORDER
#undef X

#define X(order, ns, ol) \
    void butter_bs_##ol##_init(butter_bs_##ol##_t *f, float fc1, float fc2, float fs);
FOR_EACH_BUTTER_BP_ORDER
#undef X

/* 各阶 update / reset */

#define X(order, ns, ol) \
    static inline float butter_lp_##ol##_update(butter_lp_##ol##_t *f, float input) \
    { \
        if (!f->valid) return input; \
        return biquad_cascade_update(f->sections, ns, input); \
    } \
    static inline void butter_lp_##ol##_reset(butter_lp_##ol##_t *f, float equilibrium) \
    { \
        if (!f->valid) return; \
        biquad_cascade_reset(f->sections, ns, equilibrium); \
    }
FOR_EACH_BUTTER_LP_ORDER
#undef X

#define X(order, ns, ol) \
    static inline float butter_hp_##ol##_update(butter_hp_##ol##_t *f, float input) \
    { \
        if (!f->valid) return input; \
        return biquad_cascade_update(f->sections, ns, input); \
    } \
    static inline void butter_hp_##ol##_reset(butter_hp_##ol##_t *f, float equilibrium) \
    { \
        if (!f->valid) return; \
        biquad_cascade_reset(f->sections, ns, equilibrium); \
    }
FOR_EACH_BUTTER_LP_ORDER
#undef X

#define X(order, ns, ol) \
    static inline float butter_bp_##ol##_update(butter_bp_##ol##_t *f, float input) \
    { \
        if (!f->valid) return input; \
        return biquad_cascade_update(f->sections, ns, input); \
    } \
    static inline void butter_bp_##ol##_reset(butter_bp_##ol##_t *f, float equilibrium) \
    { \
        if (!f->valid) return; \
        biquad_cascade_reset(f->sections, ns, equilibrium); \
    }
FOR_EACH_BUTTER_BP_ORDER
#undef X

#define X(order, ns, ol) \
    static inline float butter_bs_##ol##_update(butter_bs_##ol##_t *f, float input) \
    { \
        if (!f->valid) return input; \
        return biquad_cascade_update(f->sections, ns, input); \
    } \
    static inline void butter_bs_##ol##_reset(butter_bs_##ol##_t *f, float equilibrium) \
    { \
        if (!f->valid) return; \
        biquad_cascade_reset(f->sections, ns, equilibrium); \
    }
FOR_EACH_BUTTER_BP_ORDER
#undef X

#ifdef __cplusplus
}
#endif

#endif
