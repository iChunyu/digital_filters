#ifndef BUTTER_FILTER_H_
#define BUTTER_FILTER_H_

#include <stdint.h>
#include "biquad_filter.h"
#include "filter_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ── 所有静态 Butterworth 结构体的公共前缀 ───────────────────────── */

#define BUTTER_FIELDS \
    uint8_t  valid;        /* 1 = init 成功 */                         \
    uint8_t  type;         /* filter_type_e：LOWPASS、HIGHPASS 等 */   \
    uint8_t  order;        /* 滤波器阶数 N */                          \
    uint8_t  num_sections; /* 活跃 biquad 节数 */                      \
    float    fc1;          /* 截止频率 / 下带边（Hz） */               \
    float    fc2;          /* 上带边（Hz，LP/HP 为 0） */              \
    float    fs;           /* 采样频率（Hz） */

/* ── 阶数表（X-macro）────────────────────────────────────────────── */
/* order 阶数, sections_for_lp_hp LP/HP 节数, ordinal_label 序数标签 */

#define FOR_EACH_BUTTER_LP_ORDER \
    X(1, 1, 1st) \
    X(2, 1, 2nd) \
    X(3, 2, 3rd) \
    X(4, 2, 4th) \
    X(5, 3, 5th) \
    X(6, 3, 6th) \
    X(7, 4, 7th) \
    X(8, 4, 8th)

/* order 阶数, sections_for_bp_bs BP/BS 节数, ordinal_label 序数标签 */

#define FOR_EACH_BUTTER_BP_ORDER \
    X(1, 1, 1st) \
    X(2, 2, 2nd) \
    X(3, 3, 3rd) \
    X(4, 4, 4th) \
    X(5, 5, 5th) \
    X(6, 6, 6th) \
    X(7, 7, 7th) \
    X(8, 8, 8th)

/* ── 各阶结构体 typedef ──────────────────────────────────────────── */

/* 低通 */
#define X(order, ns, ol) \
    typedef struct { BUTTER_FIELDS biquad_filter_t sections[ns]; } butter_lp_##ol##_t;
FOR_EACH_BUTTER_LP_ORDER
#undef X

/* 高通 */
#define X(order, ns, ol) \
    typedef struct { BUTTER_FIELDS biquad_filter_t sections[ns]; } butter_hp_##ol##_t;
FOR_EACH_BUTTER_LP_ORDER
#undef X

/* 带通 */
#define X(order, ns, ol) \
    typedef struct { BUTTER_FIELDS biquad_filter_t sections[ns]; } butter_bp_##ol##_t;
FOR_EACH_BUTTER_BP_ORDER
#undef X

/* 带阻 */
#define X(order, ns, ol) \
    typedef struct { BUTTER_FIELDS biquad_filter_t sections[ns]; } butter_bs_##ol##_t;
FOR_EACH_BUTTER_BP_ORDER
#undef X

/* ── 各阶 init 声明 ──────────────────────────────────────────────── */

/* 低通：init(f, fc, fs) */
#define X(order, ns, ol) \
    void butter_lp_##ol##_init(butter_lp_##ol##_t *f, float fc, float fs);
FOR_EACH_BUTTER_LP_ORDER
#undef X

/* 高通：init(f, fc, fs) */
#define X(order, ns, ol) \
    void butter_hp_##ol##_init(butter_hp_##ol##_t *f, float fc, float fs);
FOR_EACH_BUTTER_LP_ORDER
#undef X

/* 带通：init(f, fc1, fc2, fs) */
#define X(order, ns, ol) \
    void butter_bp_##ol##_init(butter_bp_##ol##_t *f, float fc1, float fc2, float fs);
FOR_EACH_BUTTER_BP_ORDER
#undef X

/* 带阻：init(f, fc1, fc2, fs) */
#define X(order, ns, ol) \
    void butter_bs_##ol##_init(butter_bs_##ol##_t *f, float fc1, float fc2, float fs);
FOR_EACH_BUTTER_BP_ORDER
#undef X

/* ── 各阶 update / reset（static inline — MCU 热路径）─────────────── */

/**
 * @brief 处理一个样本，流经静态分配的 Butterworth 滤波器。
 *
 * 函数命名遵循 butter_{lp,hp,bp,bs}_{1st..8th}_update。
 * 滤波器无效（!valid）时原样返回 @p input（直通）。
 *
 * 以 static inline 定义，节数以编译期字面量传入：每样本路径
 * 编译为 biquad 循环，零函数调用、零运行时节数装载、零链接符号
 * （96 个几乎相同的 out-of-line 副本曾占 ~12.8 KB flash，
 * 每样本还多一次调用开销）。
 *
 * @param[in,out] f      滤波器结构体指针。
 * @param[in]     input  当前输入样本。
 * @return               滤波输出。
 */

/**
 * @brief 把静态分配的 Butterworth 滤波器复位到稳态。
 *
 * 函数命名遵循 butter_{lp,hp,bp,bs}_{1st..8th}_reset。
 * 滤波器无效（!valid）时为空操作。
 *
 * @param[in,out] f           滤波器结构体指针。
 * @param[in]     equilibrium  稳态常值输入。
 */

/* 低通 */
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

/* 高通 */
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

/* 带通 */
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

/* 带阻 */
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

#endif /* BUTTER_FILTER_H_ */
