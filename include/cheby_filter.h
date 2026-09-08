#ifndef CHEBY_FILTER_H_
#define CHEBY_FILTER_H_

#include <stdint.h>
#include "biquad_filter.h"
#include "filter_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ── 所有静态 Chebyshev 结构体的公共前缀 ─────────────────────────── */

#define CHEBY_FIELDS \
    uint8_t  valid;        /* 1 = init 成功 */                         \
    uint8_t  type;         /* filter_type_e：LOWPASS、HIGHPASS 等 */   \
    uint8_t  order;        /* 滤波器阶数 N */                          \
    uint8_t  num_sections; /* 活跃 biquad 节数 */                      \
    float    fc1;          /* 截止频率 / 下带边（Hz） */               \
    float    fc2;          /* 上带边（Hz，LP/HP 为 0） */              \
    float    fs;           /* 采样频率（Hz） */                        \
    float    ripple_db;    /* 通带纹波（Type I）或阻带衰减（Type II） */

/* ── 阶数表（X-macro）────────────────────────────────────────────── */
/* order, sections_for_lp_hp, ordinal_label */

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

/* ── 各阶结构体 typedef ──────────────────────────────────────────── */

/* Chebyshev I — 低通 */
#define X(order, ns, ol) \
    typedef struct { CHEBY_FIELDS biquad_filter_t sections[ns]; } cheby1_lp_##ol##_t;
FOR_EACH_CHEBY_LP_ORDER
#undef X

/* Chebyshev I — 高通 */
#define X(order, ns, ol) \
    typedef struct { CHEBY_FIELDS biquad_filter_t sections[ns]; } cheby1_hp_##ol##_t;
FOR_EACH_CHEBY_LP_ORDER
#undef X

/* Chebyshev I — 带通 */
#define X(order, ns, ol) \
    typedef struct { CHEBY_FIELDS biquad_filter_t sections[ns]; } cheby1_bp_##ol##_t;
FOR_EACH_CHEBY_BP_ORDER
#undef X

/* Chebyshev I — 带阻 */
#define X(order, ns, ol) \
    typedef struct { CHEBY_FIELDS biquad_filter_t sections[ns]; } cheby1_bs_##ol##_t;
FOR_EACH_CHEBY_BP_ORDER
#undef X

/* Chebyshev II — 低通 */
#define X(order, ns, ol) \
    typedef struct { CHEBY_FIELDS biquad_filter_t sections[ns]; } cheby2_lp_##ol##_t;
FOR_EACH_CHEBY_LP_ORDER
#undef X

/* Chebyshev II — 高通 */
#define X(order, ns, ol) \
    typedef struct { CHEBY_FIELDS biquad_filter_t sections[ns]; } cheby2_hp_##ol##_t;
FOR_EACH_CHEBY_LP_ORDER
#undef X

/* Chebyshev II — 带通 */
#define X(order, ns, ol) \
    typedef struct { CHEBY_FIELDS biquad_filter_t sections[ns]; } cheby2_bp_##ol##_t;
FOR_EACH_CHEBY_BP_ORDER
#undef X

/* Chebyshev II — 带阻 */
#define X(order, ns, ol) \
    typedef struct { CHEBY_FIELDS biquad_filter_t sections[ns]; } cheby2_bs_##ol##_t;
FOR_EACH_CHEBY_BP_ORDER
#undef X

/* ── 各阶 init 声明 ──────────────────────────────────────────────── */

/* Chebyshev I — 低通 */
#define X(order, ns, ol) \
    void cheby1_lp_##ol##_init(cheby1_lp_##ol##_t *f, float fc, float fs, float ripple_db);
FOR_EACH_CHEBY_LP_ORDER
#undef X

/* Chebyshev I — 高通 */
#define X(order, ns, ol) \
    void cheby1_hp_##ol##_init(cheby1_hp_##ol##_t *f, float fc, float fs, float ripple_db);
FOR_EACH_CHEBY_LP_ORDER
#undef X

/* Chebyshev I — 带通 */
#define X(order, ns, ol) \
    void cheby1_bp_##ol##_init(cheby1_bp_##ol##_t *f, float fc1, float fc2, float fs, float ripple_db);
FOR_EACH_CHEBY_BP_ORDER
#undef X

/* Chebyshev I — 带阻 */
#define X(order, ns, ol) \
    void cheby1_bs_##ol##_init(cheby1_bs_##ol##_t *f, float fc1, float fc2, float fs, float ripple_db);
FOR_EACH_CHEBY_BP_ORDER
#undef X

/* Chebyshev II — 低通 */
#define X(order, ns, ol) \
    void cheby2_lp_##ol##_init(cheby2_lp_##ol##_t *f, float fc, float fs, float ripple_db);
FOR_EACH_CHEBY_LP_ORDER
#undef X

/* Chebyshev II — 高通 */
#define X(order, ns, ol) \
    void cheby2_hp_##ol##_init(cheby2_hp_##ol##_t *f, float fc, float fs, float ripple_db);
FOR_EACH_CHEBY_LP_ORDER
#undef X

/* Chebyshev II — 带通 */
#define X(order, ns, ol) \
    void cheby2_bp_##ol##_init(cheby2_bp_##ol##_t *f, float fc1, float fc2, float fs, float ripple_db);
FOR_EACH_CHEBY_BP_ORDER
#undef X

/* Chebyshev II — 带阻 */
#define X(order, ns, ol) \
    void cheby2_bs_##ol##_init(cheby2_bs_##ol##_t *f, float fc1, float fc2, float fs, float ripple_db);
FOR_EACH_CHEBY_BP_ORDER
#undef X

/* ── 各阶 update / reset（static inline — MCU 热路径）─────────────── */

/**
 * @brief 处理一个样本，流经静态分配的 Chebyshev 滤波器。
 *
 * 函数命名遵循 cheby{1,2}_{lp,hp,bp,bs}_{1st..8th}_update。
 * 滤波器无效（!valid）时原样返回 @p input（直通）。
 *
 * 以 static inline 定义，节数以编译期字面量传入：每样本路径
 * 编译为 biquad 循环，零函数调用、零运行时节数装载、零链接符号。
 *
 * @param[in,out] f      滤波器结构体指针。
 * @param[in]     input  当前输入样本。
 * @return               滤波输出。
 */

/**
 * @brief 把静态分配的 Chebyshev 滤波器复位到稳态。
 *
 * 函数命名遵循 cheby{1,2}_{lp,hp,bp,bs}_{1st..8th}_reset。
 * 滤波器无效（!valid）时为空操作。
 *
 * @param[in,out] f           滤波器结构体指针。
 * @param[in]     equilibrium  稳态常值输入。
 */

/* Chebyshev I — 低通 */
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

/* Chebyshev I — 高通 */
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

/* Chebyshev I — 带通 */
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

/* Chebyshev I — 带阻 */
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

/* Chebyshev II — 低通 */
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

/* Chebyshev II — 高通 */
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

/* Chebyshev II — 带通 */
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

/* Chebyshev II — 带阻 */
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

#endif /* CHEBY_FILTER_H_ */
