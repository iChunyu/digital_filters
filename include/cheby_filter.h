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
/* order 阶数, sections_for_lp_hp LP/HP 节数, ordinal_label 序数标签 */

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
/*
 * 函数命名遵循 cheby{1,2}_{lp,hp,bp,bs}_{1st..8th}_init，
 * 签名末尾多一个 @c ripple_db（Chebyshev I 是通带纹波，Chebyshev II 是阻带纹波）。
 *
 * 参数范围：@c ripple_db > 0，@c 0 < fc < fs/2（BP/BS 另需 @c fc1 < fc2 < fs/2）。
 *
 * Fail-closed 契约：任一参数越界 / 设计管线失败 → @c valid == 0 且
 * @c num_sections == 0（两者都被显式清零，即使原值非法）；此后
 * @c _update 直通返回输入、@c _reset 为空操作。**调用方必须检查 @c valid**。
 *
 * 失败后 @c sections[] 的内容不保证，不得直接使用。
 * init 成功后内部状态为零（不是稳态）：需要带直流的信号从零相位开始时，
 * 先调 @c _reset(equilibrium)。
 *
 * 近 DC / 近 Nyquist / 极窄带 / 超宽带设计被拒是预期行为（f32 系数量化
 * 撑不住边缘增益），不是错误；边界非单调，以 @c valid 为准。
 * 偶数阶 Chebyshev I 的 DC 增益是 @c 10^(-ripple_db/20)（非 1），
 * Chebyshev II 的阻带边缘增益同理——级联增益窗口对这类边缘预期值用 ±0.25。
 * 标定与回归档位见 test/test_cheby.c。
 */

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
/*
 * 下面两个 doxygen 块描述的是**由 X-macro 生成的一整族函数**，具体函数体在
 * 本文件末尾的 `#define X` 块里。函数命名：
 *   cheby{1,2}_{lp,hp,bp,bs}_{1st..8th}_update(f, input)  → 返回滤波输出
 *   cheby{1,2}_{lp,hp,bp,bs}_{1st..8th}_reset(f, equilibrium)
 *
 * 两者都以 static inline 定义，节数以编译期字面量传入：每样本路径编译为
 * biquad 循环，零函数调用、零运行时节数装载、零链接符号。
 *
 * 无效滤波器（@c !valid）：update 原样返回输入（直通），reset 为空操作。
 *
 * 不可重入：直接读写 @p f 内的状态。同一滤波器对象被 ISR 与主循环共享时必须
 * 自行关中断或双缓冲。@c _update 刻意不做 NaN/Inf 防护（每样本分支开销），
 * 非有限输入会毒化状态直至 @c _reset。
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
