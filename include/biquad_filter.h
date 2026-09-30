/**
 * @file    biquad_filter.h
 * @brief   二阶 IIR 滤波器（直接 II 型／规范型）。
 *
 * @f[
 *   H(z) = \frac{b_0 + b_1 z^{-1} + b_2 z^{-2}}
 *                {1   + a_1 z^{-1} + a_2 z^{-2}}
 * @f]
 *
 * 直接 II 型（规范型）只用 3 个状态变量 @p w[0..2]，是二阶节最省内存的形式。
 */

#ifndef BIQUAD_FILTER_H_
#define BIQUAD_FILTER_H_

#include <stdint.h>

/*
 * fail-closed 的前提：isfinite 不能被编译器折叠成恒真。`-ffast-math` 隐含
 * `-ffinite-math-only`，编译器据此假设不存在 NaN/Inf，于是全部有限性闸静默
 * 失效（退化设计被当作 valid=1 部署出去）；重结合还会破坏 sum3f 与 Dekker
 * 分裂的误差无损性。实测细节与修复口径见 AGENTS.md 的「禁止 -ffast-math」段。
 *
 * 本头是六个公开头的唯一公共顶点（其余五族全部包含它），守卫放这一处即覆盖
 * 全部 translation unit。
 */
#if defined(__FAST_MATH__)
#error "本库要求 IEEE-754 严格语义（fail-closed 依赖 isfinite）：请移除 -ffast-math"
#endif
#if defined(__FINITE_MATH_ONLY__) && (__FINITE_MATH_ONLY__ > 0)
#error "本库要求 IEEE-754 严格语义（fail-closed 依赖 isfinite）：请移除 -ffinite-math-only"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Biquad 滤波器对象（直接 II 型）。
 *
 * 系数按归一化存储，@p den_z[0] 恒为 1.0。状态向量 @p w[] 存放
 * 规范型实现的中间值：
 *
 * @verbatim
 *   x[n] --->(+) --------> w[n] -----[b0]--->(+)---> y[n]
 *             ^-           |                  ^
 *             |           z^-1                |
 *             |            |                  |
 *             +----[a1]----w[n-1]----[b1]-----+
 *             |            |                  |
 *             |           z^-1                |
 *             |            |                  |
 *             +----[a2]----w[n-2]----[b2]-----+
 * @endverbatim
 */
typedef struct {
    float num_z[3]; /**< 分子系数 b0, b1, b2 */
    float den_z[3]; /**< 分母系数 1.0, a1, a2（已归一化） */
    float w[3];     /**< 状态变量 w[n], w[n-1], w[n-2] */
} biquad_filter_t;

/**
 * @brief 用双线性变换把 s 域 biquad 系数映射到 z 域。
 *
 * 代换 @f$ s = 2 f_s \frac{1 - z^{-1}}{1 + z^{-1}} @f$，作用于
 * @f$ H(s) = \frac{B_0 + B_1 s + B_2 s^2}{A_0 + A_1 s + A_2 s^2} @f$；
 * 输出的 @p den_z 未归一化（@p den_z[0] 未必为 1）。
 *
 * @param[out] num_z  得到的 z 域分子   [b0, b1, b2]。
 * @param[out] den_z  得到的 z 域分母   [a0, a1, a2]。
 * @param[in]  num_s  s 域分子   [B0, B1, B2]。
 * @param[in]  den_s  s 域分母   [A0, A1, A2]。
 * @param[in]  fs     采样频率（Hz）。
 */
void biquad_c2d_bilinear(float num_z[3], float den_z[3], const float num_s[3],
                        const float den_s[3], float fs);

/**
 * @brief 把 biquad 滤波器置为单位直通（恒等）。
 *
 * 初始化后 @f$ H(z) = 1 @f$——输出逐样本等于输入，无滤波作用、
 * 内部状态为零。
 *
 * @param[out] filter  滤波器对象指针。
 */
void biquad_filter_set_empty(biquad_filter_t *filter);

/**
 * @brief 用给定的 z 域系数初始化 biquad 滤波器。
 *
 * 系数归一化使 @p den_z[0] 变为 1.0。以下任一情形成立时，滤波器被静默替换为
 * 单位直通（恒等）并返回 0——宁可直通也不要发散：
 * - @p den_z[0] 为零或非有限；
 * - 任一系数非有限；
 * - 极点落在单位圆外（不稳定，Jury 三条件）；
 * - 任一极点半径 > 0.99995（1 − r < 5e-5，会振铃 ≥ 10⁴ 样本）。
 *
 * @param[out] filter  滤波器对象指针。
 * @param[in]  num_z   z 域分子系数   [b0, b1, b2]。
 * @param[in]  den_z   z 域分母系数   [a0, a1, a2]。
 * @return             部署成功返回 1，降级为恒等返回 0。
 */
uint8_t biquad_filter_init(biquad_filter_t *filter, const float num_z[3],
                           const float den_z[3]);

/**
 * @brief 处理一个输入样本并返回滤波输出。
 *
 * static inline（header-only）：输出计算与状态更新融合，@p w[] 只装载一次，
 * 消除 MCU 上每节的函数调用开销。
 *
 * @note 热路径刻意不做 NaN/Inf 防护（每样本分支开销）：非有限输入毒化状态直至
 *       reset，调用方应在源头清洗。**不可重入**：同一对象被 ISR 与主循环共享时
 *       需自行关中断或双缓冲（或让 ISR 只置标志、主循环统一推状态）。
 *
 * @param[in,out] filter  滤波器对象指针。
 * @param[in]     input   当前输入样本。
 * @return                滤波输出 @f$ y[n] @f$。
 */
static inline float biquad_filter_update(biquad_filter_t *filter, float input)
{
    /* 移入新样本前先快照当前状态（这两个值将成为 w[n-1]、w[n-2]）。 */
    const float w1 = filter->w[0];
    const float w2 = filter->w[1];

    /* 缓存系数——避免每次访问都通过指针重新装载。 */
    const float a1 = filter->den_z[1];
    const float a2 = filter->den_z[2];
    const float b0 = filter->num_z[0];
    const float b1 = filter->num_z[1];
    const float b2 = filter->num_z[2];

    /* 计算新状态：w[n] = x[n] − a1·w[n-1] − a2·w[n-2] */
    const float w0 = input - a1 * w1 - a2 * w2;

    /* 状态写回内存。 */
    filter->w[2] = w2;
    filter->w[1] = w1;
    filter->w[0] = w0;

    /* 输出：y[n] = b0·w[n] + b1·w[n-1] + b2·w[n-2] */
    return b0 * w0 + b1 * w1 + b2 * w2;
}

/**
 * @brief 返回当前输出（不推进状态）。
 *
 * @param[in] filter  滤波器对象指针。
 * @return            当前输出 @f$ y[n] @f$。
 */
float biquad_filter_get_output(const biquad_filter_t *filter);

/**
 * @brief 由状态反推当前输入。
 *
 * 调试或级联时有用——从内部 @p w[] 状态与分母系数恢复 @f$ x[n] @f$。
 *
 * @param[in] filter  滤波器对象指针。
 * @return            反推的输入 @f$ x[n] @f$。
 */
float biquad_filter_get_input(const biquad_filter_t *filter);

/**
 * @brief 把滤波器复位到稳态平衡。
 *
 * 计算状态向量 @p w[]，使常值输入 @p equilibrium 产生相同的常值
 * 输出（DC 增益匹配）。
 *
 * @note 若 @f$ 1 + a_1 + a_2 = 0 @f$（如纯积分器），稳态无定义，
 *       此时状态强制清零。非有限的 @p equilibrium 同样把状态强制
 *       清零，而不是用 NaN 毒化。@f$ 1 + a_1 + a_2 @f$ 小到让
 *       @f$ w_{ss} = x/(1+a_1+a_2) @f$ 溢出成非有限值时，同样清零。
 *
 * @param[in,out] filter       滤波器对象指针。
 * @param[in]     equilibrium  稳态常值输入。
 */
void biquad_filter_reset(biquad_filter_t *filter, float equilibrium);

/**
 * @brief 处理一个样本，流经 biquad 节级联。
 *
 * 高阶滤波器族 per-order @p _update 函数共享的核心。以 static
 * inline 定义，每样本路径零函数调用；@p num_sections 必须是
 * 编译期字面量（各滤波器族传入 X-macro 节数），循环边界在
 * 编译期折叠。
 *
 * @param[in,out] sections       biquad 节级联。
 * @param[in]     num_sections   节数（编译期字面量）。
 * @param[in]     input          当前输入样本。
 * @return                       滤波输出样本。
 */
static inline float biquad_cascade_update(biquad_filter_t *sections,
                                          uint8_t num_sections, float input)
{
    float x = input;
    for (uint8_t i = 0; i < num_sections; i++) {
        x = biquad_filter_update(&sections[i], x);
    }
    return x;
}

/**
 * @brief 把 biquad 节级联复位到稳态。
 *
 * 高阶滤波器族 per-order @p _reset 函数共享的核心。每节的稳态
 * 以前一节稳态输出为平衡值逐节计算。
 *
 * @param[in,out] sections       biquad 节级联。
 * @param[in]     num_sections   节数（编译期字面量）。
 * @param[in]     equilibrium    稳态常值输入。
 */
static inline void biquad_cascade_reset(biquad_filter_t *sections,
                                        uint8_t num_sections, float equilibrium)
{
    float x = equilibrium;
    for (uint8_t i = 0; i < num_sections; i++) {
        biquad_filter_reset(&sections[i], x);
        x = biquad_filter_get_output(&sections[i]);
    }
}


#ifdef __cplusplus
}
#endif

#endif /* BIQUAD_FILTER_H_ */
