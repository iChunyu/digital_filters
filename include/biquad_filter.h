/**
 * @file    biquad_filter.h
 * @brief   二阶 IIR 滤波器（直接 II 型／规范型）。
 *
 * 实现离散传递函数：
 * @f[
 *   H(z) = \frac{b_0 + b_1 z^{-1} + b_2 z^{-2}}
 *                {1   + a_1 z^{-1} + a_2 z^{-2}}
 * @f]
 *
 * 规范型直接 II 型仅用 3 个状态变量（@p w[0..2]），
 * 是二阶节最省内存的实现形式。
 */

#ifndef BIQUAD_FILTER_H_
#define BIQUAD_FILTER_H_

#include <stdint.h>

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
 * 把连续传递函数
 * @f[
 *   H(s) = \frac{B_0 + B_1 s + B_2 s^2}
 *                {A_0 + A_1 s + A_2 s^2}
 * @f]
 * 通过代换 @f$ s = 2 f_s \frac{1 - z^{-1}}{1 + z^{-1}} @f$
 * 映射到离散等价形式
 * @f[
 *   H(z) = \frac{b_0 + b_1 z^{-1} + b_2 z^{-2}}
 *                {a_0 + a_1 z^{-1} + a_2 z^{-2}}
 * @f]
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
 * 系数归一化使 @p den_z[0] 变为 1.0。以下任一情形成立时，滤波器
 * 被静默替换为单位直通（恒等）并返回 0——宁可直通也不要发散：
 * - @p den_z[0] 为零或非有限；
 * - 任一系数非有限；
 * - 极点落在单位圆外（不稳定，Jury 三条件）；
 * - 任一极点半径超过 0.99995（1 − r < 5e-5，振铃 ≥ 10⁴ 样本：
 *   共轭对看 |a2| > 0.9999，一阶节看 |a1| > 0.99995，不等实根对
 *   用通用求根公式）。
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
 * 内部状态前进一个时间步。
 *
 * 以 static inline（header-only）定义，消除资源受限 MCU 上每节
 * 函数调用开销；输出计算与状态更新融合，@p w[] 只装载一次。
 *
 * @note 非有限输入（NaN/Inf）会毒化状态向量，此后每个输出都为
 *       非有限值，直到 reset 恢复。这里刻意不做防护——MCU 热路径
 *       上每样本的分支开销不可接受。接入不可信或传感器数据的
 *       调用方应在源头清洗。
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
 *       清零，而不是用 NaN 毒化。
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
