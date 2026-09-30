/**
 * @file    notch_filter.h
 * @brief   二阶陷波（Notch）滤波器。
 *
 * 实现连续传递函数
 * @f[
 *   H(s) = \frac{s^2 + 2\xi g \omega_0 s + \omega_0^2}
 *                {s^2 + 2\xi   \omega_0 s + \omega_0^2}
 * @f]
 * @f$ \omega_0 = 2\pi f_0 @f$ 为陷波中心角频率，@f$ g \in (0, 1) @f$ 为陷波
 * 深度，@f$ \xi > 0 @f$ 为带宽因子。在 @f$ s = j\omega_0 @f$ 处
 * @f$ |H| = g @f$；直流与高频增益均为 1。
 *
 * 零点与极点的阻尼系数分别为 @f$ \xi g @f$ 与 @f$ \xi @f$：零点贴近 jω 轴
 * 挖出陷波，极点把它托回 0 dB——陷波的尖锐度全部来自零点。
 *
 * @f$ g = 1 @f$ 时分子与分母逐项相等，@f$ H \equiv 1 @f$——数学上就是直通，
 * 因此按**直通部署**（@p valid = 0），而不是烧一个高 Q 极点结构去实现恒等。
 *
 * 设计流程：预畸变 @f$ \omega_0 @f$ → 双线性变换 → 直接写入单个
 * @ref biquad_filter_t。二阶节只有一对一阶分子/分母，没有配对问题，因此
 * 不走 @ref design_filter / @ref zpk2sos 那条零极点管线。
 *
 * @note 预畸只保证**中心频率与深度**精确。双线性变换对整个频率轴非线性，
 *       数字域的陷波宽度 = @f$ 2\xi g f_0 @f$ 再乘一个随 @f$ f_0/f_s @f$
 *       变化的扭曲因子；@p xi 始终按模拟原型的带宽因子解释，不随 @p fs
 *       重标定。
 */

#ifndef NOTCH_FILTER_H_
#define NOTCH_FILTER_H_

#include <stdint.h>
#include "biquad_filter.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief 二阶陷波滤波器对象。
 *
 * 固定二阶、单节，没有 Butterworth / Chebyshev 那套阶数与类型族。
 */
typedef struct {
    uint8_t  valid;        /**< 1 = init 成功；0 = 部署直通（此时 num_sections 必为 0） */
    uint8_t  num_sections; /**< 活跃 biquad 节数（0 或 1） */
    float    f0;           /**< 陷波中心频率（Hz） */
    float    xi;           /**< 带宽因子 */
    float    g;            /**< 陷波深度 (0, 1)；接近 1 视为直通 */
    float    fs;           /**< 采样频率（Hz） */
    biquad_filter_t sections[1]; /**< 单个二阶节 */
} notch_filter_t;

/**
 * @brief 设计并初始化为指定参数的陷波滤波器。
 *
 * 依次通过三道闸，任一失败即置 @p valid = 0、@p num_sections = 0（两者都被
 * 显式清零），调用方按直通处理（见 @ref notch_update）：
 *
 * 1. **参数闸**：@p fs > 0，@f$ 0 < f_0 < f_s/2 @f$，@p xi > 0，
 *    @f$ 0 < g < 0.9999 @f$。判据写成 `!(x > 0)` 以拦下 NaN。
 * 2. **f32 设计参数范围闸**：@f$ \xi g (\omega_0/K)^2 \ge 5\times10^{-8} @f$
 *    （@f$ K = 2f_s @f$）。低于此值 f32 系数承载不住陷波深度，按参数直接拒，
 *    不放行已知会失真的设计。判据、标定表与范围界定见 `test/test_notch.c`
 *    的「数值闸边界」段。
 * 3. **逐级 fail-closed**：@ref biquad_filter_init 的系数/稳定性/极点半径
 *    检查，与 @ref check_cascade_gains 的 DC/Nyquist 增益窗口。
 *
 * @param[out] f   滤波器对象指针。
 * @param[in]  f0  陷波中心频率（Hz）。
 * @param[in]  xi  带宽因子（@f$ \xi > 0 @f$）。
 * @param[in]  g   陷波深度（@f$ 0 < g < 0.9999 @f$；再浅即按直通处理）。
 * @param[in]  fs  采样频率（Hz）。
 */
void notch_init(notch_filter_t *f, float f0, float xi, float g, float fs);

/**
 * @brief 处理一个输入样本并返回滤波输出。
 *
 * static inline（header-only），每样本路径零函数调用，与其余各族 `_update`
 * 同一约定。滤波器无效时原样返回 @p input（直通）。
 *
 * @note 与 biquad 层一致，热路径无 NaN/Inf 防护；非有限输入会毒化状态直至
 *       reset，调用方在源头清洗。不可重入：ISR 与主循环共享同一对象时需关
 *       中断或双缓冲。
 *
 * @param[in,out] f      滤波器对象指针。
 * @param[in]     input  当前输入样本。
 * @return               滤波输出。
 */
static inline float notch_update(notch_filter_t *f, float input)
{
    if (!f->valid) return input;
    return biquad_cascade_update(f->sections, 1, input);
}

/**
 * @brief 把陷波滤波器复位到稳态。
 *
 * 滤波器无效（@p !valid）时为空操作。
 *
 * @note 陷波 DC 增益为 1，稳态状态量 @f$ w_{ss} = e/(1+a_1+a_2) @f$ 很大
 *       （@f$ f_0 = 20 @f$ Hz、@f$ f_s = 48 @f$ kHz 时约
 *       @f$ 1.5 \times 10^5 e @f$）：复位后最初几个输出在 f32 下带绝对误差，
 *       随后收敛。
 *
 * @param[in,out] f           滤波器对象指针。
 * @param[in]     equilibrium 稳态常值输入。
 */
static inline void notch_reset(notch_filter_t *f, float equilibrium)
{
    if (!f->valid) return;
    biquad_cascade_reset(f->sections, 1, equilibrium);
}

#ifdef __cplusplus
}
#endif

#endif /* NOTCH_FILTER_H_ */
