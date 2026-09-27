/**
 * @file    notch_filter.h
 * @brief   二阶陷波（Notch）滤波器。
 *
 * 实现连续传递函数
 * @f[
 *   H(s) = \frac{s^2 + 2\xi g \omega_0 s + \omega_0^2}
 *                {s^2 + 2\xi   \omega_0 s + \omega_0^2}
 * @f]
 * 其中 @f$ \omega_0 = 2\pi f_0 @f$ 为陷波中心角频率，@f$ g \in (0, 1) @f$
 * 为陷波深度，@f$ \xi > 0 @f$ 为带宽因子。在 @f$ s = j\omega_0 @f$ 处
 * @f$ |H| = g @f$；直流与高频增益均为 1。
 *
 * @f$ g = 1 @f$ 时分子与分母逐项相等，@f$ H \equiv 1 @f$——数学上就是直通。
 * 这种设计按**直通部署**（@p valid = 0）而不是烧一个节去实现恒等：
 * 一边是什么都不做，一边是把高 Q 极点结构的状态推到 @f$ 10^5 @f$ 量级
 * 再精确相消回单位增益，后者纯属给自己找数值麻烦。
 *
 * 零点与极点是一对共轭复根，阻尼系数分别为 @f$ \xi g @f$ 与 @f$ \xi @f$：
 * 零点贴近 jω 轴挖出陷波，极点把它托回 0 dB。零点的 Q = 1/(2ξg)，
 * 比极点的 Q = 1/(2ξ) 高 1/g 倍 —— 陷波的尖锐度全部来自零点。
 *
 * 设计流程：预畸变 @f$ \omega_0 @f$ → 双线性变换 → 直接写入单个
 * @ref biquad_filter_t。二阶滤波器只有一对一阶分子/分母，没有配对问题，
 * 因此不走 @ref design_filter / @ref zpk2sos 那条零极点管线。
 *
 * @note 预畸只保证**中心频率与深度**精确。双线性变换对整个频率轴是非线性
 *       的，数字域的陷波宽度 = @f$ 2\xi g f_0 @f$ 再乘以一个随
 *       @f$ f_0/f_s @f$ 变化的扭曲因子：@f$ f_0/f_s \lesssim 0.02 @f$ 时
 *       误差 < 1%，@f$ f_0/f_s = 0.2 @f$ 时窄约 25%，再高则完全走样。
 *       @p xi 始终按模拟原型的带宽因子解释，不随 f_s 重标定。
 *
 * @note reset 的稳态量：陷波器 DC 增益为 1，故
 *       @f$ w_{ss} = e / (1 + a_1 + a_2) @f$ 是一个很大的数
 *       （@f$ f_0 = 20 @f$ Hz、@f$ f_s = 48 @f$ kHz 时约 @f$ 1.5 \times 10^5 e @f$）。
 *       这是正确的——状态必须足够大才能相消出单位输出——但复位后最初几个
 *       输出在 f32 下带约 1% 量级的绝对误差，随后收敛到稳态。
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
 * 阶数恒定（二阶，单节），因此没有 Butterworth / Chebyshev 那套
 * X-macro 阶数与类型族——就一个结构体、一个 init。
 */
typedef struct {
    uint8_t  valid;        /**< 1 = init 成功；0 = 部署直通 */
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
 * 依次通过三道闸，任一失败即置 @p valid = 0、@p num_sections = 0，
 * 调用方按直通处理（见 @ref notch_update）：
 *
 * 1. **参数闸**：@p fs > 0，@f$ 0 < f_0 < f_s/2 @f$，@p xi > 0，
 *    @f$ 0 < g < 0.9999 @f$。判据写成 `!(x > 0)` 而非 `x <= 0`——NaN 与
 *    任何数比较恒假，前者拦得住 NaN，后者会放它过去。
 *    @f$ g \ge 0.9999 @f$（陷波浅于 0.001 dB，等同不陷波）同样按直通处理，
 *    详见下文对 @f$ g = 1 @f$ 的说明。
 * 2. **f32 设计参数范围闸**：
 *    @f$ \xi \cdot g \cdot (\omega_0/K)^2 \ge 5 \times 10^{-8} @f$，
 *    其中 @f$ K = 2 f_s @f$。低于此值时 f32 系数承载不住陷波深度，谷底会
 *    静默变浅、中心偏移，而 @p valid 仍为 1。闸内跨全部实测参数组合深度
 *    误差 < 20%，闸外可达 12 倍以上。详见 `src/notch_filter.c` 文件头。
 * 3. **逐级 fail-closed**：@ref biquad_filter_init 的系数/稳定性/极点
 *    半径检查，以及 @ref check_cascade_gains 的 DC/Nyquist 增益窗口。
 *
 * @note **范围界定**：本库不负责修复、补偿或消除 f32 系数量化误差——那是
 *       数值类型的固有代价，属于调用方的选择范围。第 2 道闸不试图修好什么，
 *       它只拒绝已知会**静默失真**的参数范围，性质与 @ref biquad_filter_init
 *       里那条极点半径裕量闸（@f$ 1 - r < 5 \times 10^{-5} @f$）相同。
 *       若调用方接受深度失真，可以绕过本闸直接用 @ref biquad_c2d_bilinear
 *       自己造系数。
 *
 * @note 这道闸是**预测式**的：缺陷无法在 f32 下事后验证。在被拒的
 *       参数上求 @f$ |H(j\omega_0)| @f$ 本身就是那个病态运算，算出来的是
 *       舍入噪声；而 DC/Nyquist 增益无论设计多坏都精确是 1/1
 *       （分子分母在 @f$ z = \pm 1 @f$ 处逐项相消，是结构保证）。
 *
 * @note 实际可用区间是数值闸与 @ref biquad_filter_init 的极点半径裕量
 *       （@f$ 1 - r < 5 \times 10^{-5} @f$，对应
 *       @f$ \xi \cdot (f_0/f_s) > 7.96 \times 10^{-6} @f$）的**交集**。
 *       绝大多数参数下数值闸更紧；只有 @f$ \xi < 0.0125\,g @f$ 时极点裕量才先
 *       触发（两条边界在 @f$ \xi = 0.0125\,g @f$、
 *       @f$ \tan(\pi f_0/f_s) = 2 \times 10^{-3}/g @f$ 处相交）。两者都
 *       fail-closed，只是被拒的原因不同。
 *
 *       实用刻度：@f$ \xi = 0.05 @f$、@f$ g = 0.1 @f$ 时数值闸要求
 *       @f$ f_0/f_s \ge 1.007 \times 10^{-3} @f$，即 48 kHz 下
 *       @f$ f_0 \ge 48 @f$ Hz、8 kHz 下 @f$ f_0 \ge 8 @f$ Hz。
 *       工频陷波（50/60 Hz @ 48 kHz）在这个区间内。
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
 * 以 static inline（header-only）定义，每样本路径零函数调用——与
 * Butterworth / Chebyshev 族的 `_update` 同一约定。
 *
 * 滤波器无效（@p !valid）时原样返回 @p input（直通）。
 *
 * @note 与 biquad 层一致，热路径无 NaN/Inf 防护；非有限输入会毒化
 *       状态直至 reset。调用方在源头清洗。
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
