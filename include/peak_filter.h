/**
 * @file    peak_filter.h
 * @brief   二阶峰值（Peak / boost）滤波器。
 *
 * 实现连续传递函数
 * @f[
 *   H(s) = \frac{s^2 + 2\xi     \omega_0 s + \omega_0^2}
 *                {s^2 + 2\xi/g   \omega_0 s + \omega_0^2}
 * @f]
 * 其中 @f$ \omega_0 = 2\pi f_0 @f$ 为峰值中心角频率，@f$ g > 1 @f$
 * 为**线性峰值增益**。在 @f$ s = j\omega_0 @f$ 处 @f$ |H| = g @f$；
 * 直流与高频增益均为 1。
 *
 * @note g 的语义与 @ref notch_filter_t 严格互为倒数。两者满足
 * @f[
 *   H_{\text{peak}}(\xi, g) \equiv \frac{1}{H_{\text{notch}}(\xi, 1/g)}
 * @f]
 * ——**逐点**成立（模拟域如此，双线性之后依然如此，因为双线性是点对点
 * 映射）。所以本文件里所有数值分析都用一个统一参数：**等效陷波深度**
 * @f$ d = 1/g \in (0, 1) @f$。g 越大，等效陷波越深、越难被 f32 承载。
 *
 * 零点阻尼 @f$ \xi @f$（宽）、极点阻尼 @f$ \xi/g @f$（窄），
 * 因此峰值的高度与尖锐度全部来自**极点**——正好与陷波器相反
 * （陷波器的尖锐度来自零点）。@f$ g \to \infty @f$ 时极点阻尼
 * 趋近 0，这是本族唯一会真正发散的方向，靠 @ref biquad_filter_init
 * 的极点半径裕量闸 fail-closed 拦住。
 *
 * 设计流程：预畸变 @f$ \omega_0 @f$ → 双线性变换 → 直接写入单个
 * @ref biquad_filter_t。与陷波器同为二阶、无配对问题，不走
 * @ref design_filter / @ref zpk2sos。
 *
 * @note 预畸只保证**中心频率与增益**精确（与陷波器同一条语义）。
 *       双线性对整个频率轴非线性，数字域的峰宽 = @f$ 2\xi f_0 @f$
 *       （@f$ \sqrt{g} @f$ 电平处）再乘一个随 @f$ f_0/f_s @f$ 变化的
 *       扭曲因子；@p xi 始终按模拟原型的带宽因子解释。
 *
 * @note 复位语义与陷波器不同：峰值滤波器 DC 增益为 1，故
 *       @f$ w_{ss} = e/(1+a_1+a_2) @f$ 是**很大**的数
 *       （@f$ f_0 = 20 @f$ Hz、@f$ f_s = 48 @f$ kHz、@f$ g = 10 @f$ 时
 *       约 @f$ 1.5 \times 10^7 e @f$）。参数合法时这是对的——状态必须
 *       足够大才能相消出正确的输出；超出数值包络时极点半径闸会先拒。
 */

#ifndef PEAK_FILTER_H_
#define PEAK_FILTER_H_

#include <stdint.h>
#include "biquad_filter.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief 二阶峰值滤波器对象。
 *
 * 与 @ref notch_filter_t 同构（阶数恒定、单节），字段语义一一对应，
 * 只有 @p g 的解释不同：这里是 @f$ g > 1 @f$ 的线性增益。
 */
typedef struct {
    uint8_t  valid;        /**< 1 = init 成功；0 = 部署直通 */
    uint8_t  num_sections; /**< 活跃 biquad 节数（0 或 1） */
    float    f0;           /**< 峰值中心频率（Hz） */
    float    xi;           /**< 带宽因子 */
    float    g;            /**< 线性峰值增益（> 1）；接近 1 视为直通 */
    float    fs;           /**< 采样频率（Hz） */
    biquad_filter_t sections[1]; /**< 单个二阶节 */
} peak_filter_t;

/**
 * @brief 设计并初始化为指定参数的峰值滤波器。
 *
 * 依次通过三道闸，任一失败即置 @p valid = 0、@p num_sections = 0，
 * 调用方按直通处理（见 @ref peak_update）：
 *
 * 1. **参数闸**：@p fs > 0，@f$ 0 < f_0 < f_s/2 @f$，@p xi > 0，
 *    @f$ g > 1 @f$。判据写成 `!(x > 0)` 而非 `x <= 0`——NaN 与任何数
 *    比较恒假，前者拦得住 NaN，后者会放它过去。
 *    @f$ g \le 1/0.9999 \approx 1.0001 @f$（峰值浅于 0.001 dB，等同不
 *    boost）同样按直通处理，与 @ref notch_init 的 @f$ g \ge 0.9999 @f$
 *    边界互为镜像。
 * 2. **f32 设计参数范围闸**：@f$ \xi \cdot d \cdot (\omega_0/K)^2
 *    \ge 5 \times 10^{-8} @f$，其中 @f$ d = 1/g @f$、@f$ K = 2 f_s @f$。
 *    **必须用 d 而不是 g**：本滤波器的病态程度由等效陷波深度决定，
 *    用 g 判会让闸门随增益被放大约 @f$ g^2 @f$ 倍——g = 10 时
 *    48 kHz 下 20 Hz 的 boost 会"过闸"，而它的镜像陷波在同一点是被
 *    拒的。详见 `src/peak_filter.c` 文件头。
 * 3. **逐级 fail-closed**：@ref biquad_filter_init 的系数/稳定性/极点
 *    半径检查，以及 @ref check_cascade_gains 的 DC/Nyquist 增益窗口。
 *    对峰值滤波器而言，**极点半径闸通常是更紧的那一道**：极点阻尼
 *    为 @f$ \xi/g @f$，故 @f$ 1 - r \approx 2\xi\tan(\pi f_0/f_s)/g @f$，
 *    要它 @f$ \ge 5\times10^{-5} @f$ 即
 *    @f$ g \le 4\times10^{4}\,\xi\tan(\pi f_0/f_s) @f$。
 *    实用刻度（@f$ \xi = 0.05 @f$）：@f$ f_0/f_s = 0.05 @f$ 时
 *    g 上限约 317（50 dB）；@f$ f_0/f_s = 2\times10^{-3} @f$ 时
 *    上限只有约 12.6（22 dB）。
 *
 * @note **范围界定**：本库不负责修复、补偿或消除 f32 系数量化误差。
 *       第 2 道闸只拒绝已知会**静默失真**的参数范围（峰高变矮、
 *       中心偏移，而 @p valid 仍为 1），性质与 @ref notch_init 的
 *       同名闸门完全相同——两者共用同一个标定常数 5e-8。
 *
 * @note 实测标定（f0 = 20 Hz、ξ = 0.05、g = 10，f32 系数 + f32 DF-II
 *       状态递推，扫 f0 的 ±1% 邻域取最坏）：fs = 1 kHz 误差 < 3e-5、
 *       8 kHz 为 2.8e-4、48 kHz 时峰高在 7.11 ~ 11.06 之间（理想 10，
 *       最坏 29%）。镜像核对：同点 @f$ g = 0.1 @f$ 的陷波深度在
 *       0.126 ~ 0.166 之间（理想 0.1，最坏 66%）——两者就是互为倒数
 *       的同一个病态。48 kHz/20 Hz 这一档被极点半径闸拒绝（@f$ 1-r @f$
 *       = 1.31e-5 < 5e-5）。
 *
 * @note 数字域的可用区间是数值闸、极点半径闸两者的**交集**，两者都
 *       fail-closed。@f$ \xi = 0.05 @f$、@f$ g = 10 @f$（@f$ d = 0.1 @f$）
 *       时数值闸要求 @f$ f_0/f_s \ge 1.007\times10^{-3} @f$（与陷波器
 *       同值，因为用的是同一个 d），而半径闸要求
 *       @f$ f_0/f_s \ge 1.59\times10^{-3} @f$——对峰值滤波器，
 *       半径闸先触发。
 *
 * @param[out] f   滤波器对象指针。
 * @param[in]  f0  峰值中心频率（Hz）。
 * @param[in]  xi  带宽因子（@f$ \xi > 0 @f$）。
 * @param[in]  g   线性峰值增益（@f$ g > 1 @f$；@f$ g \le 1.0001 @f$
 *                 按直通处理）。
 * @param[in]  fs  采样频率（Hz）。
 */
void peak_init(peak_filter_t *f, float f0, float xi, float g, float fs);

/**
 * @brief 处理一个输入样本并返回滤波输出。
 *
 * 以 static inline（header-only）定义，每样本路径零函数调用——与
 * 其余各族 `_update` 同一约定。
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
static inline float peak_update(peak_filter_t *f, float input)
{
    if (!f->valid) return input;
    return biquad_cascade_update(f->sections, 1, input);
}

/**
 * @brief 把峰值滤波器复位到稳态。
 *
 * 滤波器无效（@p !valid）时为空操作。
 *
 * @param[in,out] f           滤波器对象指针。
 * @param[in]     equilibrium 稳态常值输入。
 */
static inline void peak_reset(peak_filter_t *f, float equilibrium)
{
    if (!f->valid) return;
    biquad_cascade_reset(f->sections, 1, equilibrium);
}

#ifdef __cplusplus
}
#endif

#endif /* PEAK_FILTER_H_ */
