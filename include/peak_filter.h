/**
 * @file    peak_filter.h
 * @brief   二阶峰值（Peak / boost）滤波器。
 *
 * 实现连续传递函数
 * @f[
 *   H(s) = \frac{s^2 + 2\xi     \omega_0 s + \omega_0^2}
 *                {s^2 + 2\xi/g   \omega_0 s + \omega_0^2}
 * @f]
 * @f$ \omega_0 = 2\pi f_0 @f$ 为峰值中心角频率，@f$ g > 1 @f$ 为**线性峰值
 * 增益**。在 @f$ s = j\omega_0 @f$ 处 @f$ |H| = g @f$；直流与高频增益均为 1。
 *
 * @note g 的语义与 @ref notch_filter_t 严格互为倒数：两者满足
 *       @f$ H_{\text{peak}}(\xi, g) \equiv 1/H_{\text{notch}}(\xi, 1/g) @f$，
 *       **逐点**成立（模拟域如此，双线性之后依然如此）。因此本族全部数值
 *       分析都用一个统一参数：**等效陷波深度** @f$ d = 1/g \in (0, 1) @f$——
 *       g 越大，等效陷波越深、越难被 f32 承载。
 *
 * 零点阻尼 @f$ \xi @f$（宽）、极点阻尼 @f$ \xi/g @f$（窄），因此峰值的高度与
 * 尖锐度全部来自**极点**——正好与陷波器相反。@f$ g \to \infty @f$ 时极点阻尼
 * 趋近 0，这是本族唯一会真正发散的方向，靠 @ref biquad_filter_init 的极点
 * 半径裕量闸 fail-closed 拦住。
 *
 * 设计流程：预畸变 @f$ \omega_0 @f$ → 双线性变换 → 直接写入单个
 * @ref biquad_filter_t。与陷波器同为二阶、无配对问题，不走
 * @ref design_filter / @ref zpk2sos。
 *
 * @note 预畸只保证**中心频率与增益**精确（与陷波器同一条语义）。数字域的
 *       峰宽 = @f$ 2\xi f_0/\sqrt{g} @f$（@f$ \sqrt{g} @f$ 电平处的全宽）再乘
 *       一个随 @f$ f_0/f_s @f$ 变化的扭曲因子；@p xi 始终按模拟原型的带宽
 *       因子解释，不随 @p fs 重标定。
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
 * 与 @ref notch_filter_t 同构（阶数恒定、单节），字段语义一一对应，只有
 * @p g 的解释不同：这里是 @f$ g > 1 @f$ 的线性增益。
 */
typedef struct {
    uint8_t  valid;        /**< 1 = init 成功；0 = 部署直通（此时 num_sections 必为 0） */
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
 * 依次通过三道闸，任一失败即置 @p valid = 0、@p num_sections = 0（两者都被
 * 显式清零），调用方按直通处理（见 @ref peak_update）：
 *
 * 1. **参数闸**：@p fs > 0，@f$ 0 < f_0 < f_s/2 @f$，@p xi > 0，@f$ g > 1 @f$。
 *    判据写成 `!(x > 0)` 以拦下 NaN。@f$ g \le 1/0.9999 @f$（浅于 0.001 dB）
 *    按直通处理，与 @ref notch_init 的边界互为镜像。
 * 2. **f32 设计参数范围闸**：@f$ \xi d (\omega_0/K)^2 \ge 5\times10^{-8} @f$
 *    （@f$ d = 1/g @f$）。**必须用 d 而不是 g**：病态程度由等效陷波深度决定，
 *    用 g 会让闸门随增益被放大约 @f$ g^2 @f$ 倍。
 * 3. **逐级 fail-closed**：@ref biquad_filter_init 的三道闸与
 *    @ref check_cascade_gains 的增益窗口。本族**极点半径闸通常更紧**，实用上限
 *    约 @f$ g \le 4\times10^{4}\xi\tan(\pi f_0/f_s) @f$。
 *
 * 判据推导、两条闸的先后关系与标定表见 `test/test_peak.c` 的「数值闸边界」段。
 * `PEAK_GATE` 与 `NOTCH_GATE` 取相同值但各自独立定义，改一个需同步另一个。
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
 * @note 与陷波器同理：DC 增益为 1 要求稳态状态量在 @f$ 10^5 e @f$ 量级，
 *       复位后最初几个输出在 f32 下带绝对误差，随后收敛。
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
