#ifndef FILTER_UTILS_H_
#define FILTER_UTILS_H_

#include <math.h>
#include <stdint.h>
#include "biquad_filter.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    FILTER_LOWPASS,
    FILTER_HIGHPASS,
    FILTER_BANDPASS,
    FILTER_BANDSTOP
} filter_type_e;

typedef struct {
    float re;
    float im;
} complex_t;

/**
 * @brief 把数字截止频率预畸变到模拟等效频率。
 *
 * 采用双线性变换的频率映射：
 *   f_analog = fs/pi * tan(pi * f_digital / fs)
 *
 * @param fd  目标数字截止频率（Hz），须满足 0 < fd < fs/2。
 * @param fs  采样频率（Hz）。
 * @return    等效模拟截止频率（Hz）。
 */
float prewarp(float fd, float fs);

/**
 * @brief 模拟低通频率变换：零极点按 wc 缩放。
 *
 * @param[in,out] poles  复数极点数组（np 个）。
 * @param[in]     np     极点数。
 * @param[in,out] zeros  复数零点数组（nz 个）。
 * @param[in]     nz     零点数。
 * @param[in]     wc     模拟截止频率（rad/s，即 2·π·prewarp）。
 */
void analog_lp_transform(complex_t *poles, uint8_t np,
                         complex_t *zeros, uint8_t nz, float wc);

/**
 * @brief 模拟高通频率变换：零极点取倒数并缩放。
 *
 * 每个极点 p 变为 wc/p，每个零点 z 变为 wc/z。
 * 调用方必须在本函数之后追加 np − nz 个原点零点 (0, 0)，
 * 以对应原型中的无穷远零点。
 *
 * @param[in,out] poles  复数极点数组（np 个）。
 * @param[in]     np     极点数。
 * @param[in,out] zeros  复数零点数组（nz 个）。
 * @param[in]     nz     零点数。
 * @param[in]     wc     模拟截止频率（rad/s）。
 */
void analog_hp_transform(complex_t *poles, uint8_t np,
                         complex_t *zeros, uint8_t nz, float wc);

/**
 * @brief 双线性变换：把 s 域零极点原地映射到 z 域。
 *
 * 对每个元素应用 z = (2·fs + s) / (2·fs − s)。
 * 调用方须在之后追加 (np − nz) 个 z = −1 的零点，以对应多余的极点。
 *
 * @param[in,out] zp  复数极点或零点数组（n 个，原地修改）。
 * @param[in]     n   元素个数。
 * @param[in]     fs  采样频率（Hz）。
 */
void bilinear_transform(complex_t *zp, uint8_t n, float fs);

/**
 * @brief 模拟带通频率变换。
 *
 * 代换 s → (s² + ω₀²) / (ξ·s)，阶数翻倍。每个极点与零点经二次
 * 公式一分为二。函数追加 (np_old − nz_old) 个原点零点 (0, 0)，
 * 对应原型的无穷远零点。
 *
 * @param[in,out] poles  极点数组（容量 ≥ 2·*np，原地变换）。
 * @param[in,out] np     输入：原型极点数；输出：2 × 输入。
 * @param[in,out] zeros  零点数组（容量 ≥ *np + *nz，原地变换）。
 * @param[in,out] nz     输入：原型零点数；输出：*np_old + *nz_old。
 * @param[in]     w0     模拟中心频率 ω₀ = √(ω₁·ω₂)（rad/s）。
 * @param[in]     xi     带宽 ξ = ω₂ − ω₁（rad/s）。
 */
void analog_bp_transform(complex_t *poles, uint8_t *np,
                         complex_t *zeros, uint8_t *nz,
                         float w0, float xi);

/**
 * @brief 模拟带阻频率变换。
 *
 * 代换 s → ξ·s / (s² + ω₀²)，阶数翻倍。函数追加
 * 2·(np_old − nz_old) 个 ±jω₀ 处零点，对应原型的无穷远零点。
 *
 * @param[in,out] poles  极点数组（容量 ≥ 2·*np，原地变换）。
 * @param[in,out] np     输入：原型极点数；输出：2 × 输入。
 * @param[in,out] zeros  零点数组（容量 ≥ 2·*np，原地变换）。
 * @param[in,out] nz     输入：原型零点数；输出：2 × *np_old。
 * @param[in]     w0     模拟中心频率 ω₀ = √(ω₁·ω₂)（rad/s）。
 * @param[in]     xi     带宽 ξ = ω₂ − ω₁（rad/s）。
 */
void analog_bs_transform(complex_t *poles, uint8_t *np,
                         complex_t *zeros, uint8_t *nz,
                         float w0, float xi);

/**
 * @brief 计算模拟 LP→HP 或 LP→BS 变换的增益调整。
 *
 * k = k · Re(∏(−z) / ∏(−p))，在原型（变换前）的零极点上计算。
 * nz = 0 时空积 ∏(−z) 取 1。
 *
 * @param k   当前系统增益。
 * @param z   原型零点（变换前），nz 个。
 * @param nz  原型零点数。
 * @param p   原型极点（变换前），np 个。
 * @param np  原型极点数。
 * @return    调整后的增益。
 */
float zpk_hp_bs_gain(float k, const complex_t *z, uint8_t nz,
                     const complex_t *p, uint8_t np);

/**
 * @brief 计算双线性变换的增益调整。
 *
 * k_z = k · Re(∏(K − z) / ∏(K − p))，其中 K = 2·fs，z、p 为
 * 双线性代换**之前**的模拟域（s 平面）零点与极点。
 *
 * @param k   当前系统增益。
 * @param z   模拟域零点，nz 个。
 * @param nz  模拟域零点数。
 * @param p   模拟域极点，np 个。
 * @param np  模拟域极点数。
 * @param K   双线性常数 = 2·fs。
 * @return    调整后的增益。
 */
float bilinear_zpk_gain(float k, const complex_t *z, uint8_t nz,
                         const complex_t *p, uint8_t np, float K);

/**
 * @brief 计算 LP/BP 频率变换加双线性变换的合并增益调整。
 *
 * k = k · s^degree · Re(∏(K − z) / ∏(K − p))，K = 2·fs。
 *
 * s 因子与 (K − p) 除法交错折叠，避免任何中间量溢出 f32——
 * 单独计算 s^degree 在近 Nyquist 设计（wc^8 ≈ FLT_MAX）上先溢出，
 * 之后同量级的 ∏(K − p) 因子来不及抵消。
 *
 * @param k       当前系统增益。
 * @param s       每阶缩放因子：LP 为 wc（rad/s），BP 为带宽 ξ（rad/s）。
 * @param degree  原型相对阶数 = np − nz。
 * @param z       模拟域零点（频率变换后），nz 个。
 * @param nz      模拟域零点数。
 * @param p       模拟域极点（频率变换后），np 个。
 * @param np      模拟域极点数。
 * @param K       双线性常数 = 2·fs。
 * @return        调整后的增益。
 */
float bilinear_zpk_gain_scaled(float k, float s, uint8_t degree,
                               const complex_t *z, uint8_t nz,
                               const complex_t *p, uint8_t np, float K);

/**
 * @brief 把 z 域零极点数组转换为二阶节（SOS）系数。
 *
 * 用"最不利极点优先"算法配对极点与最近的零点，生成 biquad 系数。
 * 每节分子增益为 1；系统总增益 @p k 只施加到第一节分子
 * （与 scipy zpk2sos 约定一致）。
 *
 * 直接在调用方的可变零极点数组上工作，不做内部拷贝（省 ~256 字节
 * 栈——zpk2sos 在 init 调用链深处被调用，MCU 上意义显著）；
 * 数组只读不改，归属跟踪用内部 used[] 位图。
 *
 * Fail-closed：零极点集不平衡、任一元素未配对（如合成共轭）、
 * 或各节根无法复现输入零极点多重集（跨对误认领）时返回 0。
 * 调用方必须把 0 视为"设计失败 → 部署直通"。
 *
 * @param[in]  zeros  z 域零点数组（n 个）。
 * @param[in]  poles  z 域极点数组（n 个）。
 * @param[in]  n      极点数（必须等于零点数）。
 * @param[out] sos    SOS 矩阵，ceil(n/2) 行，每行 [b0,b1,b2, 1,a1,a2]。
 *                    调用方须分配 ceil(n/2) 行。
 * @param[in]  k      系统总增益，施加到 sos[0] 分子。
 * @return            SOS 节数 = ceil(n/2)，失败返回 0。
 */
uint8_t zpk2sos(complex_t *zeros, complex_t *poles, uint8_t n,
                float (*sos)[6], float k);

/**
 * @brief 共享 IIR 设计管线：模拟原型 → 频率变换 → 增益折叠 →
 *        双线性 → zpk2sos → 部署。
 *
 * Butterworth 与 Chebyshev 族共用（管线相同、原型不同）：
 * Butterworth 传 ROM 极点表 + k=1、nz=0；Chebyshev 运行时计算
 * 原型，传自己的 k 与有限零点。管线只此一份，每个 fail-closed
 * 咽喉点（k 有限非零、节数上限、逐节 init、增益窗口）恰好
 * 存在一次。
 *
 * Fail-closed：任一咽喉点失败返回 0（调用方部署直通）。
 *
 * @param[out] sections      输出 biquad 数组。
 * @param[in]  max_sections  sections 容量。
 * @param[in]  type          filter_type_e。
 * @param[in]  wc1           预畸变后的下截止频率 rad/s（2π·prewarp(fc, fs)）。
 * @param[in]  wc2           BP/BS 预畸变上截止频率 rad/s（LP/HP 不用）。
 * @param[in]  fs            采样频率（Hz）。
 * @param[in]  k             原型增益。
 * @param[in]  proto_poles   原型极点，np 个（拷贝进内部数组）。
 * @param[in]  np            原型极点数。
 * @param[in]  proto_zeros   原型有限零点，nz 个（可为 NULL）。
 * @param[in]  nz            原型有限零点数。
 * @return                   部署的节数，0 表示失败。
 */
uint8_t design_filter(biquad_filter_t *sections, uint8_t max_sections,
                      uint8_t type,
                      float wc1, float wc2, float fs,
                      float k,
                      const complex_t *proto_poles, uint8_t np,
                      const complex_t *proto_zeros, uint8_t nz);

/**
 * @brief 解析级联 DC 与 Nyquist 增益校验。
 *
 * 错但稳定的零极点配对与增益缩放错误能通过所有逐节检查（每节都
 * 有限且 Jury 稳定），却毁掉响应形状——已确认的缺陷实测 DC 增益
 * 97.9、带阻上 350× 谐振（本应 ≈ 1）。H(0) 与 H(π) 是精确有理
 * 求值（无采样、无三角，init 时 O(节数) 成本）。
 *
 * 窗口标定（在部署后的 f32 系数上实测）：接受的设计实现误差
 * ≤ ~3.4e-2（最坏：fc ≈ 6 Hz@48 kHz 窄带，f32 系数量化真实体现）；
 * 已确认的配对缺陷偏差 ≥ 0.45。±0.1 / ±0.25 窗口在两侧均保持
 * ≥ 3× 分离。
 *
 * @param[in] sections      已部署的 biquad 级联。
 * @param[in] num_sections  已部署的节数。
 * @param[in] dc_exp        期望的级联 DC 增益（0、1 或纹波边缘）。
 * @param[in] ny_exp        期望的级联 Nyquist 增益（0、1 或纹波边缘）。
 * @return                  两个增益均在容差内返回 1。
 */
uint8_t check_cascade_gains(const biquad_filter_t *sections, uint8_t num_sections,
                            float dc_exp, float ny_exp);

#ifdef __cplusplus
}
#endif

#endif /* FILTER_UTILS_H_ */
