#include "cheby_filter.h"
#include <math.h>
#include <stddef.h>

/* ================================================================== */
/*  Chebyshev 原型（运行时计算——依赖纹波）                            */
/* ================================================================== */

/**
 * @brief 运行时计算 Chebyshev I 原型极点。
 *
 * 极点角度 θ_k = π·(2k + N + 1) / (2N)，半径由 ε 经
 * μ = asinh(1/ε) / N 得到：
 *   p_k = sinh μ · cos θ_k + j · cosh μ · sin θ_k
 *
 * @param[out] poles    输出极点数组（n 个）。
 * @param[in]  n        阶数。
 * @param[in]  epsilon  通带纹波参数 ε = √(10^(rp/10) − 1)。
 */
static void cheby1_proto(complex_t *poles, uint8_t n, float epsilon)
{
    float inv_eps = 1.0f / epsilon;
    float mu = logf(inv_eps + sqrtf(inv_eps * inv_eps + 1.0f)) / (float)n;
    float sinh_mu = sinhf(mu);
    float cosh_mu = coshf(mu);

    for (uint8_t k = 0; k < n; k++) {
        float theta = (float)(2 * (k + 1) + n - 1)
                    / (2.0f * (float)n) * (float)M_PI;
        poles[k].re = sinh_mu * cosf(theta);
        poles[k].im = cosh_mu * sinf(theta);
    }
}

/**
 * @brief 运行时计算 Chebyshev II 原型极点与零点。
 *
 * 极点为 Chebyshev I 极点的倒数（conj(p1) / |p1|²）；
 * 零点 z_k = j / sin(θ_k)。
 *
 * @note 跳过 sin(θ) ≈ 0 的 θ（奇数阶该零点位于 s = ∞，不输出）。
 *       合法零点在 θ = π ± π/N 或更远处，|sin(θ)| ≥ sin(π/8)
 *       ≈ 0.38，比阈值高三个数量级。阈值只需吞掉 θ = π_float 处
 *       sinf 的求值误差：glibc 给 ~8.7e-8，软浮点 libm 差 2 ulp
 *       就会输出 j·1e7 的幻影零点，静默重塑响应（实测级联
 *       DC 增益 1.0 → 0.9998）。1e-4 给 libm 误差留三个数量级
 *       余量，同时远低于任何合法零点。
 *
 * @param[out] poles    输出极点数组（n 个）。
 * @param[out] zeros    输出零点数组（最多 n 个）。
 * @param[in]  n        阶数。
 * @param[in]  epsilon  阻带纹波参数 ε = 1/√(10^(rs/10) − 1)。
 * @return              有限零点个数（偶数阶 n，奇数阶 n−1）。
 */
static uint8_t cheby2_proto(complex_t *poles, complex_t *zeros, uint8_t n,
                             float epsilon)
{
    float inv_eps = 1.0f / epsilon;
    float mu = logf(inv_eps + sqrtf(inv_eps * inv_eps + 1.0f)) / (float)n;
    float sinh_mu = sinhf(mu);
    float cosh_mu = coshf(mu);
    uint8_t nz = 0;

    for (uint8_t k = 0; k < n; k++) {
        float theta = (float)(2 * (k + 1) + n - 1)
                    / (2.0f * (float)n) * (float)M_PI;

        /* Pole: 1 / (Chebyshev I pole) */
        float den = sinh_mu * sinh_mu * cosf(theta) * cosf(theta)
                  + cosh_mu * cosh_mu * sinf(theta) * sinf(theta);
        poles[k].re =  sinh_mu * cosf(theta) / den;
        poles[k].im = -cosh_mu * sinf(theta) / den;

        /* 零点跳过阈值 1e-4 的选取论证见函数 doxygen @note。 */
        float s = sinf(theta);
        if (fabsf(s) > 1e-4f) {
            zeros[nz].re = 0.0f;
            zeros[nz].im = 1.0f / s;
            nz++;
        }
    }
    return nz;
}

/* ================================================================== */
/*  各类型 init 辅助函数                                               */
/* ================================================================== */

/*
 * 各族 / 类型 / 阶数奇偶的期望 DC/Nyquist 增益（与 scipy 核对）：
 *   cheby1：通带边缘为 10^(−rp/20)（偶数阶），奇数阶为 1。
 *   cheby2：阻带边缘为 10^(−rs/20)（偶数阶），奇数阶为 0。
 */

/**
 * @brief Chebyshev I 纹波边缘增益：偶数阶 10^(−rp/20)，奇数阶 1。
 *
 * @param order      阶数。
 * @param ripple_db  通带纹波（dB）。
 * @return           边缘增益。
 */
static float cheby1_edge_gain(uint8_t order, float ripple_db)
{
    return (order % 2 == 0) ? powf(10.0f, -ripple_db / 20.0f) : 1.0f;
}

/**
 * @brief Chebyshev II 阻带边缘增益：偶数阶 10^(−rs/20)，奇数阶 0。
 *
 * @param order      阶数。
 * @param ripple_db  阻带衰减（dB）。
 * @return           边缘增益。
 */
static float cheby2_edge_gain(uint8_t order, float ripple_db)
{
    return (order % 2 == 0) ? powf(10.0f, -ripple_db / 20.0f) : 0.0f;
}

/* ── Chebyshev I ──────────────────────────────────────────────────── */

/**
 * @brief Chebyshev I 低通设计：原型 → 管线 → 增益校验。
 *
 * 参数越界（order 不在 1..8、fc 不在 (0, fs/2)、ripple_db ≤ 0）
 * 或校验失败返回 0（调用方部署直通）。
 * 期望级联增益：DC 1（奇数阶）或 10^(−rp/20)（偶数阶），Nyquist 0。
 *
 * @param[out] sections      biquad 节数组。
 * @param[in]  max_sections  sections 容量。
 * @param[in]  order         原型阶数。
 * @param[in]  fc            截止频率（Hz）。
 * @param[in]  fs            采样频率（Hz）。
 * @param[in]  ripple_db     通带纹波（dB）。
 * @return                   部署的节数，0 表示失败。
 */
static uint8_t cheby1_lp_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order, float fc, float fs,
                                      float ripple_db)
{
    if (order == 0 || order > 8 || fc <= 0.0f || fc >= fs * 0.5f || ripple_db <= 0.0f)
        return 0;

    float epsilon = sqrtf(powf(10.0f, ripple_db / 10.0f) - 1.0f);
    float wc = 2.0f * (float)M_PI * prewarp(fc, fs);

    complex_t poles[8];
    cheby1_proto(poles, order, epsilon);

    float k = 1.0f / zpk_hp_bs_gain(1.0f, NULL, 0, poles, order);
    if (order % 2 == 0) k /= sqrtf(1.0f + epsilon * epsilon);

    uint8_t n = design_filter(sections, max_sections,
                              FILTER_LOWPASS,
                              wc, 0.0f, fs, k,
                              poles, order, NULL, 0);
    if (n == 0) return 0;
    /* cheby1 LP: DC gain 1 (odd) / 10^(−rp/20) (even), Nyquist 0 */
    if (!check_cascade_gains(sections, n,
                             cheby1_edge_gain(order, ripple_db), 0.0f)) return 0;
    return n;
}

/**
 * @brief Chebyshev I 高通设计：原型 → 管线 → 增益校验。
 *
 * 失败返回 0（调用方部署直通）。
 * 期望级联增益：DC 0，Nyquist 1（奇数阶）或 10^(−rp/20)（偶数阶）。
 *
 * @param[out] sections      biquad 节数组。
 * @param[in]  max_sections  sections 容量。
 * @param[in]  order         原型阶数。
 * @param[in]  fc            截止频率（Hz）。
 * @param[in]  fs            采样频率（Hz）。
 * @param[in]  ripple_db     通带纹波（dB）。
 * @return                   部署的节数，0 表示失败。
 */
static uint8_t cheby1_hp_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order, float fc, float fs,
                                      float ripple_db)
{
    if (order == 0 || order > 8 || fc <= 0.0f || fc >= fs * 0.5f || ripple_db <= 0.0f)
        return 0;

    float epsilon = sqrtf(powf(10.0f, ripple_db / 10.0f) - 1.0f);
    float wc = 2.0f * (float)M_PI * prewarp(fc, fs);

    complex_t poles[8];
    cheby1_proto(poles, order, epsilon);

    float k = 1.0f / zpk_hp_bs_gain(1.0f, NULL, 0, poles, order);
    if (order % 2 == 0) k /= sqrtf(1.0f + epsilon * epsilon);

    uint8_t n = design_filter(sections, max_sections,
                              FILTER_HIGHPASS,
                              wc, 0.0f, fs, k,
                              poles, order, NULL, 0);
    if (n == 0) return 0;
    /* cheby1 HP: DC 0, Nyquist gain 1 (odd) / 10^(−rp/20) (even) */
    if (!check_cascade_gains(sections, n, 0.0f,
                             cheby1_edge_gain(order, ripple_db))) return 0;
    return n;
}

/**
 * @brief Chebyshev I 带通设计：原型 → 管线 → 增益校验。
 *
 * 失败返回 0（调用方部署直通）。
 * 期望级联增益：DC 与 Nyquist 均为 0。
 *
 * @param[out] sections      biquad 节数组。
 * @param[in]  max_sections  sections 容量。
 * @param[in]  order         原型阶数。
 * @param[in]  fc1           下带边（Hz）。
 * @param[in]  fc2           上带边（Hz）。
 * @param[in]  fs            采样频率（Hz）。
 * @param[in]  ripple_db     通带纹波（dB）。
 * @return                   部署的节数，0 表示失败。
 */
static uint8_t cheby1_bp_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order,
                                      float fc1, float fc2, float fs,
                                      float ripple_db)
{
    if (order == 0 || order > 8 || fc1 <= 0.0f || fc1 >= fs * 0.5f || ripple_db <= 0.0f)
        return 0;
    if (fc2 <= fc1 || fc2 >= fs * 0.5f) return 0;

    float epsilon = sqrtf(powf(10.0f, ripple_db / 10.0f) - 1.0f);
    float wc1 = 2.0f * (float)M_PI * prewarp(fc1, fs);
    float wc2 = 2.0f * (float)M_PI * prewarp(fc2, fs);

    complex_t poles[8];
    cheby1_proto(poles, order, epsilon);

    float k = 1.0f / zpk_hp_bs_gain(1.0f, NULL, 0, poles, order);
    if (order % 2 == 0) k /= sqrtf(1.0f + epsilon * epsilon);

    uint8_t n = design_filter(sections, max_sections,
                              FILTER_BANDPASS,
                              wc1, wc2, fs, k,
                              poles, order, NULL, 0);
    if (n == 0) return 0;
    /* cheby1 BP: DC and Nyquist gains 0 */
    if (!check_cascade_gains(sections, n, 0.0f, 0.0f)) return 0;
    return n;
}

/**
 * @brief Chebyshev I 带阻设计：原型 → 管线 → 增益校验。
 *
 * 失败返回 0（调用方部署直通）。
 * 期望级联增益：DC 与 Nyquist 均为 1（奇数阶）或 10^(−rp/20)
 * （偶数阶）。
 *
 * @param[out] sections      biquad 节数组。
 * @param[in]  max_sections  sections 容量。
 * @param[in]  order         原型阶数。
 * @param[in]  fc1           下带边（Hz）。
 * @param[in]  fc2           上带边（Hz）。
 * @param[in]  fs            采样频率（Hz）。
 * @param[in]  ripple_db     通带纹波（dB）。
 * @return                   部署的节数，0 表示失败。
 */
static uint8_t cheby1_bs_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order,
                                      float fc1, float fc2, float fs,
                                      float ripple_db)
{
    if (order == 0 || order > 8 || fc1 <= 0.0f || fc1 >= fs * 0.5f || ripple_db <= 0.0f)
        return 0;
    if (fc2 <= fc1 || fc2 >= fs * 0.5f) return 0;

    float epsilon = sqrtf(powf(10.0f, ripple_db / 10.0f) - 1.0f);
    float wc1 = 2.0f * (float)M_PI * prewarp(fc1, fs);
    float wc2 = 2.0f * (float)M_PI * prewarp(fc2, fs);

    complex_t poles[8];
    cheby1_proto(poles, order, epsilon);

    float k = 1.0f / zpk_hp_bs_gain(1.0f, NULL, 0, poles, order);
    if (order % 2 == 0) k /= sqrtf(1.0f + epsilon * epsilon);

    uint8_t n = design_filter(sections, max_sections,
                              FILTER_BANDSTOP,
                              wc1, wc2, fs, k,
                              poles, order, NULL, 0);
    if (n == 0) return 0;
    /* cheby1 BS: DC and Nyquist gains 1 (odd) / 10^(−rp/20) (even) */
    if (!check_cascade_gains(sections, n,
                             cheby1_edge_gain(order, ripple_db),
                             cheby1_edge_gain(order, ripple_db))) return 0;
    return n;
}

/* ── Chebyshev II ─────────────────────────────────────────────────── */

/**
 * @brief Chebyshev II 低通设计：原型 → 管线 → 增益校验。
 *
 * 失败返回 0（调用方部署直通）。
 * 期望级联增益：DC 1，Nyquist 0（奇数阶）或 10^(−rs/20)（偶数阶）。
 *
 * @param[out] sections      biquad 节数组。
 * @param[in]  max_sections  sections 容量。
 * @param[in]  order         原型阶数。
 * @param[in]  fc            截止频率（Hz）。
 * @param[in]  fs            采样频率（Hz）。
 * @param[in]  ripple_db     阻带衰减（dB）。
 * @return                   部署的节数，0 表示失败。
 */
static uint8_t cheby2_lp_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order, float fc, float fs,
                                      float ripple_db)
{
    if (order == 0 || order > 8 || fc <= 0.0f || fc >= fs * 0.5f || ripple_db <= 0.0f)
        return 0;

    float epsilon = 1.0f / sqrtf(powf(10.0f, ripple_db / 10.0f) - 1.0f);
    float wc = 2.0f * (float)M_PI * prewarp(fc, fs);

    complex_t poles[8], zeros[8];
    uint8_t nz = cheby2_proto(poles, zeros, order, epsilon);

    float k = 1.0f / zpk_hp_bs_gain(1.0f, zeros, nz, poles, order);

    uint8_t n = design_filter(sections, max_sections,
                              FILTER_LOWPASS,
                              wc, 0.0f, fs, k,
                              poles, order, zeros, nz);
    if (n == 0) return 0;
    /* cheby2 LP: DC gain 1, Nyquist 0 (odd) / 10^(−rs/20) (even) */
    if (!check_cascade_gains(sections, n, 1.0f,
                             cheby2_edge_gain(order, ripple_db))) return 0;
    return n;
}

/**
 * @brief Chebyshev II 高通设计：原型 → 管线 → 增益校验。
 *
 * 失败返回 0（调用方部署直通）。
 * 期望级联增益：DC 0（奇数阶）或 10^(−rs/20)（偶数阶），Nyquist 1。
 *
 * @param[out] sections      biquad 节数组。
 * @param[in]  max_sections  sections 容量。
 * @param[in]  order         原型阶数。
 * @param[in]  fc            截止频率（Hz）。
 * @param[in]  fs            采样频率（Hz）。
 * @param[in]  ripple_db     阻带衰减（dB）。
 * @return                   部署的节数，0 表示失败。
 */
static uint8_t cheby2_hp_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order, float fc, float fs,
                                      float ripple_db)
{
    if (order == 0 || order > 8 || fc <= 0.0f || fc >= fs * 0.5f || ripple_db <= 0.0f)
        return 0;

    float epsilon = 1.0f / sqrtf(powf(10.0f, ripple_db / 10.0f) - 1.0f);
    float wc = 2.0f * (float)M_PI * prewarp(fc, fs);

    complex_t poles[8], zeros[8];
    uint8_t nz = cheby2_proto(poles, zeros, order, epsilon);

    float k = 1.0f / zpk_hp_bs_gain(1.0f, zeros, nz, poles, order);

    uint8_t n = design_filter(sections, max_sections,
                              FILTER_HIGHPASS,
                              wc, 0.0f, fs, k,
                              poles, order, zeros, nz);
    if (n == 0) return 0;
    /* cheby2 HP: DC 0 (odd) / 10^(−rs/20) (even), Nyquist gain 1 */
    if (!check_cascade_gains(sections, n,
                             cheby2_edge_gain(order, ripple_db), 1.0f)) return 0;
    return n;
}

/**
 * @brief Chebyshev II 带通设计：原型 → 管线 → 增益校验。
 *
 * 失败返回 0（调用方部署直通）。
 * 期望级联增益：DC 与 Nyquist 均为 0（奇数阶）或 10^(−rs/20)
 * （偶数阶）。
 *
 * @param[out] sections      biquad 节数组。
 * @param[in]  max_sections  sections 容量。
 * @param[in]  order         原型阶数。
 * @param[in]  fc1           下带边（Hz）。
 * @param[in]  fc2           上带边（Hz）。
 * @param[in]  fs            采样频率（Hz）。
 * @param[in]  ripple_db     阻带衰减（dB）。
 * @return                   部署的节数，0 表示失败。
 */
static uint8_t cheby2_bp_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order,
                                      float fc1, float fc2, float fs,
                                      float ripple_db)
{
    if (order == 0 || order > 8 || fc1 <= 0.0f || fc1 >= fs * 0.5f || ripple_db <= 0.0f)
        return 0;
    if (fc2 <= fc1 || fc2 >= fs * 0.5f) return 0;

    float epsilon = 1.0f / sqrtf(powf(10.0f, ripple_db / 10.0f) - 1.0f);
    float wc1 = 2.0f * (float)M_PI * prewarp(fc1, fs);
    float wc2 = 2.0f * (float)M_PI * prewarp(fc2, fs);

    complex_t poles[8], zeros[8];
    uint8_t nz = cheby2_proto(poles, zeros, order, epsilon);

    float k = 1.0f / zpk_hp_bs_gain(1.0f, zeros, nz, poles, order);

    uint8_t n = design_filter(sections, max_sections,
                              FILTER_BANDPASS,
                              wc1, wc2, fs, k,
                              poles, order, zeros, nz);
    if (n == 0) return 0;
    /* cheby2 BP: DC and Nyquist gains 0 (odd) / 10^(−rs/20) (even) */
    if (!check_cascade_gains(sections, n,
                             cheby2_edge_gain(order, ripple_db),
                             cheby2_edge_gain(order, ripple_db))) return 0;
    return n;
}

/**
 * @brief Chebyshev II 带阻设计：原型 → 管线 → 增益校验。
 *
 * 失败返回 0（调用方部署直通）。
 * 期望级联增益：DC 与 Nyquist 均为 1。
 *
 * @param[out] sections      biquad 节数组。
 * @param[in]  max_sections  sections 容量。
 * @param[in]  order         原型阶数。
 * @param[in]  fc1           下带边（Hz）。
 * @param[in]  fc2           上带边（Hz）。
 * @param[in]  fs            采样频率（Hz）。
 * @param[in]  ripple_db     阻带衰减（dB）。
 * @return                   部署的节数，0 表示失败。
 */
static uint8_t cheby2_bs_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order,
                                      float fc1, float fc2, float fs,
                                      float ripple_db)
{
    if (order == 0 || order > 8 || fc1 <= 0.0f || fc1 >= fs * 0.5f || ripple_db <= 0.0f)
        return 0;
    if (fc2 <= fc1 || fc2 >= fs * 0.5f) return 0;

    float epsilon = 1.0f / sqrtf(powf(10.0f, ripple_db / 10.0f) - 1.0f);
    float wc1 = 2.0f * (float)M_PI * prewarp(fc1, fs);
    float wc2 = 2.0f * (float)M_PI * prewarp(fc2, fs);

    complex_t poles[8], zeros[8];
    uint8_t nz = cheby2_proto(poles, zeros, order, epsilon);

    float k = 1.0f / zpk_hp_bs_gain(1.0f, zeros, nz, poles, order);

    uint8_t n = design_filter(sections, max_sections,
                              FILTER_BANDSTOP,
                              wc1, wc2, fs, k,
                              poles, order, zeros, nz);
    if (n == 0) return 0;
    /* cheby2 BS: DC and Nyquist gains 1 */
    if (!check_cascade_gains(sections, n, 1.0f, 1.0f)) return 0;
    return n;
}

/* ================================================================== */
/*  各阶 init 函数（宏生成）                                           */
/* ================================================================== */

/*
 * 宏生成的各阶 init 共享上述辅助函数的语义（详见头文件声明）：
 * 各自填写元数据字段，经对应类型的辅助函数校验，并按返回值
 * 设置 num_sections 与 valid。
 */

/* Chebyshev I — 低通 */
#define X(ord, ns, ol) \
    void cheby1_lp_##ol##_init(cheby1_lp_##ol##_t *f, float fc, float fs, float ripple_db) { \
        f->type = FILTER_LOWPASS; \
        f->order = ord; \
        f->fc1 = fc; \
        f->fc2 = 0.0f; \
        f->fs = fs; \
        f->ripple_db = ripple_db; \
        f->valid = 0; \
        uint8_t n = cheby1_lp_init(f->sections, ns, ord, fc, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_CHEBY_LP_ORDER
#undef X

/* Chebyshev I — 高通 */
#define X(ord, ns, ol) \
    void cheby1_hp_##ol##_init(cheby1_hp_##ol##_t *f, float fc, float fs, float ripple_db) { \
        f->type = FILTER_HIGHPASS; \
        f->order = ord; \
        f->fc1 = fc; \
        f->fc2 = 0.0f; \
        f->fs = fs; \
        f->ripple_db = ripple_db; \
        f->valid = 0; \
        uint8_t n = cheby1_hp_init(f->sections, ns, ord, fc, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_CHEBY_LP_ORDER
#undef X

/* Chebyshev I — 带通 */
#define X(ord, ns, ol) \
    void cheby1_bp_##ol##_init(cheby1_bp_##ol##_t *f, float fc1, float fc2, float fs, float ripple_db) { \
        f->type = FILTER_BANDPASS; \
        f->order = ord; \
        f->fc1 = fc1; \
        f->fc2 = fc2; \
        f->fs = fs; \
        f->ripple_db = ripple_db; \
        f->valid = 0; \
        uint8_t n = cheby1_bp_init(f->sections, ns, ord, fc1, fc2, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_CHEBY_BP_ORDER
#undef X

/* Chebyshev I — 带阻 */
#define X(ord, ns, ol) \
    void cheby1_bs_##ol##_init(cheby1_bs_##ol##_t *f, float fc1, float fc2, float fs, float ripple_db) { \
        f->type = FILTER_BANDSTOP; \
        f->order = ord; \
        f->fc1 = fc1; \
        f->fc2 = fc2; \
        f->fs = fs; \
        f->ripple_db = ripple_db; \
        f->valid = 0; \
        uint8_t n = cheby1_bs_init(f->sections, ns, ord, fc1, fc2, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_CHEBY_BP_ORDER
#undef X

/* Chebyshev II — 低通 */
#define X(ord, ns, ol) \
    void cheby2_lp_##ol##_init(cheby2_lp_##ol##_t *f, float fc, float fs, float ripple_db) { \
        f->type = FILTER_LOWPASS; \
        f->order = ord; \
        f->fc1 = fc; \
        f->fc2 = 0.0f; \
        f->fs = fs; \
        f->ripple_db = ripple_db; \
        f->valid = 0; \
        uint8_t n = cheby2_lp_init(f->sections, ns, ord, fc, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_CHEBY_LP_ORDER
#undef X

/* Chebyshev II — 高通 */
#define X(ord, ns, ol) \
    void cheby2_hp_##ol##_init(cheby2_hp_##ol##_t *f, float fc, float fs, float ripple_db) { \
        f->type = FILTER_HIGHPASS; \
        f->order = ord; \
        f->fc1 = fc; \
        f->fc2 = 0.0f; \
        f->fs = fs; \
        f->ripple_db = ripple_db; \
        f->valid = 0; \
        uint8_t n = cheby2_hp_init(f->sections, ns, ord, fc, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_CHEBY_LP_ORDER
#undef X

/* Chebyshev II — 带通 */
#define X(ord, ns, ol) \
    void cheby2_bp_##ol##_init(cheby2_bp_##ol##_t *f, float fc1, float fc2, float fs, float ripple_db) { \
        f->type = FILTER_BANDPASS; \
        f->order = ord; \
        f->fc1 = fc1; \
        f->fc2 = fc2; \
        f->fs = fs; \
        f->ripple_db = ripple_db; \
        f->valid = 0; \
        uint8_t n = cheby2_bp_init(f->sections, ns, ord, fc1, fc2, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_CHEBY_BP_ORDER
#undef X

/* Chebyshev II — 带阻 */
#define X(ord, ns, ol) \
    void cheby2_bs_##ol##_init(cheby2_bs_##ol##_t *f, float fc1, float fc2, float fs, float ripple_db) { \
        f->type = FILTER_BANDSTOP; \
        f->order = ord; \
        f->fc1 = fc1; \
        f->fc2 = fc2; \
        f->fs = fs; \
        f->ripple_db = ripple_db; \
        f->valid = 0; \
        uint8_t n = cheby2_bs_init(f->sections, ns, ord, fc1, fc2, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_CHEBY_BP_ORDER
#undef X
