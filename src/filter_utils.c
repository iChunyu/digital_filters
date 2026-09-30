#include "filter_utils.h"
#include <math.h>
#include <string.h>

/* ------------------------------------------------------------------ */

float prewarp(float fd, float fs)
{
    /*
     * Fail-closed：区间 (0, fs/2) 之外双线性映射无定义（过了 Nyquist，tan()
     * 变号，根本不存在调用方能用的模拟频率）。写成取反比较，NaN 的 fd 或 fs
     * 也会 fail-closed 而不是漏到 tanf()；返回的 NaN 由下游 design_filter 的
     * 有限性闸兜住，错误截止频率得到一个直通，而不是悄悄错掉的滤波器。
     */
    if (!(fd > 0.0f) || !(fd < fs * 0.5f)) {
        return NAN;
    }
    return fs / ((float)M_PI) * tanf((float)M_PI * fd / fs);
}

/* ------------------------------------------------------------------ */

/**
 * @brief 复数数组整体乘以标量。
 *
 * @param[in,out] zp  复数数组。
 * @param[in]     n   元素个数。
 * @param[in]     s   缩放系数。
 */
static void scale_complex(complex_t *zp, uint8_t n, float s)
{
    for (uint8_t i = 0; i < n; i++) {
        zp[i].re *= s;
        zp[i].im *= s;
    }
}

void analog_lp_transform(complex_t *poles, uint8_t np,
                         complex_t *zeros, uint8_t nz, float wc)
{
    (void)np;
    (void)nz;
    scale_complex(poles, np, wc);
    if (nz > 0) {
        scale_complex(zeros, nz, wc);
    }
}

/* ------------------------------------------------------------------ */

void analog_hp_transform(complex_t *poles, uint8_t np,
                         complex_t *zeros, uint8_t nz, float wc)
{
    for (uint8_t i = 0; i < np; i++) {
        float den = poles[i].re * poles[i].re + poles[i].im * poles[i].im;
        /* p = wc / p_proto = wc * conj(p_proto) / |p_proto|^2 */
        /* 1 / (a + jb) = (a - jb) / (a^2 + b^2) */
        poles[i].re =  wc * poles[i].re / den;
        poles[i].im = -wc * poles[i].im / den;
    }
    for (uint8_t i = 0; i < nz; i++) {
        float den = zeros[i].re * zeros[i].re + zeros[i].im * zeros[i].im;
        zeros[i].re =  wc * zeros[i].re / den;
        zeros[i].im = -wc * zeros[i].im / den;
    }
}

/* ------------------------------------------------------------------ */

/**
 * @brief 复数平方根：计算 sqrt(re + j·im)。
 *
 * 返回主支（实部非负）。
 *
 * @note 幅值采用缩放计算：近 Nyquist 的 BP/BS 判别式 |re|,|im| ~ 1.5e10，
 *       裸平方在任一分量超过 √FLT_MAX ≈ 1.8e19 时溢出，而管线要撑到 prewarp
 *       的 tanf() 饱和处。详细推导与档位见 test/test_butter.c 的近 Nyquist sweep 段。
 *
 * @param re  实部。
 * @param im  虚部。
 * @return    平方根（主支）。
 */
static complex_t c_sqrt(float re, float im)
{
    complex_t r;
    /* 按最大分量缩放，避免近 Nyquist 判别式的平方溢出（见上方 @note）。 */
    float m = fmaxf(fabsf(re), fabsf(im));
    if (m < 1e-20f) {
        r.re = 0.0f;
        r.im = 0.0f;
        return r;
    }
    /* 处理负实部情形，避免 (mag + re)/2 中的灾难性消减。 */
    if (re < 0.0f && fabsf(im) < 1e-12f * fabsf(re)) {
        r.re = 0.0f;
        r.im = sqrtf(-re);
        if (im < 0.0f) r.im = -r.im;
        return r;
    }
    float sr = re / m;
    float si = im / m;
    float mag = m * sqrtf(sr * sr + si * si);
    float s1 = sqrtf(0.5f * (mag + re));
    float s2 = sqrtf(0.5f * (mag - re));
    r.re = s1;
    r.im = (im >= 0.0f) ? s2 : -s2;
    return r;
}

/* ------------------------------------------------------------------ */

static void c_mul(float *rr, float *ri, float ar, float ai, float br, float bi);
static void c_div(float *rr, float *ri, float ar, float ai, float br, float bi);

/**
 * @brief 求 s² + B·s + C = 0（B、C 为复数）的稳定根。
 *
 * root1 = (−B − √(B²−4C))/2，root2 = C / root1。直接公式 (−B + √…)/2 在
 * |B|² ≫ 4|C| 时灾难性消减（宽带 BP/BS 的实原型极点）；小根改用 C / root1
 * 绕开。标定与回归档位见 test/test_butter.c 的「回归：含近实极点对的宽带 BP」段。
 *
 * @param[in]  B_re  一次项系数 B 的实部。
 * @param[in]  B_im  一次项系数 B 的虚部。
 * @param[in]  C_re  常数项系数 C 的实部。
 * @param[in]  C_im  常数项系数 C 的虚部。
 * @param[out] r1    大根（−B − √disc）/2。
 * @param[out] r2    小根 C / r1（r1 退化时退化为直接公式）。
 */
static void stable_roots(float B_re, float B_im, float C_re, float C_im,
                         complex_t *r1, complex_t *r2)
{
    float d_re, d_im;
    c_mul(&d_re, &d_im, B_re, B_im, B_re, B_im);
    d_re -= 4.0f * C_re;
    d_im -= 4.0f * C_im;
    complex_t sd = c_sqrt(d_re, d_im);

    r1->re = 0.5f * (-B_re - sd.re);
    r1->im = 0.5f * (-B_im - sd.im);

    if (r1->re * r1->re + r1->im * r1->im < 1e-20f) {
        /* 大根退化——退回直接公式。 */
        r2->re = 0.5f * (-B_re + sd.re);
        r2->im = 0.5f * (-B_im + sd.im);
    } else {
        c_div(&r2->re, &r2->im, C_re, C_im, r1->re, r1->im);
    }
}

void analog_bp_transform(complex_t *poles, uint8_t *np,
                         complex_t *zeros, uint8_t *nz,
                         float w0, float xi)
{
    uint8_t np_old = *np;
    uint8_t nz_old = *nz;
    float w0_sq = w0 * w0;

    /* 极点倒序遍历变换，避免覆盖尚未读取的元素。
       s → (s² + w0²) / (xi·s)  ⇒  s² − xi·p·s + w0² = 0。 */
    for (int i = np_old - 1; i >= 0; i--) {
        stable_roots(-xi * poles[i].re, -xi * poles[i].im,
                     w0_sq, 0.0f,
                     &poles[2 * i], &poles[2 * i + 1]);
    }
    *np = 2 * np_old;

    /* 有限零点用同一公式变换。 */
    for (int i = nz_old - 1; i >= 0; i--) {
        stable_roots(-xi * zeros[i].re, -xi * zeros[i].im,
                     w0_sq, 0.0f,
                     &zeros[2 * i], &zeros[2 * i + 1]);
    }

    /* 为原型中的无穷远零点补 (np_old − nz_old) 个原点零点。
       双线性变换后它们映射到 z = +1。 */
    for (uint8_t i = 0; i < np_old - nz_old; i++) {
        zeros[2 * nz_old + i].re = 0.0f;
        zeros[2 * nz_old + i].im = 0.0f;
    }
    *nz = nz_old + np_old;  /* = 2·nz_old + (np_old − nz_old) */
}

/* ------------------------------------------------------------------ */

void analog_bs_transform(complex_t *poles, uint8_t *np,
                         complex_t *zeros, uint8_t *nz,
                         float w0, float xi)
{
    uint8_t np_old = *np;
    uint8_t nz_old = *nz;
    float w0_sq = w0 * w0;

    /* 极点倒序遍历变换。p → ξ·s / (s² + ω₀²)
       ⇒  s² − (ξ/p)·s + ω₀² = 0。 */
    for (int i = np_old - 1; i >= 0; i--) {
        float pr = poles[i].re;
        float pi = poles[i].im;
        float mag2 = pr * pr + pi * pi;
        if (mag2 < 1e-20f) {
            /* 退化情形：原样保留（稳定原型不应出现）。 */
            poles[2 * i].re     = poles[i].re;
            poles[2 * i].im     = poles[i].im;
            poles[2 * i + 1].re = poles[i].re;
            poles[2 * i + 1].im = poles[i].im;
            continue;
        }

        float b_re, b_im;
        c_div(&b_re, &b_im, -xi, 0.0f, pr, pi); /* B = −ξ / p */
        stable_roots(b_re, b_im, w0_sq, 0.0f,
                     &poles[2 * i], &poles[2 * i + 1]);
    }
    *np = 2 * np_old;

    /* 有限零点用同一公式变换。 */
    for (int i = nz_old - 1; i >= 0; i--) {
        float zr = zeros[i].re;
        float zi = zeros[i].im;
        float mag2 = zr * zr + zi * zi;
        if (mag2 < 1e-20f) {
            zeros[2 * i].re     = zeros[i].re;
            zeros[2 * i].im     = zeros[i].im;
            zeros[2 * i + 1].re = zeros[i].re;
            zeros[2 * i + 1].im = zeros[i].im;
            continue;
        }

        float b_re, b_im;
        c_div(&b_re, &b_im, -xi, 0.0f, zr, zi); /* B = −ξ / z */
        stable_roots(b_re, b_im, w0_sq, 0.0f,
                     &zeros[2 * i], &zeros[2 * i + 1]);
    }

    /* 为原型中的无穷远零点补 2·(np_old − nz_old) 个 ±jω₀ 零点。 */
    for (uint8_t i = 0; i < np_old - nz_old; i++) {
        zeros[2 * nz_old + 2 * i].re     =  0.0f;
        zeros[2 * nz_old + 2 * i].im     =  w0;
        zeros[2 * nz_old + 2 * i + 1].re =  0.0f;
        zeros[2 * nz_old + 2 * i + 1].im = -w0;
    }
    *nz = 2 * np_old;
}

/* ------------------------------------------------------------------ */

void bilinear_transform(complex_t *zp, uint8_t n, float fs)
{
    float K = 2.0f * fs;
    float K2 = K * K;

    for (uint8_t i = 0; i < n; i++) {
        float s_re = zp[i].re;
        float s_im = zp[i].im;

        /* z = (K + s) / (K - s) = (K + s)(K - conj(s)) / |K − s|² */
        float den = K2 - 2.0f * K * s_re + s_re * s_re + s_im * s_im;

        zp[i].re = (K2 - s_re * s_re - s_im * s_im) / den;
        zp[i].im = (2.0f * K * s_im) / den;
    }
}

/* ------------------------------------------------------------------ */

/* 两点间欧氏距离 |a − b|。 */
static inline float c_dist(const complex_t *a, const complex_t *b)
{
    float dr = a->re - b->re;
    float di = a->im - b->im;
    return sqrtf(dr * dr + di * di);
}

/* 相对容差实数判定：|im| ≤ eps·|z|（与 scipy _cplxreal 一致）。 */
static inline int is_real(const complex_t *a, float eps)
{
    return fabsf(a->im) <= eps * sqrtf(a->re * a->re + a->im * a->im);
}

/**
 * @brief 由一对极点与一对零点计算 biquad 系数。
 *
 * 极点 p1,p2 与零点 z1,z2 各自要么同为实数、要么互为共轭。
 *
 * @param[in]  p1  第一个极点。
 * @param[in]  p2  第二个极点。
 * @param[in]  z1  第一个零点。
 * @param[in]  z2  第二个零点。
 * @param[out] b   分子系数 [b1, b2]（b0 恒为 1）。
 * @param[out] a   分母系数 [a1, a2]（a0 恒为 1）。
 */
static void make_biquad(const complex_t *p1, const complex_t *p2,
                        const complex_t *z1, const complex_t *z2,
                        float *b, float *a)
{
    /* a1 = −(p1 + p2)，a2 = p1 · p2 */
    a[0] = -(p1->re + p2->re);          /* a1 */
    a[1] = p1->re * p2->re - p1->im * p2->im;  /* a2 = Re(p1·p2)      */
    /* （共轭时即 |p1|²；实数时 im = 0，同样成立）                     */

    /* b0 = 1，b1 = −(z1 + z2)，b2 = z1 · z2 */
    b[0] = 1.0f;
    b[1] = -(z1->re + z2->re);
    b[2] = z1->re * z2->re - z1->im * z2->im;
}

/* 复数乘法 r = a · b，实虚部分开传出。 */
static void c_mul(float *rr, float *ri, float ar, float ai, float br, float bi)
{
    *rr = ar * br - ai * bi;
    *ri = ar * bi + ai * br;
}

/* 复数除法 r = a / b；b 为零时结果置 0（防御路径，正常流程不会发生）。 */
static void c_div(float *rr, float *ri, float ar, float ai, float br, float bi)
{
    float den = br * br + bi * bi;
    if (den == 0.0f) { *rr = 0.0f; *ri = 0.0f; return; }
    *rr = (ar * br + ai * bi) / den;
    *ri = (ai * br - ar * bi) / den;
}

/* ================================================================== */
/*  ZPK 增益辅助（对齐 scipy 的 zpk 增益追踪）                         */
/* ================================================================== */

float zpk_hp_bs_gain(float k, const complex_t *z, uint8_t nz,
                     const complex_t *p, uint8_t np)
{
    /* 增量计算 k · ∏(−z_i) / ∏(−p_i)，避免中间量溢出。 */
    float rr = k, ri = 0.0f;
    uint8_t i;
    for (i = 0; i < nz && i < np; i++) {
        /* 乘 (−z[i])，除 (−p[i])。 */
        c_mul(&rr, &ri, rr, ri, -z[i].re, -z[i].im);
        c_div(&rr, &ri, rr, ri, -p[i].re, -p[i].im);
    }
    for (; i < nz; i++)
        c_mul(&rr, &ri, rr, ri, -z[i].re, -z[i].im);
    for (; i < np; i++)
        c_div(&rr, &ri, rr, ri, -p[i].re, -p[i].im);
    return rr; /* 共轭对下虚部相消 */
}

float bilinear_zpk_gain(float k, const complex_t *z, uint8_t nz,
                         const complex_t *p, uint8_t np, float K)
{
    /* 增量计算 k · ∏(K−z_i) / ∏(K−p_i)，避免中间量溢出。 */
    float rr = k, ri = 0.0f;
    uint8_t i;
    for (i = 0; i < nz && i < np; i++) {
        c_mul(&rr, &ri, rr, ri, K - z[i].re, -z[i].im);
        c_div(&rr, &ri, rr, ri, K - p[i].re, -p[i].im);
    }
    for (; i < nz; i++)
        c_mul(&rr, &ri, rr, ri, K - z[i].re, -z[i].im);
    for (; i < np; i++)
        c_div(&rr, &ri, rr, ri, K - p[i].re, -p[i].im);
    return rr;
}

float bilinear_zpk_gain_scaled(float k, float s, uint8_t degree,
                               const complex_t *z, uint8_t nz,
                               const complex_t *p, uint8_t np, float K)
{
    /* k · s^degree · ∏(K−z_i) / ∏(K−p_i)，按交错顺序算让运行乘积始终有界：
       单独先算 s^degree 会在近 Nyquist 的 LP/BP 上溢出 f32（wc^8 ≈ FLT_MAX），
       而中间量一旦溢出，后面同量级的 ∏(K−p) 因子再也抵消不回来。 */
    float rr = k, ri = 0.0f;
    uint8_t i_s = 0, i_z = 0, i_p = 0;

    while (i_s < degree || i_z < nz || i_p < np) {
        if (i_s < degree && i_p < np) {
            c_mul(&rr, &ri, rr, ri, s, 0.0f);
            c_div(&rr, &ri, rr, ri, K - p[i_p].re, -p[i_p].im);
            i_s++;
            i_p++;
        } else if (i_z < nz && i_p < np) {
            c_mul(&rr, &ri, rr, ri, K - z[i_z].re, -z[i_z].im);
            c_div(&rr, &ri, rr, ri, K - p[i_p].re, -p[i_p].im);
            i_z++;
            i_p++;
        } else if (i_s < degree) {
            c_mul(&rr, &ri, rr, ri, s, 0.0f);
            i_s++;
        } else if (i_z < nz) {
            c_mul(&rr, &ri, rr, ri, K - z[i_z].re, -z[i_z].im);
            i_z++;
        } else {
            c_div(&rr, &ri, rr, ri, K - p[i_p].re, -p[i_p].im);
            i_p++;
        }
    }
    return rr; /* 共轭对下虚部相消 */
}

/**
 * @brief 找离单位圆最近的未用极点（幅值最大者）。
 *
 * 调用前提：至少存在一个未用元素。
 *
 * @param[in] poles  极点数组。
 * @param[in] used   已用位图。
 * @param[in] n      数组长度。
 * @return           最不利极点的下标。
 */
static uint8_t find_worst_pole(const complex_t *poles, const uint8_t *used, uint8_t n)
{
    uint8_t idx = 0;
    float best = -1.0f;
    for (uint8_t i = 0; i < n; i++) {
        if (used[i]) continue;
        float mag = poles[i].re * poles[i].re + poles[i].im * poles[i].im;
        if (mag > best) {
            best = mag;
            idx = i;
        }
    }
    return idx;
}

/**
 * @brief 找未用元素中离 target 最近者的下标。
 *
 * 调用前提：至少存在一个未用元素。
 *
 * @param[in] arr     数组。
 * @param[in] used    已用位图。
 * @param[in] n       数组长度。
 * @param[in] target  目标点。
 * @return            最近元素的下标。
 */
static uint8_t find_nearest(const complex_t *arr, const uint8_t *used, uint8_t n,
                            const complex_t *target)
{
    uint8_t idx = 0;
    float best = 1e30f;
    for (uint8_t i = 0; i < n; i++) {
        if (used[i]) continue;
        float d = c_dist(&arr[i], target);
        if (d < best) {
            best = d;
            idx = i;
        }
    }
    return idx;
}

/**
 * @brief 同 find_nearest，但限定实数（want_real=1）或复数元素。
 *
 * @param[in]  arr       数组。
 * @param[in]  used      已用位图。
 * @param[in]  n         数组长度。
 * @param[in]  target    目标点。
 * @param[in]  want_real 1 只考虑实数元素，0 只考虑复数元素。
 * @param[in]  eps       实/复分类容差。
 * @param[out] best      最近距离；无匹配元素时为 1e30f。
 * @return               下标（无匹配时以 *best == 1e30f 判定，
 *                       返回值无意义）。
 */
static uint8_t find_nearest_typed(const complex_t *arr, const uint8_t *used,
                                  uint8_t n, const complex_t *target,
                                  int want_real, float eps, float *best)
{
    uint8_t idx = 0;
    *best = 1e30f;
    for (uint8_t i = 0; i < n; i++) {
        if (used[i]) continue;
        if (is_real(&arr[i], eps) != want_real) continue;
        float d = c_dist(&arr[i], target);
        if (d < *best) {
            *best = d;
            idx = i;
        }
    }
    return idx;
}

/**
 * @brief 统计未用元素中具有指定实/复性质的数量。
 *
 * @param[in] used      已用位图。
 * @param[in] n         数组长度。
 * @param[in] arr       数组（用于分类判定）。
 * @param[in] want_real 1 统计实数元素，0 统计复数元素。
 * @param[in] eps       实/复分类容差。
 * @return               符合条件的未用元素个数。
 */
static uint8_t count_used(const uint8_t *used, uint8_t n,
                          const complex_t *arr, int want_real, float eps)
{
    uint8_t cnt = 0;
    for (uint8_t i = 0; i < n; i++) {
        if (!used[i] && is_real(&arr[i], eps) == want_real) {
            cnt++;
        }
    }
    return cnt;
}

/**
 * @brief 在未用元素中认领 arr[idx] 的共轭并标记已用。
 *
 * 匹配为相对容差：|candidate − conj(target)| ≤ eps·|target|。
 * 失败时合成共轭（不消耗任何元素——调用方的全认领不变量
 * 随后 fail-closed）。
 *
 * @note 盒内取**最近**未用元素，不是第一个：否则高 Q BP/BS 簇里会偷走别对
 *       的伴侣，部署出两节完全重复、丢失一对。回归档位见 test/test_cheby.c。
 *
 * @param[in]     arr  数组。
 * @param[in,out] used 已用位图（命中时置位）。
 * @param[in]     n    数组长度。
 * @param[in]     idx  目标元素下标（调用方须已置 used[idx]）。
 * @param[out]    out  共轭元素（命中为对称平均对，未命中为合成值）。
 * @param[in]     eps  相对匹配容差。
 */
static void claim_conjugate(const complex_t *arr, uint8_t *used, uint8_t n,
                            uint8_t idx, complex_t *out, float eps)
{
    float tr = arr[idx].re;
    float ti = arr[idx].im;
    float mag = sqrtf(tr * tr + ti * ti);
    float thresh = eps * (mag > 1e-12f ? mag : 1.0f);
    uint8_t best_i = n;
    float best_d = 1e30f;
    for (uint8_t i = 0; i < n; i++) {
        if (used[i]) continue;
        if (fabsf(arr[i].re - tr) > thresh) continue;
        if (fabsf(arr[i].im + ti) > thresh) continue;
        float dr = arr[i].re - tr;
        float di = arr[i].im + ti;
        float d = sqrtf(dr * dr + di * di);
        if (d < best_d) {
            best_d = d;
            best_i = i;
        }
    }
    if (best_i != n) {
        used[best_i] = 1;
        /* 对两点取平均（对齐 scipy 的 _cplxreal）。
           arr[best_i] 的虚部符号与 arr[idx] 相反。
           返回 arr[best_i] 与 conj(arr[idx]) 的平均，使调用方拿到
           主元素的真共轭。 */
        out->re = 0.5f * (tr + arr[best_i].re);
        out->im = 0.5f * (arr[best_i].im - ti);
        return;
    }
    /* 最后手段：数值不匹配——直接合成共轭。 */
    out->re =  tr;
    out->im = -ti;
}

/**
 * @brief 求 z² + c1·z + c2 的根。
 *
 * c2 == 0 时该节实际是一阶（实极点/零点对），只返回实根 −c1，忽略补位产生的
 * z = 0 根——输入零极点多重集里没有它（见下方的 match_roots）。
 *
 * @param[in]  c1    一次项系数。
 * @param[in]  c2    常数项系数。
 * @param[out] roots 最多 2 个根。
 * @return           根个数（1 或 2）。
 */
static uint8_t poly_roots(float c1, float c2, complex_t roots[2])
{
    if (c2 == 0.0f) {
        roots[0].re = -c1;
        roots[0].im = 0.0f;
        return 1;
    }
    float disc = c1 * c1 - 4.0f * c2;
    if (disc >= 0.0f) {
        float s = sqrtf(disc);
        roots[0].re = 0.5f * (-c1 + s);
        roots[0].im = 0.0f;
        roots[1].re = 0.5f * (-c1 - s);
        roots[1].im = 0.0f;
    } else {
        float im = 0.5f * sqrtf(-disc);
        roots[0].re = -0.5f * c1;
        roots[0].im =  im;
        roots[1].re = -0.5f * c1;
        roots[1].im = -im;
    }
    return 2;
}

/**
 * @brief 多重集匹配：roots 的每个根须命中 pool 中一个不同的
 *        未匹配元素（欧氏容差）。
 *
 * 紧对成员允许互换命中——无害；拒绝的是复现不出任何输入元素的根。
 *
 * @param[in]     roots    待匹配根数组。
 * @param[in]     nr       根个数。
 * @param[in]     pool     输入池。
 * @param[in]     n        池大小。
 * @param[in,out] matched  匹配位图（命中时置位）。
 * @param[in]     tol      欧氏容差。
 * @return                 全部命中返回 1。
 */
static int match_roots(const complex_t *roots, uint8_t nr,
                       const complex_t *pool, uint8_t n,
                       uint8_t *matched, float tol)
{
    for (uint8_t r = 0; r < nr; r++) {
        uint8_t bi = n;
        float best = tol;
        for (uint8_t i = 0; i < n; i++) {
            if (matched[i]) continue;
            float d = c_dist(&roots[r], &pool[i]);
            if (d < best) {
                best = d;
                bi = i;
            }
        }
        if (bi == n) return 0;
        matched[bi] = 1;
    }
    return 1;
}

/* zpk2sos 支持的最大零极点对数。
 * 8 阶原型 → BP/BS 翻倍到 16。 */
#define ZPK2SOS_MAX_N 16

uint8_t zpk2sos(const complex_t *zeros, const complex_t *poles, uint8_t n,
                float (*sos)[6], float k)
{
    if (n == 0 || n > ZPK2SOS_MAX_N) return 0;

    uint8_t used_p[ZPK2SOS_MAX_N];
    uint8_t used_z[ZPK2SOS_MAX_N];
    memset(used_p, 0, n * sizeof(uint8_t));
    memset(used_z, 0, n * sizeof(uint8_t));

    uint8_t n_p = n;
    uint8_t n_z = n;
    uint8_t section = 0;
    uint8_t max_sections = (n + 1) / 2;
    /* eps_class（实/复分类）与 eps_claim（共轭认领盒）是两个不同的量，不可合并：
       - 真·实零极点的虚部严格为 0.0f，而真·共轭对 |im| ≥ ~1e-5。用 1e-3 做分类
         容差会把宽带 BP/BS 的合法近实共轭对压平成"实"，与另一个对的成员跨对
         配对，造出极点恰好在 z = 1 的节，整个设计随之失败。
       - 认领盒必须吞下双线性变换在单位圆附近的 f32 舍入（可把共轭对劈开 ~2e-4）。 */
    const float eps_class = 1e-5f;
    const float eps_claim = 1e-3f;

    while (n_p > 0) {
        /* 安全上限：不超出已分配的行数，
           也不把计数器减过零。 */
        if (section >= max_sections) break;
        if (n_z == 0) break;

        /* 1. 挑出剩余极点中最不利的（|p| 最大）。 */
        uint8_t p1_i = find_worst_pole(poles, used_p, n);
        complex_t p1 = poles[p1_i];
        used_p[p1_i] = 1;
        n_p--;

        complex_t p2, z1, z2;

        if (is_real(&p1, eps_class)) {
            /* p1 是实数——尝试与另一个实极点配对。 */
            if (count_used(used_p, n, poles, 1, eps_class) > 0) {
                float best;
                uint8_t p2_i = find_nearest_typed(poles, used_p, n, &p1,
                                                  1, eps_class, &best);
                p2 = poles[p2_i];
                used_p[p2_i] = 1;
                n_p--;
            } else {
                /* 没有更多实极点 → 用实零点构成一阶节。
                   若实零点也不剩，说明零极点集不平衡；此时 fail-closed，
                   而不是把一个复零点压成实数、把它的共轭丢出级联。 */
                float best;
                uint8_t z1_i = find_nearest_typed(zeros, used_z, n, &p1,
                                                  1, eps_class, &best);
                if (best == 1e30f) return 0;
                z1 = zeros[z1_i];
                used_z[z1_i] = 1;
                n_z--;

                float b[3] = {1.0f, -z1.re, 0.0f};
                float a[2] = {-p1.re, 0.0f};

                sos[section][0] = b[0];
                sos[section][1] = b[1];
                sos[section][2] = b[2];
                sos[section][3] = 1.0f;
                sos[section][4] = a[0];
                sos[section][5] = a[1];
                section++;
                continue;
            }
        } else {
            /* p1 是复数——与其共轭配对。 */
            claim_conjugate(poles, used_p, n, p1_i, &p2, eps_claim);
            n_p--;
        }

        /* 2. 为这对极点匹配两个零点。 */
        uint8_t z1_i = find_nearest(zeros, used_z, n, &p1);

        if (is_real(&zeros[z1_i], eps_class)) {
            if (count_used(used_z, n, zeros, 1, eps_class) > 1) {
                /* 有两个实零点可用——取 z1 与最近的另一个实零点。 */
                z1 = zeros[z1_i];
                used_z[z1_i] = 1;
                n_z--;

                float best;
                uint8_t z2_i = find_nearest_typed(zeros, used_z, n, &p1,
                                                  1, eps_class, &best);
                z2 = zeros[z2_i];
                used_z[z2_i] = 1;
                n_z--;
            } else {
                /* 只剩一个实零点——留给后面的一阶节，
                   这里改用一对复零点。 */
                uint8_t best_i = z1_i;
                float best_d = 1e30f;
                for (uint8_t i = 0; i < n; i++) {
                    if (used_z[i]) continue;
                    if (is_real(&zeros[i], eps_class)) continue;
                    float d = c_dist(&zeros[i], &p1);
                    if (d < best_d) { best_d = d; best_i = i; }
                }
                z1 = zeros[best_i];
                used_z[best_i] = 1;
                n_z--;

                claim_conjugate(zeros, used_z, n, best_i, &z2, eps_claim);
                n_z--;
            }
        } else {
            /* z1 是复数——取它与其共轭。 */
            z1 = zeros[z1_i];
            used_z[z1_i] = 1;
            n_z--;

            claim_conjugate(zeros, used_z, n, z1_i, &z2, eps_claim);
            n_z--;
        }

        /* 3. 计算 biquad 系数。 */
        float b[3], a[2];
        make_biquad(&p1, &p2, &z1, &z2, b, a);

        sos[section][0] = b[0];
        sos[section][1] = b[1];
        sos[section][2] = b[2];
        sos[section][3] = 1.0f;
        sos[section][4] = a[0];
        sos[section][5] = a[1];
        section++;
    }

    /* 3b. 不变量：全部零极点都必须被认领。留下未认领的空隙说明 n_p/n_z 记账与
       数组脱节、有零极点被悄悄丢出级联（回归见 test/test_cheby.c）。 */
    for (uint8_t i = 0; i < n; i++) {
        if (!used_p[i] || !used_z[i]) return 0;
    }

    /* 3c. 校验每节的根能否复现输入的零极点多重集。防的是跨对误认领——
       used[] 全置位了，但配对配错了元素。 */
    {
        uint8_t matched[ZPK2SOS_MAX_N];
        float tol = 1e-3f;

        memset(matched, 0, sizeof(matched));
        for (uint8_t s = 0; s < section; s++) {
            complex_t roots[2];
            uint8_t nr = poly_roots(sos[s][4], sos[s][5], roots);
            if (!match_roots(roots, nr, poles, n, matched, tol)) return 0;
        }
        memset(matched, 0, sizeof(matched));
        for (uint8_t s = 0; s < section; s++) {
            complex_t roots[2];
            uint8_t nr = poly_roots(sos[s][1], sos[s][2], roots); /* b0 == 1 */
            if (!match_roots(roots, nr, zeros, n, matched, tol)) return 0;
        }
    }

    /* 4. 反转节顺序：最慢的极点排最后 → 最快的排最前。 */
    for (uint8_t i = 0; i < section / 2; i++) {
        for (uint8_t j = 0; j < 6; j++) {
            float tmp = sos[i][j];
            sos[i][j] = sos[section - 1 - i][j];
            sos[section - 1 - i][j] = tmp;
        }
    }

    /* 5. 把系统总增益施加到第一节的分子。 */
    sos[0][0] *= k;
    sos[0][1] *= k;
    sos[0][2] *= k;

    return section;
}

/* ================================================================== */
/*  共享设计管线（Butterworth / Chebyshev）                            */
/* ================================================================== */

/* BP/BS 变换后的最大原型阶数：2 × 8 = 16 个零极点，ceil(16/2) = 8 节。 */

uint8_t design_filter(biquad_filter_t *sections, uint8_t max_sections,
                      uint8_t type,
                      float wc1, float wc2, float fs,
                      float k,
                      const complex_t *proto_poles, uint8_t np,
                      const complex_t *proto_zeros, uint8_t nz)
{
    /*
     * Fail-closed 输入边界闸：poles/zeros 只容纳 ZPK2SOS_MAX_N 个元素、degree 是
     * uint8_t，np 超出包络会越界写栈，nz > np 会让相对阶数下溢。
     *
     * BP/BS 另需再紧一半（上限 ZPK2SOS_MAX_N/2 = 8 阶）：这两种类型的变换把每个
     * 原型零极点一分为二原地写回（容量要求 2·np，见 analog_bp_transform 文档）。
     * 族 API 恰好卡在 8 阶，但本函数是导出符号，边界得自己守——越界发生在随后的
     * BP/BS 分支，而唯一能拦住它的 ns > max_sections 在越界之后。档位见
     * test/test_butter.c 的「回归：design_filter 输入边界闸」段。
     */
    if (np == 0 || np > ZPK2SOS_MAX_N || nz > np) return 0;
    if ((type == FILTER_BANDPASS || type == FILTER_BANDSTOP)
        && np > ZPK2SOS_MAX_N / 2) return 0;

    complex_t poles[ZPK2SOS_MAX_N];
    complex_t zeros[ZPK2SOS_MAX_N];
    uint8_t degree = np - nz; /* 原型相对阶数 */

    memcpy(poles, proto_poles, (size_t)np * sizeof(complex_t));
    if (nz > 0) {
        /* nz > 0 却没有零点数据：不能继续。下面整条管线会在**未初始化**的栈数组
           上做频率变换、增益折叠与零极点配对（读未初始化值是 UB，而且很可能
           部署出一个看似有效的垃圾设计）。fail-closed 拒绝，而不是静默跳过
           memcpy 接着跑。 */
        if (proto_zeros == NULL) return 0;
        memcpy(zeros, proto_zeros, (size_t)nz * sizeof(complex_t));
    }

    /* LP/BP 增益缩放：wc^degree（LP）或 xi^degree（BP），折叠进下面的
       双线性增益（见 bilinear_zpk_gain_scaled）。 */
    float gs = 0.0f;

    /* 1. 模拟频率变换 */
    switch (type) {
    case FILTER_LOWPASS:
        analog_lp_transform(poles, np, zeros, nz, wc1);
        gs = wc1;
        break;
    case FILTER_HIGHPASS:
        k = zpk_hp_bs_gain(k, zeros, nz, poles, np);
        analog_hp_transform(poles, np, zeros, nz, wc1);
        for (uint8_t i = nz; i < np; i++) {
            zeros[i].re = 0.0f;
            zeros[i].im = 0.0f;
        }
        nz = np;
        break;
    case FILTER_BANDPASS: {
        float w0 = sqrtf(wc1 * wc2);
        float xi = wc2 - wc1;
        analog_bp_transform(poles, &np, zeros, &nz, w0, xi);
        gs = xi;
        break;
    }
    case FILTER_BANDSTOP: {
        float w0 = sqrtf(wc1 * wc2);
        float xi = wc2 - wc1;
        k = zpk_hp_bs_gain(k, zeros, nz, poles, np);
        analog_bs_transform(poles, &np, zeros, &nz, w0, xi);
        break;
    }
    default:
        return 0;
    }

    /* 2. 双线性增益（在 s 域零极点上算，须在双线性变换覆写它们之前）。
       k == 0 会部署出分子全零的"静音"滤波器——那是退化结果，
       不是有效设计。 */
    if (gs != 0.0f) {
        k = bilinear_zpk_gain_scaled(k, gs, degree, zeros, nz, poles, np,
                                     2.0f * fs);
    } else {
        k = bilinear_zpk_gain(k, zeros, nz, poles, np, 2.0f * fs);
    }
    if (!isfinite(k) || k == 0.0f) return 0;

    /* 3. 双线性变换：s → z */
    bilinear_transform(poles, np, fs);
    bilinear_transform(zeros, nz, fs);

    /* 4. 补零：原型在 s=∞ 的零点 → z = −1（BS 不做）。 */
    if (type != FILTER_BANDSTOP) {
        for (uint8_t i = nz; i < np; i++) {
            zeros[i].re = -1.0f;
            zeros[i].im =  0.0f;
        }
        nz = np;
    }

    /* 5. 零极点配对 → SOS 系数。 */
    uint8_t ns = (np + 1) / 2;
    if (ns > max_sections) return 0;

    float sos[ZPK2SOS_MAX_N / 2][6];
    uint8_t n_sections = zpk2sos(zeros, poles, np, sos, k);
    if (n_sections != ns) return 0;

    /* 6. 部署到 biquad 节。 */
    for (uint8_t i = 0; i < n_sections; i++) {
        float num[3] = {sos[i][0], sos[i][1], sos[i][2]};
        float den[3] = {sos[i][3], sos[i][4], sos[i][5]};
        if (!biquad_filter_init(&sections[i], num, den)) return 0;
    }

    return n_sections;
}

uint8_t check_cascade_gains(const biquad_filter_t *sections,
                            uint8_t num_sections,
                            float dc_exp, float ny_exp)
{
    /* 空级联两端乘积都是 1.0，带阻类的 (1, 1) 期望会被空洞通过。0 节不是
       设计——拒绝，而不是盖章放行。 */
    if (num_sections == 0) return 0;

    float h0 = 1.0f, hn = 1.0f;
    for (uint8_t i = 0; i < num_sections; i++) {
        const biquad_filter_t *b = &sections[i];
        h0 *= (b->num_z[0] + b->num_z[1] + b->num_z[2])
            / (1.0f + b->den_z[1] + b->den_z[2]);
        hn *= (b->num_z[0] - b->num_z[1] + b->num_z[2])
            / (1.0f - b->den_z[1] + b->den_z[2]);
    }
    /* 精确结构增益（0 或 1，由 z = ±1 处的零点保证）容差 ±0.1；纹波边缘增益
       （cheby1/2 偶数阶的 10^(−rp/rs/20)）落在带边响应最陡处，放宽到 ±0.25。
       窗口标定与配对缺陷偏差见 test/test_cheby.c 的 sweep 段。 */
    float tol0 = (dc_exp == 0.0f || dc_exp == 1.0f) ? 0.1f : 0.25f;
    float toln = (ny_exp == 0.0f || ny_exp == 1.0f) ? 0.1f : 0.25f;
    return fabsf(h0 - dc_exp) <= tol0 && fabsf(hn - ny_exp) <= toln;
}
