/**
 * @file    peak_filter.c
 * @brief   二阶峰值滤波器（Peak / boost）的设计。
 */

#include "peak_filter.h"
#include "filter_utils.h"
#include <math.h>

/* ================================================================== */
/*  原型选择：分母带 g 的 H1                                           */
/* ================================================================== */
/*
 * H1 = (s² + 2ξω₀s + ω₀²)/(s² + 2ξω₀/g·s + ω₀²) 与陷波原型
 * H2 = (s² + 2ξgω₀s + ω₀²)/(s² + 2ξω₀s + ω₀²) 满足 H1(ξ,g) ≡ 1/H2(ξ,1/g)，
 * 这是本族与 notch 互为倒数、boost 与 cut 带宽语义对称的来源。
 * 对偶不变量在 test/test_peak.c 里逐系数校验。
 */

/* ================================================================== */
/*  f32 设计参数范围闸                                                  */
/* ================================================================== */
/*
 * 判据 Q = ξ·d·(ω0/K)² ≥ 5e-8（d = 1/g）。必须用 d 而不是 g：由上面的对偶，
 * "g 的 boost 有多难做"与"d 的陷波有多难做"是同一个问题，写成 g 会把闸门随
 * 增益放大约 g² 倍——同一个病态 boost 被放行，它的镜像陷波却被拒。
 *
 * 低于此值时峰高与中心都不再可信（实测偏差双符号，目标 g = 10 时两侧出现过
 * 7.113 与 11.06），本函数直接拒（valid = 0）而不放行。标定表、两条闸的先后
 * 关系、以及半径闸刻度 g ≤ 4e4·ξ·tan(πf0/fs) 的来历见 test/test_peak.c 的
 * 「数值闸边界」段。
 */
#define PEAK_GATE 5.0e-8f

/* d ≥ 0.9999（g ≤ 1.0001）时峰值浅于 0.001 dB，等同不 boost，按直通部署。 */
#define PEAK_D_MAX 0.9999f

void peak_init(peak_filter_t *f, float f0, float xi, float g, float fs)
{
    f->f0 = f0;
    f->xi = xi;
    f->g = g;
    f->fs = fs;
    f->valid = 0;
    f->num_sections = 0; /* valid=0 时绝不能留下这个垃圾值 */

    /* 参数闸一律写 !(x > 0) 而非 x <= 0：后者与 NaN 比较恒假，会放它过去。 */
    if (!(fs > 0.0f)) return;
    if (!(f0 > 0.0f) || !(f0 < fs * 0.5f)) return;
    if (!(xi > 0.0f)) return;
    if (!(g > 1.0f)) return;

    /* 等效陷波深度：本族数值分析的统一参数。 */
    float d = 1.0f / g;

    /* g = +Inf → d = 0，不在此拦，由下面的数值闸以 Q = 0 拒掉。 */
    if (!(d < PEAK_D_MAX)) return; /* 浅到等同不 boost → 直通 */

    /* 预畸：ω0 必须精确落到 f0 上，否则峰中心偏移。prewarp 对越界输入返回
       NaN，上面的闸已挡掉越界，这里只兜有限性。 */
    float w0 = 2.0f * (float)M_PI * prewarp(f0, fs);
    float K = 2.0f * fs;
    if (!isfinite(w0)) return;

    /* f32 设计参数范围闸——见文件头。判据用 d 而不是 g。 */
    float r = w0 / K;
    if (!(xi * d * r * r >= PEAK_GATE)) return;

    /* 直写 biquad：二阶节没有配对问题，不走 design_filter/zpk2sos。
       num_s/den_s 按升幂存放（[s^0, s^1, s^2]），与 biquad_c2d_bilinear 一致。
       分母写成 2·xi·d·w0（而非 2·xi·w0/g）是为了与 notch 的分子逐位同形，
       使两者的 s 域系数向量互为精确交换——校验见 test/test_peak.c。 */
    const float num_s[3] = { w0 * w0, 2.0f * xi * w0, 1.0f };
    const float den_s[3] = { w0 * w0, 2.0f * xi * d * w0, 1.0f };

    float num_z[3];
    float den_z[3];
    biquad_c2d_bilinear(num_z, den_z, num_s, den_s, fs);

    /* 三道闸都在 init 内，失败即直通。对本族这是最常触发的一道（阻尼 ξ·d）。 */
    if (!biquad_filter_init(&f->sections[0], num_z, den_z)) return;

    /* DC/Nyquist 恒为 1 是结构保证；保留此闸与另两族同构，并兜粗粒度算错。 */
    if (!check_cascade_gains(f->sections, 1, 1.0f, 1.0f)) return;

    f->num_sections = 1;
    f->valid = 1;
}
