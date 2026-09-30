/**
 * @file    notch_filter.c
 * @brief   二阶陷波滤波器的设计。
 */

#include "notch_filter.h"
#include "filter_utils.h"
#include <math.h>

/* ================================================================== */
/*  f32 设计参数范围闸                                                  */
/* ================================================================== */
/*
 * 判据 Q = ξ·g·(ω0/K)² ≥ 5e-8（K = 2fs）。ω0/K = tan(π·f0/fs) 是预畸本就要
 * 算的量，判据不引入新的 libm 调用。
 *
 * 低于此值时谷底静默变浅、中心偏移，而下游察觉不到：在被拒的参数上求
 * |H(jω0)| 本身就是那个病态运算，而 DC/Nyquist 增益无论设计多坏都精确是 1
 * （z = ±1 处逐项相消，结构保证）。所以只能按参数判，且必须在设计之前判。
 * 标定表（邻域最坏值）、"不能采孤立点"的告诫与范围界定见 test/test_notch.c
 * 的「数值闸边界」段。
 */
#define NOTCH_GATE 5.0e-8f

/* g ≥ 0.9999 时陷波浅于 0.001 dB（g = 1 即 H ≡ 1），按直通部署比烧一个
   高 Q 节划算。 */
#define NOTCH_G_MAX 0.9999f

void notch_init(notch_filter_t *f, float f0, float xi, float g, float fs)
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
    if (!(g > 0.0f)) return;
    if (g >= NOTCH_G_MAX) return; /* 浅到等同不陷波 → 直通 */

    /* 预畸：ω0 必须精确落到 f0 上，否则陷波中心偏移。prewarp 对越界输入返回
       NaN，上面的闸已挡掉越界，这里只兜有限性。 */
    float w0 = 2.0f * (float)M_PI * prewarp(f0, fs);
    float K = 2.0f * fs;
    if (!isfinite(w0)) return;

    /* f32 设计参数范围闸——见文件头。 */
    float r = w0 / K;
    if (!(xi * g * r * r >= NOTCH_GATE)) return;

    /* 直写 biquad：二阶节没有配对问题，不走 design_filter/zpk2sos。
       num_s/den_s 按升幂存放（[s^0, s^1, s^2]），与 biquad_c2d_bilinear 一致。 */
    const float num_s[3] = { w0 * w0, 2.0f * xi * g * w0, 1.0f };
    const float den_s[3] = { w0 * w0, 2.0f * xi * w0, 1.0f };

    float num_z[3];
    float den_z[3];
    biquad_c2d_bilinear(num_z, den_z, num_s, den_s, fs);

    /* 系数有限性、Jury 稳定性、极点半径裕量全在这一步，失败即直通。 */
    if (!biquad_filter_init(&f->sections[0], num_z, den_z)) return;

    /* DC/Nyquist 恒为 1 是结构保证；保留此闸与另两族同构，并兜粗粒度算错。 */
    if (!check_cascade_gains(f->sections, 1, 1.0f, 1.0f)) return;

    f->num_sections = 1;
    f->valid = 1;
}
