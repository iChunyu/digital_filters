/* 对偶、统一参数 d = 1/g 与闸门依据见 docs/peak_filter.md。 */

#include "peak_filter.h"
#include "filter_utils.h"
#include <math.h>

/* 判据用 d = 1/g 而非 g：写成 g 会让阈值随增益放大约 g² 倍，放行病态 boost 而拒掉其镜像陷波。 */
#define PEAK_GATE 5.0e-8f

#define PEAK_D_MAX 0.9999f

void peak_init(peak_filter_t *f, float f0, float xi, float g, float fs)
{
    f->f0 = f0;
    f->xi = xi;
    f->g = g;
    f->fs = fs;
    f->valid = 0;
    f->num_sections = 0;

    if (!(fs > 0.0f)) return;
    if (!(f0 > 0.0f) || !(f0 < fs * 0.5f)) return;
    if (!(xi > 0.0f)) return;
    if (!(g > 1.0f)) return;

    float d = 1.0f / g;

    /* g = +Inf → d = 0，由数值闸以 Q = 0 拒掉；通常由极点半径闸先触发。 */
    if (!(d < PEAK_D_MAX)) return;

    float w0 = 2.0f * (float)M_PI * prewarp(f0, fs);
    float K = 2.0f * fs;
    if (!isfinite(w0)) return;

    float r = w0 / K;
    if (!(xi * d * r * r >= PEAK_GATE)) return;

    /* 分母写 2·xi·d·w0（而非 2·xi·w0/g）是为与 notch 分子逐位同形，使 s 域系数向量互为精确交换。 */
    const float num_s[3] = { w0 * w0, 2.0f * xi * w0, 1.0f };
    const float den_s[3] = { w0 * w0, 2.0f * xi * d * w0, 1.0f };

    float num_z[3];
    float den_z[3];
    biquad_c2d_bilinear(num_z, den_z, num_s, den_s, fs);

    if (!biquad_filter_init(&f->sections[0], num_z, den_z)) return;

    if (!check_cascade_gains(f->sections, 1, 1.0f, 1.0f)) return;

    f->num_sections = 1;
    f->valid = 1;
}
