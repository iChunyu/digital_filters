/* 设计流程与闸门推导见 docs/notch_filter.md。 */

#include "notch_filter.h"
#include "filter_utils.h"
#include <math.h>

/* 预测式闸：低于此值陷波深度承载不住，按参数拒而不放行。标定见 tests/test_notch.c。 */
#define NOTCH_GATE 5.0e-8f

/* g ≥ 0.9999 时数学上恒等于直通，不值得烧一个高 Q 节。 */
#define NOTCH_G_MAX 0.9999f

void notch_init(notch_filter_t *f, float f0, float xi, float g, float fs)
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
    if (!(g > 0.0f)) return;
    if (g >= NOTCH_G_MAX) return;

    float w0 = 2.0f * (float)M_PI * prewarp(f0, fs);
    float K = 2.0f * fs;
    if (!isfinite(w0)) return;

    float r = w0 / K;
    if (!(xi * g * r * r >= NOTCH_GATE)) return;

    const float num_s[3] = { w0 * w0, 2.0f * xi * g * w0, 1.0f };
    const float den_s[3] = { w0 * w0, 2.0f * xi * w0, 1.0f };

    float num_z[3];
    float den_z[3];
    biquad_c2d_bilinear(num_z, den_z, num_s, den_s, fs);

    if (!biquad_filter_init(&f->sections[0], num_z, den_z)) return;

    if (!check_cascade_gains(f->sections, 1, 1.0f, 1.0f)) return;

    f->num_sections = 1;
    f->valid = 1;
}
