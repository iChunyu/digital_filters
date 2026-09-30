/* 实现与闸门依据见 docs/biquad_filter.md。 */

#include "biquad_filter.h"
#include <math.h>

/* 三项补偿求和。收益在 reset 的分母：真余量可小到 ~3e-8，裸 f32 会舍成 0；Jury 闸同口径但判决不变。 */
static float sum3f(float x, float y, float z)
{
    float s = x;
    float c = 0.0f;
    float t = s + y;
    c += (fabsf(s) >= fabsf(y)) ? (s - t) + y : (y - t) + s;
    s = t;
    t = s + z;
    c += (fabsf(s) >= fabsf(z)) ? (s - t) + z : (z - t) + s;
    s = t;
    return s + c;
}

static void biquad_zero_state(biquad_filter_t *filter)
{
    filter->w[0] = 0.0f;
    filter->w[1] = 0.0f;
    filter->w[2] = 0.0f;
}

void biquad_filter_set_empty(biquad_filter_t *filter)
{
    filter->num_z[0] = 1.0f;
    filter->num_z[1] = 0.0f;
    filter->num_z[2] = 0.0f;

    filter->den_z[0] = 1.0f;
    filter->den_z[1] = 0.0f;
    filter->den_z[2] = 0.0f;

    filter->w[0] = 0.0f;
    filter->w[1] = 0.0f;
    filter->w[2] = 0.0f;
}

void biquad_c2d_bilinear(float num_z[3], float den_z[3], const float num_s[3],
                        const float den_s[3], float fs)
{
    float K = 2.0f * fs;
    float K2 = K * K;

    num_z[0] = num_s[0] + num_s[1] * K + num_s[2] * K2;
    num_z[1] = 2.0f * num_s[0] - 2.0f * num_s[2] * K2;
    num_z[2] = num_s[0] - num_s[1] * K + num_s[2] * K2;

    den_z[0] = den_s[0] + den_s[1] * K + den_s[2] * K2;
    den_z[1] = 2.0f * den_s[0] - 2.0f * den_s[2] * K2;
    den_z[2] = den_s[0] - den_s[1] * K + den_s[2] * K2;
}

uint8_t biquad_filter_init(biquad_filter_t *filter, const float num_z[3],
                           const float den_z[3])
{
    /* den_z[0] 零/非有限 → 直通。 */
    /* 有限但极大的 den_z[0]（如 1e30）会归一化出近静音滤波器并返回 1，管道内由增益窗口兜住。 */
    if (den_z[0] == 0.0f || !isfinite(den_z[0])) {
        biquad_filter_set_empty(filter);
        return 0;
    }

    float inv = 1.0f / den_z[0];

    filter->num_z[0] = num_z[0] * inv;
    filter->num_z[1] = num_z[1] * inv;
    filter->num_z[2] = num_z[2] * inv;

    filter->den_z[0] = 1.0f;
    filter->den_z[1] = den_z[1] * inv;
    filter->den_z[2] = den_z[2] * inv;

    if (!(isfinite(filter->num_z[0]) && isfinite(filter->num_z[1])
          && isfinite(filter->num_z[2]) && isfinite(filter->den_z[1])
          && isfinite(filter->den_z[2]))) {
        biquad_filter_set_empty(filter);
        return 0;
    }

    /* Jury 三条件。教科书契约，独立于下面的半径裕量闸（后者蕴含它）。 */
    float a1 = filter->den_z[1];
    float a2 = filter->den_z[2];

    if (!(a2 > -1.0f && a2 < 1.0f
          && sum3f(1.0f, a1, a2) > 0.0f
          && sum3f(1.0f, -a1, a2) > 0.0f)) {
        biquad_filter_set_empty(filter);
        return 0;
    }

    /* r > 0.99995 拒（f32 下会振铃数百万样本）。按半径判而非按 a2：a2 = r1·r2 对实根对有盲区。判别式用 Dekker 分裂，宽带 BP/BS 的近实极点对那里有 ~5e-7 噪声。 */
    float p = a1 * 4097.0f;
    float hi = p - (p - a1);
    float lo = a1 - hi;
    float sq = a1 * a1;
    float err = ((hi * hi - sq) + 2.0f * hi * lo) + lo * lo;
    float disc = (sq - 4.0f * a2) + err;
    int reject;
    if (disc < 0.0f) {
        reject = a2 > 0.9999f;
    } else {
        /* r_max > 0.99995 ⟺ √disc > 1.9999 − |a1|，两边平方以免 sqrtf。 */
        float rhs = 1.9999f - fabsf(a1);
        reject = (rhs < 0.0f) || (disc > rhs * rhs);
    }
    if (reject) {
        biquad_filter_set_empty(filter);
        return 0;
    }

    filter->w[0] = 0.0f;
    filter->w[1] = 0.0f;
    filter->w[2] = 0.0f;
    return 1;
}

float biquad_filter_get_output(const biquad_filter_t *filter)
{
    return filter->num_z[0] * filter->w[0]
         + filter->num_z[1] * filter->w[1]
         + filter->num_z[2] * filter->w[2];
}

float biquad_filter_get_input(const biquad_filter_t *filter)
{
    return filter->w[0]
         + filter->den_z[1] * filter->w[1]
         + filter->den_z[2] * filter->w[2];
}

void biquad_filter_reset(biquad_filter_t *filter, float equilibrium)
{
    if (!isfinite(equilibrium)) {
        biquad_zero_state(filter);
        return;
    }

    /* 稳态 w_ss = x/(1+a1+a2)；分母为 0/非有限或 w_ss 溢出 → 清零，永不产生 inf/NaN。 */
    float denom = sum3f(1.0f, filter->den_z[1], filter->den_z[2]);
    if (denom == 0.0f || !isfinite(denom)) {
        biquad_zero_state(filter);
        return;
    }
    float w_ss = equilibrium / denom;
    if (!isfinite(w_ss)) {
        biquad_zero_state(filter);
        return;
    }
    filter->w[0] = w_ss;
    filter->w[1] = w_ss;
    filter->w[2] = w_ss;
}
