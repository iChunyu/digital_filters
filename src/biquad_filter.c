/**
 * @file    biquad_filter.c
 * @brief   直接 II 型 biquad 滤波器的实现。
 */

#include "biquad_filter.h"
#include <math.h>

/**
 * @brief 浮点绝对值，不依赖 libm，也不依赖编译器内建。
 *
 * 与 `fabsf` 在所有输入上判定等价：有限值与 ±inf 逐位一致；-0.0f
 * 返回 -0.0f（与 +0.0f 比较相等，本文件的用法全是幅值比较，不受
 * 影响）；NaN 返回 NaN，而 `NaN >= x` 恒假，判定结果与 `fabsf(NaN)`
 * 相同。
 *
 * 存在的理由：`fabsf` 在 GCC/Clang 上默认是内建并被内联，但
 * `-fno-builtin` 下退化成真实的 libm 符号。本库对外承诺 biquad 层
 * 零 libm 依赖（裸机不链 libm 也能用），这个承诺不该建立在编译器
 * 的内建行为上。
 *
 * @param x  输入值。
 * @return   |x|。
 */
static inline float abs_f(float x)
{
    return (x < 0.0f) ? -x : x;
}

/**
 * @brief 三项补偿求和（TwoSum 式）。
 *
 * 裸 f32 计算 1 + a1 + a2 在真余量为 1.0 的若干 ulp 时会恰好消成
 * 0.0f（窄带设计：a1 ≈ −2，a2 ≈ 1 − ε）——把极点远离单位圆的
 * 稳定滤波器误拒，同样也会误算 reset 的分母。把每次加法的
 * 舍入误差折回即可在 f32 代价下恢复余量。
 *
 * @param x  第一加数。
 * @param y  第二加数。
 * @param z  第三加数。
 * @return   补偿后的和 x + y + z。
 */
static float sum3f(float x, float y, float z)
{
    float s = x;
    float c = 0.0f;
    float t = s + y;
    c += (abs_f(s) >= abs_f(y)) ? (s - t) + y : (y - t) + s;
    s = t;
    t = s + z;
    c += (abs_f(s) >= abs_f(z)) ? (s - t) + z : (z - t) + s;
    s = t;
    return s + c;
}

/**
 * @brief 状态向量清零。
 *
 * @param[out] filter  滤波器对象指针。
 */
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
    /* 拒绝为零/非有限的首项分母。零会导致除零；±Inf 会让 inv = 0，
       悄悄部署出分子全零的"静音"滤波器——直通是更安全的兜底。 */
    if (den_z[0] == 0.0f || !isfinite(den_z[0])) {
        biquad_filter_set_empty(filter);
        return 0;
    }

    /* 归一化，使 den_z[0] == 1.0 */
    float inv = 1.0f / den_z[0];

    filter->num_z[0] = num_z[0] * inv;
    filter->num_z[1] = num_z[1] * inv;
    filter->num_z[2] = num_z[2] * inv;

    filter->den_z[0] = 1.0f;
    filter->den_z[1] = den_z[1] * inv;
    filter->den_z[2] = den_z[2] * inv;

    /*
     * 拒绝非有限系数（例如设计管线里的增益溢出）。否则 NaN/Inf 分子会
     * 通过下面的 Jury 检查——它只看分母——然后毒化输出。
     */
    if (!(isfinite(filter->num_z[0]) && isfinite(filter->num_z[1])
          && isfinite(filter->num_z[2]) && isfinite(filter->den_z[1])
          && isfinite(filter->den_z[2]))) {
        biquad_filter_set_empty(filter);
        return 0;
    }

    /*
     * 稳定性检查——二阶系统的三个 Jury 条件：
     *   |a2| < 1
     *   1 + a1 + a2 > 0
     *   1 - a1 + a2 > 0
     * 任一条件不满足即不稳定；退回单位直通。
     * 两个 Jury 和用补偿求和（见 sum3f）：对真余量为 1.0 的若干 ulp 的
     * 窄带设计，裸 f32 求值会恰好消成 0.0f，把稳定滤波器误拒。
     * 反过来，补偿和**恰好**为 0.0f 说明 f32 系数把极点放到了单位圆上
     * （量化塌缩到 z = ±1）——那种确实该拒，且下面的裕量检查分辨不出
     * 它的半径。
     */
    float a1 = filter->den_z[1];
    float a2 = filter->den_z[2];

    if (!(a2 > -1.0f && a2 < 1.0f
          && sum3f(1.0f, a1, a2) > 0.0f
          && sum3f(1.0f, -a1, a2) > 0.0f)) {
        biquad_filter_set_empty(filter);
        return 0;
    }

    /*
     * 稳定性裕量。上面的 Jury 条件允许极点任意贴近单位圆；f32 下
     * a2 = 1 − 2⁻²⁴ 照样通过，滤波器会振铃数百万个样本（每样本收缩
     * ~6e-8）。拒绝任何半径 r > 0.99995（1 − r < 5e-5）的极点。
     *
     * 判定按极点半径本身做，不按 a2：a2 = r1·r2 对非等实根对的主导极点
     * 有盲区（r1 ≈ 1、r2 ≈ 0.85 → a2 ≈ 0.85 轻松过关），而单侧的
     * a2 > 0.9999 检查会漏掉异号实根对（a2 ≈ −0.99995）。
     * 共轭对 r² = a2（于是 r > 0.99995 ⟺ a2 > 0.9999）；实根
     * r_max = (|a1| + √(a1² − 4·a2)) / 2，经过改写使其不需要 sqrtf——
     * biquad init 路径必须保持可被不链 libm 的裸机固件调用。该公式
     * 同样覆盖一阶节（a2 = 0 → r_max = |a1|）。
     *
     * a1² − 4·a2 用 Dekker 分裂补偿计算：宽带 BP/BS 设计的近实极点对
     * 落在距单位圆 ~1e-4 处，那里裸 f32 判别式带 ~5e-7 噪声，而真值
     * （如 −4·(虚部)² ≈ −1e-7）并不更大——共轭/实根分支之间的符号翻转
     * 会把裕量判定变成对合法设计的抛硬币。
     */
    float p = a1 * 4097.0f;          /* Dekker 分裂：a1 = hi + lo */
    float hi = p - (p - a1);
    float lo = a1 - hi;
    float sq = a1 * a1;
    float err = ((hi * hi - sq) + 2.0f * hi * lo) + lo * lo;
    float disc = (sq - 4.0f * a2) + err;
    int reject;
    if (disc < 0.0f) {
        reject = a2 > 0.9999f;                 /* 共轭对：a2 = r² */
    } else {
        /* r_max > 0.99995  ⟺  √disc > 1.9999 − |a1|（两边平方）。 */
        float rhs = 1.9999f - abs_f(a1);
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
    /*
     * 非有限的 equilibrium 会用 NaN 毒化状态向量，并传到此后每一个输出。
     * 退回零状态（单位直通自身的稳态）。
     */
    if (!isfinite(equilibrium)) {
        biquad_zero_state(filter);
        return;
    }

    /*
     * 稳态：x 为常数时 w[0]=w[1]=w[2]=w_ss。
     * 由状态方程：
     *   x = w_ss + a1*w_ss + a2*w_ss  =>  w_ss = x / (1 + a1 + a2)
     *
     * 这是公开结构体上的公开 API——系数不必先经过 biquad_filter_init()。
     * 补偿求和避免分母的 f32 求值在 init 通过的窄带滤波器上恰好消成
     * 0.0f；下面的守卫实现 @note 里的契约：1 + a1 + a2 == 0（如纯积分器）
     * 无稳态 → 强制清零；分母小到非规格化会让 w_ss 溢出为 inf、进而用
     * NaN 毒化此后每次 update，所以在那里也强制清零。
     */
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
