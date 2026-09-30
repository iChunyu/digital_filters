/**
 * @file    biquad_filter.c
 * @brief   直接 II 型 biquad 滤波器的实现。
 */

#include "biquad_filter.h"
#include <math.h>

#include "biquad_filter.h"
#include <math.h>

/**
 * @brief 三项补偿求和（TwoSum 式）。
 *
 * 补偿的收益在 reset，不在 Jury 闸：`biquad_filter_reset` 拿 1 + a1 + a2 做
 * 分母，裸 f32 求和可能把真值 ~3e-8 的余量舍成 0.0f，于是复位把状态清零，
 * 而按 f32 系数的真稳态该是 equilibrium / 3e-8。折回每次加法的舍入误差即可
 * 在 f32 代价下拿回真余量。
 *
 * 注意 Jury 那条的真余量要么落在 Sterbenz 精确区（|a1| ∈ [0.5, 2]），要么
 * 对应极点落在 z ≈ 1 而被下面的半径裕量闸正当拒掉，所以补偿从不改变 init
 * 的判决——此处用它与 reset 同口径，也为将来单独放宽半径闸时仍按真余量判。
 * 实测与标定见 test/test_biquad.c 的「sum3f 补偿」段。
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
    c += (fabsf(s) >= fabsf(y)) ? (s - t) + y : (y - t) + s;
    s = t;
    t = s + z;
    c += (fabsf(s) >= fabsf(z)) ? (s - t) + z : (z - t) + s;
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
       悄悄部署出分子全零的"静音"滤波器——直通是更安全的兜底。
       注意本闸只拦 0/非有限：den_z[0] = 1e30f 这种**有限但极大**的值会归一化
       出 b0 ≈ 1e-30 的近静音滤波器并返回 1。管线内由 check_cascade_gains 的
       DC/Nyquist 窗口兜住；直接调用本函数的调用方需自行确认增益量级。 */
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

    /* 拒绝非有限系数（如管线里的增益溢出）：NaN/Inf 分子会通过下面的 Jury
       检查（它只看分母）然后毒化输出。 */
    if (!(isfinite(filter->num_z[0]) && isfinite(filter->num_z[1])
          && isfinite(filter->num_z[2]) && isfinite(filter->den_z[1])
          && isfinite(filter->den_z[2]))) {
        biquad_filter_set_empty(filter);
        return 0;
    }

    /*
     * 稳定性检查——二阶系统的三个 Jury 条件：
     *   |a2| < 1、1 + a1 + a2 > 0、1 - a1 + a2 > 0
     * 任一不满足即不稳定，退回单位直通。判定用补偿和（见 sum3f）：本闸独立
     * 承担教科书稳定性契约，不依赖下面半径裕量闸的存在——后者虽然蕴含了这
     * 三条，但两层判据写法的意图不同（见 test/test_biquad.c 的档位）。
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
     * a2 = 1 − 2⁻²⁴ 照样通过，滤波器会振铃数百万个样本（每样本收缩 ~6e-8）。
     * 拒绝任何极点半径 r > 0.99995（1 − r < 5e-5）。判定按半径本身做、不按
     * a2：a2 = r1·r2 对主导极点不是 a2 的实根对有盲区，单侧的 a2 > 0.9999
     * 检查也会漏掉异号实根对。共轭对看 a2 = r²；实根用
     * r_max = (|a1| + √(a1² − 4a2))/2，两边平方改写成不需要 sqrtf 的形式——
     * 这条 init 路径要保持可在中断上下文中廉价执行。判别式用 Dekker 分裂
     * 补偿：宽带 BP/BS 的近实极点对落在距单位圆 ~1e-4 处，裸 f32 判别式在那
     * 里带 ~5e-7 噪声，足以对合法设计抛硬币。边界档位见 test/test_biquad.c。
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
    /*
     * 非有限的 equilibrium 会用 NaN 毒化状态向量，并传到此后每一个输出。
     * 退回零状态（单位直通自身的稳态）。
     */
    if (!isfinite(equilibrium)) {
        biquad_zero_state(filter);
        return;
    }

    /*
     * 稳态：x 为常数时 w[0]=w[1]=w[2]=w_ss，由状态方程
     *   x = w_ss + a1·w_ss + a2·w_ss  =>  w_ss = x / (1 + a1 + a2)。
     *
     * 这是公开结构体上的公开 API，系数不必先过 biquad_filter_init()。补偿求和
     * 避免分母在窄带滤波器上恰好消成 0；下面的守卫实现 @note 里的契约：
     * 1 + a1 + a2 == 0（如纯积分器）无稳态 → 清零，分母小到让 w_ss 溢出为
     * inf、进而毒化此后每次 update 也清零。
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
