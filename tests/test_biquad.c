#include "biquad_filter.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures = 0;

#define CHECK(cond, msg)                                         \
    do {                                                         \
        if (!(cond)) {                                           \
            failures++;                                          \
            fprintf(stderr, "FAIL: %s\n", msg);                 \
        }                                                        \
    } while (0)

#define CLOSE(a, b, eps) (fabsf((a) - (b)) <= (eps))

int main(void)
{
    /* ── 单位直通 / 透传 ────────────────────────────────────────── */

    biquad_filter_t id;
    biquad_filter_set_empty(&id);

    float y = biquad_filter_update(&id, 0.5f);
    CHECK(y == 0.5f, "identity DC");

    y = biquad_filter_update(&id, -0.25f);
    CHECK(y == -0.25f, "identity negative");

    /* ── 已知滤波器：二阶低通 ───────────────────────────────────── */
    /* s 域：H(s) = 1 / (s² + √2·s + 1)（Butterworth，fc = 1 Hz） */
    float num_s[3] = {1.0f, 0.0f, 0.0f};
    float den_s[3] = {1.0f, 1.41421356f, 1.0f};

    float num_z[3], den_z[3];
    biquad_c2d_bilinear(num_z, den_z, num_s, den_s, 10.0f);

    biquad_filter_t lpf;
    biquad_filter_init(&lpf, num_z, den_z);

    /* DC 增益应为 1.0 */
    biquad_filter_reset(&lpf, 1.0f);
    y = biquad_filter_update(&lpf, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "LPF DC gain ~ 1");

    /* ── 零系数分母的降级 ───────────────────────────────────────── */

    float zero_den[3] = {0.0f, 1.0f, 1.0f};
    biquad_filter_init(&lpf, num_z, zero_den);
    y = biquad_filter_update(&lpf, 0.5f);
    CHECK(y == 0.5f, "zero den_z[0] falls back to identity");

    /* ── 非有限分子 → 降级为单位直通 ────────────────────────────── */

    float nan_num[3] = {NAN, 1.0f, 1.0f};
    biquad_filter_t nf;
    uint8_t rc = biquad_filter_init(&nf, nan_num, num_z);
    CHECK(rc == 0, "NaN numerator rejected");
    y = biquad_filter_update(&nf, 0.75f);
    CHECK(y == 0.75f, "NaN numerator → identity passthrough");

    /* ── 首项分母为 Inf → 直通，而不是静音 ─────────────────────────── */

    float inf_den[3] = {INFINITY, 0.0f, 0.0f};
    rc = biquad_filter_init(&nf, num_z, inf_den);
    CHECK(rc == 0, "den_z[0] = Inf rejected");
    y = biquad_filter_update(&nf, 0.5f);
    CHECK(y == 0.5f, "den_z[0] = Inf → identity passthrough");

    /* ── 稳定性裕量：贴近单位圆的极点被拒 ─────────────────────────────
       阈值 r > 0.99995：f32 下 a2 = 1 − 2⁻²⁴ 照样通过 Jury，但滤波器会振铃数百
       万个样本（每样本收缩 ~6e-8），故按极点半径本身拒绝。─────── */

    float marg_num[3] = {1.0f, 0.0f, 0.0f};

    float marg_den_bad[3] = {1.0f, 0.0f, 0.99999f};   /* a2 略小于 1 */
    rc = biquad_filter_init(&nf, marg_num, marg_den_bad);
    CHECK(rc == 0, "a2 = 0.99999 rejected (margin)");

    float marg_den_ok[3] = {1.0f, 0.0f, 0.9998f};     /* 在裕量之内 */
    rc = biquad_filter_init(&nf, marg_num, marg_den_ok);
    CHECK(rc == 1, "a2 = 0.9998 accepted");

    float marg1_bad[3] = {1.0f, 0.99999f, 0.0f};      /* 一阶 |a1| ~ 1 */
    rc = biquad_filter_init(&nf, marg_num, marg1_bad);
    CHECK(rc == 0, "1st-order |a1| = 0.99999 rejected (margin)");

    /* ── 非有限 equilibrium 的 reset → 状态清零 ─────────────────────── */

    biquad_filter_t fr;
    biquad_filter_init(&fr, num_z, den_z);
    biquad_filter_reset(&fr, NAN);
    CHECK(fr.w[0] == 0.0f && fr.w[1] == 0.0f && fr.w[2] == 0.0f,
          "reset(NaN) zeroes the state");
    y = biquad_filter_get_output(&fr);
    CHECK(y == 0.0f, "reset(NaN) → zero output");

    /* ── get_output / get_input 一致性 ──────────────────────────── */

    float num[3] = {0.2f, 0.4f, 0.2f};
    float den[3] = {1.0f, 0.0f, 0.0f};
    biquad_filter_t f;

    /* ── get_output / get_input 一致性（非退化系数）────────────────────
       旧档位用 den = {1,0,0} 且只推一次，于是 w[1] = w[2] = 0，两个函数里的
       a·w 与 b·w 多项式项永不参与——把 get_input 简化成只返回 w[0]、
       get_output 简化成 b0·w[0] 也不会变红（实测：删掉 b1·w1 + b2·w2
       两个变异体都能存活）。下面推三次让 w[1]、w[2] 都非零，三项真参与。 */
    float gnum[3] = {0.2f, 0.4f, 0.2f};
    float gden[3] = {1.0f, -0.3f, -0.2f};
    CHECK(biquad_filter_init(&f, gnum, gden) == 1, "get_* fixture deploys");

    (void)biquad_filter_update(&f, 0.3f);
    (void)biquad_filter_update(&f, 0.7f);
    float out = biquad_filter_update(&f, 0.5f);
    CHECK(f.w[1] != 0.0f && f.w[2] != 0.0f,
          "get_* fixture keeps w[1] and w[2] nonzero (non-degenerate)");
    CHECK(CLOSE(biquad_filter_get_input(&f), 0.5f, 1e-6f),
          "get_input reconstructs x[n] with nonzero a1/a2 and w[1]");
    CHECK(CLOSE(biquad_filter_get_output(&f), out, 1e-6f),
          "get_output matches the last update() with nonzero b1/b2");

    /* ── 裕量对称性：异号实根对同样被拒 ──────────────────────────────────
       单侧的 a2 > 0.9999 检查会漏掉 a2 ≈ −0.99995（实根 ±0.99997），
       而基于乘积的检查对非等实根对的主导极点有盲区。 ──────────────────── */

    float marg_den_neg[3] = {1.0f, 0.0f, -0.99999f};   /* 实根对 ±0.999995 */
    rc = biquad_filter_init(&nf, marg_num, marg_den_neg);
    CHECK(rc == 0, "a2 = -0.99999 rejected (real pair, symmetric margin)");

    float dom_num[3] = {1.0f, 0.0f, 0.0f};
    /* 极点 0.999999 与 0.85：a2 = 0.85 能通过基于乘积的检查，
       且 Jury 和 1 + a1 + a2 = 1e-6 为正——只有对主导极点的
       半径检查才能拒掉它。 */
    float dom_den[3] = {1.0f, -1.849999f, 0.85f};
    rc = biquad_filter_init(&nf, dom_num, dom_den);
    CHECK(rc == 0, "dominant pole 0.999999 rejected (product a2=0.85 blind spot)");

    /* ── 裕量阈值两侧的紧档位 ──────────────────────────────────────
       上面的档位只把常数夹在 (3e-5, 2.5e-4) 这个很宽的区间里；这里用复根
       （a2 = r²）与实根（主导极点 = |a1|）两条路径各自把 0.9999f / 1.9999f
       两个常数钉到 ±2e-5。────────────────────────────────────── */

    float mth_c_ok[3]  = {1.0f, 0.0f, 0.99988f};    /* 复根：r = 0.999940 */
    CHECK(biquad_filter_init(&nf, marg_num, mth_c_ok) == 1,
          "conjugate pair r = 0.99994 accepted (margin upper side)");

    float mth_c_bad[3] = {1.0f, 0.0f, 0.99992f};    /* 复根：r = 0.999960 */
    CHECK(biquad_filter_init(&nf, marg_num, mth_c_bad) == 0,
          "conjugate pair r = 0.99996 rejected (margin lower side)");

    float mth_r_ok[3]  = {1.0f, -0.99994f, 0.0f};   /* 实根：主导极点 0.999940 */
    CHECK(biquad_filter_init(&nf, marg_num, mth_r_ok) == 1,
          "real root 0.99994 accepted (margin upper side)");

    float mth_r_bad[3] = {1.0f, -0.99996f, 0.0f};   /* 实根：主导极点 0.999960 */
    CHECK(biquad_filter_init(&nf, marg_num, mth_r_bad) == 0,
          "real root 0.99996 rejected (margin lower side)");

    /* ── 不稳定极点在单位圆外必须被拒（教科书 Jury 契约）───────────
       当前标定下这三条 Jury 判据被上面的半径裕量闸蕴含（要让 Jury 单独变红
       需要一个极点落在单位圆外，而那种半径必被裕量闸拒），所以这里测的是
       **外部契约**：任何不稳定极点都不得被部署。────────────── */

    static const float unst[4][3] = {
        {1.0f, -2.05f,  1.05f},   /* 实根 1.0 与 1.05 */
        {1.0f,  0.0f,   1.0001f}, /* 复根半径 1.00005 */
        {1.0f,  0.0f,  -1.0001f}, /* 实根对 ±1.00005 */
        {1.0f, -2.0f,   1.0f},    /* 重根 z = 1（1 + a1 + a2 == 0） */
    };
    for (int i = 0; i < 4; i++) {
        CHECK(biquad_filter_init(&nf, marg_num, unst[i]) == 0,
              "poles outside the unit circle rejected");
        y = biquad_filter_update(&nf, 1.0f);
        CHECK(y == 1.0f, "rejected unstable design → identity passthrough");
    }

    /* ── 纯积分器（1 + a1 + a2 == 0）的 reset ────────────────────────────
       @note 契约：稳态不存在，状态必须强制清零——绝不出现 inf/NaN。
       biquad_filter_reset 是公开结构体上的公开 API，系数不必先过 init。 ── */

    biquad_filter_t integ;
    biquad_filter_set_empty(&integ);
    integ.den_z[1] = -1.0f;    /* H(z) = 1 / (1 − z⁻¹)：极点在 z = 1 */
    biquad_filter_reset(&integ, 1.0f);
    CHECK(integ.w[0] == 0.0f && integ.w[1] == 0.0f && integ.w[2] == 0.0f,
          "reset(integrator) forces zero state (not inf)");
    y = biquad_filter_update(&integ, 0.5f);
    CHECK(isfinite(y), "integrator update stays finite after guarded reset");

    /* ── reset 的 w_ss 溢出分支 → 状态清零 ────────────────────────────
       分母非 0 但极小、equilibrium 又极大时 w_ss 会溢出成 inf，该分支此前
       无档位覆盖（只要分母不过 0，@note 里那三道守卫就只剩这一道可练）。 */
    biquad_filter_t ovf;
    biquad_filter_set_empty(&ovf);
    ovf.den_z[1] = -0.999999f;           /* 1 + a1 + a2 = 1e-6，不过 0 */
    biquad_filter_reset(&ovf, 3.4e38f);  /* FLT_MAX 量级 → w_ss 溢出 */
    CHECK(ovf.w[0] == 0.0f && ovf.w[1] == 0.0f && ovf.w[2] == 0.0f,
          "reset(w_ss overflow) zeroes the state");
    y = biquad_filter_update(&ovf, 1.0f);
    CHECK(isfinite(y), "update after w_ss-overflow reset stays finite");

    /* ── sum3f 补偿：为什么有、以及为什么只在 reset 上可观测 ──────────────
       sum3f 把每次加法的舍入误差折回，服务于 1 + a1 + a2 的两条路径：
       biquad_filter_init 的 Jury 闸，与 biquad_filter_reset 的分母。

       标定（实测，三种手段互相印证）：

       1. |a1| ∈ [0.5, 2] 落在 Sterbenz 引理的精确区，两次加法都无舍入。在
          a1 ≈ −2、a2 ≈ 1 − ε（源码注释原先举的例子）上系统搜 68008 组 f32
          系数对，裸和与补偿和**逐位相同**——补偿在那里是 no-op。
       2. 真正会让裸和塌成 0 的区域是 |a1| < 0.5 且 a2 ≈ −1 − a1（搜 52006
          组找到 5 例）。但那种系数的极点落在 z ≈ 1，init 的极点半径裕量闸
          会正当拒绝它，所以补偿**不改变 init 的判决**。
       3. 全套件插桩：sum3f 被调 13304 次，补偿与裸和的判决差异 **0 次**。
          把 sum3f 整体换成裸求和的变异体，除本档位外全部测试仍绿。

       结论：sum3f 的收益只在 reset——让复位用**按 f32 系数的真稳态**，而不是
       把 1e-8 量级的真余量误判成 0 后清零。Jury 闸里保留补偿只是与 reset 同
       口径。以下档位是它为唯一的回归点。 */

    biquad_filter_t s3;
    biquad_filter_set_empty(&s3);
    s3.den_z[1] = -0.49975f;
    s3.den_z[2] = -0.50025f;   /* 裸 f32 和 == 0.0f；按 f32 输入的真值 +2.98e-8 */
    biquad_filter_reset(&s3, 1.0f);
    CHECK(s3.w[0] > 1e6f,
          "sum3f: reset uses the true residual, not a 0.0f collapse");
    CHECK(s3.w[0] == s3.w[1] && s3.w[1] == s3.w[2],
          "sum3f: reset sets all three state variables");

    /* ── 非有限输入毒化状态（文档化的热路径语义）────────────────────────
       每样本 update 刻意不做 NaN 防护（MCU 热路径上的分支开销）；
       NaN 输入会一路传播，直到 reset 恢复。 ──────────────────────────── */

    biquad_filter_t ns;
    biquad_filter_init(&ns, num, den);
    biquad_filter_reset(&ns, 0.0f);
    y = biquad_filter_update(&ns, 1.0f);
    CHECK(CLOSE(y, 0.2f, 1e-6f), "NaN-poison test: sane state before");
    y = biquad_filter_update(&ns, NAN);
    CHECK(isnan(y), "NaN input → NaN output (no hot-path guard)");
    y = biquad_filter_update(&ns, 1.0f);
    CHECK(isnan(y), "NaN input poisons subsequent samples");
    biquad_filter_reset(&ns, 0.0f);
    y = biquad_filter_update(&ns, 1.0f);
    CHECK(CLOSE(y, 0.2f, 1e-6f), "reset recovers from NaN-poisoned state");

    /* ── 汇总 ─────────────────────────────────────────────────────── */

    if (failures) {
        fprintf(stderr, "%d test(s) FAILED.\n", failures);
        return EXIT_FAILURE;
    }

    printf("All tests passed.\n");
    return EXIT_SUCCESS;
}
