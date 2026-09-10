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

    /* ── 稳定性裕量：贴近单位圆的极点被拒 ───────────────────────────── */

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
    biquad_filter_init(&f, num, den);

    float x = 0.3f;
    biquad_filter_update(&f, x);
    float out = biquad_filter_get_output(&f);
    float in  = biquad_filter_get_input(&f);
    CHECK(CLOSE(in, x, 1e-6f), "get_input reconstructs x[n]");
    CHECK(CLOSE(out, 0.2f * x, 1e-6f), "get_output matches y[n] = b0·x");

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
