#include "notch_filter.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

static int failures = 0;

#define CHECK(cond, msg)                                         \
    do {                                                         \
        if (!(cond)) {                                           \
            failures++;                                          \
            fprintf(stderr, "FAIL: %s\n", msg);                 \
        }                                                        \
    } while (0)

#define CLOSE(a, b, eps) (fabsf((a) - (b)) <= (eps))

/**
 * @brief 喂 freq 正弦，测稳态增益 |H(freq)|。
 *
 * 两处都必须小心：
 *
 * 1. **用相位累加器**而不是 sinf(2π·f·n/fs)。后者在长序列 + 高 f0/fs 下
 *    参数涨到千弧度，f32 的相位量化会把实际频率拽偏，而陷波只有 2ξg·f0
 *    宽，偏一点就掉出谷底。累加器的相位漂移在 2e4 样本上约 4e-4 rad。
 *
 * 2. **返回 RMS 比值而不是输出峰值**。正弦在低采样密度下峰值被低估：
 *    fs=1kHz、f=200Hz 只有 5 点/周期，采样峰值只有 0.951；f=500Hz 正好是
 *    Nyquist，sinf(πn) ≡ 1e-7，激励整个退化。用峰值会把这种采样伪影记成
 *    "通带增益 0.951"，凭空造出 5% 的假偏差。RMS 比值 √(Σy²/Σx²) 对采样
 *    密度免疫（两者同比例缩放），直接消掉这个伪影。
 *
 * 无效滤波器返回 1.0f——直通部署必须让衰减检查失败，而不是用 0.0f
 * 空过检查（与 butter/cheby 测试同一约定）。
 *
 * @param[in] f      滤波器对象指针。
 * @param[in] freq   测试正弦频率（Hz）。
 * @param[in] steps  总采样点数。
 * @return           后四分之一样本上的 |H(freq)|。
 */
static float measure_notch_gain(notch_filter_t *f, float freq, int steps)
{
    if (!f->valid) return 1.0f;

    notch_reset(f, 0.0f);

    const float two_pi = 2.0f * (float)M_PI;
    const float dphase = two_pi * freq / f->fs;
    float phase = 0.0f;
    /* float 累加（库的约定是零 double）。这里量级完全够：窗长约 3000 样本、
       Σx² ≈ 1500，f32 累加的相对误差 ~√N·eps ≈ 3e-6，比 0.01 的判据低四个
       数量级；且 x、y 同窗同比例缩放，误差在比值里还会进一步对消。 */
    float sum_x2 = 0.0f, sum_y2 = 0.0f;

    for (int n = 0; n < steps; n++) {
        float x = sinf(phase);
        float y = notch_update(f, x);
        phase += dphase;
        if (phase >= two_pi) phase -= two_pi;
        if (n >= steps - steps / 4) {
            sum_x2 += x * x;
            sum_y2 += y * y;
        }
    }
    if (sum_x2 <= 1e-6f) return 1.0f; /* 激励退化（如 Nyquist），无法测量 */
    return sqrtf(sum_y2 / sum_x2);
}

/**
 * @brief 由节系数解析计算 DC 增益 H(1)。
 *
 * @param[in] f  滤波器对象指针。
 * @return       H(1)；无效滤波器返回 NAN，让调用方的窗口检查失败。
 */
static float notch_dc_gain(const notch_filter_t *f)
{
    if (!f->valid) return NAN;
    const biquad_filter_t *s = &f->sections[0];
    return (s->num_z[0] + s->num_z[1] + s->num_z[2])
         / (1.0f + s->den_z[1] + s->den_z[2]);
}

/**
 * @brief 由节系数解析计算 Nyquist 增益 H(−1)。
 *
 * @param[in] f  滤波器对象指针。
 * @return       H(−1)；无效滤波器返回 NAN。
 */
static float notch_nyq_gain(const notch_filter_t *f)
{
    if (!f->valid) return NAN;
    const biquad_filter_t *s = &f->sections[0];
    return (s->num_z[0] - s->num_z[1] + s->num_z[2])
         / (1.0f - s->den_z[1] + s->den_z[2]);
}

/* 深度测量需要的样本数：陷波零点的衰减时间常数 τ ≈ 1/(2ξg·π·f0/fs)
   个样本。f0/fs = 0.05、ξg = 0.005 时 τ ≈ 640，2e4 样本≈31τ，后四分之一
   早已稳态。 */
#define DEPTH_STEPS 20000

int main(void)
{
    /* ── 基本设计：fs=1kHz, f0=50Hz, xi=0.05, g=0.1 ───────────────────── */

    notch_filter_t n;
    n.valid = 0; /* 压掉宏里未初始化的告警 */
    notch_init(&n, 50.0f, 0.05f, 0.1f, 1000.0f);
    CHECK(n.valid == 1, "basic: init valid");
    CHECK(n.num_sections == 1, "basic: 1 section");
    CHECK(CLOSE(n.f0, 50.0f, 1e-6f), "basic: f0 stored");
    CHECK(CLOSE(n.xi, 0.05f, 1e-6f), "basic: xi stored");
    CHECK(CLOSE(n.g, 0.1f, 1e-6f), "basic: g stored");
    CHECK(CLOSE(n.fs, 1000.0f, 1e-6f), "basic: fs stored");

    /* 陷波深度 ≈ g */
    float y = measure_notch_gain(&n, 50.0f, DEPTH_STEPS);
    CHECK(CLOSE(y, 0.1f, 0.01f), "basic: |H(f0)| ~ g within 10%");

    /* 通带：直流与 Nyquist 精确为 1（结构保证：b1 == a1 且
       b0-1 = -(b2-a2)，在 z = ±1 处分子分母逐项相消） */
    CHECK(CLOSE(notch_dc_gain(&n), 1.0f, 1e-4f), "basic: DC gain = 1");
    CHECK(CLOSE(notch_nyq_gain(&n), 1.0f, 1e-4f), "basic: Nyquist gain = 1");

    /* 它确实是"陷波"而不是宽带衰减：f0 的 1/2 处应≈直通 */
    y = measure_notch_gain(&n, 25.0f, DEPTH_STEPS);
    CHECK(CLOSE(y, 1.0f, 0.02f), "basic: |H(f0/2)| ~ 1 (narrow, not broadband)");

    /* ── 深度扫掠 ─────────────────────────────────────────────────────── */

    const float g_grid[] = { 0.5f, 0.2f, 0.1f, 0.05f, 0.02f, 0.01f };
    for (unsigned i = 0; i < sizeof(g_grid) / sizeof(g_grid[0]); i++) {
        notch_filter_t gn;
        gn.valid = 0;
        notch_init(&gn, 50.0f, 0.05f, g_grid[i], 1000.0f);
        CHECK(gn.valid == 1, "g sweep: valid");
        if (!gn.valid) continue;
        float d = measure_notch_gain(&gn, 50.0f, DEPTH_STEPS);
        CHECK(CLOSE(d, g_grid[i], 0.05f * g_grid[i] + 0.002f),
              "g sweep: depth matches g within 5%");
    }

    /* ── 频率扫掠：全落在闸内，且必须全部精确 —— 这一段直接对着
           "失效在 0.002 Hz 尺度上剧烈震荡" 那个病态来守 ────────────────── */

    for (float f0 = 50.0f; f0 <= 400.0f; f0 += 7.0f) {
        notch_filter_t fn;
        fn.valid = 0;
        notch_init(&fn, f0, 0.05f, 0.1f, 1000.0f);
        CHECK(fn.valid == 1, "f sweep: valid");
        if (!fn.valid) continue;
        float d = measure_notch_gain(&fn, f0, DEPTH_STEPS);
        CHECK(CLOSE(d, 0.1f, 0.01f), "f sweep: depth within 10%");
        CHECK(CLOSE(notch_dc_gain(&fn), 1.0f, 1e-4f), "f sweep: DC = 1");
        CHECK(CLOSE(notch_nyq_gain(&fn), 1.0f, 1e-4f), "f sweep: Nyquist = 1");
    }

    /* ── xi 扫掠 ──────────────────────────────────────────────────────── */

    const float xi_grid[] = { 0.3f, 0.1f, 0.05f, 0.02f, 0.01f };
    for (unsigned i = 0; i < sizeof(xi_grid) / sizeof(xi_grid[0]); i++) {
        notch_filter_t xn;
        xn.valid = 0;
        notch_init(&xn, 100.0f, xi_grid[i], 0.1f, 1000.0f);
        CHECK(xn.valid == 1, "xi sweep: valid");
        if (!xn.valid) continue;
        float d = measure_notch_gain(&xn, 100.0f, DEPTH_STEPS);
        CHECK(CLOSE(d, 0.1f, 0.012f), "xi sweep: depth within 12%");
    }

    /* ── 数值闸边界：Q = xi*g*(w0/K)^2 跨过 5e-8 ──────────────────────── */

    /* f0=1.2Hz @ fs=1kHz：Q = 0.05*0.1*(2π*1.2/2000)^2 = 7.11e-8 > 5e-8 */
    notch_filter_t gpass;
    gpass.valid = 0;
    notch_init(&gpass, 1.2f, 0.05f, 0.1f, 1000.0f);
    CHECK(gpass.valid == 1, "gate: Q=7.11e-8 passes");

    /* f0=0.9Hz @ fs=1kHz：Q = 4.00e-8 < 5e-8。
       注意此处极点半径裕量 (1-r ≈ 2.8e-4) 远离 5e-5 阈值，拒的是数值闸。 */
    notch_filter_t grej;
    grej.valid = 1;
    grej.num_sections = 9;
    notch_init(&grej, 0.9f, 0.05f, 0.1f, 1000.0f);
    CHECK(grej.valid == 0, "gate: Q=4.00e-8 rejected");
    CHECK(grej.num_sections == 0, "gate: rejected → num_sections cleared");
    /* 被拒的滤波器必须直通 */
    CHECK(measure_notch_gain(&grej, 0.9f, 64) == 1.0f,
          "gate: rejected filter passes through");

    /*
     * 反直觉但有物理意义：Q 与 g 成正比，所以同样的 f0 下**浅**陷波
     * （g 大）反而更容易过闸。浅陷波的分子分母差异大、需要相消的位数少，
     * f32 承载得住；深陷波才是难的那一端。
     */
    notch_filter_t gshallow, gdeep;
    gshallow.valid = 0;
    gdeep.valid = 1;
    gdeep.num_sections = 9;
    notch_init(&gshallow, 1.0f, 0.05f, 0.9f, 1000.0f); /* Q=4.44e-7 > 5e-8 */
    notch_init(&gdeep, 1.0f, 0.05f, 0.1f, 1000.0f);    /* Q=4.93e-8 < 5e-8 */
    CHECK(gshallow.valid == 1, "gate: g=0.9 at f0/fs=1e-3 passes");
    CHECK(gdeep.valid == 0, "gate: g=0.1 at f0/fs=1e-3 rejected");
    CHECK(gdeep.num_sections == 0, "gate: rejected → num_sections cleared");

    /* ── g 上限：浅到等同不陷波 → 直通 ────────────────────────────────── */

    notch_filter_t g1;
    g1.valid = 1;
    g1.num_sections = 9;
    notch_init(&g1, 50.0f, 0.05f, 1.0f, 1000.0f);
    CHECK(g1.valid == 0, "g=1 → passthrough (valid=0)");
    CHECK(g1.num_sections == 0, "g=1 → num_sections cleared");

    notch_filter_t gmax;
    gmax.valid = 1;
    gmax.num_sections = 9;
    notch_init(&gmax, 50.0f, 0.05f, 0.9999f, 1000.0f);
    CHECK(gmax.valid == 0, "g=0.9999 → passthrough");

    notch_filter_t gok;
    gok.valid = 0;
    notch_init(&gok, 50.0f, 0.05f, 0.999f, 1000.0f);
    CHECK(gok.valid == 1, "g=0.999 → valid (0.0087 dB notch)");

    /* ── 参数闸：非法输入一律 valid=0 且 num_sections 清零 ─────────────── */

    struct { float f0, xi, g, fs; const char *msg; } bad[] = {
        {  0.0f,  0.05f, 0.1f,   1000.0f, "f0 = 0"          },
        { -50.0f, 0.05f, 0.1f,   1000.0f, "f0 < 0"          },
        { 500.0f, 0.05f, 0.1f,   1000.0f, "f0 = fs/2"       },
        { 600.0f, 0.05f, 0.1f,   1000.0f, "f0 > fs/2"       },
        { 50.0f,  0.0f,  0.1f,   1000.0f, "xi = 0"          },
        { 50.0f, -0.05f, 0.1f,   1000.0f, "xi < 0"          },
        { 50.0f,  0.05f, 0.0f,   1000.0f, "g = 0"           },
        { 50.0f,  0.05f, -0.1f,  1000.0f, "g < 0"           },
        { 50.0f,  0.05f, 1.5f,   1000.0f, "g > 1"           },
        { 50.0f,  0.05f, 0.1f,   0.0f,    "fs = 0"          },
        { 50.0f,  0.05f, 0.1f,  -1000.0f, "fs < 0"          },
        { NAN,    0.05f, 0.1f,   1000.0f, "f0 = NaN"        },
        { 50.0f,  NAN,   0.1f,   1000.0f, "xi = NaN"        },
        { 50.0f,  0.05f, NAN,    1000.0f, "g = NaN"         },
        { 50.0f,  0.05f, 0.1f,   NAN,     "fs = NaN"        },
        { INFINITY, 0.05f, 0.1f, 1000.0f, "f0 = +Inf"       },
        { 50.0f,  0.05f, 0.1f,   INFINITY, "fs = +Inf"      },
    };
    for (unsigned i = 0; i < sizeof(bad) / sizeof(bad[0]); i++) {
        notch_filter_t bf;
        bf.valid = 1;
        bf.num_sections = 9;
        notch_init(&bf, bad[i].f0, bad[i].xi, bad[i].g, bad[i].fs);
        CHECK(bf.valid == 0, bad[i].msg);
        CHECK(bf.num_sections == 0, bad[i].msg);
        /* 无效 → 原样直通 */
        CHECK(notch_update(&bf, 0.375f) == 0.375f, bad[i].msg);
        notch_reset(&bf, 1.0f); /* 空操作，不得崩溃 */
    }

    /* ── reset：DC 稳态 ───────────────────────────────────────────────── */

    notch_filter_t rn;
    rn.valid = 0;
    notch_init(&rn, 50.0f, 0.05f, 0.1f, 1000.0f);
    CHECK(rn.valid == 1, "reset: init valid");
    /* 直流增益为 1，故 reset(0.5) 后紧接着的输出应≈0.5 */
    notch_reset(&rn, 0.5f);
    y = notch_update(&rn, 0.5f);
    CHECK(CLOSE(y, 0.5f, 5e-3f), "reset(0.5) → steady output 0.5");

    notch_reset(&rn, 0.0f);
    y = notch_update(&rn, 0.0f);
    CHECK(CLOSE(y, 0.0f, 1e-4f), "reset(0) → output 0");

    /* 非有限 equilibrium 不得毒化状态 */
    notch_reset(&rn, NAN);
    y = notch_update(&rn, 0.25f);
    CHECK(!isnan(y), "reset(NaN) → state not poisoned");

    /* ── 汇总 ─────────────────────────────────────────────────────────── */

    if (failures) {
        fprintf(stderr, "%d test(s) FAILED.\n", failures);
        return EXIT_FAILURE;
    }
    printf("All tests passed.\n");
    return EXIT_SUCCESS;
}
