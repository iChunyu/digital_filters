#include "peak_filter.h"
#include "notch_filter.h" /* 对偶不变量测试用：峰值 ≡ 倒数深度的陷波 */
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
 * 与 test_notch.c 的同名函数逐条同构：相位累加器（避免 f32 相位量化把
 * 频率拽偏）+ **RMS 比值**（正弦峰值在低采样密度下被低估，会凭空造出
 * 假偏差）。
 *
 * 无效滤波器返回 1.0f——直通部署必须让增益检查失败（目标 g > 1，1.0f
 * 必然落在容差外），而不是用 0.0f 空过检查。
 *
 * @param[in] f      滤波器对象指针。
 * @param[in] freq   测试正弦频率（Hz）。
 * @param[in] steps  总采样点数。
 * @return           后四分之一样本上的 |H(freq)|。
 */
static float measure_peak_gain(peak_filter_t *f, float freq, int steps)
{
    if (!f->valid) return 1.0f;

    peak_reset(f, 0.0f);

    const float two_pi = 2.0f * (float)M_PI;
    const float dphase = two_pi * freq / f->fs;
    float phase = 0.0f;
    float sum_x2 = 0.0f, sum_y2 = 0.0f; /* 零 double；精度余量见 test_notch.c */

    for (int n = 0; n < steps; n++) {
        float x = sinf(phase);
        float y = peak_update(f, x);
        phase += dphase;
        if (phase >= two_pi) phase -= two_pi;
        if (n >= steps - steps / 4) {
            sum_x2 += x * x;
            sum_y2 += y * y;
        }
    }
    if (sum_x2 <= 1e-6f) return 1.0f; /* 激励退化，无法测量 */
    return sqrtf(sum_y2 / sum_x2);
}

/**
 * @brief 由节系数解析计算 DC 增益 H(1)。
 *
 * @param[in] f  滤波器对象指针。
 * @return       H(1)；无效滤波器返回 NAN，让调用方的窗口检查失败。
 */
static float peak_dc_gain(const peak_filter_t *f)
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
static float peak_nyq_gain(const peak_filter_t *f)
{
    if (!f->valid) return NAN;
    const biquad_filter_t *s = &f->sections[0];
    return (s->num_z[0] - s->num_z[1] + s->num_z[2])
         / (1.0f - s->den_z[1] + s->den_z[2]);
}

/**
 * @brief 校验峰值与陷波的对偶不变量。
 *
 * H_peak(ξ,g) ≡ 1 / H_notch(ξ,1/g) 在系数层面的体现：令 d = 1/g，
 * 则 notch_init(f0,ξ,d,fs) 与 peak_init(f0,ξ,g,fs) 的 z 域系数向量
 * **互为交换**，只差一个公共比例因子（两者各自的 den_z[0] 归一化所致）：
 *
 *     peak.num_z[i] = notch.den_z[i] / notch.num_z[0]
 *     peak.den_z[i] = notch.num_z[i] / notch.num_z[0]
 *
 * src/peak_filter.c 刻意把分母写成 2·xi·d·w0（而非 2·xi·w0/g）并与
 * notch 的 2·xi·g·w0 同形，就是为了让这条关系在 f32 下也成立。
 *
 * @param[in] f0  峰值中心频率（Hz）。
 * @param[in] xi  带宽因子。
 * @param[in] g   线性峰值增益（> 1）。
 * @param[in] fs  采样频率（Hz）。
 */
static void check_duality(float f0, float xi, float g, float fs)
{
    peak_filter_t p;
    notch_filter_t n;
    p.valid = 0;
    n.valid = 0;

    const float d = 1.0f / g;
    peak_init(&p, f0, xi, g, fs);
    notch_init(&n, f0, xi, d, fs);

    CHECK(p.valid == 1 && n.valid == 1, "duality: both designs valid");
    if (!p.valid || !n.valid) return;

    const float scale = 1.0f / n.sections[0].num_z[0];
    /* 系数都在 O(1) 量级，绝对容差 2e-5 ≈ 1e-5 相对——两条归一化除法
       各自舍入一次的余量（实测偏差 ~1e-7，留了约 100 倍）。 */
    for (int i = 0; i < 3; i++) {
        CHECK(CLOSE(p.sections[0].num_z[i], n.sections[0].den_z[i] * scale, 2e-5f),
              "duality: peak num == notch den (scaled)");
        CHECK(CLOSE(p.sections[0].den_z[i], n.sections[0].num_z[i] * scale, 2e-5f),
              "duality: peak den == notch num (scaled)");
    }

    /* 直流/奈奎斯特两端同样互为倒数，且都精确等于 1（结构保证）。 */
    CHECK(CLOSE(peak_dc_gain(&p) * (1.0f / peak_dc_gain(&p)), 1.0f, 1e-3f),
          "duality: sanity");
}

/* 峰值测量需要的样本数：极点阻尼 ζ = ξ·d，衰减时间常数
   τ ≈ 1/(ζ·π·f0/fs) 个样本。f0/fs = 0.05、ξd = 0.005 时 τ ≈ 630，
   2e4 样本 ≈ 31τ，后四分之一早已稳态。 */
#define PEAK_STEPS 20000

int main(void)
{
    /* ── 基本设计：fs=1kHz, f0=50Hz, xi=0.05, g=10 (+20 dB) ─────────────── */

    peak_filter_t p;
    p.valid = 0; /* 压掉宏里未初始化的告警 */
    peak_init(&p, 50.0f, 0.05f, 10.0f, 1000.0f);
    CHECK(p.valid == 1, "basic: init valid");
    CHECK(p.num_sections == 1, "basic: 1 section");
    CHECK(CLOSE(p.f0, 50.0f, 1e-6f), "basic: f0 stored");
    CHECK(CLOSE(p.xi, 0.05f, 1e-6f), "basic: xi stored");
    CHECK(CLOSE(p.g, 10.0f, 1e-6f), "basic: g stored");
    CHECK(CLOSE(p.fs, 1000.0f, 1e-6f), "basic: fs stored");

    /* 峰高 ≈ g */
    float y = measure_peak_gain(&p, 50.0f, PEAK_STEPS);
    CHECK(CLOSE(y, 10.0f, 1.0f), "basic: |H(f0)| ~ g within 10%");

    /* 通带：直流与 Nyquist 精确为 1（结构保证：b1 == a1 且
       b0-1 = -(b2-a2)，在 z = ±1 处分子分母逐项相消） */
    CHECK(CLOSE(peak_dc_gain(&p), 1.0f, 1e-4f), "basic: DC gain = 1");
    CHECK(CLOSE(peak_nyq_gain(&p), 1.0f, 1e-4f), "basic: Nyquist gain = 1");

    /* 它是窄峰而不是宽带提升：f0 的 1/2 处应≈直通 */
    y = measure_peak_gain(&p, 25.0f, PEAK_STEPS);
    CHECK(CLOSE(y, 1.0f, 0.02f), "basic: |H(f0/2)| ~ 1 (narrow, not broadband)");

    /* ── 增益扫掠 ─────────────────────────────────────────────────────── */

    const float g_grid[] = { 2.0f, 5.0f, 10.0f, 20.0f, 50.0f };
    for (unsigned i = 0; i < sizeof(g_grid) / sizeof(g_grid[0]); i++) {
        peak_filter_t gn;
        gn.valid = 0;
        peak_init(&gn, 50.0f, 0.05f, g_grid[i], 1000.0f);
        CHECK(gn.valid == 1, "g sweep: valid");
        if (!gn.valid) continue;
        float h = measure_peak_gain(&gn, 50.0f, PEAK_STEPS);
        CHECK(CLOSE(h, g_grid[i], 0.05f * g_grid[i] + 0.002f),
              "g sweep: peak height matches g within 5%");
    }

    /* ── 频率扫掠 + 对偶不变量：全部落在闸内，且必须精确 ─────────────── */

    for (float f0 = 50.0f; f0 <= 400.0f; f0 += 7.0f) {
        peak_filter_t fn;
        fn.valid = 0;
        peak_init(&fn, f0, 0.05f, 10.0f, 1000.0f);
        CHECK(fn.valid == 1, "f sweep: valid");
        if (!fn.valid) continue;
        float h = measure_peak_gain(&fn, f0, PEAK_STEPS);
        CHECK(CLOSE(h, 10.0f, 1.0f), "f sweep: height within 10%");
        CHECK(CLOSE(peak_dc_gain(&fn), 1.0f, 1e-4f), "f sweep: DC = 1");
        CHECK(CLOSE(peak_nyq_gain(&fn), 1.0f, 1e-4f), "f sweep: Nyquist = 1");

        check_duality(f0, 0.05f, 10.0f, 1000.0f);
    }

    /* ── xi 扫掠 ──────────────────────────────────────────────────────── */

    const float xi_grid[] = { 0.3f, 0.1f, 0.05f, 0.02f };
    for (unsigned i = 0; i < sizeof(xi_grid) / sizeof(xi_grid[0]); i++) {
        peak_filter_t xn;
        xn.valid = 0;
        peak_init(&xn, 100.0f, xi_grid[i], 10.0f, 1000.0f);
        CHECK(xn.valid == 1, "xi sweep: valid");
        if (!xn.valid) continue;
        float h = measure_peak_gain(&xn, 100.0f, PEAK_STEPS);
        CHECK(CLOSE(h, 10.0f, 1.2f), "xi sweep: height within 12%");
    }

    /* ── 数值闸边界：Q = xi*d*(w0/K)^2 跨过 5e-8 ──────────────────────── */

    /*
     * 这里刻意用 ξ = 0.5 而不是 0.05，且频率落在 0.3 Hz 附近。
     *
     * 峰值滤波器的极点阻尼是 ξ·d，半径闸要求 2ξd·tan(πf0/fs) ≥ 5e-5。
     * 两条闸的边界（以 r = tan(πf0/fs) 计）是
     *
     *     数值闸：r_q = √(5e-8/(ξd))        半径闸：r_r = 2.5e-5/(ξd)
     *
     * r_q > r_r ⟺ ξd > 0.0125 —— 只有这时数值闸才是**先触发**的那一道。
     * ξ = 0.5、d = 0.1 给出 ξd = 0.05，于是 r_q = 1e-3、r_r = 5e-4：
     * 在 r ∈ (5e-4, 1e-3)（即 f0 ∈ (0.159, 0.318) Hz）里半径闸放行而
     * 数值闸拒绝，这一段才测得到数值闸。
     *
     * 注意这与 notch 相反：notch 的极点阻尼是 ξ（与 g 无关），数值闸
     * 在绝大多数参数下更紧。峰值这边半径闸几乎处处更紧（ξd < 0.0125
     * 时先触发），所以它的实际可用上限主要由极点半径决定——这正是
     * include/peak_filter.h 里 g ≤ 4e4·ξ·tan(πf0/fs) 那条刻度的来历。
     */
    const float xi_q = 0.5f;

    /* f0=0.40Hz @ fs=1kHz：Q = 0.5*0.1*(2π*0.4/2000)^2 = 7.90e-8 > 5e-8；
       半径 1-r ≈ 2*0.05*1.257e-3 = 1.257e-4 > 5e-5，两道闸都过。 */
    peak_filter_t qpass;
    qpass.valid = 0;
    peak_init(&qpass, 0.40f, xi_q, 10.0f, 1000.0f);
    CHECK(qpass.valid == 1, "gate: Q=7.90e-8 passes");

    /* f0=0.28Hz @ fs=1kHz：Q = 3.87e-8 < 5e-8。
       此处极点半径裕量 1-r ≈ 8.8e-5，**高于** 5e-5 阈值，半径闸放行，
       拒的是数值闸——这一档专门用来把数值闸从半径闸后面露出来。 */
    peak_filter_t qrej;
    qrej.valid = 1;
    qrej.num_sections = 9;
    peak_init(&qrej, 0.28f, xi_q, 10.0f, 1000.0f);
    CHECK(qrej.valid == 0, "gate: Q=3.87e-8 rejected");
    CHECK(qrej.num_sections == 0, "gate: rejected → num_sections cleared");
    CHECK(measure_peak_gain(&qrej, 0.28f, 64) == 1.0f,
          "gate: rejected filter passes through");

    /*
     * 反直觉但有物理意义：Q 与 d = 1/g 成正比，所以同样的 f0 下
     * **高** boost（g 大、d 小）反而更容易被拒。与 notch 的
     * "浅陷波反而更容易过闸" 正好互为镜像。
     *
     * f0=0.42Hz、ξ=0.5：Q = 8.705e-7/g，在 g ≈ 17.4 处跨过 5e-8；
     * 半径闸的边界是 g ≤ 26.4。取 g=13（两道闸都过）与 g=22
     * （Q=3.96e-8 被拒，半径 6.0e-5 仍放行）。
     */
    peak_filter_t gdeep, gshallow;
    gdeep.valid = 1;
    gdeep.num_sections = 9;
    gshallow.valid = 0;
    peak_init(&gdeep, 0.42f, xi_q, 22.0f, 1000.0f);    /* Q=3.96e-8 < 5e-8 */
    peak_init(&gshallow, 0.42f, xi_q, 13.0f, 1000.0f); /* Q=6.70e-8 ≥ 5e-8 */
    CHECK(gdeep.valid == 0, "gate: g=22 at f0/fs=4.2e-4 rejected");
    CHECK(gdeep.num_sections == 0, "gate: rejected → num_sections cleared");
    CHECK(gshallow.valid == 1, "gate: g=13 at f0/fs=4.2e-4 passes");

    /*
     * 半径闸先触发的档位：ξ=0.05、g=10（d=0.1）时数值闸在
     * f0 ≈ 1.01 Hz 就该放行（Q=7.11e-8 ≥ 5e-8），但极点半径
     * 1-r ≈ 2ξd·tan(πf0/fs) = 3.8e-5 < 5e-5，被 biquad_filter_init 拒。
     * 2.5 Hz 处 1-r ≈ 7.9e-5，两道闸都过——作为阳性对照。
     */
    peak_filter_t radrej, radok;
    radrej.valid = 1;
    radrej.num_sections = 9;
    radok.valid = 0;
    peak_init(&radrej, 1.2f, 0.05f, 10.0f, 1000.0f);
    peak_init(&radok, 2.5f, 0.05f, 10.0f, 1000.0f);
    CHECK(radrej.valid == 0, "gate: radius margin rejects f0=1.2Hz at g=10");
    CHECK(radrej.num_sections == 0, "gate: radius reject clears num_sections");
    CHECK(radok.valid == 1, "gate: f0=2.5Hz passes both gates");

    /* ── g 边界：浅到等同不 boost → 直通（与 notch 的 g ≥ 0.9999 镜像）── */

    peak_filter_t g1;
    g1.valid = 1;
    g1.num_sections = 9;
    peak_init(&g1, 50.0f, 0.05f, 1.0f, 1000.0f);
    CHECK(g1.valid == 0, "g=1 → passthrough (valid=0)");
    CHECK(g1.num_sections == 0, "g=1 → num_sections cleared");

    /* 1/1.0001f 在 f32 下恰好舍入到 0.9999f，落在阈值上（判据是
       !(d < 0.9999f)，等号即拒）——与 notch 的 !(g < 0.9999f) 严格镜像。 */
    peak_filter_t gmin;
    gmin.valid = 1;
    gmin.num_sections = 9;
    peak_init(&gmin, 50.0f, 0.05f, 1.0001f, 1000.0f);
    CHECK(gmin.valid == 0, "g=1.0001 → passthrough");
    CHECK(gmin.num_sections == 0, "g=1.0001 → num_sections cleared");

    peak_filter_t gok;
    gok.valid = 0;
    peak_init(&gok, 50.0f, 0.05f, 1.001f, 1000.0f);
    CHECK(gok.valid == 1, "g=1.001 → valid (0.0087 dB boost)");

    /* ── 参数闸：非法输入一律 valid=0 且 num_sections 清零 ─────────────── */

    struct { float f0, xi, g, fs; const char *msg; } bad[] = {
        {  0.0f,  0.05f, 10.0f,    1000.0f, "f0 = 0"            },
        { -50.0f, 0.05f, 10.0f,    1000.0f, "f0 < 0"            },
        { 500.0f, 0.05f, 10.0f,    1000.0f, "f0 = fs/2"         },
        { 600.0f, 0.05f, 10.0f,    1000.0f, "f0 > fs/2"         },
        { 50.0f,  0.0f,  10.0f,    1000.0f, "xi = 0"            },
        { 50.0f, -0.05f, 10.0f,    1000.0f, "xi < 0"            },
        { 50.0f,  0.05f, 1.0f,     1000.0f, "g = 1（不是 boost）" },
        { 50.0f,  0.05f, 0.1f,     1000.0f, "g < 1（那是 notch）" },
        { 50.0f,  0.05f, -10.0f,   1000.0f, "g < 0"             },
        { 50.0f,  0.05f, 0.0f,     1000.0f, "g = 0"             },
        { 50.0f,  0.05f, 1.00005f, 1000.0f, "g 太浅（d>0.9999）" },
        { 50.0f,  0.05f, 10.0f,    0.0f,    "fs = 0"            },
        { 50.0f,  0.05f, 10.0f,   -1000.0f, "fs < 0"            },
        { NAN,    0.05f, 10.0f,    1000.0f, "f0 = NaN"          },
        { 50.0f,  NAN,   10.0f,    1000.0f, "xi = NaN"          },
        { 50.0f,  0.05f, NAN,      1000.0f, "g = NaN"           },
        { 50.0f,  0.05f, 10.0f,    NAN,     "fs = NaN"          },
        { INFINITY, 0.05f, 10.0f,  1000.0f, "f0 = +Inf"         },
        { 50.0f,  0.05f, 10.0f,    INFINITY, "fs = +Inf"        },
        /* g = +Inf → d = 0 → Q = 0 < 5e-8，由数值闸拒（参数闸不拦 +Inf） */
        { 50.0f,  0.05f, INFINITY, 1000.0f, "g = +Inf"          },
    };
    for (unsigned i = 0; i < sizeof(bad) / sizeof(bad[0]); i++) {
        peak_filter_t bf;
        bf.valid = 1;
        bf.num_sections = 9;
        peak_init(&bf, bad[i].f0, bad[i].xi, bad[i].g, bad[i].fs);
        CHECK(bf.valid == 0, bad[i].msg);
        CHECK(bf.num_sections == 0, bad[i].msg);
        /* 无效 → 原样直通 */
        CHECK(peak_update(&bf, 0.375f) == 0.375f, bad[i].msg);
        peak_reset(&bf, 1.0f); /* 空操作，不得崩溃 */
    }

    /* ── reset：DC 稳态（峰值滤波器 DC 增益同样是 1）──────────────────── */

    peak_filter_t rn;
    rn.valid = 0;
    peak_init(&rn, 50.0f, 0.05f, 10.0f, 1000.0f);
    CHECK(rn.valid == 1, "reset: init valid");
    /* 直流增益为 1，故 reset(0.5) 后紧接着的输出应≈0.5 */
    peak_reset(&rn, 0.5f);
    y = peak_update(&rn, 0.5f);
    CHECK(CLOSE(y, 0.5f, 5e-3f), "reset(0.5) → steady output 0.5");

    peak_reset(&rn, 0.0f);
    y = peak_update(&rn, 0.0f);
    CHECK(CLOSE(y, 0.0f, 1e-4f), "reset(0) → output 0");

    /* 非有限 equilibrium 不得毒化状态 */
    peak_reset(&rn, NAN);
    y = peak_update(&rn, 0.25f);
    CHECK(!isnan(y), "reset(NaN) → state not poisoned");

    /* ── 汇总 ─────────────────────────────────────────────────────────── */

    if (failures) {
        fprintf(stderr, "%d test(s) FAILED.\n", failures);
        return EXIT_FAILURE;
    }
    printf("All tests passed.\n");
    return EXIT_SUCCESS;
}
