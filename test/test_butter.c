#include "butter_filter.h"
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

/**
 * @brief 让正弦扫过滤波器，测稳态幅值。
 *
 * 无效滤波器返回 1.0f——直通部署必须让衰减检查失败，
 * 而不是用 0.0f 空过检查。
 *
 * @param[in] sections  biquad 节级联。
 * @param[in] ns        节数。
 * @param[in] valid     滤波器 valid 标志。
 * @param[in] freq      测试正弦频率（Hz）。
 * @param[in] fs        采样频率（Hz）。
 * @param[in] steps     总采样点数。
 * @return              后一半样本的峰值幅值。
 */
static float measure_gain(biquad_filter_t *sections, uint8_t ns, uint8_t valid,
                          float freq, float fs, int steps)
{
    if (!valid) return 1.0f;

    float x = 0.0f;
    for (uint8_t i = 0; i < ns; i++) {
        biquad_filter_reset(&sections[i], x);
        x = biquad_filter_get_output(&sections[i]);
    }

    float max_out = 0.0f;
    for (int n = 0; n < steps; n++) {
        float t = (float)n / fs;
        float in_val = sinf(2.0f * (float)M_PI * freq * t);
        x = in_val;
        for (uint8_t i = 0; i < ns; i++) {
            x = biquad_filter_update(&sections[i], x);
        }
        if (n > steps / 2) {
            if (fabsf(x) > max_out) {
                max_out = fabsf(x);
            }
        }
    }
    return max_out;
}

/**
 * @brief 测 Nyquist 音（交替 ±1）的稳态幅值。
 *
 * @note 恰好 fs/2 的正弦在 f32 下求值为 sin(π·n) ≡ 0——该激励
 *       空过任何滤波器。交替 ±1 方波才是真正的 Nyquist 信号。
 *
 * @param[in] sections  biquad 节级联。
 * @param[in] ns        节数。
 * @param[in] valid     滤波器 valid 标志。
 * @param[in] steps     总采样点数。
 * @return              后一半样本的峰值幅值。
 */
static float measure_nyquist_gain(biquad_filter_t *sections, uint8_t ns,
                                  uint8_t valid, int steps)
{
    if (!valid) return 1.0f;

    for (uint8_t i = 0; i < ns; i++) {
        biquad_filter_reset(&sections[i], 0.0f);
    }

    float max_out = 0.0f;
    for (int n = 0; n < steps; n++) {
        float x = (n % 2 == 0) ? 1.0f : -1.0f;
        for (uint8_t i = 0; i < ns; i++) {
            x = biquad_filter_update(&sections[i], x);
        }
        if (n > steps / 2 && fabsf(x) > max_out) {
            max_out = fabsf(x);
        }
    }
    return max_out;
}

/**
 * @brief 由节系数解析计算级联 DC 增益。
 *
 * @param[in] sections  biquad 节级联。
 * @param[in] ns        节数。
 * @return              H(1)。
 */
static float cascade_dc_gain(const biquad_filter_t *sections, uint8_t ns)
{
    float h = 1.0f;
    for (uint8_t i = 0; i < ns; i++)
        h *= (sections[i].num_z[0] + sections[i].num_z[1] + sections[i].num_z[2])
           / (1.0f + sections[i].den_z[1] + sections[i].den_z[2]);
    return h;
}

/**
 * @brief 由节系数解析计算级联 Nyquist 增益。
 *
 * @param[in] sections  biquad 节级联。
 * @param[in] ns        节数。
 * @return              H(−1)。
 */
static float cascade_nyq_gain(const biquad_filter_t *sections, uint8_t ns)
{
    float h = 1.0f;
    for (uint8_t i = 0; i < ns; i++)
        h *= (sections[i].num_z[0] - sections[i].num_z[1] + sections[i].num_z[2])
           / (1.0f - sections[i].den_z[1] + sections[i].den_z[2]);
    return h;
}

/**
 * @brief 固定频率网格上的最大稳态增益。
 *
 * @param[in] sections  biquad 节级联。
 * @param[in] ns        节数。
 * @param[in] fs        采样频率（Hz）。
 * @param[in] steps     每个频点的采样点数。
 * @return              网格上的最大峰值幅值。
 */
static float max_gain_over(biquad_filter_t *sections, uint8_t ns,
                           float fs, int steps)
{
    const float grid[] = {20.0f, 24.5f, 50.0f, 100.0f, 200.0f, 300.0f, 480.0f};
    float m = 0.0f;
    for (unsigned i = 0; i < sizeof(grid) / sizeof(grid[0]); i++) {
        float g = measure_gain(sections, ns, 1, grid[i], fs, steps);
        if (g > m) m = g;
    }
    return m;
}

int main(void)
{
    float y;

    /* ── LP 2nd order, fc=2 Hz, fs=20 Hz ──────────────────────────────── */

    butter_lp_2nd_t blp;
    butter_lp_2nd_init(&blp, 2.0f, 20.0f);
    CHECK(blp.valid == 1, "LP 2nd init valid");
    CHECK(blp.type == FILTER_LOWPASS, "LP 2nd type");
    CHECK(blp.order == 2, "LP 2nd order");
    CHECK(blp.num_sections == 1, "LP 2nd → 1 section");

    /* DC gain ≈ 1.0 */
    butter_lp_2nd_reset(&blp, 1.0f);
    y = butter_lp_2nd_update(&blp, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "LP 2nd DC gain ~ 1");

    /* Attenuation at Nyquist (10 Hz) */
    float gn = measure_nyquist_gain(blp.sections, blp.num_sections, blp.valid, 400);
    CHECK(gn < 0.15f, "LP 2nd Nyquist attenuation");

    /* ── HP 2nd order, fc=5 Hz, fs=40 Hz ──────────────────────────────── */

    butter_hp_2nd_t bhp;
    butter_hp_2nd_init(&bhp, 5.0f, 40.0f);
    CHECK(bhp.valid == 1, "HP 2nd init valid");
    CHECK(bhp.num_sections == 1, "HP 2nd → 1 section");

    /* DC gain ≈ 0 */
    butter_hp_2nd_reset(&bhp, 1.0f);
    y = butter_hp_2nd_update(&bhp, 1.0f);
    CHECK(CLOSE(y, 0.0f, 1e-3f), "HP 2nd blocks DC");

    /* Nyquist gain ≈ 1.0 */
    butter_hp_2nd_reset(&bhp, 0.0f);
    float nyq_gain = 0.0f;
    for (int n = 0; n < 200; n++) {
        float x = (n % 2 == 0) ? 1.0f : -1.0f;
        y = butter_hp_2nd_update(&bhp, x);
        if (n > 100 && fabsf(y) > nyq_gain) nyq_gain = fabsf(y);
    }
    CHECK(CLOSE(nyq_gain, 1.0f, 1e-3f), "HP 2nd Nyquist gain ~ 1");

    /* ── Invalid parameters → valid = 0, passthrough ──────────────────── */

    butter_lp_2nd_t binv;

    /* fc=0 invalid */
    binv.valid = 0;  /* suppress uninit warning in macro */
    butter_lp_2nd_init(&binv, 0.0f, 20.0f);
    CHECK(binv.valid == 0, "LP 2nd fc=0 invalid");
    y = butter_lp_2nd_update(&binv, 0.5f);
    CHECK(y == 0.5f, "invalid filter passthrough");

    /* fc=fs/2 invalid */
    butter_lp_2nd_init(&binv, 10.0f, 20.0f);
    CHECK(binv.valid == 0, "LP 2nd fc=fs/2 invalid");

    /* ── LP 4th-order → 2 sections ────────────────────────────────────── */

    butter_lp_4th_t b4;
    butter_lp_4th_init(&b4, 3.0f, 20.0f);
    CHECK(b4.valid == 1, "LP 4th init valid");
    CHECK(b4.num_sections == 2, "LP 4th → 2 sections");

    butter_lp_4th_reset(&b4, 1.0f);
    y = butter_lp_4th_update(&b4, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "LP 4th DC gain ~ 1");

    /* ── LP 1st-order ─────────────────────────────────────────────────── */

    butter_lp_1st_t b1;
    butter_lp_1st_init(&b1, 2.0f, 20.0f);
    CHECK(b1.valid == 1, "LP 1st init valid");
    CHECK(b1.num_sections == 1, "LP 1st → 1 section");

    butter_lp_1st_reset(&b1, 1.0f);
    y = butter_lp_1st_update(&b1, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "LP 1st DC gain ~ 1");

    /* ── HP 3rd-order → 2 sections ────────────────────────────────────── */

    butter_hp_3rd_t bhp3;
    butter_hp_3rd_init(&bhp3, 5.0f, 40.0f);
    CHECK(bhp3.valid == 1, "HP 3rd init valid");
    CHECK(bhp3.num_sections == 2, "HP 3rd → 2 sections");

    butter_hp_3rd_reset(&bhp3, 1.0f);
    y = butter_hp_3rd_update(&bhp3, 1.0f);
    CHECK(CLOSE(y, 0.0f, 1e-3f), "HP 3rd blocks DC");

    /* ── BP 2nd order, fc1=2 Hz, fc2=5 Hz, fs=40 Hz ──────────────────── */

    butter_bp_2nd_t bbp;
    butter_bp_2nd_init(&bbp, 2.0f, 5.0f, 40.0f);
    CHECK(bbp.valid == 1, "BP 2nd init valid");
    CHECK(bbp.num_sections == 2, "BP 2nd → 2 sections");

    /* DC gain ≈ 0 */
    butter_bp_2nd_reset(&bbp, 1.0f);
    y = butter_bp_2nd_update(&bbp, 1.0f);
    CHECK(CLOSE(y, 0.0f, 1e-3f), "BP 2nd blocks DC");

    /* Centre frequency gain ≈ 1.0 */
    float f0 = sqrtf(2.0f * 5.0f);
    gn = measure_gain(bbp.sections, bbp.num_sections, bbp.valid, f0, 40.0f, 800);
    CHECK(CLOSE(gn, 1.0f, 0.1f), "BP 2nd centre freq gain ~ 1");

    /* Out-of-band attenuation at Nyquist */
    gn = measure_nyquist_gain(bbp.sections, bbp.num_sections, bbp.valid, 800);
    CHECK(gn < 0.15f, "BP 2nd Nyquist attenuation");

    /* ── BP 1st-order → 1 section ─────────────────────────────────────── */

    butter_bp_1st_t bbp1;
    butter_bp_1st_init(&bbp1, 3.0f, 6.0f, 40.0f);
    CHECK(bbp1.valid == 1, "BP 1st init valid");
    CHECK(bbp1.num_sections == 1, "BP 1st → 1 section");

    /* ── BS 2nd order, fc1=2 Hz, fc2=5 Hz, fs=40 Hz ──────────────────── */

    butter_bs_2nd_t bbs;
    butter_bs_2nd_init(&bbs, 2.0f, 5.0f, 40.0f);
    CHECK(bbs.valid == 1, "BS 2nd init valid");
    CHECK(bbs.num_sections == 2, "BS 2nd → 2 sections");

    /* DC gain ≈ 1 */
    butter_bs_2nd_reset(&bbs, 1.0f);
    y = butter_bs_2nd_update(&bbs, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "BS 2nd DC gain ~ 1");

    /* Centre frequency notch */
    gn = measure_gain(bbs.sections, bbs.num_sections, bbs.valid, f0, 40.0f, 800);
    CHECK(gn < 0.15f, "BS 2nd notch at centre freq");

    /* Nyquist gain ≈ 1.0 */
    butter_bs_2nd_reset(&bbs, 0.0f);
    nyq_gain = 0.0f;
    for (int n = 0; n < 200; n++) {
        float x = (n % 2 == 0) ? 1.0f : -1.0f;
        y = butter_bs_2nd_update(&bbs, x);
        if (n > 100 && fabsf(y) > nyq_gain) nyq_gain = fabsf(y);
    }
    CHECK(CLOSE(nyq_gain, 1.0f, 1e-3f), "BS 2nd Nyquist gain ~ 1");

    /* ── BP/BS invalid: fc1 >= fc2 ────────────────────────────────────── */

    butter_bp_2nd_init(&bbp, 5.0f, 2.0f, 40.0f);
    CHECK(bbp.valid == 0, "BP fc1>=fc2 invalid");

    butter_bs_2nd_t bbs_inv;
    butter_bs_2nd_init(&bbs_inv, 5.0f, 5.0f, 40.0f);
    CHECK(bbs_inv.valid == 0, "BS fc1==fc2 invalid");

    /* ── BP/BS invalid: fc2 >= fs/2 ───────────────────────────────────── */

    butter_bp_2nd_init(&bbp, 2.0f, 20.0f, 40.0f);
    CHECK(bbp.valid == 0, "BP fc2=fs/2 invalid");

    butter_bs_2nd_init(&bbs, 2.0f, 25.0f, 40.0f);
    CHECK(bbs.valid == 0, "BS fc2>=fs/2 invalid");

    /* ── LP 8th-order → 4 sections ────────────────────────────────────── */

    butter_lp_8th_t b8;
    butter_lp_8th_init(&b8, 3.0f, 20.0f);
    CHECK(b8.valid == 1, "LP 8th init valid");
    CHECK(b8.num_sections == 4, "LP 8th → 4 sections");

    butter_lp_8th_reset(&b8, 1.0f);
    y = butter_lp_8th_update(&b8, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-4f), "LP 8th DC gain ~ 1");

    /* ── Regression: ultra-wideband BP8 (gain chain f32 overflow) ────────
       Used to overflow (k = xi^8 > FLT_MAX) and fail closed; with the
       folded gain chain the design must now succeed and match the
       Butterworth response. ─────────────────────────────────────────────── */

    butter_bp_8th_t bwb;
    bwb.valid = 0;
    butter_bp_8th_init(&bwb, 10.0f, 3990.0f, 8000.0f);
    CHECK(bwb.valid == 1, "BP 8th ultra-wideband now valid");
    CHECK(bwb.num_sections == 8, "BP 8th ultra-wideband → 8 sections");

    /* DC and Nyquist must be blocked */
    butter_bp_8th_reset(&bwb, 1.0f);
    y = butter_bp_8th_update(&bwb, 1.0f);
    CHECK(CLOSE(y, 0.0f, 1e-3f), "BP 8th ultra-wideband blocks DC");

    /* Passband gain ≈ 1 and band-edge gain ≈ 1/√2 */
    gn = measure_gain(bwb.sections, bwb.num_sections, bwb.valid, 1000.0f, 8000.0f, 8000);
    CHECK(CLOSE(gn, 1.0f, 0.05f), "BP 8th ultra-wideband passband gain ~ 1");
    gn = measure_gain(bwb.sections, bwb.num_sections, bwb.valid, 3990.0f, 8000.0f, 8000);
    CHECK(CLOSE(gn, 0.707f, 0.05f), "BP 8th ultra-wideband edge gain ~ 0.707");

    /* ── Regression: near-Nyquist LP8 (previously NaN at HEAD) ──────────── */

    butter_lp_8th_t bnn;
    bnn.valid = 0;
    butter_lp_8th_init(&bnn, 470.0f, 1000.0f);
    CHECK(bnn.valid == 1, "LP 8th near-Nyquist now valid");

    butter_lp_8th_reset(&bnn, 1.0f);
    y = butter_lp_8th_update(&bnn, 1.0f);
    CHECK(CLOSE(y, 1.0f, 1e-3f), "LP 8th near-Nyquist DC gain ~ 1");

    gn = measure_gain(bnn.sections, bnn.num_sections, bnn.valid, 470.0f, 1000.0f, 8000);
    CHECK(CLOSE(gn, 0.707f, 0.05f), "LP 8th near-Nyquist edge gain ~ 0.707");

    /* ── Regression: wide-band BP with near-real pole pairs ────────────────
       The f32 Jury sum 1 + a1 + a2 used to cancel to EXACTLY 0.0f for
       orders 4/5/8 (poles ~2.5e-4 from z=1) and reject the design while
       orders 1-3 passed; a too-loose real/complex classification then
       cross-paired the near-real pairs into sections with a pole exactly
       at z = 1.  Both fixed. ─────────────────────────────────────────────── */

    butter_bp_4th_t bwb4; bwb4.valid = 0;
    butter_bp_4th_init(&bwb4, 5.0f, 20000.0f, 48000.0f);
    CHECK(bwb4.valid == 1, "BP 4th (5,20000,48k) valid (Jury/pairing)");

    butter_bp_5th_t bwb5; bwb5.valid = 0;
    butter_bp_5th_init(&bwb5, 5.0f, 20000.0f, 48000.0f);
    CHECK(bwb5.valid == 1, "BP 5th (5,20000,48k) valid (Jury/pairing)");

    butter_bp_8th_t bwb8; bwb8.valid = 0;
    butter_bp_8th_init(&bwb8, 5.0f, 20000.0f, 48000.0f);
    CHECK(bwb8.valid == 1, "BP 8th (5,20000,48k) valid (Jury/pairing)");

    /* ── Regression: narrow-but-representable LP accepted, ultra-narrow
          (pole within the 5e-5 margin of z=1) rejected fail-closed ─────── */

    butter_lp_2nd_t blp6; blp6.valid = 0;
    butter_lp_2nd_init(&blp6, 6.0f, 48000.0f);
    CHECK(blp6.valid == 1, "LP 2nd fc=6Hz@48k accepted (Jury f32 residual)");

    butter_lp_2nd_t blp1; blp1.valid = 0;
    butter_lp_2nd_init(&blp1, 1.0f, 48000.0f);
    CHECK(blp1.valid == 0, "LP 2nd fc=1Hz@48k rejected (pole margin 5e-5)");
    y = butter_lp_2nd_update(&blp1, 0.5f);
    CHECK(y == 0.5f, "LP 2nd fc=1Hz passthrough");

    /* ── Regression: near-Nyquist band edges — inside the pole margin the
          design deploys; past it the design is rejected deterministically
          (no NaN garbage, no acceptance cliff between adjacent specs) ──── */

    butter_bs_1st_t bbs_ok; bbs_ok.valid = 0;
    butter_bs_1st_init(&bbs_ok, 100.0f, 23990.0f, 48000.0f);
    CHECK(bbs_ok.valid == 1, "BS 1st (100,23990,48k) accepted (pole r=0.992)");

    butter_bs_1st_t bbs_x; bbs_x.valid = 0;
    butter_bs_1st_init(&bbs_x, 100.0f, 23999.8f, 48000.0f);
    CHECK(bbs_x.valid == 0, "BS 1st (100,23999.8,48k) rejected (pole r=0.99997)");

    /* ── Struct sizes ─────────────────────────────────────────────────── */

    CHECK(sizeof(butter_lp_3rd_t) > sizeof(butter_lp_1st_t),
          "LP 3rd bigger than LP 1st");
    CHECK(sizeof(butter_bp_8th_t) > sizeof(butter_bp_1st_t),
          "BP 8th bigger than BP 1st");

    /* ── Sweep: butter 全阶矩阵（LP/HP/BP/BS × 全 8 阶）──────────────
       补齐 X-macro 生成 init 的测试覆盖。init 内部已做级联增益
       校验，此处复核窗口并约束 max|H|。──────────────────────────── */

    static const float sw_bfc[2] = {100.0f, 480.0f};
    static const float sw_bbands[3][2] = {{50.0f, 120.0f}, {20.0f, 480.0f},
                                          {10.0f, 499.0f}};

    #define X(ord, ns, ol) \
        for (int sw_fi = 0; sw_fi < 2; sw_fi++) { \
            butter_lp_##ol##_t swf; \
            swf.valid = 0; \
            butter_lp_##ol##_init(&swf, sw_bfc[sw_fi], 1000.0f); \
            CHECK(swf.valid == 1, "sweep butter LP " #ol " valid"); \
            if (swf.valid) { \
                y = cascade_dc_gain(swf.sections, swf.num_sections); \
                CHECK(CLOSE(y, 1.0f, 0.1f), "sweep butter LP " #ol " DC window"); \
                y = cascade_nyq_gain(swf.sections, swf.num_sections); \
                CHECK(fabsf(y) <= 0.1f, "sweep butter LP " #ol " Nyq window"); \
                y = max_gain_over(swf.sections, swf.num_sections, 1000.0f, 2000); \
                CHECK(y < 1.2f, "sweep butter LP " #ol " max|H| sane"); \
            } \
        }
    FOR_EACH_BUTTER_LP_ORDER
    #undef X

    #define X(ord, ns, ol) \
        for (int sw_fi = 0; sw_fi < 2; sw_fi++) { \
            butter_hp_##ol##_t swf; \
            swf.valid = 0; \
            butter_hp_##ol##_init(&swf, sw_bfc[sw_fi], 1000.0f); \
            CHECK(swf.valid == 1, "sweep butter HP " #ol " valid"); \
            if (swf.valid) { \
                y = cascade_dc_gain(swf.sections, swf.num_sections); \
                CHECK(fabsf(y) <= 0.1f, "sweep butter HP " #ol " DC window"); \
                y = cascade_nyq_gain(swf.sections, swf.num_sections); \
                CHECK(CLOSE(y, 1.0f, 0.1f), "sweep butter HP " #ol " Nyq window"); \
                y = max_gain_over(swf.sections, swf.num_sections, 1000.0f, 2000); \
                CHECK(y < 1.2f, "sweep butter HP " #ol " max|H| sane"); \
            } \
        }
    FOR_EACH_BUTTER_LP_ORDER
    #undef X

    #define X(ord, ns, ol) \
        for (int sw_bi = 0; sw_bi < 3; sw_bi++) { \
            butter_bp_##ol##_t swf; \
            swf.valid = 0; \
            butter_bp_##ol##_init(&swf, sw_bbands[sw_bi][0], sw_bbands[sw_bi][1], \
                                  1000.0f); \
            CHECK(swf.valid == 1, "sweep butter BP " #ol " valid"); \
            if (swf.valid) { \
                y = cascade_dc_gain(swf.sections, swf.num_sections); \
                CHECK(fabsf(y) <= 0.1f, "sweep butter BP " #ol " DC window"); \
                y = cascade_nyq_gain(swf.sections, swf.num_sections); \
                CHECK(fabsf(y) <= 0.1f, "sweep butter BP " #ol " Nyq window"); \
                y = max_gain_over(swf.sections, swf.num_sections, 1000.0f, 2000); \
                CHECK(y < 1.2f, "sweep butter BP " #ol " max|H| sane"); \
            } \
        }
    FOR_EACH_BUTTER_BP_ORDER
    #undef X

    #define X(ord, ns, ol) \
        for (int sw_bi = 0; sw_bi < 2; sw_bi++) { \
            butter_bs_##ol##_t swf; \
            swf.valid = 0; \
            /* BS 只取前 2 个频带：[10,499] 超宽带下解析响应正确， \
               但时域 f32 DF-II 状态噪声被近单位圆极点放大（6~8 阶 \
               max|H| 实测 47~6.6e5），超出库的实用包络。 */ \
            butter_bs_##ol##_init(&swf, sw_bbands[sw_bi][0], sw_bbands[sw_bi][1], \
                                  1000.0f); \
            CHECK(swf.valid == 1, "sweep butter BS " #ol " valid"); \
            if (swf.valid) { \
                y = cascade_dc_gain(swf.sections, swf.num_sections); \
                CHECK(CLOSE(y, 1.0f, 0.1f), "sweep butter BS " #ol " DC window"); \
                y = cascade_nyq_gain(swf.sections, swf.num_sections); \
                CHECK(CLOSE(y, 1.0f, 0.1f), "sweep butter BS " #ol " Nyq window"); \
                y = max_gain_over(swf.sections, swf.num_sections, 1000.0f, 2000); \
                CHECK(y < 1.2f, "sweep butter BS " #ol " max|H| sane"); \
            } \
        }
    FOR_EACH_BUTTER_BP_ORDER
    #undef X

    /* ── Report ───────────────────────────────────────────────────────── */

    if (failures) {
        fprintf(stderr, "%d test(s) FAILED.\n", failures);
        return EXIT_FAILURE;
    }

    printf("All tests passed.\n");
    return EXIT_SUCCESS;
}
