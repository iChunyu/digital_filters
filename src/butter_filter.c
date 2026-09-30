#include "butter_filter.h"
#include <math.h>
#include <stddef.h>

/* ================================================================== */
/*  Butterworth 原型极点预计算表（1–8 阶）                             */
/* ================================================================== */
/*
 * 公式：θ_k = π · (2k + N + 1) / (2N)
 *       p_k = cos(θ_k) − j·sin(θ_k)，k = 0..N−1
 *
 * 所有极点位于左半平面（实部为负）。
 * 索引：butter_proto[order − 1][pole_index]。
 */

static const complex_t butter_proto[8][8] = {
    /* N=1 */ {{-1.0000000000f, -0.0000000000f}},
    /* N=2 */ {{-0.7071067812f, -0.7071067812f},
               {-0.7071067812f,  0.7071067812f}},
    /* N=3 */ {{-0.5000000000f, -0.8660254038f},
               {-1.0000000000f, -0.0000000000f},
               {-0.5000000000f,  0.8660254038f}},
    /* N=4 */ {{-0.3826834324f, -0.9238795325f},
               {-0.9238795325f, -0.3826834324f},
               {-0.9238795325f,  0.3826834324f},
               {-0.3826834324f,  0.9238795325f}},
    /* N=5 */ {{-0.3090169944f, -0.9510565163f},
               {-0.8090169944f, -0.5877852523f},
               {-1.0000000000f, -0.0000000000f},
               {-0.8090169944f,  0.5877852523f},
               {-0.3090169944f,  0.9510565163f}},
    /* N=6 */ {{-0.2588190451f, -0.9659258263f},
               {-0.7071067812f, -0.7071067812f},
               {-0.9659258263f, -0.2588190451f},
               {-0.9659258263f,  0.2588190451f},
               {-0.7071067812f,  0.7071067812f},
               {-0.2588190451f,  0.9659258263f}},
    /* N=7 */ {{-0.2225209340f, -0.9749279122f},
               {-0.6234898019f, -0.7818314825f},
               {-0.9009688679f, -0.4338837391f},
               {-1.0000000000f, -0.0000000000f},
               {-0.9009688679f,  0.4338837391f},
               {-0.6234898019f,  0.7818314825f},
               {-0.2225209340f,  0.9749279122f}},
    /* N=8 */ {{-0.1950903220f, -0.9807852804f},
               {-0.5555702330f, -0.8314696123f},
               {-0.8314696123f, -0.5555702330f},
               {-0.9807852804f, -0.1950903220f},
               {-0.9807852804f,  0.1950903220f},
               {-0.8314696123f,  0.5555702330f},
               {-0.5555702330f,  0.8314696123f},
               {-0.1950903220f,  0.9807852804f}},
};


/* ================================================================== */
/*  各类型设计辅助（static；公开 API 见 include/butter_filter.h）        */
/* ================================================================== */
/*
 * 签名：butter_{lp,hp,bp,bs}_init(sections, max_sections, order, ...)
 * 流程：ROM 原型表 → 共享管线 design_filter → 级联增益校验。
 *
 * 参数越界（order 不在 1..8，或频率不满足 0 < fc < fs/2、fc1 < fc2 < fs/2）或任一
 * 校验失败都返回 0，由 X-macro 宏体负责把 valid / num_sections 清零。
 *
 * 期望的 DC/Nyquist 增益按类型分别是 (1,0) (0,1) (0,0) (1,1)，与该类型的结构零点
 * 位置一一对应，全是精确值——所以 butter 只用 ±0.1 窗口，从不使用纹波档 ±0.25。
 * 标定与回归档位见 test/test_butter.c。
 */

/* 低通：期望 DC 1、Nyquist 0 */
static uint8_t butter_lp_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order, float fc, float fs)
{
    /* 入口闸一律写 !(x > 0) 而不是 x <= 0：后者与 NaN 比较恒假，会把 NaN 放行
       （NaN 会一路传到 prewarp，靠下游 design_filter 的有限性闸兜成 fail-closed，
       防线就只剩一层）。本文件 3 处入口闸同此写法。 */
    if (order == 0 || order > 8 || !(fc > 0.0f) || !(fc < fs * 0.5f))
        return 0;

    float wc = 2.0f * (float)M_PI * prewarp(fc, fs);
    uint8_t n = design_filter(sections, max_sections,
                              FILTER_LOWPASS,
                              wc, 0.0f, fs, 1.0f,
                              butter_proto[order - 1], order, NULL, 0);
    if (n == 0) return 0;
    /* LP：DC 增益 1，Nyquist 增益 0 */
    if (!check_cascade_gains(sections, n, 1.0f, 0.0f)) return 0;
    return n;
}

/* 高通：期望 DC 0、Nyquist 1 */
static uint8_t butter_hp_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order, float fc, float fs)
{
    if (order == 0 || order > 8 || !(fc > 0.0f) || !(fc < fs * 0.5f))
        return 0;

    float wc = 2.0f * (float)M_PI * prewarp(fc, fs);
    uint8_t n = design_filter(sections, max_sections,
                              FILTER_HIGHPASS,
                              wc, 0.0f, fs, 1.0f,
                              butter_proto[order - 1], order, NULL, 0);
    if (n == 0) return 0;
    /* HP：DC 增益 0，Nyquist 增益 1 */
    if (!check_cascade_gains(sections, n, 0.0f, 1.0f)) return 0;
    return n;
}

/* 带通：期望 DC 0、Nyquist 0 */
static uint8_t butter_bp_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order,
                                      float fc1, float fc2, float fs)
{
    if (order == 0 || order > 8 || !(fc1 > 0.0f) || !(fc1 < fs * 0.5f))
        return 0;
    if (!(fc2 > fc1) || !(fc2 < fs * 0.5f)) return 0;

    float wc1 = 2.0f * (float)M_PI * prewarp(fc1, fs);
    float wc2 = 2.0f * (float)M_PI * prewarp(fc2, fs);

    uint8_t n = design_filter(sections, max_sections,
                              FILTER_BANDPASS,
                              wc1, wc2, fs, 1.0f,
                              butter_proto[order - 1], order, NULL, 0);
    if (n == 0) return 0;
    /* BP：DC 与 Nyquist 增益均为 0 */
    if (!check_cascade_gains(sections, n, 0.0f, 0.0f)) return 0;
    return n;
}

/* 带阻：期望 DC 1、Nyquist 1 */
static uint8_t butter_bs_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order,
                                      float fc1, float fc2, float fs)
{
    if (order == 0 || order > 8 || !(fc1 > 0.0f) || !(fc1 < fs * 0.5f))
        return 0;
    if (!(fc2 > fc1) || !(fc2 < fs * 0.5f)) return 0;

    float wc1 = 2.0f * (float)M_PI * prewarp(fc1, fs);
    float wc2 = 2.0f * (float)M_PI * prewarp(fc2, fs);

    uint8_t n = design_filter(sections, max_sections,
                              FILTER_BANDSTOP,
                              wc1, wc2, fs, 1.0f,
                              butter_proto[order - 1], order, NULL, 0);
    if (n == 0) return 0;
    /* BS：DC 与 Nyquist 增益均为 1 */
    if (!check_cascade_gains(sections, n, 1.0f, 1.0f)) return 0;
    return n;
}

/* ================================================================== */
/*  各阶 init 函数（宏生成）                                           */
/* ================================================================== */

/*
 * LP init:  butter_lp_Nth_init(f, fc, fs)
 * HP init:  butter_hp_Nth_init(f, fc, fs)
 * BP init:  butter_bp_Nth_init(f, fc1, fc2, fs)
 * BS init:  butter_bs_Nth_init(f, fc1, fc2, fs)
 *
 * 宏生成的各阶 init 共享上述辅助函数的语义（详见头文件声明）：
 * 各自填写元数据字段，经对应类型的辅助函数校验，并按返回值
 * 设置 num_sections 与 valid。
 */

/* 低通 */
#define X(ord, ns, ol) \
    void butter_lp_##ol##_init(butter_lp_##ol##_t *f, float fc, float fs) { \
        f->type = FILTER_LOWPASS; \
        f->order = ord; \
        f->fc1 = fc; \
        f->fc2 = 0.0f; \
        f->fs = fs; \
        f->valid = 0; \
        f->num_sections = 0; /* valid=0 时绝不能留下这个垃圾值 */ \
        uint8_t n = butter_lp_init(f->sections, ns, ord, fc, fs); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_BUTTER_LP_ORDER
#undef X

/* 高通 */
#define X(ord, ns, ol) \
    void butter_hp_##ol##_init(butter_hp_##ol##_t *f, float fc, float fs) { \
        f->type = FILTER_HIGHPASS; \
        f->order = ord; \
        f->fc1 = fc; \
        f->fc2 = 0.0f; \
        f->fs = fs; \
        f->valid = 0; \
        f->num_sections = 0; /* valid=0 时绝不能留下这个垃圾值 */ \
        uint8_t n = butter_hp_init(f->sections, ns, ord, fc, fs); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_BUTTER_LP_ORDER
#undef X

/* 带通 */
#define X(ord, ns, ol) \
    void butter_bp_##ol##_init(butter_bp_##ol##_t *f, float fc1, float fc2, float fs) { \
        f->type = FILTER_BANDPASS; \
        f->order = ord; \
        f->fc1 = fc1; \
        f->fc2 = fc2; \
        f->fs = fs; \
        f->valid = 0; \
        f->num_sections = 0; /* valid=0 时绝不能留下这个垃圾值 */ \
        uint8_t n = butter_bp_init(f->sections, ns, ord, fc1, fc2, fs); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_BUTTER_BP_ORDER
#undef X

/* 带阻 */
#define X(ord, ns, ol) \
    void butter_bs_##ol##_init(butter_bs_##ol##_t *f, float fc1, float fc2, float fs) { \
        f->type = FILTER_BANDSTOP; \
        f->order = ord; \
        f->fc1 = fc1; \
        f->fc2 = fc2; \
        f->fs = fs; \
        f->valid = 0; \
        f->num_sections = 0; /* valid=0 时绝不能留下这个垃圾值 */ \
        uint8_t n = butter_bs_init(f->sections, ns, ord, fc1, fc2, fs); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_BUTTER_BP_ORDER
#undef X
