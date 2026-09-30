/* 直接 II 型二阶 IIR 节：3 个状态变量 w[]。设计与闸门依据见 docs/biquad_filter.md。 */

#ifndef BIQUAD_FILTER_H_
#define BIQUAD_FILTER_H_

#include <stdint.h>

/* 有限性闸依赖 isfinite 不被折叠成恒真；-ffast-math 会静默废掉全部闸门。本头是六个公开头的唯一公共顶点，守卫放这一处即覆盖全部 TU。 */
#if defined(__FAST_MATH__)
#error "本库要求 IEEE-754 严格语义（fail-closed 依赖 isfinite）：请移除 -ffast-math"
#endif
#if defined(__FINITE_MATH_ONLY__) && (__FINITE_MATH_ONLY__ > 0)
#error "本库要求 IEEE-754 严格语义（fail-closed 依赖 isfinite）：请移除 -ffinite-math-only"
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    float num_z[3];
    float den_z[3]; /* 分母 1.0, a1, a2（已归一化） */
    float w[3];     /* 状态 w[n], w[n-1], w[n-2] */
} biquad_filter_t;

/* 输出的 den_z 未归一化（den_z[0] 未必为 1）。 */
void biquad_c2d_bilinear(float num_z[3], float den_z[3], const float num_s[3],
                        const float den_s[3], float fs);

void biquad_filter_set_empty(biquad_filter_t *filter);

uint8_t biquad_filter_init(biquad_filter_t *filter, const float num_z[3],
                           const float den_z[3]);

static inline float biquad_filter_update(biquad_filter_t *filter, float input)
{
    const float w1 = filter->w[0];
    const float w2 = filter->w[1];

    const float a1 = filter->den_z[1];
    const float a2 = filter->den_z[2];
    const float b0 = filter->num_z[0];
    const float b1 = filter->num_z[1];
    const float b2 = filter->num_z[2];

    const float w0 = input - a1 * w1 - a2 * w2;

    filter->w[2] = w2;
    filter->w[1] = w1;
    filter->w[0] = w0;

    return b0 * w0 + b1 * w1 + b2 * w2;
}

float biquad_filter_get_output(const biquad_filter_t *filter);

float biquad_filter_get_input(const biquad_filter_t *filter);

/* 公开 API，系数不必先过 init：分母为 0/非有限或 w_ss 溢出时清零状态，永不产生 inf/NaN。 */
void biquad_filter_reset(biquad_filter_t *filter, float equilibrium);

static inline float biquad_cascade_update(biquad_filter_t *sections,
                                          uint8_t num_sections, float input)
{
    float x = input;
    for (uint8_t i = 0; i < num_sections; i++) {
        x = biquad_filter_update(&sections[i], x);
    }
    return x;
}

static inline void biquad_cascade_reset(biquad_filter_t *sections,
                                        uint8_t num_sections, float equilibrium)
{
    float x = equilibrium;
    for (uint8_t i = 0; i < num_sections; i++) {
        biquad_filter_reset(&sections[i], x);
        x = biquad_filter_get_output(&sections[i]);
    }
}


#ifdef __cplusplus
}
#endif

#endif
