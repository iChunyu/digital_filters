/* 二阶陷波：H(s) 零点阻尼 ξg、极点阻尼 ξ，|H(jω0)| = g，DC/Nyquist 增益 1。
   依据见 docs/notch_filter.md。 */

#ifndef NOTCH_FILTER_H_
#define NOTCH_FILTER_H_

#include <stdint.h>
#include "biquad_filter.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint8_t  valid;
    uint8_t  num_sections;
    float    f0;
    float    xi;
    float    g;
    float    fs;
    biquad_filter_t sections[1];
} notch_filter_t;

/* 三道闸（参数 / f32 范围 / 逐级）任一失败 → valid=0 且 num_sections=0，_update 直通。 */
void notch_init(notch_filter_t *f, float f0, float xi, float g, float fs);

static inline float notch_update(notch_filter_t *f, float input)
{
    if (!f->valid) return input;
    return biquad_cascade_update(f->sections, 1, input);
}

static inline void notch_reset(notch_filter_t *f, float equilibrium)
{
    if (!f->valid) return;
    biquad_cascade_reset(f->sections, 1, equilibrium);
}

#ifdef __cplusplus
}
#endif

#endif
