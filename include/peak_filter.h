/* 二阶峰值：极点阻尼 ξ/g，峰值来自极点，|H(jω0)| = g。
   依据见 docs/peak_filter.md。 */

#ifndef PEAK_FILTER_H_
#define PEAK_FILTER_H_

#include <stdint.h>
#include "biquad_filter.h"

#ifdef __cplusplus
extern "C" {
#endif

/* 与 notch 严格互为倒数：H_peak(ξ,g) ≡ 1/H_notch(ξ,1/g)，逐点成立。 */
typedef struct {
    uint8_t  valid;
    uint8_t  num_sections;
    float    f0;
    float    xi;
    float    g;
    float    fs;
    biquad_filter_t sections[1];
} peak_filter_t;

/* 三道闸任一失败 → valid=0 且 num_sections=0，_update 直通。 */
void peak_init(peak_filter_t *f, float f0, float xi, float g, float fs);

static inline float peak_update(peak_filter_t *f, float input)
{
    if (!f->valid) return input;
    return biquad_cascade_update(f->sections, 1, input);
}

static inline void peak_reset(peak_filter_t *f, float equilibrium)
{
    if (!f->valid) return;
    biquad_cascade_reset(f->sections, 1, equilibrium);
}

#ifdef __cplusplus
}
#endif

#endif
