/* 原型设计管线：频率变换、zpk2sos、增益折叠与增益校验。依据见 docs/filter_utils.md。 */

#ifndef FILTER_UTILS_H_
#define FILTER_UTILS_H_

#include <math.h>
#include <stdint.h>
#include "biquad_filter.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    FILTER_LOWPASS,
    FILTER_HIGHPASS,
    FILTER_BANDPASS,
    FILTER_BANDSTOP
} filter_type_e;

typedef struct {
    float re;
    float im;
} complex_t;

/* 区间 (0, fs/2) 之外返回 NaN（双线性映射在 Nyquist 以外无定义），由下游有限性闸拦下。 */
float prewarp(float fd, float fs);

void analog_lp_transform(complex_t *poles, uint8_t np,
                         complex_t *zeros, uint8_t nz, float wc);

void analog_hp_transform(complex_t *poles, uint8_t np,
                         complex_t *zeros, uint8_t nz, float wc);

void bilinear_transform(complex_t *zp, uint8_t n, float fs);

void analog_bp_transform(complex_t *poles, uint8_t *np,
                         complex_t *zeros, uint8_t *nz,
                         float w0, float xi);

void analog_bs_transform(complex_t *poles, uint8_t *np,
                         complex_t *zeros, uint8_t *nz,
                         float w0, float xi);

float zpk_hp_bs_gain(float k, const complex_t *z, uint8_t nz,
                     const complex_t *p, uint8_t np);

float bilinear_zpk_gain(float k, const complex_t *z, uint8_t nz,
                         const complex_t *p, uint8_t np, float K);

/* s^degree 必须与 ∏(K−p) 交错折叠：单独算 s^degree 会在近 Nyquist 的 wc^8 ≈ FLT_MAX 处先溢出。 */
float bilinear_zpk_gain_scaled(float k, float s, uint8_t degree,
                               const complex_t *z, uint8_t nz,
                               const complex_t *p, uint8_t np, float K);

/* fail-closed：全部零极点被认领 + 每节根与输入多重集一致，任一失败返回 0 → 部署直通。 */
uint8_t zpk2sos(const complex_t *zeros, const complex_t *poles, uint8_t n,
                float (*sos)[6], float k);

/* np==0 || np>16 || nz>np → 0；BP/BS 把每个零极点一分为二原地写回同一 16 元素数组，故这两种类型原型上限 8 阶；proto_zeros==NULL 仅当 nz==0 合法。 */
uint8_t design_filter(biquad_filter_t *sections, uint8_t max_sections,
                      uint8_t type,
                      float wc1, float wc2, float fs,
                      float k,
                      const complex_t *proto_poles, uint8_t np,
                      const complex_t *proto_zeros, uint8_t nz);

/* num_sections==0 返回 0：空级联两端增益都是 1.0，否则带阻的 (1,1) 期望会被空洞通过。 */
uint8_t check_cascade_gains(const biquad_filter_t *sections, uint8_t num_sections,
                            float dc_exp, float ny_exp);

#ifdef __cplusplus
}
#endif

#endif
