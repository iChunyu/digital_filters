/* 原型计算与边缘增益，依据见 docs/cheby_filter.md。 */

#include "cheby_filter.h"
#include <math.h>
#include <stddef.h>

/* Chebyshev I 原型极点：角度同 Butterworth，半径按 μ 的 sinh / cosh 缩放成椭圆。 */
static void cheby1_proto(complex_t *poles, uint8_t n, float epsilon)
{
    float inv_eps = 1.0f / epsilon;
    float mu = logf(inv_eps + sqrtf(inv_eps * inv_eps + 1.0f)) / (float)n;
    float sinh_mu = sinhf(mu);
    float cosh_mu = coshf(mu);

    for (uint8_t k = 0; k < n; k++) {
        float theta = (float)(2 * (k + 1) + n - 1)
                    / (2.0f * (float)n) * (float)M_PI;
        poles[k].re = sinh_mu * cosf(theta);
        poles[k].im = cosh_mu * sinf(theta);
    }
}

/* Chebyshev II 原型：极点取 Chebyshev I 的倒数，零点 z_k = j/sin θ_k。 */
static uint8_t cheby2_proto(complex_t *poles, complex_t *zeros, uint8_t n,
                             float epsilon)
{
    float inv_eps = 1.0f / epsilon;
    float mu = logf(inv_eps + sqrtf(inv_eps * inv_eps + 1.0f)) / (float)n;
    float sinh_mu = sinhf(mu);
    float cosh_mu = coshf(mu);
    uint8_t nz = 0;

    for (uint8_t k = 0; k < n; k++) {
        float theta = (float)(2 * (k + 1) + n - 1)
                    / (2.0f * (float)n) * (float)M_PI;

        float den = sinh_mu * sinh_mu * cosf(theta) * cosf(theta)
                  + cosh_mu * cosh_mu * sinf(theta) * sinf(theta);
        poles[k].re =  sinh_mu * cosf(theta) / den;
        poles[k].im = -cosh_mu * sinf(theta) / den;

        /* 阈值 1e-4 只需吞 θ = π 处 sinf 的求值误差；n ≤ 8 时最近合法 θ 的 |sin θ| ≥ sin(π/16) ≈ 0.195，约 1950× 余量，兼兜软浮点的幻影零点。 */
        float s = sinf(theta);
        if (fabsf(s) > 1e-4f) {
            zeros[nz].re = 0.0f;
            zeros[nz].im = 1.0f / s;
            nz++;
        }
    }
    return nz;
}


/* 奇偶阶与类型的边缘增益不同，见 docs/cheby_filter.md；纹波边缘走 ±0.25 窗口。 */

static float cheby1_edge_gain(uint8_t order, float ripple_db)
{
    return (order % 2 == 0) ? powf(10.0f, -ripple_db / 20.0f) : 1.0f;
}

static float cheby2_edge_gain(uint8_t order, float ripple_db)
{
    return (order % 2 == 0) ? powf(10.0f, -ripple_db / 20.0f) : 0.0f;
}


static uint8_t cheby1_lp_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order, float fc, float fs,
                                      float ripple_db)
{
    if (order == 0 || order > 8 || !(fc > 0.0f) || !(fc < fs * 0.5f)
        || !(ripple_db > 0.0f))
        return 0;

    float epsilon = sqrtf(powf(10.0f, ripple_db / 10.0f) - 1.0f);
    float wc = 2.0f * (float)M_PI * prewarp(fc, fs);

    complex_t poles[8];
    cheby1_proto(poles, order, epsilon);

    float k = 1.0f / zpk_hp_bs_gain(1.0f, NULL, 0, poles, order);
    if (order % 2 == 0) k /= sqrtf(1.0f + epsilon * epsilon);

    uint8_t n = design_filter(sections, max_sections,
                              FILTER_LOWPASS,
                              wc, 0.0f, fs, k,
                              poles, order, NULL, 0);
    if (n == 0) return 0;
    if (!check_cascade_gains(sections, n,
                             cheby1_edge_gain(order, ripple_db), 0.0f)) return 0;
    return n;
}

static uint8_t cheby1_hp_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order, float fc, float fs,
                                      float ripple_db)
{
    if (order == 0 || order > 8 || !(fc > 0.0f) || !(fc < fs * 0.5f)
        || !(ripple_db > 0.0f))
        return 0;

    float epsilon = sqrtf(powf(10.0f, ripple_db / 10.0f) - 1.0f);
    float wc = 2.0f * (float)M_PI * prewarp(fc, fs);

    complex_t poles[8];
    cheby1_proto(poles, order, epsilon);

    float k = 1.0f / zpk_hp_bs_gain(1.0f, NULL, 0, poles, order);
    if (order % 2 == 0) k /= sqrtf(1.0f + epsilon * epsilon);

    uint8_t n = design_filter(sections, max_sections,
                              FILTER_HIGHPASS,
                              wc, 0.0f, fs, k,
                              poles, order, NULL, 0);
    if (n == 0) return 0;
    if (!check_cascade_gains(sections, n, 0.0f,
                             cheby1_edge_gain(order, ripple_db))) return 0;
    return n;
}

static uint8_t cheby1_bp_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order,
                                      float fc1, float fc2, float fs,
                                      float ripple_db)
{
    if (order == 0 || order > 8 || !(fc1 > 0.0f) || !(fc1 < fs * 0.5f)
        || !(ripple_db > 0.0f))
        return 0;
    if (!(fc2 > fc1) || !(fc2 < fs * 0.5f)) return 0;

    float epsilon = sqrtf(powf(10.0f, ripple_db / 10.0f) - 1.0f);
    float wc1 = 2.0f * (float)M_PI * prewarp(fc1, fs);
    float wc2 = 2.0f * (float)M_PI * prewarp(fc2, fs);

    complex_t poles[8];
    cheby1_proto(poles, order, epsilon);

    float k = 1.0f / zpk_hp_bs_gain(1.0f, NULL, 0, poles, order);
    if (order % 2 == 0) k /= sqrtf(1.0f + epsilon * epsilon);

    uint8_t n = design_filter(sections, max_sections,
                              FILTER_BANDPASS,
                              wc1, wc2, fs, k,
                              poles, order, NULL, 0);
    if (n == 0) return 0;
    if (!check_cascade_gains(sections, n, 0.0f, 0.0f)) return 0;
    return n;
}

static uint8_t cheby1_bs_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order,
                                      float fc1, float fc2, float fs,
                                      float ripple_db)
{
    if (order == 0 || order > 8 || !(fc1 > 0.0f) || !(fc1 < fs * 0.5f)
        || !(ripple_db > 0.0f))
        return 0;
    if (!(fc2 > fc1) || !(fc2 < fs * 0.5f)) return 0;

    float epsilon = sqrtf(powf(10.0f, ripple_db / 10.0f) - 1.0f);
    float wc1 = 2.0f * (float)M_PI * prewarp(fc1, fs);
    float wc2 = 2.0f * (float)M_PI * prewarp(fc2, fs);

    complex_t poles[8];
    cheby1_proto(poles, order, epsilon);

    float k = 1.0f / zpk_hp_bs_gain(1.0f, NULL, 0, poles, order);
    if (order % 2 == 0) k /= sqrtf(1.0f + epsilon * epsilon);

    uint8_t n = design_filter(sections, max_sections,
                              FILTER_BANDSTOP,
                              wc1, wc2, fs, k,
                              poles, order, NULL, 0);
    if (n == 0) return 0;
    if (!check_cascade_gains(sections, n,
                             cheby1_edge_gain(order, ripple_db),
                             cheby1_edge_gain(order, ripple_db))) return 0;
    return n;
}

static uint8_t cheby2_lp_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order, float fc, float fs,
                                      float ripple_db)
{
    if (order == 0 || order > 8 || !(fc > 0.0f) || !(fc < fs * 0.5f)
        || !(ripple_db > 0.0f))
        return 0;

    float epsilon = 1.0f / sqrtf(powf(10.0f, ripple_db / 10.0f) - 1.0f);
    float wc = 2.0f * (float)M_PI * prewarp(fc, fs);

    complex_t poles[8], zeros[8];
    uint8_t nz = cheby2_proto(poles, zeros, order, epsilon);

    float k = 1.0f / zpk_hp_bs_gain(1.0f, zeros, nz, poles, order);

    uint8_t n = design_filter(sections, max_sections,
                              FILTER_LOWPASS,
                              wc, 0.0f, fs, k,
                              poles, order, zeros, nz);
    if (n == 0) return 0;
    if (!check_cascade_gains(sections, n, 1.0f,
                             cheby2_edge_gain(order, ripple_db))) return 0;
    return n;
}

static uint8_t cheby2_hp_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order, float fc, float fs,
                                      float ripple_db)
{
    if (order == 0 || order > 8 || !(fc > 0.0f) || !(fc < fs * 0.5f)
        || !(ripple_db > 0.0f))
        return 0;

    float epsilon = 1.0f / sqrtf(powf(10.0f, ripple_db / 10.0f) - 1.0f);
    float wc = 2.0f * (float)M_PI * prewarp(fc, fs);

    complex_t poles[8], zeros[8];
    uint8_t nz = cheby2_proto(poles, zeros, order, epsilon);

    float k = 1.0f / zpk_hp_bs_gain(1.0f, zeros, nz, poles, order);

    uint8_t n = design_filter(sections, max_sections,
                              FILTER_HIGHPASS,
                              wc, 0.0f, fs, k,
                              poles, order, zeros, nz);
    if (n == 0) return 0;
    if (!check_cascade_gains(sections, n,
                             cheby2_edge_gain(order, ripple_db), 1.0f)) return 0;
    return n;
}

static uint8_t cheby2_bp_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order,
                                      float fc1, float fc2, float fs,
                                      float ripple_db)
{
    if (order == 0 || order > 8 || !(fc1 > 0.0f) || !(fc1 < fs * 0.5f)
        || !(ripple_db > 0.0f))
        return 0;
    if (!(fc2 > fc1) || !(fc2 < fs * 0.5f)) return 0;

    float epsilon = 1.0f / sqrtf(powf(10.0f, ripple_db / 10.0f) - 1.0f);
    float wc1 = 2.0f * (float)M_PI * prewarp(fc1, fs);
    float wc2 = 2.0f * (float)M_PI * prewarp(fc2, fs);

    complex_t poles[8], zeros[8];
    uint8_t nz = cheby2_proto(poles, zeros, order, epsilon);

    float k = 1.0f / zpk_hp_bs_gain(1.0f, zeros, nz, poles, order);

    uint8_t n = design_filter(sections, max_sections,
                              FILTER_BANDPASS,
                              wc1, wc2, fs, k,
                              poles, order, zeros, nz);
    if (n == 0) return 0;
    if (!check_cascade_gains(sections, n,
                             cheby2_edge_gain(order, ripple_db),
                             cheby2_edge_gain(order, ripple_db))) return 0;
    return n;
}

static uint8_t cheby2_bs_init(biquad_filter_t *sections,
                                      uint8_t max_sections,
                                      uint8_t order,
                                      float fc1, float fc2, float fs,
                                      float ripple_db)
{
    if (order == 0 || order > 8 || !(fc1 > 0.0f) || !(fc1 < fs * 0.5f)
        || !(ripple_db > 0.0f))
        return 0;
    if (!(fc2 > fc1) || !(fc2 < fs * 0.5f)) return 0;

    float epsilon = 1.0f / sqrtf(powf(10.0f, ripple_db / 10.0f) - 1.0f);
    float wc1 = 2.0f * (float)M_PI * prewarp(fc1, fs);
    float wc2 = 2.0f * (float)M_PI * prewarp(fc2, fs);

    complex_t poles[8], zeros[8];
    uint8_t nz = cheby2_proto(poles, zeros, order, epsilon);

    float k = 1.0f / zpk_hp_bs_gain(1.0f, zeros, nz, poles, order);

    uint8_t n = design_filter(sections, max_sections,
                              FILTER_BANDSTOP,
                              wc1, wc2, fs, k,
                              poles, order, zeros, nz);
    if (n == 0) return 0;
    if (!check_cascade_gains(sections, n, 1.0f, 1.0f)) return 0;
    return n;
}


#define X(ord, ns, ol) \
    void cheby1_lp_##ol##_init(cheby1_lp_##ol##_t *f, float fc, float fs, float ripple_db) { \
        f->type = FILTER_LOWPASS; \
        f->order = ord; \
        f->fc1 = fc; \
        f->fc2 = 0.0f; \
        f->fs = fs; \
        f->ripple_db = ripple_db; \
        f->valid = 0; \
        f->num_sections = 0; \
        uint8_t n = cheby1_lp_init(f->sections, ns, ord, fc, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_CHEBY_LP_ORDER
#undef X

#define X(ord, ns, ol) \
    void cheby1_hp_##ol##_init(cheby1_hp_##ol##_t *f, float fc, float fs, float ripple_db) { \
        f->type = FILTER_HIGHPASS; \
        f->order = ord; \
        f->fc1 = fc; \
        f->fc2 = 0.0f; \
        f->fs = fs; \
        f->ripple_db = ripple_db; \
        f->valid = 0; \
        f->num_sections = 0; \
        uint8_t n = cheby1_hp_init(f->sections, ns, ord, fc, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_CHEBY_LP_ORDER
#undef X

#define X(ord, ns, ol) \
    void cheby1_bp_##ol##_init(cheby1_bp_##ol##_t *f, float fc1, float fc2, float fs, float ripple_db) { \
        f->type = FILTER_BANDPASS; \
        f->order = ord; \
        f->fc1 = fc1; \
        f->fc2 = fc2; \
        f->fs = fs; \
        f->ripple_db = ripple_db; \
        f->valid = 0; \
        f->num_sections = 0; \
        uint8_t n = cheby1_bp_init(f->sections, ns, ord, fc1, fc2, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_CHEBY_BP_ORDER
#undef X

#define X(ord, ns, ol) \
    void cheby1_bs_##ol##_init(cheby1_bs_##ol##_t *f, float fc1, float fc2, float fs, float ripple_db) { \
        f->type = FILTER_BANDSTOP; \
        f->order = ord; \
        f->fc1 = fc1; \
        f->fc2 = fc2; \
        f->fs = fs; \
        f->ripple_db = ripple_db; \
        f->valid = 0; \
        f->num_sections = 0; \
        uint8_t n = cheby1_bs_init(f->sections, ns, ord, fc1, fc2, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_CHEBY_BP_ORDER
#undef X

#define X(ord, ns, ol) \
    void cheby2_lp_##ol##_init(cheby2_lp_##ol##_t *f, float fc, float fs, float ripple_db) { \
        f->type = FILTER_LOWPASS; \
        f->order = ord; \
        f->fc1 = fc; \
        f->fc2 = 0.0f; \
        f->fs = fs; \
        f->ripple_db = ripple_db; \
        f->valid = 0; \
        f->num_sections = 0; \
        uint8_t n = cheby2_lp_init(f->sections, ns, ord, fc, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_CHEBY_LP_ORDER
#undef X

#define X(ord, ns, ol) \
    void cheby2_hp_##ol##_init(cheby2_hp_##ol##_t *f, float fc, float fs, float ripple_db) { \
        f->type = FILTER_HIGHPASS; \
        f->order = ord; \
        f->fc1 = fc; \
        f->fc2 = 0.0f; \
        f->fs = fs; \
        f->ripple_db = ripple_db; \
        f->valid = 0; \
        f->num_sections = 0; \
        uint8_t n = cheby2_hp_init(f->sections, ns, ord, fc, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_CHEBY_LP_ORDER
#undef X

#define X(ord, ns, ol) \
    void cheby2_bp_##ol##_init(cheby2_bp_##ol##_t *f, float fc1, float fc2, float fs, float ripple_db) { \
        f->type = FILTER_BANDPASS; \
        f->order = ord; \
        f->fc1 = fc1; \
        f->fc2 = fc2; \
        f->fs = fs; \
        f->ripple_db = ripple_db; \
        f->valid = 0; \
        f->num_sections = 0; \
        uint8_t n = cheby2_bp_init(f->sections, ns, ord, fc1, fc2, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_CHEBY_BP_ORDER
#undef X

#define X(ord, ns, ol) \
    void cheby2_bs_##ol##_init(cheby2_bs_##ol##_t *f, float fc1, float fc2, float fs, float ripple_db) { \
        f->type = FILTER_BANDSTOP; \
        f->order = ord; \
        f->fc1 = fc1; \
        f->fc2 = fc2; \
        f->fs = fs; \
        f->ripple_db = ripple_db; \
        f->valid = 0; \
        f->num_sections = 0; \
        uint8_t n = cheby2_bs_init(f->sections, ns, ord, fc1, fc2, fs, ripple_db); \
        if (n == 0) return; \
        f->num_sections = n; \
        f->valid = 1; \
    }
FOR_EACH_CHEBY_BP_ORDER
#undef X
