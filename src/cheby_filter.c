#include "cheby_filter.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ================================================================== */
/*  Chebyshev Type I                                                  */
/* ================================================================== */

static void cheby1_prototype(complex_t *poles, uint8_t n, float epsilon)
{
    float mu = asinhf(1.0f / epsilon) / n;
    float sinh_mu = sinhf(mu);
    float cosh_mu = coshf(mu);

    for (uint8_t k = 0; k < n; k++) {
        float theta = (float)(2 * (k + 1) + n - 1) / (2.0f * n) * (float)M_PI;
        poles[k].re = sinh_mu * cosf(theta);
        poles[k].im = cosh_mu * sinf(theta);
    }
}

void cheby1_init(cheby1_t *c, uint8_t type, uint8_t order,
                 float fc1, float fc2, float fs, float ripple_db)
{
    c->valid     = 0;
    c->type      = type;
    c->order     = order;
    c->fc1       = fc1;
    c->fc2       = fc2;
    c->fs        = fs;
    c->ripple_db = ripple_db;
    c->sos       = NULL;

    if (order == 0 || fc1 <= 0.0f || fc1 >= fs * 0.5f || ripple_db <= 0.0f)
        return;
    if (type == FILTER_BANDPASS || type == FILTER_BANDSTOP) {
        if (fc2 <= fc1 || fc2 >= fs * 0.5f) return;
    }

    float epsilon = sqrtf(powf(10.0f, ripple_db / 10.0f) - 1.0f);
    float wc1 = 2.0f * (float)M_PI * prewarp(fc1, fs);
    float w0_analog, xi;
    if (type == FILTER_BANDPASS || type == FILTER_BANDSTOP) {
        float wc2 = 2.0f * (float)M_PI * prewarp(fc2, fs);
        w0_analog = sqrtf(wc1 * wc2);
        xi        = wc2 - wc1;
    } else {
        w0_analog = wc1;
        xi        = 0.0f;
    }

    /* Chebyshev I prototype — N poles, no zeros */
    complex_t poles[64];
    complex_t zeros[64];
    uint8_t np = order;
    uint8_t nz = 0;

    cheby1_prototype(poles, order, epsilon);

    /* Prototype gain: k = ∏(−p), divided by √(1+ε²) for even N. */
    float k = 1.0f / zpk_hp_bs_gain(1.0f, NULL, 0, poles, order);
    if (order % 2 == 0) {
        k /= sqrtf(1.0f + epsilon * epsilon);
    }

    /* Analog frequency transform */
    switch (type) {
    case FILTER_LOWPASS:
        analog_lp_transform(poles, np, zeros, nz, w0_analog);
        k = zpk_lp_gain(k, w0_analog, order);
        break;
    case FILTER_HIGHPASS:
        k = zpk_hp_bs_gain(k, zeros, nz, poles, np);
        analog_hp_transform(poles, np, zeros, nz, w0_analog);
        for (uint8_t i = nz; i < np; i++) {
            zeros[i].re = 0.0f;
            zeros[i].im = 0.0f;
        }
        nz = np;
        break;
    case FILTER_BANDPASS:
        analog_bp_transform(poles, &np, zeros, &nz, w0_analog, xi);
        k = zpk_bp_gain(k, xi, order);
        break;
    case FILTER_BANDSTOP:
        k = zpk_hp_bs_gain(k, zeros, nz, poles, np);
        analog_bs_transform(poles, &np, zeros, &nz, w0_analog, xi);
        break;
    default:
        return;
    }

    /* Bilinear gain (on s-domain zp, before the substitution). */
    {
        complex_t z_analog[64], p_analog[64];
        memcpy(z_analog, zeros, (size_t)nz * sizeof(complex_t));
        memcpy(p_analog, poles, (size_t)np * sizeof(complex_t));
        k = bilinear_zpk_gain(k, z_analog, nz, p_analog, np, 2.0f * fs);
    }

    /* Bilinear transform */
    bilinear_transform(poles, np, fs);
    bilinear_transform(zeros, nz, fs);

    /* Zero-pad: prototype zeros at s=∞ → z = -1 (not for BS). */
    if (type != FILTER_BANDSTOP) {
        for (uint8_t i = nz; i < np; i++) {
            zeros[i].re = -1.0f;
            zeros[i].im =  0.0f;
        }
        nz = np;
    }

    /* ZPK → SOS → deploy */
    uint8_t ns = (np + 1) / 2;
    float (*sos)[6] = (float (*)[6])malloc((size_t)ns * 6 * sizeof(float));
    if (sos == NULL) return;

    uint8_t n_sections = zpk2sos(zeros, poles, np, sos, k);

    c->sos = (sos_filter_t *)malloc(sizeof(sos_filter_t));
    if (c->sos == NULL) { free(sos); return; }

    sos_filter_init(c->sos, n_sections);
    for (uint8_t i = 0; i < n_sections; i++) {
        float num[3] = {sos[i][0], sos[i][1], sos[i][2]};
        float den[3] = {sos[i][3], sos[i][4], sos[i][5]};
        sos_filter_set_section(c->sos, i, num, den);
    }

    free(sos);
    c->valid = 1;
}

void cheby1_destroy(cheby1_t *c)
{
    if (c->sos != NULL) {
        sos_filter_destroy(c->sos);
        free(c->sos);
        c->sos = NULL;
    }
    c->valid = 0;
}

float cheby1_update(cheby1_t *c, float input)
{
    if (!c->valid || c->sos == NULL) return input;
    return sos_filter_update(c->sos, input);
}

void cheby1_reset(cheby1_t *c, float equilibrium)
{
    if (!c->valid || c->sos == NULL) return;
    sos_filter_reset(c->sos, equilibrium);
}

/* ================================================================== */
/*  Chebyshev Type II                                                 */
/* ================================================================== */

/**
 * @return Number of finite zeros placed (N for even order, N-1 for odd).
 */
static uint8_t cheby2_prototype(complex_t *poles, complex_t *zeros, uint8_t n,
                                float epsilon)
{
    float mu = asinhf(1.0f / epsilon) / n;
    float sinh_mu = sinhf(mu);
    float cosh_mu = coshf(mu);
    uint8_t nz = 0;

    for (uint8_t k = 0; k < n; k++) {
        float theta = (float)(2 * (k + 1) + n - 1) / (2.0f * n) * (float)M_PI;

        /* Pole: inverse of Chebyshev I pole */
        float den = sinh_mu * sinh_mu * cosf(theta) * cosf(theta)
                  + cosh_mu * cosh_mu * sinf(theta) * sinf(theta);
        poles[k].re =  sinh_mu * cosf(theta) / den;
        poles[k].im = -cosh_mu * sinf(theta) / den;

        /* Zero: j / sin(theta).  When sin(theta) ≈ 0 (N odd) this zero is
           at infinity — skip it; the caller will add a z=-1 zero later. */
        float s = sinf(theta);
        if (fabsf(s) > 1e-7f) {
            zeros[nz].re = 0.0f;
            zeros[nz].im = 1.0f / s;
            nz++;
        }
    }
    return nz;
}

void cheby2_init(cheby2_t *c, uint8_t type, uint8_t order,
                 float fc1, float fc2, float fs, float ripple_db)
{
    c->valid     = 0;
    c->type      = type;
    c->order     = order;
    c->fc1       = fc1;
    c->fc2       = fc2;
    c->fs        = fs;
    c->ripple_db = ripple_db;
    c->sos       = NULL;

    if (order == 0 || fc1 <= 0.0f || fc1 >= fs * 0.5f || ripple_db <= 0.0f)
        return;
    if (type == FILTER_BANDPASS || type == FILTER_BANDSTOP) {
        if (fc2 <= fc1 || fc2 >= fs * 0.5f) return;
    }

    /* ε for Chebyshev II: stopband gain = ε / sqrt(1+ε²) = 10^(-dB/20) */
    float epsilon = 1.0f / sqrtf(powf(10.0f, ripple_db / 10.0f) - 1.0f);
    float wc1 = 2.0f * (float)M_PI * prewarp(fc1, fs);
    float w0_analog, xi;
    if (type == FILTER_BANDPASS || type == FILTER_BANDSTOP) {
        float wc2 = 2.0f * (float)M_PI * prewarp(fc2, fs);
        w0_analog = sqrtf(wc1 * wc2);
        xi        = wc2 - wc1;
    } else {
        w0_analog = wc1;
        xi        = 0.0f;
    }

    /* Chebyshev II prototype — N poles, N (even) or N-1 (odd) finite zeros */
    complex_t poles[64];
    complex_t zeros[64];
    uint8_t np = order;
    uint8_t nz = cheby2_prototype(poles, zeros, order, epsilon);
    uint8_t degree = np - nz; /* 0 for even N, 1 for odd N */

    /* Prototype gain: k = ∏(−p) / ∏(−z) = 1 / (∏(−z)/∏(−p)). */
    float k = 1.0f / zpk_hp_bs_gain(1.0f, zeros, nz, poles, np);

    /* Analog frequency transform */
    switch (type) {
    case FILTER_LOWPASS:
        analog_lp_transform(poles, np, zeros, nz, w0_analog);
        k = zpk_lp_gain(k, w0_analog, degree);
        break;
    case FILTER_HIGHPASS:
        k = zpk_hp_bs_gain(k, zeros, nz, poles, np);
        analog_hp_transform(poles, np, zeros, nz, w0_analog);
        for (uint8_t i = nz; i < np; i++) {
            zeros[i].re = 0.0f;
            zeros[i].im = 0.0f;
        }
        nz = np;
        break;
    case FILTER_BANDPASS:
        analog_bp_transform(poles, &np, zeros, &nz, w0_analog, xi);
        k = zpk_bp_gain(k, xi, degree);
        break;
    case FILTER_BANDSTOP:
        k = zpk_hp_bs_gain(k, zeros, nz, poles, np);
        analog_bs_transform(poles, &np, zeros, &nz, w0_analog, xi);
        break;
    default:
        return;
    }

    /* Bilinear gain (on s-domain zp, before the substitution). */
    {
        complex_t z_analog[64], p_analog[64];
        memcpy(z_analog, zeros, (size_t)nz * sizeof(complex_t));
        memcpy(p_analog, poles, (size_t)np * sizeof(complex_t));
        k = bilinear_zpk_gain(k, z_analog, nz, p_analog, np, 2.0f * fs);
    }

    /* Bilinear transform */
    bilinear_transform(poles, np, fs);
    bilinear_transform(zeros, nz, fs);

    /* Zero-pad: prototype zeros at s=∞ → z = -1 (not for BS). */
    if (type != FILTER_BANDSTOP) {
        for (uint8_t i = nz; i < np; i++) {
            zeros[i].re = -1.0f;
            zeros[i].im =  0.0f;
        }
        nz = np;
    }

    /* ZPK → SOS → deploy */
    uint8_t ns = (np + 1) / 2;
    float (*sos)[6] = (float (*)[6])malloc((size_t)ns * 6 * sizeof(float));
    if (sos == NULL) return;

    uint8_t n_sections = zpk2sos(zeros, poles, np, sos, k);

    c->sos = (sos_filter_t *)malloc(sizeof(sos_filter_t));
    if (c->sos == NULL) { free(sos); return; }

    sos_filter_init(c->sos, n_sections);
    for (uint8_t i = 0; i < n_sections; i++) {
        float num[3] = {sos[i][0], sos[i][1], sos[i][2]};
        float den[3] = {sos[i][3], sos[i][4], sos[i][5]};
        sos_filter_set_section(c->sos, i, num, den);
    }

    free(sos);
    c->valid = 1;
}

void cheby2_destroy(cheby2_t *c)
{
    if (c->sos != NULL) {
        sos_filter_destroy(c->sos);
        free(c->sos);
        c->sos = NULL;
    }
    c->valid = 0;
}

float cheby2_update(cheby2_t *c, float input)
{
    if (!c->valid || c->sos == NULL) return input;
    return sos_filter_update(c->sos, input);
}

void cheby2_reset(cheby2_t *c, float equilibrium)
{
    if (!c->valid || c->sos == NULL) return;
    sos_filter_reset(c->sos, equilibrium);
}
