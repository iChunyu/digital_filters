#include "cheby_filter.h"
#include <math.h>
#include <stdlib.h>

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
    if (type == FILTER_LOWPASS || type == FILTER_HIGHPASS)
        (void)fc2;

    float epsilon = sqrtf(powf(10.0f, ripple_db / 10.0f) - 1.0f);
    float wc = 2.0f * (float)M_PI * prewarp(fc1, fs);

    /* 1. Chebyshev I prototype — N poles, no zeros */
    complex_t poles[32];
    complex_t zeros[32];
    uint8_t np = order;
    uint8_t nz = 0;

    cheby1_prototype(poles, order, epsilon);

    /* 2. Analog frequency transform */
    switch (type) {
    case FILTER_LOWPASS:
        analog_lp_transform(poles, np, zeros, nz, wc);
        break;
    case FILTER_HIGHPASS:
        analog_hp_transform(poles, np, zeros, nz, wc);
        for (uint8_t i = nz; i < np; i++) {
            zeros[i].re = 0.0f;
            zeros[i].im = 0.0f;
        }
        nz = np;
        break;
    default:
        return;
    }

    /* 3. Bilinear transform */
    bilinear_transform(poles, np, fs);
    bilinear_transform(zeros, nz, fs);

    /* 4. Zero-pad: prototype zeros at s=∞ → z = -1 */
    for (uint8_t i = nz; i < np; i++) {
        zeros[i].re = -1.0f;
        zeros[i].im =  0.0f;
    }
    nz = np;

    /* 5. ZPK → SOS → deploy */
    uint8_t ns = (np + 1) / 2;
    float (*sos)[6] = (float (*)[6])malloc((size_t)ns * 6 * sizeof(float));
    if (sos == NULL) return;

    uint8_t n_sections = zpk2sos(zeros, poles, np, sos, type);

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

float cheby1_reset(cheby1_t *c, float equilibrium)
{
    if (!c->valid || c->sos == NULL) return equilibrium;
    sos_filter_reset(c->sos, equilibrium);
    return sos_filter_get_output(c->sos);
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
    if (type == FILTER_LOWPASS || type == FILTER_HIGHPASS)
        (void)fc2;

    /* ε for Chebyshev II: stopband gain = ε / sqrt(1+ε²) = 10^(-dB/20) */
    float epsilon = 1.0f / sqrtf(powf(10.0f, ripple_db / 10.0f) - 1.0f);
    float wc = 2.0f * (float)M_PI * prewarp(fc1, fs);

    /* 1. Chebyshev II prototype — N poles, N (even) or N-1 (odd) finite zeros */
    complex_t poles[32];
    complex_t zeros[32];
    uint8_t np = order;
    uint8_t nz = cheby2_prototype(poles, zeros, order, epsilon);

    /* 2. Analog frequency transform */
    switch (type) {
    case FILTER_LOWPASS:
        analog_lp_transform(poles, np, zeros, nz, wc);
        break;
    case FILTER_HIGHPASS:
        analog_hp_transform(poles, np, zeros, nz, wc);
        for (uint8_t i = nz; i < np; i++) {
            zeros[i].re = 0.0f;
            zeros[i].im = 0.0f;
        }
        nz = np;
        break;
    default:
        return;
    }

    /* 3. Bilinear transform */
    bilinear_transform(poles, np, fs);
    bilinear_transform(zeros, nz, fs);

    /* 4. Zero-pad: prototype zeros at s=∞ → z = -1 */
    for (uint8_t i = nz; i < np; i++) {
        zeros[i].re = -1.0f;
        zeros[i].im =  0.0f;
    }
    nz = np;

    /* 5. ZPK → SOS → deploy */
    uint8_t ns = (np + 1) / 2;
    float (*sos)[6] = (float (*)[6])malloc((size_t)ns * 6 * sizeof(float));
    if (sos == NULL) return;

    uint8_t n_sections = zpk2sos(zeros, poles, np, sos, type);

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

float cheby2_reset(cheby2_t *c, float equilibrium)
{
    if (!c->valid || c->sos == NULL) return equilibrium;
    sos_filter_reset(c->sos, equilibrium);
    return sos_filter_get_output(c->sos);
}
