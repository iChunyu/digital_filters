#include "butter_filter.h"
#include <math.h>
#include <stdlib.h>

static void butter_prototype(complex_t *poles, uint8_t n)
{
    for (uint8_t k = 0; k < n; k++) {
        float theta = (float)(2 * (k + 1) + n - 1) / (2.0f * n) * (float)M_PI;
        poles[k].re = cosf(theta);
        poles[k].im = -sinf(theta);
    }
}

void butter_init(butter_t *b, uint8_t type, uint8_t order,
                 float fc1, float fc2, float fs)
{
    b->valid = 0;
    b->type  = type;
    b->order = order;
    b->fc1   = fc1;
    b->fc2   = fc2;
    b->fs    = fs;
    b->sos   = NULL;

    /* Validate */
    if (order == 0 || fc1 <= 0.0f || fc1 >= fs * 0.5f) return;
    if (type == FILTER_BANDPASS || type == FILTER_BANDSTOP) {
        if (fc2 <= fc1 || fc2 >= fs * 0.5f) return;
    }

    /* Digital normalisation frequency (rad/sample). */
    float w0_norm;
    switch (type) {
    case FILTER_LOWPASS:  w0_norm = 0.0f;            break;
    case FILTER_HIGHPASS: w0_norm = (float)M_PI;     break;
    case FILTER_BANDPASS: w0_norm = 2.0f * (float)M_PI
                                  * sqrtf(fc1 * fc2) / fs; break;
    case FILTER_BANDSTOP: w0_norm = 0.0f;            break;
    default: return;
    }

    /* 1. Pre-warp digital cutoff(s) → analog radian frequency */
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

    /* 2. Butterworth prototype — N poles, no zeros */
    complex_t poles[64];
    complex_t zeros[64];
    uint8_t np = order;
    uint8_t nz = 0;

    butter_prototype(poles, order);

    /* 3. Analog frequency transform */
    switch (type) {
    case FILTER_LOWPASS:
        analog_lp_transform(poles, np, zeros, nz, w0_analog);
        break;
    case FILTER_HIGHPASS:
        analog_hp_transform(poles, np, zeros, nz, w0_analog);
        for (uint8_t i = nz; i < np; i++) {
            zeros[i].re = 0.0f;
            zeros[i].im = 0.0f;
        }
        nz = np;
        break;
    case FILTER_BANDPASS:
        analog_bp_transform(poles, &np, zeros, &nz, w0_analog, xi);
        break;
    case FILTER_BANDSTOP:
        analog_bs_transform(poles, &np, zeros, &nz, w0_analog, xi);
        break;
    default:
        return;
    }

    /* 4. Bilinear transform: s → z */
    bilinear_transform(poles, np, fs);
    bilinear_transform(zeros, nz, fs);

    /* 5. Add zeros at z = -1 for excess poles (prototype zeros at s = ∞).
          For BS, all zeros are already finite (none at ∞). */
    if (type != FILTER_BANDSTOP) {
        for (uint8_t i = nz; i < np; i++) {
            zeros[i].re = -1.0f;
            zeros[i].im =  0.0f;
        }
        nz = np;
    }

    /* 6. Pair poles and zeros → SOS coefficients. */
    uint8_t ns = (np + 1) / 2; /* ceil(N/2) sections */
    float (*sos)[6] = (float (*)[6])malloc((size_t)ns * 6 * sizeof(float));
    if (sos == NULL) return;

    uint8_t n_sections = zpk2sos(zeros, poles, np, sos, w0_norm);

    /* 7. Deploy to sos_filter_t. */
    b->sos = (sos_filter_t *)malloc(sizeof(sos_filter_t));
    if (b->sos == NULL) {
        free(sos);
        return;
    }

    sos_filter_init(b->sos, n_sections);
    for (uint8_t i = 0; i < n_sections; i++) {
        float num[3] = {sos[i][0], sos[i][1], sos[i][2]};
        float den[3] = {sos[i][3], sos[i][4], sos[i][5]};
        sos_filter_set_section(b->sos, i, num, den);
    }

    free(sos);
    b->valid = 1;
}

void butter_destroy(butter_t *b)
{
    if (b->sos != NULL) {
        sos_filter_destroy(b->sos);
        free(b->sos);
        b->sos = NULL;
    }
    b->valid = 0;
}

float butter_update(butter_t *b, float input)
{
    if (!b->valid || b->sos == NULL) {
        return input;
    }
    return sos_filter_update(b->sos, input);
}

float butter_reset(butter_t *b, float equilibrium)
{
    if (!b->valid || b->sos == NULL) {
        return equilibrium;
    }
    sos_filter_reset(b->sos, equilibrium);
    return sos_filter_get_output(b->sos);
}
