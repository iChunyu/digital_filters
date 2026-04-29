#include "filter_utils.h"
#include <math.h>
#include <string.h>

/* ------------------------------------------------------------------ */

float prewarp(float fd, float fs)
{
    if (fd <= 0.0f || fd >= fs * 0.5f) {
        return fd; /* let caller clamp — degenerate case returns as-is */
    }
    return fs / ((float)M_PI) * tanf((float)M_PI * fd / fs);
}

/* ------------------------------------------------------------------ */

static void scale_complex(complex_t *zp, uint8_t n, float s)
{
    for (uint8_t i = 0; i < n; i++) {
        zp[i].re *= s;
        zp[i].im *= s;
    }
}

void analog_lp_transform(complex_t *poles, uint8_t np,
                         complex_t *zeros, uint8_t nz, float wc)
{
    (void)np;
    (void)nz;
    scale_complex(poles, np, wc);
    if (nz > 0) {
        scale_complex(zeros, nz, wc);
    }
}

/* ------------------------------------------------------------------ */

void analog_hp_transform(complex_t *poles, uint8_t np,
                         complex_t *zeros, uint8_t nz, float wc)
{
    for (uint8_t i = 0; i < np; i++) {
        float den = poles[i].re * poles[i].re + poles[i].im * poles[i].im;
        /* p = wc / p_proto = wc * conj(p_proto) / |p_proto|^2 */
        /* 1 / (a + jb) = (a - jb) / (a^2 + b^2) */
        poles[i].re =  wc * poles[i].re / den;
        poles[i].im = -wc * poles[i].im / den;
    }
    for (uint8_t i = 0; i < nz; i++) {
        float den = zeros[i].re * zeros[i].re + zeros[i].im * zeros[i].im;
        zeros[i].re =  wc * zeros[i].re / den;
        zeros[i].im = -wc * zeros[i].im / den;
    }
}

/* ------------------------------------------------------------------ */

void bilinear_transform(complex_t *zp, uint8_t n, float fs)
{
    float K = 2.0f * fs;
    float K2 = K * K;

    for (uint8_t i = 0; i < n; i++) {
        float s_re = zp[i].re;
        float s_im = zp[i].im;

        /* z = (K + s) / (K - s) = (K + s)(K - conj(s)) / |K - s|^2 */
        float den = K2 - 2.0f * K * s_re + s_re * s_re + s_im * s_im;

        zp[i].re = (K2 - s_re * s_re - s_im * s_im) / den;
        zp[i].im = (2.0f * K * s_im) / den;
    }
}

/* ------------------------------------------------------------------ */

static inline float c_abs(const complex_t *a)
{
    return sqrtf(a->re * a->re + a->im * a->im);
}

static inline float c_dist(const complex_t *a, const complex_t *b)
{
    float dr = a->re - b->re;
    float di = a->im - b->im;
    return sqrtf(dr * dr + di * di);
}

static inline int is_real(const complex_t *a, float eps)
{
    return fabsf(a->im) <= eps;
}

/* Remove element at idx from array of length *len (shift remaining left).
   Used with a working copy of the pole/zero lists. */
static void remove_elem(complex_t *arr, uint8_t *len, uint8_t idx)
{
    if (idx < *len - 1) {
        memmove(&arr[idx], &arr[idx + 1],
                (size_t)(*len - idx - 1) * sizeof(complex_t));
    }
    (*len)--;
}

/* Compute biquad coefficients from a pole pair and zero pair.
   Poles p1,p2 and zeros z1,z2 are either both real or conjugate pairs. */
static void make_biquad(const complex_t *p1, const complex_t *p2,
                        const complex_t *z1, const complex_t *z2,
                        float *b, float *a)
{
    /* a1 = -(p1 + p2),  a2 = p1 * p2 */
    a[0] = -(p1->re + p2->re);          /* a1 */
    a[1] = p1->re * p2->re - p1->im * p2->im;  /* a2 = Re(p1*p2)      */
    /* (when conjugates this is just |p1|^2; when real im=0 so fine)   */

    /* b0 = 1, b1 = -(z1 + z2), b2 = z1 * z2 */
    b[0] = 1.0f;
    b[1] = -(z1->re + z2->re);
    b[2] = z1->re * z2->re - z1->im * z2->im;
}

/* Normalise b coefs for unity gain at reference frequency. */
static void norm_gain(float *b, const float *a, filter_type_e type)
{
    float gain_num, gain_den;

    if (type == FILTER_HIGHPASS) {
        /* Unity gain at Nyquist (z = -1) */
        gain_num = b[0] - b[1] + b[2];
        gain_den = 1.0f - a[0] + a[1];
    } else {
        /* Unity gain at DC (z = 1) — default for LP, BP, BS */
        gain_num = b[0] + b[1] + b[2];
        gain_den = 1.0f + a[0] + a[1];
    }

    if (fabsf(gain_num) < 1e-12f || fabsf(gain_den) < 1e-12f) {
        return; /* degenerate — leave as-is */
    }

    float inv_gain = gain_den / gain_num;
    b[0] *= inv_gain;
    b[1] *= inv_gain;
    b[2] *= inv_gain;
}

/* Find index of the pole closest to the unit circle (largest magnitude). */
static uint8_t find_worst_pole(const complex_t *poles, const uint8_t *used, uint8_t n)
{
    uint8_t idx = 0;
    float best = -1.0f;
    for (uint8_t i = 0; i < n; i++) {
        if (used[i]) continue;
        float mag = poles[i].re * poles[i].re + poles[i].im * poles[i].im;
        if (mag > best) {
            best = mag;
            idx = i;
        }
    }
    return idx;
}

/* Find index of the (used) element nearest to target. */
static uint8_t find_nearest(const complex_t *arr, const uint8_t *used, uint8_t n,
                            const complex_t *target)
{
    uint8_t idx = 0;
    float best = INFINITY;
    for (uint8_t i = 0; i < n; i++) {
        if (used[i]) continue;
        float d = c_dist(&arr[i], target);
        if (d < best) {
            best = d;
            idx = i;
        }
    }
    return idx;
}

/* Count used elements with a given real/complex property. */
static uint8_t count_used(const uint8_t *used, uint8_t n,
                          const complex_t *arr, int want_real)
{
    uint8_t cnt = 0;
    for (uint8_t i = 0; i < n; i++) {
        if (!used[i] && is_real(&arr[i], 1e-7f) == want_real) {
            cnt++;
        }
    }
    return cnt;
}

/* Find the conjugate of zp[idx] and mark it used. */
static void claim_conjugate(const complex_t *arr, uint8_t *used, uint8_t n,
                            uint8_t idx, complex_t *out)
{
    float eps = 1e-5f;
    for (uint8_t i = 0; i < n; i++) {
        if (used[i]) continue;
        if (fabsf(arr[i].re - arr[idx].re) < eps
            && fabsf(arr[i].im + arr[idx].im) < eps) {
            used[i] = 1;
            *out = arr[i];
            return;
        }
    }
    /* Fallback: exact conjugate not found; synthesise it.  This can happen
       due to floating-point rounding in the bilinear transform. */
    out->re =  arr[idx].re;
    out->im = -arr[idx].im;
}

uint8_t zpk2sos(const complex_t *zeros, const complex_t *poles, uint8_t n,
                float (*sos)[6], filter_type_e type)
{
    if (n == 0) return 0;

    /* Working copies so we can mutate ownership via used[] flags. */
    uint8_t used_p[32];
    uint8_t used_z[32];
    complex_t wp[32], wz[32];

    memcpy(wp, poles, n * sizeof(complex_t));
    memcpy(wz, zeros, n * sizeof(complex_t));
    memset(used_p, 0, n * sizeof(uint8_t));
    memset(used_z, 0, n * sizeof(uint8_t));

    uint8_t n_p = n;
    uint8_t n_z = n;
    uint8_t section = 0;
    const float eps_real = 1e-7f;

    while (n_p > 0) {
        /* 1. Pick the most unfavorable (largest |p|) remaining pole. */
        uint8_t p1_i = find_worst_pole(wp, used_p, n);
        complex_t p1 = wp[p1_i];
        used_p[p1_i] = 1;
        n_p--;

        complex_t p2, z1, z2;

        if (is_real(&p1, eps_real)) {
            /* p1 is real — try to pair with another real pole. */
            if (count_used(used_p, n, wp, 1) > 0) {
                uint8_t p2_i = find_nearest(wp, used_p, n, &p1);
                p2 = wp[p2_i];
                used_p[p2_i] = 1;
                n_p--;
            } else {
                /* No more real poles → first-order section with a real zero. */
                uint8_t z1_i = find_nearest(wz, used_z, n, &p1);
                z1 = wz[z1_i];
                used_z[z1_i] = 1;
                n_z--;

                float b[3] = {1.0f, -z1.re, 0.0f};
                float a[2] = {-p1.re, 0.0f};

                norm_gain(b, a, type);

                sos[section][0] = b[0];
                sos[section][1] = b[1];
                sos[section][2] = b[2];
                sos[section][3] = 1.0f;
                sos[section][4] = a[0];
                sos[section][5] = a[1];
                section++;
                continue;
            }
        } else {
            /* p1 is complex — pair with its conjugate. */
            claim_conjugate(wp, used_p, n, p1_i, &p2);
            n_p--;
        }

        /* 2. Match two zeros to this pole pair. */
        uint8_t z1_i = find_nearest(wz, used_z, n, &p1);

        if (is_real(&wz[z1_i], eps_real)) {
            if (count_used(used_z, n, wz, 1) > 1) {
                /* Two real zeros available — use z1 and the nearest other real. */
                z1 = wz[z1_i];
                used_z[z1_i] = 1;
                n_z--;

                uint8_t z2_i = find_nearest(wz, used_z, n, &p1);
                z2 = wz[z2_i];
                used_z[z2_i] = 1;
                n_z--;
            } else {
                /* Only one real zero left — save it for a later 1st-order
                   section and use a complex pair instead. */
                uint8_t best_i = z1_i;
                float best_d = INFINITY;
                for (uint8_t i = 0; i < n; i++) {
                    if (used_z[i]) continue;
                    if (is_real(&wz[i], eps_real)) continue;
                    float d = c_dist(&wz[i], &p1);
                    if (d < best_d) { best_d = d; best_i = i; }
                }
                z1 = wz[best_i];
                used_z[best_i] = 1;
                n_z--;

                claim_conjugate(wz, used_z, n, best_i, &z2);
                n_z--;
            }
        } else {
            /* z1 is complex — use it and its conjugate. */
            z1 = wz[z1_i];
            used_z[z1_i] = 1;
            n_z--;

            claim_conjugate(wz, used_z, n, z1_i, &z2);
            n_z--;
        }

        /* 3. Compute biquad coefficients. */
        float b[3], a[2];
        make_biquad(&p1, &p2, &z1, &z2, b, a);

        norm_gain(b, a, type);

        sos[section][0] = b[0];
        sos[section][1] = b[1];
        sos[section][2] = b[2];
        sos[section][3] = 1.0f;
        sos[section][4] = a[0];
        sos[section][5] = a[1];
        section++;
    }

    /* 4. Reverse section order: slowest poles last → fastest first. */
    for (uint8_t i = 0; i < section / 2; i++) {
        for (uint8_t j = 0; j < 6; j++) {
            float tmp = sos[i][j];
            sos[i][j] = sos[section - 1 - i][j];
            sos[section - 1 - i][j] = tmp;
        }
    }

    return section;
}
