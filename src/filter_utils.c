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

/* Complex square root: sqrt(re + j·im).
   Returns the principal branch (non-negative real part). */
static complex_t c_sqrt(float re, float im)
{
    complex_t r;
    float mag = sqrtf(re * re + im * im);
    if (mag < 1e-20f) {
        r.re = 0.0f;
        r.im = 0.0f;
        return r;
    }
    /* Handle negative-real case to avoid cancellation in (mag + re)/2. */
    if (re < 0.0f && fabsf(im) < 1e-12f * fabsf(re)) {
        r.re = 0.0f;
        r.im = sqrtf(-re);
        if (im < 0.0f) r.im = -r.im;
        return r;
    }
    float sr = sqrtf(0.5f * (mag + re));
    float si = sqrtf(0.5f * (mag - re));
    r.re = sr;
    r.im = (im >= 0.0f) ? si : -si;
    return r;
}

/* ------------------------------------------------------------------ */

void analog_bp_transform(complex_t *poles, uint8_t *np,
                         complex_t *zeros, uint8_t *nz,
                         float w0, float xi)
{
    uint8_t np_old = *np;
    uint8_t nz_old = *nz;
    float w0_sq = w0 * w0;

    /* Transform poles backward to avoid clobbering. */
    for (int i = np_old - 1; i >= 0; i--) {
        float pr = poles[i].re;
        float pi = poles[i].im;

        /* disc = (p·xi)² − 4·w0² */
        float a = xi * pr;
        float b = xi * pi;
        float disc_re = a * a - b * b - 4.0f * w0_sq;
        float disc_im = 2.0f * a * b;

        complex_t sd = c_sqrt(disc_re, disc_im);

        /* p1 = (p·xi + sqrt(disc)) / 2,  p2 = (p·xi − sqrt(disc)) / 2 */
        poles[2 * i].re     = 0.5f * (a + sd.re);
        poles[2 * i].im     = 0.5f * (b + sd.im);
        poles[2 * i + 1].re = 0.5f * (a - sd.re);
        poles[2 * i + 1].im = 0.5f * (b - sd.im);
    }
    *np = 2 * np_old;

    /* Transform finite zeros with the same formula. */
    for (int i = nz_old - 1; i >= 0; i--) {
        float zr = zeros[i].re;
        float zi = zeros[i].im;

        float a = xi * zr;
        float b = xi * zi;
        float disc_re = a * a - b * b - 4.0f * w0_sq;
        float disc_im = 2.0f * a * b;

        complex_t sd = c_sqrt(disc_re, disc_im);

        zeros[2 * i].re     = 0.5f * (a + sd.re);
        zeros[2 * i].im     = 0.5f * (b + sd.im);
        zeros[2 * i + 1].re = 0.5f * (a - sd.re);
        zeros[2 * i + 1].im = 0.5f * (b - sd.im);
    }

    /* Append (np_old − nz_old) zeros at the origin for infinite prototype zeros.
       These map to z = +1 after bilinear transform. */
    for (uint8_t i = 0; i < np_old - nz_old; i++) {
        zeros[2 * nz_old + i].re = 0.0f;
        zeros[2 * nz_old + i].im = 0.0f;
    }
    *nz = nz_old + np_old;  /* = 2·nz_old + (np_old − nz_old) */
}

/* ------------------------------------------------------------------ */

void analog_bs_transform(complex_t *poles, uint8_t *np,
                         complex_t *zeros, uint8_t *nz,
                         float w0, float xi)
{
    uint8_t np_old = *np;
    uint8_t nz_old = *nz;
    float xi_sq = xi * xi;
    float w0_sq = w0 * w0;

    /* Transform poles backward.  p_k = (ξ ± √(ξ² − 4·p²·ω₀²)) / (2·p) */
    for (int i = np_old - 1; i >= 0; i--) {
        float pr = poles[i].re;
        float pi = poles[i].im;
        float mag2 = pr * pr + pi * pi;
        if (mag2 < 1e-20f) {
            /* Degenerate: leave as-is (should not happen for stable prototypes). */
            poles[2 * i].re     = poles[i].re;
            poles[2 * i].im     = poles[i].im;
            poles[2 * i + 1].re = poles[i].re;
            poles[2 * i + 1].im = poles[i].im;
            continue;
        }

        /* disc = ξ² − 4·p²·ω₀² */
        float p2_re = pr * pr - pi * pi;
        float p2_im = 2.0f * pr * pi;
        float disc_re = xi_sq - 4.0f * w0_sq * p2_re;
        float disc_im =        - 4.0f * w0_sq * p2_im;

        complex_t sd = c_sqrt(disc_re, disc_im);

        /* p1 = (ξ + √disc) / (2·p),  p2 = (ξ − √disc) / (2·p) */
        /* Division by complex p:  (num) / p = num · conj(p) / |p|² */
        float inv = 0.5f / mag2;

        float n1_re = xi + sd.re;
        float n1_im =      sd.im;
        poles[2 * i].re     = inv * (n1_re * pr + n1_im * pi);
        poles[2 * i].im     = inv * (n1_im * pr - n1_re * pi);

        float n2_re = xi - sd.re;
        float n2_im =     -sd.im;
        poles[2 * i + 1].re = inv * (n2_re * pr + n2_im * pi);
        poles[2 * i + 1].im = inv * (n2_im * pr - n2_re * pi);
    }
    *np = 2 * np_old;

    /* Transform finite zeros with the same formula. */
    for (int i = nz_old - 1; i >= 0; i--) {
        float zr = zeros[i].re;
        float zi = zeros[i].im;
        float mag2 = zr * zr + zi * zi;
        if (mag2 < 1e-20f) {
            zeros[2 * i].re     = zeros[i].re;
            zeros[2 * i].im     = zeros[i].im;
            zeros[2 * i + 1].re = zeros[i].re;
            zeros[2 * i + 1].im = zeros[i].im;
            continue;
        }

        float z2_re = zr * zr - zi * zi;
        float z2_im = 2.0f * zr * zi;
        float disc_re = xi_sq - 4.0f * w0_sq * z2_re;
        float disc_im =        - 4.0f * w0_sq * z2_im;

        complex_t sd = c_sqrt(disc_re, disc_im);

        float inv = 0.5f / mag2;

        float n1_re = xi + sd.re;
        float n1_im =      sd.im;
        zeros[2 * i].re     = inv * (n1_re * zr + n1_im * zi);
        zeros[2 * i].im     = inv * (n1_im * zr - n1_re * zi);

        float n2_re = xi - sd.re;
        float n2_im =     -sd.im;
        zeros[2 * i + 1].re = inv * (n2_re * zr + n2_im * zi);
        zeros[2 * i + 1].im = inv * (n2_im * zr - n2_re * zi);
    }

    /* Append 2·(np_old − nz_old) zeros at ±jω₀ for infinite prototype zeros. */
    for (uint8_t i = 0; i < np_old - nz_old; i++) {
        zeros[2 * nz_old + 2 * i].re     =  0.0f;
        zeros[2 * nz_old + 2 * i].im     =  w0;
        zeros[2 * nz_old + 2 * i + 1].re =  0.0f;
        zeros[2 * nz_old + 2 * i + 1].im = -w0;
    }
    *nz = 2 * np_old;
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

/* Normalise b coefs for unity gain at digital frequency w0 (rad/sample).
   w0 = 0 → DC (z=1), w0 = π → Nyquist (z=-1), otherwise general z=e^{j·w0}. */
static void norm_gain(float *b, const float *a, float w0)
{
    float gain_num, gain_den;

    if (fabsf(w0) < 1e-12f) {
        /* DC (z = 1) — LP, BS */
        gain_num = b[0] + b[1] + b[2];
        gain_den = 1.0f + a[0] + a[1];
    } else if (fabsf(w0 - (float)M_PI) < 1e-12f) {
        /* Nyquist (z = -1) — HP */
        gain_num = b[0] - b[1] + b[2];
        gain_den = 1.0f - a[0] + a[1];
    } else {
        /* General frequency z = e^{j·w0} — BP centre frequency */
        float c1 = cosf(w0);
        float s1 = sinf(w0);
        float c2 = cosf(2.0f * w0);
        float s2 = sinf(2.0f * w0);

        float num_re = b[0] + b[1] * c1 + b[2] * c2;
        float num_im =      - b[1] * s1 - b[2] * s2;
        float den_re = 1.0f + a[0] * c1 + a[1] * c2;
        float den_im =      - a[0] * s1 - a[1] * s2;

        gain_num = sqrtf(num_re * num_re + num_im * num_im);
        gain_den = sqrtf(den_re * den_re + den_im * den_im);
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
                float (*sos)[6], float w0)
{
    if (n == 0) return 0;

    /* Working copies so we can mutate ownership via used[] flags. */
    uint8_t used_p[64];
    uint8_t used_z[64];
    complex_t wp[64], wz[64];

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

                norm_gain(b, a, w0);

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

        norm_gain(b, a, w0);

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
