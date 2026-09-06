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

static void c_mul(float *rr, float *ri, float ar, float ai, float br, float bi);
static void c_div(float *rr, float *ri, float ar, float ai, float br, float bi);

/* Stable roots of s² + B·s + C = 0 (B, C complex): root1 = (−B − √(B²−4C))/2
   and root2 = C / root1.  The direct formula (−B + √…)/2 cancels
   catastrophically when |B|² ≫ 4|C| (real prototype poles of wideband
   BP/BS designs, |xi·p| ≈ √disc): the f32 result carries ~ulp(|B|)/2 of
   imaginary noise, which then breaks conjugate pairing downstream.
   Computing the small root as C / root1 avoids the cancellation entirely. */
static void stable_roots(float B_re, float B_im, float C_re, float C_im,
                         complex_t *r1, complex_t *r2)
{
    float d_re, d_im;
    c_mul(&d_re, &d_im, B_re, B_im, B_re, B_im);
    d_re -= 4.0f * C_re;
    d_im -= 4.0f * C_im;
    complex_t sd = c_sqrt(d_re, d_im);

    r1->re = 0.5f * (-B_re - sd.re);
    r1->im = 0.5f * (-B_im - sd.im);

    if (r1->re * r1->re + r1->im * r1->im < 1e-20f) {
        /* Degenerate large root — fall back to the direct formula. */
        r2->re = 0.5f * (-B_re + sd.re);
        r2->im = 0.5f * (-B_im + sd.im);
    } else {
        c_div(&r2->re, &r2->im, C_re, C_im, r1->re, r1->im);
    }
}

void analog_bp_transform(complex_t *poles, uint8_t *np,
                         complex_t *zeros, uint8_t *nz,
                         float w0, float xi)
{
    uint8_t np_old = *np;
    uint8_t nz_old = *nz;
    float w0_sq = w0 * w0;

    /* Transform poles backward to avoid clobbering.
       s → (s² + w0²) / (xi·s)  ⇒  s² − xi·p·s + w0² = 0. */
    for (int i = np_old - 1; i >= 0; i--) {
        stable_roots(-xi * poles[i].re, -xi * poles[i].im,
                     w0_sq, 0.0f,
                     &poles[2 * i], &poles[2 * i + 1]);
    }
    *np = 2 * np_old;

    /* Transform finite zeros with the same formula. */
    for (int i = nz_old - 1; i >= 0; i--) {
        stable_roots(-xi * zeros[i].re, -xi * zeros[i].im,
                     w0_sq, 0.0f,
                     &zeros[2 * i], &zeros[2 * i + 1]);
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
    float w0_sq = w0 * w0;

    /* Transform poles backward.  p → ξ·s / (s² + ω₀²)
       ⇒  s² − (ξ/p)·s + ω₀² = 0. */
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

        float b_re, b_im;
        c_div(&b_re, &b_im, -xi, 0.0f, pr, pi); /* B = −ξ / p */
        stable_roots(b_re, b_im, w0_sq, 0.0f,
                     &poles[2 * i], &poles[2 * i + 1]);
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

        float b_re, b_im;
        c_div(&b_re, &b_im, -xi, 0.0f, zr, zi); /* B = −ξ / z */
        stable_roots(b_re, b_im, w0_sq, 0.0f,
                     &zeros[2 * i], &zeros[2 * i + 1]);
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

/* Relative tolerance: |im| <= eps * |z|  (matches scipy _cplxreal). */
static inline int is_real(const complex_t *a, float eps)
{
    return fabsf(a->im) <= eps * sqrtf(a->re * a->re + a->im * a->im);
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

/* Complex multiply: r = a * b. */
static void c_mul(float *rr, float *ri, float ar, float ai, float br, float bi)
{
    *rr = ar * br - ai * bi;
    *ri = ar * bi + ai * br;
}

/* Complex divide: r = a / b  (b ≠ 0). */
static void c_div(float *rr, float *ri, float ar, float ai, float br, float bi)
{
    float den = br * br + bi * bi;
    if (den == 0.0f) { *rr = 0.0f; *ri = 0.0f; return; }
    *rr = (ar * br + ai * bi) / den;
    *ri = (ai * br - ar * bi) / den;
}

/* ================================================================== */
/*  ZPK gain helpers (mirrors scipy zpk gain tracking)                */
/* ================================================================== */

float zpk_lp_gain(float k, float wo, uint8_t degree)
{
    float result = k;
    for (uint8_t i = 0; i < degree; i++)
        result *= wo;
    return result;
}

float zpk_hp_bs_gain(float k, const complex_t *z, uint8_t nz,
                     const complex_t *p, uint8_t np)
{
    /* Compute k * ∏(−z_i) / ∏(−p_i) incrementally to avoid overflow. */
    float rr = k, ri = 0.0f;
    uint8_t i;
    for (i = 0; i < nz && i < np; i++) {
        /* Multiply by (-z[i]), divide by (-p[i]). */
        c_mul(&rr, &ri, rr, ri, -z[i].re, -z[i].im);
        c_div(&rr, &ri, rr, ri, -p[i].re, -p[i].im);
    }
    for (; i < nz; i++)
        c_mul(&rr, &ri, rr, ri, -z[i].re, -z[i].im);
    for (; i < np; i++)
        c_div(&rr, &ri, rr, ri, -p[i].re, -p[i].im);
    return rr; /* imaginary part cancels for conjugate pairs */
}

float zpk_bp_gain(float k, float bw, uint8_t degree)
{
    float result = k;
    for (uint8_t i = 0; i < degree; i++)
        result *= bw;
    return result;
}

float bilinear_zpk_gain(float k, const complex_t *z, uint8_t nz,
                         const complex_t *p, uint8_t np, float K)
{
    /* Compute k * ∏(K−z_i) / ∏(K−p_i) incrementally to avoid overflow. */
    float rr = k, ri = 0.0f;
    uint8_t i;
    for (i = 0; i < nz && i < np; i++) {
        c_mul(&rr, &ri, rr, ri, K - z[i].re, -z[i].im);
        c_div(&rr, &ri, rr, ri, K - p[i].re, -p[i].im);
    }
    for (; i < nz; i++)
        c_mul(&rr, &ri, rr, ri, K - z[i].re, -z[i].im);
    for (; i < np; i++)
        c_div(&rr, &ri, rr, ri, K - p[i].re, -p[i].im);
    return rr;
}

float bilinear_zpk_gain_scaled(float k, float s, uint8_t degree,
                               const complex_t *z, uint8_t nz,
                               const complex_t *p, uint8_t np, float K)
{
    /* k · s^degree · ∏(K−z_i) / ∏(K−p_i) in interleaved order so the
       running product stays bounded.  Computing s^degree standalone
       overflows f32 for near-Nyquist LP/BP (e.g. wc^8 ≈ FLT_MAX) even
       though the final k is small — the ∏(K−p) factors are the same size
       as s (LP: |K−p| ≈ wc) but the cancellation never happens once an
       intermediate has overflowed.  Pairing each s factor with one
       (K−p) division keeps the per-step ratio s/|K−p| bounded (~1.03 for
       wc = 63600, K = 2000). */
    float rr = k, ri = 0.0f;
    uint8_t i_s = 0, i_z = 0, i_p = 0;

    while (i_s < degree || i_z < nz || i_p < np) {
        if (i_s < degree && i_p < np) {
            c_mul(&rr, &ri, rr, ri, s, 0.0f);
            c_div(&rr, &ri, rr, ri, K - p[i_p].re, -p[i_p].im);
            i_s++;
            i_p++;
        } else if (i_z < nz && i_p < np) {
            c_mul(&rr, &ri, rr, ri, K - z[i_z].re, -z[i_z].im);
            c_div(&rr, &ri, rr, ri, K - p[i_p].re, -p[i_p].im);
            i_z++;
            i_p++;
        } else if (i_s < degree) {
            c_mul(&rr, &ri, rr, ri, s, 0.0f);
            i_s++;
        } else if (i_z < nz) {
            c_mul(&rr, &ri, rr, ri, K - z[i_z].re, -z[i_z].im);
            i_z++;
        } else {
            c_div(&rr, &ri, rr, ri, K - p[i_p].re, -p[i_p].im);
            i_p++;
        }
    }
    return rr; /* imaginary part cancels for conjugate pairs */
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
    float best = 1e30f;
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

/* Like find_nearest, restricted to real (want_real=1) or complex elements.
   Returns nearest distance via *best; 1e30f if no element matches. */
static uint8_t find_nearest_typed(const complex_t *arr, const uint8_t *used,
                                  uint8_t n, const complex_t *target,
                                  int want_real, float eps, float *best)
{
    uint8_t idx = 0;
    *best = 1e30f;
    for (uint8_t i = 0; i < n; i++) {
        if (used[i]) continue;
        if (is_real(&arr[i], eps) != want_real) continue;
        float d = c_dist(&arr[i], target);
        if (d < *best) {
            *best = d;
            idx = i;
        }
    }
    return idx;
}

/* Count unused elements with a given real/complex property. */
static uint8_t count_used(const uint8_t *used, uint8_t n,
                          const complex_t *arr, int want_real, float eps)
{
    uint8_t cnt = 0;
    for (uint8_t i = 0; i < n; i++) {
        if (!used[i] && is_real(&arr[i], eps) == want_real) {
            cnt++;
        }
    }
    return cnt;
}

/* Find the conjugate of zp[idx] among unused elements and mark it used.
   Matching tolerance is relative: |candidate - conj(target)| <= eps * |target|.
   On failure, synthesises the conjugate (no element consumed). */
static void claim_conjugate(const complex_t *arr, uint8_t *used, uint8_t n,
                            uint8_t idx, complex_t *out, float eps)
{
    float mag = sqrtf(arr[idx].re * arr[idx].re + arr[idx].im * arr[idx].im);
    float thresh = eps * (mag > 1e-12f ? mag : 1.0f);
    for (uint8_t i = 0; i < n; i++) {
        if (used[i]) continue;
        if (fabsf(arr[i].re - arr[idx].re) <= thresh
            && fabsf(arr[i].im + arr[idx].im) <= thresh) {
            used[i] = 1;
            /* Average the pair (mirrors scipy _cplxreal).
               arr[i] has opposite imag sign to arr[idx].
               Return arr[i] averaged with conj(arr[idx]) so the caller
               gets the true conjugate of the primary element. */
            out->re = 0.5f * (arr[idx].re + arr[i].re);
            out->im = 0.5f * (arr[i].im - arr[idx].im);
            return;
        }
    }
    /* Last resort: numerical mismatch — synthesise the conjugate. */
    out->re =  arr[idx].re;
    out->im = -arr[idx].im;
}

/* Roots of z² + c1·z + c2.  c2 == 0 → single finite root at −c1 (the
   second root sits at infinity and is ignored).  Returns the count. */
static uint8_t poly_roots(float c1, float c2, complex_t roots[2])
{
    if (c2 == 0.0f) {
        roots[0].re = -c1;
        roots[0].im = 0.0f;
        return 1;
    }
    float disc = c1 * c1 - 4.0f * c2;
    if (disc >= 0.0f) {
        float s = sqrtf(disc);
        roots[0].re = 0.5f * (-c1 + s);
        roots[0].im = 0.0f;
        roots[1].re = 0.5f * (-c1 - s);
        roots[1].im = 0.0f;
    } else {
        float im = 0.5f * sqrtf(-disc);
        roots[0].re = -0.5f * c1;
        roots[0].im =  im;
        roots[1].re = -0.5f * c1;
        roots[1].im = -im;
    }
    return 2;
}

/* True if every root matches a distinct unused element of pool within tol
   (Euclidean).  A multiset check: members of a tight pair may match either
   way around, which is harmless — what it rejects is a root that reproduces
   no input element at all (duplicated pair displacing a real one). */
static int match_roots(const complex_t *roots, uint8_t nr,
                       const complex_t *pool, uint8_t n,
                       uint8_t *matched, float tol)
{
    for (uint8_t r = 0; r < nr; r++) {
        uint8_t bi = n;
        float best = tol;
        for (uint8_t i = 0; i < n; i++) {
            if (matched[i]) continue;
            float d = c_dist(&roots[r], &pool[i]);
            if (d < best) {
                best = d;
                bi = i;
            }
        }
        if (bi == n) return 0;
        matched[bi] = 1;
    }
    return 1;
}

/* Maximum number of pole-zero pairs supported by zpk2sos.
 * 8th-order prototype → BP/BS doubles to 16. */
#define ZPK2SOS_MAX_N 16

uint8_t zpk2sos_impl(complex_t *zeros, complex_t *poles, uint8_t n,
                     float (*sos)[6], float k)
{
    if (n == 0 || n > ZPK2SOS_MAX_N) return 0;

    uint8_t used_p[ZPK2SOS_MAX_N];
    uint8_t used_z[ZPK2SOS_MAX_N];
    memset(used_p, 0, n * sizeof(uint8_t));
    memset(used_z, 0, n * sizeof(uint8_t));

    uint8_t n_p = n;
    uint8_t n_z = n;
    uint8_t section = 0;
    uint8_t max_sections = (n + 1) / 2;
    /* Real/complex classification and conjugate-claim tolerance.  The
       bilinear transform's f32 rounding can split a conjugate pair by
       ~2e-4 absolute (near the unit circle, K² − |s|² cancellation), so
       1e-4 was too tight: a failed claim synthesises a conjugate and the
       n_p bookkeeping silently loses the real partner.  1e-3 is still far
       below the separation between distinct pole pairs (≈ 2.0). */
    const float eps_real = 1e-3f;

    while (n_p > 0) {
        /* Safety caps: never exceed allocated rows
           or decrement counters past zero. */
        if (section >= max_sections) break;
        if (n_z == 0) break;

        /* 1. Pick the most unfavorable (largest |p|) remaining pole. */
        uint8_t p1_i = find_worst_pole(poles, used_p, n);
        complex_t p1 = poles[p1_i];
        used_p[p1_i] = 1;
        n_p--;

        complex_t p2, z1, z2;

        if (is_real(&p1, eps_real)) {
            /* p1 is real — try to pair with another real pole. */
            if (count_used(used_p, n, poles, 1, eps_real) > 0) {
                float best;
                uint8_t p2_i = find_nearest_typed(poles, used_p, n, &p1,
                                                  1, eps_real, &best);
                p2 = poles[p2_i];
                used_p[p2_i] = 1;
                n_p--;
            } else {
                /* No more real poles → first-order section with a real zero.
                   If no real zero remains, the pole/zero sets are unbalanced;
                   fail closed rather than flatten a complex zero into a real
                   one and drop its conjugate from the cascade. */
                float best;
                uint8_t z1_i = find_nearest_typed(zeros, used_z, n, &p1,
                                                  1, eps_real, &best);
                if (best == 1e30f) return 0;
                z1 = zeros[z1_i];
                used_z[z1_i] = 1;
                n_z--;

                float b[3] = {1.0f, -z1.re, 0.0f};
                float a[2] = {-p1.re, 0.0f};

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
            claim_conjugate(poles, used_p, n, p1_i, &p2, eps_real);
            n_p--;
        }

        /* 2. Match two zeros to this pole pair. */
        uint8_t z1_i = find_nearest(zeros, used_z, n, &p1);

        if (is_real(&zeros[z1_i], eps_real)) {
            if (count_used(used_z, n, zeros, 1, eps_real) > 1) {
                /* Two real zeros available — use z1 and the nearest other real. */
                z1 = zeros[z1_i];
                used_z[z1_i] = 1;
                n_z--;

                float best;
                uint8_t z2_i = find_nearest_typed(zeros, used_z, n, &p1,
                                                  1, eps_real, &best);
                z2 = zeros[z2_i];
                used_z[z2_i] = 1;
                n_z--;
            } else {
                /* Only one real zero left — save it for a later 1st-order
                   section and use a complex pair instead. */
                uint8_t best_i = z1_i;
                float best_d = 1e30f;
                for (uint8_t i = 0; i < n; i++) {
                    if (used_z[i]) continue;
                    if (is_real(&zeros[i], eps_real)) continue;
                    float d = c_dist(&zeros[i], &p1);
                    if (d < best_d) { best_d = d; best_i = i; }
                }
                z1 = zeros[best_i];
                used_z[best_i] = 1;
                n_z--;

                claim_conjugate(zeros, used_z, n, best_i, &z2, eps_real);
                n_z--;
            }
        } else {
            /* z1 is complex — use it and its conjugate. */
            z1 = zeros[z1_i];
            used_z[z1_i] = 1;
            n_z--;

            claim_conjugate(zeros, used_z, n, z1_i, &z2, eps_real);
            n_z--;
        }

        /* 3. Compute biquad coefficients. */
        float b[3], a[2];
        make_biquad(&p1, &p2, &z1, &z2, b, a);

        sos[section][0] = b[0];
        sos[section][1] = b[1];
        sos[section][2] = b[2];
        sos[section][3] = 1.0f;
        sos[section][4] = a[0];
        sos[section][5] = a[1];
        section++;
    }

    /* 3b. Invariant: every pole and zero must have been claimed.  A
       synthesised conjugate leaves an unclaimed gap, meaning the n_p/n_z
       bookkeeping diverged from the arrays and poles/zeros were silently
       dropped from the cascade (observed: a misplaced near-duplicate
       section displacing a missing low-Q pair → ~350x resonance).
       Fail closed instead. */
    for (uint8_t i = 0; i < n; i++) {
        if (!used_p[i] || !used_z[i]) return 0;
    }

    /* 3c. Verify each section's roots reproduce the input pole/zero
       multisets.  Guards against cross-pair misclaims that leave used[]
       fully set yet pair the wrong elements together. */
    {
        uint8_t matched[ZPK2SOS_MAX_N];
        float tol = 1e-3f;

        memset(matched, 0, sizeof(matched));
        for (uint8_t s = 0; s < section; s++) {
            complex_t roots[2];
            uint8_t nr = poly_roots(sos[s][4], sos[s][5], roots);
            if (!match_roots(roots, nr, poles, n, matched, tol)) return 0;
        }
        memset(matched, 0, sizeof(matched));
        for (uint8_t s = 0; s < section; s++) {
            complex_t roots[2];
            uint8_t nr = poly_roots(sos[s][1], sos[s][2], roots); /* b0 == 1 */
            if (!match_roots(roots, nr, zeros, n, matched, tol)) return 0;
        }
    }

    /* 4. Reverse section order: slowest poles last → fastest first. */
    for (uint8_t i = 0; i < section / 2; i++) {
        for (uint8_t j = 0; j < 6; j++) {
            float tmp = sos[i][j];
            sos[i][j] = sos[section - 1 - i][j];
            sos[section - 1 - i][j] = tmp;
        }
    }

    /* 5. Apply overall system gain to the first section numerator. */
    sos[0][0] *= k;
    sos[0][1] *= k;
    sos[0][2] *= k;

    return section;
}

uint8_t zpk2sos(const complex_t *zeros, const complex_t *poles, uint8_t n,
                float (*sos)[6], float k)
{
    if (n == 0 || n > ZPK2SOS_MAX_N) return 0;

    complex_t wz[ZPK2SOS_MAX_N], wp[ZPK2SOS_MAX_N];
    memcpy(wp, poles, n * sizeof(complex_t));
    memcpy(wz, zeros, n * sizeof(complex_t));

    return zpk2sos_impl(wz, wp, n, sos, k);
}
