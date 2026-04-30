#include <stdio.h>
#include <math.h>
#include "cheby_filter.h"
#include "static_cheby_filter.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

#define FS          400.0f
#define FC_LP       50.0f
#define FC_HP       20.0f
#define FC1_BP      20.0f
#define FC2_BP      50.0f
#define ORDER       7
#define CHEBY1_RIPPLE_DB 3.0f
#define CHEBY2_RIPPLE_DB 40.0f
#define NUM_SAMPLES 400
#define INPUT_FREQ  30.0f

int main(void)
{
    FILE *f = fopen("test_cheby_data.csv", "w");
    if (!f) { perror("test_cheby_data.csv"); return 1; }

    /* ── Chebyshev I dynamic ──────────────────────────────────────────── */
    cheby1_t c1_lp, c1_hp, c1_bp, c1_bs;
    fprintf(stderr, "c1_lp "); fflush(stderr);
    cheby1_init(&c1_lp, FILTER_LOWPASS,  ORDER, FC_LP,  0.0f,   FS, CHEBY1_RIPPLE_DB);
    fprintf(stderr, "v=%d ", c1_lp.valid); fflush(stderr);
    fprintf(stderr, "c1_hp "); fflush(stderr);
    cheby1_init(&c1_hp, FILTER_HIGHPASS, ORDER, FC_HP,  0.0f,   FS, CHEBY1_RIPPLE_DB);
    fprintf(stderr, "v=%d ", c1_hp.valid); fflush(stderr);
    fprintf(stderr, "c1_bp "); fflush(stderr);
    cheby1_init(&c1_bp, FILTER_BANDPASS,  ORDER, FC1_BP, FC2_BP, FS, CHEBY1_RIPPLE_DB);
    fprintf(stderr, "v=%d ", c1_bp.valid); fflush(stderr);
    fprintf(stderr, "c1_bs "); fflush(stderr);
    cheby1_init(&c1_bs, FILTER_BANDSTOP,  ORDER, FC1_BP, FC2_BP, FS, CHEBY1_RIPPLE_DB);
    fprintf(stderr, "v=%d\n", c1_bs.valid); fflush(stderr);

    /* ── Chebyshev I static ───────────────────────────────────────────── */
    cheby1_lp_7th_t sc1_lp;
    cheby1_hp_7th_t sc1_hp;
    cheby1_bp_7th_t sc1_bp;
    cheby1_bs_7th_t sc1_bs;
    fprintf(stderr, "sc1_lp "); fflush(stderr);
    cheby1_lp_7th_init(&sc1_lp, FC_LP,  FS, CHEBY1_RIPPLE_DB);
    fprintf(stderr, "v=%d ", sc1_lp.valid); fflush(stderr);
    fprintf(stderr, "sc1_hp "); fflush(stderr);
    cheby1_hp_7th_init(&sc1_hp, FC_HP,  FS, CHEBY1_RIPPLE_DB);
    fprintf(stderr, "v=%d ", sc1_hp.valid); fflush(stderr);
    fprintf(stderr, "sc1_bp "); fflush(stderr);
    cheby1_bp_7th_init(&sc1_bp, FC1_BP, FC2_BP, FS, CHEBY1_RIPPLE_DB);
    fprintf(stderr, "v=%d ", sc1_bp.valid); fflush(stderr);
    fprintf(stderr, "sc1_bs "); fflush(stderr);
    cheby1_bs_7th_init(&sc1_bs, FC1_BP, FC2_BP, FS, CHEBY1_RIPPLE_DB);
    fprintf(stderr, "v=%d\n", sc1_bs.valid); fflush(stderr);

    /* ── Chebyshev II dynamic ─────────────────────────────────────────── */
    cheby2_t c2_lp, c2_hp, c2_bp, c2_bs;
    fprintf(stderr, "c2_lp "); fflush(stderr);
    cheby2_init(&c2_lp, FILTER_LOWPASS,  ORDER, FC_LP,  0.0f,   FS, CHEBY2_RIPPLE_DB);
    fprintf(stderr, "v=%d ", c2_lp.valid); fflush(stderr);
    fprintf(stderr, "c2_hp "); fflush(stderr);
    cheby2_init(&c2_hp, FILTER_HIGHPASS, ORDER, FC_HP,  0.0f,   FS, CHEBY2_RIPPLE_DB);
    fprintf(stderr, "v=%d ", c2_hp.valid); fflush(stderr);
    fprintf(stderr, "c2_bp "); fflush(stderr);
    cheby2_init(&c2_bp, FILTER_BANDPASS,  ORDER, FC1_BP, FC2_BP, FS, CHEBY2_RIPPLE_DB);
    fprintf(stderr, "v=%d ", c2_bp.valid); fflush(stderr);
    fprintf(stderr, "c2_bs "); fflush(stderr);
    cheby2_init(&c2_bs, FILTER_BANDSTOP,  ORDER, FC1_BP, FC2_BP, FS, CHEBY2_RIPPLE_DB);
    fprintf(stderr, "v=%d\n", c2_bs.valid); fflush(stderr);

    /* ── Chebyshev II static ──────────────────────────────────────────── */
    cheby2_lp_7th_t sc2_lp;
    cheby2_hp_7th_t sc2_hp;
    cheby2_bp_7th_t sc2_bp;
    cheby2_bs_7th_t sc2_bs;
    fprintf(stderr, "sc2_lp "); fflush(stderr);
    cheby2_lp_7th_init(&sc2_lp, FC_LP,  FS, CHEBY2_RIPPLE_DB);
    fprintf(stderr, "v=%d ", sc2_lp.valid); fflush(stderr);
    fprintf(stderr, "sc2_hp "); fflush(stderr);
    cheby2_hp_7th_init(&sc2_hp, FC_HP,  FS, CHEBY2_RIPPLE_DB);
    fprintf(stderr, "v=%d ", sc2_hp.valid); fflush(stderr);
    fprintf(stderr, "sc2_bp "); fflush(stderr);
    cheby2_bp_7th_init(&sc2_bp, FC1_BP, FC2_BP, FS, CHEBY2_RIPPLE_DB);
    fprintf(stderr, "v=%d ", sc2_bp.valid); fflush(stderr);
    fprintf(stderr, "sc2_bs "); fflush(stderr);
    cheby2_bs_7th_init(&sc2_bs, FC1_BP, FC2_BP, FS, CHEBY2_RIPPLE_DB);
    fprintf(stderr, "v=%d\n", sc2_bs.valid); fflush(stderr);

    /* ── Reset with first input ───────────────────────────────────────── */
    fprintf(stderr, "reset...\n"); fflush(stderr);
    float first = cosf(0.0f);

    fprintf(stderr, "r c1_lp\n"); fflush(stderr); cheby1_reset(&c1_lp, first);
    fprintf(stderr, "r c1_hp\n"); fflush(stderr); cheby1_reset(&c1_hp, first);
    fprintf(stderr, "r c1_bp\n"); fflush(stderr); cheby1_reset(&c1_bp, first);
    fprintf(stderr, "r c1_bs\n"); fflush(stderr); cheby1_reset(&c1_bs, first);

    fprintf(stderr, "r sc1\n"); fflush(stderr);
    static_cheby_reset((static_cheby_t *)&sc1_lp, first);
    static_cheby_reset((static_cheby_t *)&sc1_hp, first);
    static_cheby_reset((static_cheby_t *)&sc1_bp, first);
    static_cheby_reset((static_cheby_t *)&sc1_bs, first);

    fprintf(stderr, "r c2\n"); fflush(stderr);
    cheby2_reset(&c2_lp, first);  cheby2_reset(&c2_hp, first);
    cheby2_reset(&c2_bp, first);  cheby2_reset(&c2_bs, first);

    fprintf(stderr, "r sc2\n"); fflush(stderr);
    static_cheby_reset((static_cheby_t *)&sc2_lp, first);
    static_cheby_reset((static_cheby_t *)&sc2_hp, first);
    static_cheby_reset((static_cheby_t *)&sc2_bp, first);
    static_cheby_reset((static_cheby_t *)&sc2_bs, first);

    fprintf(stderr, "r done\n"); fflush(stderr);

    fprintf(stderr, "reset done, header...\n"); fflush(stderr);

    /* ── CSV header ───────────────────────────────────────────────────── */
    fprintf(f, "timestamp,input,"
            "cheby1_lp,cheby1_hp,cheby1_bp,cheby1_bs,"
            "cheby2_lp,cheby2_hp,cheby2_bp,cheby2_bs,"
            "static_cheby1_lp,static_cheby1_hp,static_cheby1_bp,static_cheby1_bs,"
            "static_cheby2_lp,static_cheby2_hp,static_cheby2_bp,static_cheby2_bs\n");

    /* ── Main loop ────────────────────────────────────────────────────── */
    fprintf(stderr, "header done, loop...\n"); fflush(stderr);
    for (int n = 0; n < NUM_SAMPLES; n++) {
        float t = (float)n / FS;
        float x = cosf(2.0f * (float)M_PI * INPUT_FREQ * t);

        fprintf(f, "%.6f,%.6f", t, x);

        /* Chebyshev I dynamic */
        fprintf(f, ",%.6f,%.6f,%.6f,%.6f",
                cheby1_update(&c1_lp, x), cheby1_update(&c1_hp, x),
                cheby1_update(&c1_bp, x), cheby1_update(&c1_bs, x));

        /* Chebyshev II dynamic */
        fprintf(f, ",%.6f,%.6f,%.6f,%.6f",
                cheby2_update(&c2_lp, x), cheby2_update(&c2_hp, x),
                cheby2_update(&c2_bp, x), cheby2_update(&c2_bs, x));

        /* Chebyshev I static */
        fprintf(f, ",%.6f,%.6f,%.6f,%.6f",
                static_cheby_update((static_cheby_t *)&sc1_lp, x),
                static_cheby_update((static_cheby_t *)&sc1_hp, x),
                static_cheby_update((static_cheby_t *)&sc1_bp, x),
                static_cheby_update((static_cheby_t *)&sc1_bs, x));

        /* Chebyshev II static */
        fprintf(f, ",%.6f,%.6f,%.6f,%.6f",
                static_cheby_update((static_cheby_t *)&sc2_lp, x),
                static_cheby_update((static_cheby_t *)&sc2_hp, x),
                static_cheby_update((static_cheby_t *)&sc2_bp, x),
                static_cheby_update((static_cheby_t *)&sc2_bs, x));

        fprintf(f, "\n");
    }

    fclose(f);
    cheby1_destroy(&c1_lp); cheby1_destroy(&c1_hp);
    cheby1_destroy(&c1_bp); cheby1_destroy(&c1_bs);
    cheby2_destroy(&c2_lp); cheby2_destroy(&c2_hp);
    cheby2_destroy(&c2_bp); cheby2_destroy(&c2_bs);
    return 0;
}
