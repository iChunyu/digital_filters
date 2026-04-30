#include <stdio.h>
#include <math.h>
#include "butter_filter.h"
#include "static_butter_filter.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

#define FS       400.0f
#define FC_LP    50.0f
#define FC_HP    20.0f
#define FC1_BP   20.0f
#define FC2_BP   50.0f
#define ORDER    7
#define NUM_SAMPLES 400
#define INPUT_FREQ  30.0f

int main(void)
{
    FILE *f = fopen("test_butter_data.csv", "w");
    if (!f) { perror("test_butter_data.csv"); return 1; }

    /* ── Dynamic API ──────────────────────────────────────────────────── */
    butter_t b_lp, b_hp, b_bp, b_bs;
    butter_init(&b_lp, FILTER_LOWPASS,  ORDER, FC_LP,  0.0f,   FS);
    butter_init(&b_hp, FILTER_HIGHPASS, ORDER, FC_HP,  0.0f,   FS);
    butter_init(&b_bp, FILTER_BANDPASS,  ORDER, FC1_BP, FC2_BP, FS);
    butter_init(&b_bs, FILTER_BANDSTOP,  ORDER, FC1_BP, FC2_BP, FS);

    /* ── Static API ───────────────────────────────────────────────────── */
    butter_lp_7th_t sb_lp;
    butter_hp_7th_t sb_hp;
    butter_bp_7th_t sb_bp;
    butter_bs_7th_t sb_bs;
    butter_lp_7th_init(&sb_lp, FC_LP,  FS);
    butter_hp_7th_init(&sb_hp, FC_HP,  FS);
    butter_bp_7th_init(&sb_bp, FC1_BP, FC2_BP, FS);
    butter_bs_7th_init(&sb_bs, FC1_BP, FC2_BP, FS);

    /* ── Reset with first input ───────────────────────────────────────── */
    float first = cosf(0.0f);

    butter_reset(&b_lp, first);  butter_reset(&b_hp, first);
    butter_reset(&b_bp, first);  butter_reset(&b_bs, first);
    static_butter_reset((static_butter_t *)&sb_lp, first);
    static_butter_reset((static_butter_t *)&sb_hp, first);
    static_butter_reset((static_butter_t *)&sb_bp, first);
    static_butter_reset((static_butter_t *)&sb_bs, first);

    /* ── CSV header ───────────────────────────────────────────────────── */
    fprintf(f, "timestamp,input,"
            "butter_lp,butter_hp,butter_bp,butter_bs,"
            "static_butter_lp,static_butter_hp,static_butter_bp,static_butter_bs\n");

    /* ── Main loop ────────────────────────────────────────────────────── */
    for (int n = 0; n < NUM_SAMPLES; n++) {
        float t = (float)n / FS;
        float x = cosf(2.0f * (float)M_PI * INPUT_FREQ * t);

        fprintf(f, "%.6f,%.6f", t, x);

        fprintf(f, ",%.6f,%.6f,%.6f,%.6f",
                butter_update(&b_lp, x), butter_update(&b_hp, x),
                butter_update(&b_bp, x), butter_update(&b_bs, x));

        fprintf(f, ",%.6f,%.6f,%.6f,%.6f",
                static_butter_update((static_butter_t *)&sb_lp, x),
                static_butter_update((static_butter_t *)&sb_hp, x),
                static_butter_update((static_butter_t *)&sb_bp, x),
                static_butter_update((static_butter_t *)&sb_bs, x));

        fprintf(f, "\n");
    }

    fclose(f);
    butter_destroy(&b_lp);  butter_destroy(&b_hp);
    butter_destroy(&b_bp);  butter_destroy(&b_bs);
    return 0;
}
