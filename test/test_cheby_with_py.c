#include <stdio.h>
#include <math.h>
#include "cheby_filter.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

#define FS               400.0f
#define FC_LP            50.0f
#define FC_HP            20.0f
#define FC1_BP           20.0f
#define FC2_BP           50.0f
#define ORDER            7
#define CHEBY1_RIPPLE_DB 3.0f
#define CHEBY2_RIPPLE_DB 40.0f
#define NUM_SAMPLES      2000
#define INPUT_FREQ       30.0f

int main(void)
{
    FILE *f = fopen("test_cheby_data.csv", "w");
    if (!f) { perror("test_cheby_data.csv"); return 1; }

    /* Chebyshev I */
    cheby1_lp_7th_t c1_lp;
    cheby1_hp_7th_t c1_hp;
    cheby1_bp_7th_t c1_bp;
    cheby1_bs_7th_t c1_bs;
    cheby1_lp_7th_init(&c1_lp, FC_LP,  FS, CHEBY1_RIPPLE_DB);
    cheby1_hp_7th_init(&c1_hp, FC_HP,  FS, CHEBY1_RIPPLE_DB);
    cheby1_bp_7th_init(&c1_bp, FC1_BP, FC2_BP, FS, CHEBY1_RIPPLE_DB);
    cheby1_bs_7th_init(&c1_bs, FC1_BP, FC2_BP, FS, CHEBY1_RIPPLE_DB);

    /* Chebyshev II */
    cheby2_lp_7th_t c2_lp;
    cheby2_hp_7th_t c2_hp;
    cheby2_bp_7th_t c2_bp;
    cheby2_bs_7th_t c2_bs;
    cheby2_lp_7th_init(&c2_lp, FC_LP,  FS, CHEBY2_RIPPLE_DB);
    cheby2_hp_7th_init(&c2_hp, FC_HP,  FS, CHEBY2_RIPPLE_DB);
    cheby2_bp_7th_init(&c2_bp, FC1_BP, FC2_BP, FS, CHEBY2_RIPPLE_DB);
    cheby2_bs_7th_init(&c2_bs, FC1_BP, FC2_BP, FS, CHEBY2_RIPPLE_DB);

    float first = cosf(0.0f);
    cheby1_lp_7th_reset(&c1_lp, first); cheby1_hp_7th_reset(&c1_hp, first);
    cheby1_bp_7th_reset(&c1_bp, first); cheby1_bs_7th_reset(&c1_bs, first);
    cheby2_lp_7th_reset(&c2_lp, first); cheby2_hp_7th_reset(&c2_hp, first);
    cheby2_bp_7th_reset(&c2_bp, first); cheby2_bs_7th_reset(&c2_bs, first);

    fprintf(f, "timestamp,input,"
            "cheby1_lp,cheby1_hp,cheby1_bp,cheby1_bs,"
            "cheby2_lp,cheby2_hp,cheby2_bp,cheby2_bs\n");

    for (int n = 0; n < NUM_SAMPLES; n++) {
        float t = (float)n / FS;
        float x = cosf(2.0f * (float)M_PI * INPUT_FREQ * t);
        fprintf(f, "%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                t, x,
                cheby1_lp_7th_update(&c1_lp, x),
                cheby1_hp_7th_update(&c1_hp, x),
                cheby1_bp_7th_update(&c1_bp, x),
                cheby1_bs_7th_update(&c1_bs, x),
                cheby2_lp_7th_update(&c2_lp, x),
                cheby2_hp_7th_update(&c2_hp, x),
                cheby2_bp_7th_update(&c2_bp, x),
                cheby2_bs_7th_update(&c2_bs, x));
    }

    fclose(f);
    return 0;
}
