#include <stdio.h>
#include <math.h>
#include "butter_filter.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

#define FS          400.0f
#define FC_LP       50.0f
#define FC_HP       20.0f
#define FC1_BP      20.0f
#define FC2_BP      50.0f
#define ORDER       7
#define NUM_SAMPLES 2000
#define INPUT_FREQ  30.0f

int main(void)
{
    FILE *f = fopen("test_butter_data.csv", "w");
    if (!f) { perror("test_butter_data.csv"); return 1; }

    butter_lp_7th_t b_lp;
    butter_hp_7th_t b_hp;
    butter_bp_7th_t b_bp;
    butter_bs_7th_t b_bs;
    butter_lp_7th_init(&b_lp, FC_LP,  FS);
    butter_hp_7th_init(&b_hp, FC_HP,  FS);
    butter_bp_7th_init(&b_bp, FC1_BP, FC2_BP, FS);
    butter_bs_7th_init(&b_bs, FC1_BP, FC2_BP, FS);

    float first = cosf(0.0f);
    butter_lp_7th_reset(&b_lp, first);
    butter_hp_7th_reset(&b_hp, first);
    butter_bp_7th_reset(&b_bp, first);
    butter_bs_7th_reset(&b_bs, first);

    fprintf(f, "timestamp,input,butter_lp,butter_hp,butter_bp,butter_bs\n");

    for (int n = 0; n < NUM_SAMPLES; n++) {
        float t = (float)n / FS;
        float x = cosf(2.0f * (float)M_PI * INPUT_FREQ * t);
        fprintf(f, "%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                t, x,
                butter_lp_7th_update(&b_lp, x),
                butter_hp_7th_update(&b_hp, x),
                butter_bp_7th_update(&b_bp, x),
                butter_bs_7th_update(&b_bs, x));
    }

    fclose(f);
    return 0;
}
