/**
 * @file    test_notch_with_py.c
 * @brief   生成陷波滤波器的频响 CSV，供 Python 侧 f64 参考实现对拍。
 *
 * 与 butter/cheby 的 CSV 生成器不同，这里扫的是**频响**而不是单频时间序列：
 * 陷波器的全部信息都在 "谷底多深、多宽、位置在哪" 上，单频点看不出所以然。
 *
 * 输出列：config,f0,xi,g,freq,measured
 */

#include <stdio.h>
#include <math.h>
#include "notch_filter.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

#define FS      1000.0f
#define STEPS   12000

/* 相对 f0 的频率网格：谷底附近加密，两侧放缓 */
static const float FREQ_REL[] = {
    0.10f, 0.20f, 0.30f, 0.40f, 0.50f, 0.60f, 0.70f, 0.80f,
    0.90f, 0.95f, 0.98f, 0.99f, 1.00f, 1.01f, 1.02f, 1.05f,
    1.10f, 1.20f, 1.30f, 1.50f, 1.70f, 2.00f, 2.50f, 3.00f, 4.00f
};

static const float CFG[][3] = {
    /* f0,   xi,    g     */
    {  50.0f, 0.05f, 0.10f },
    { 100.0f, 0.02f, 0.05f },
    { 200.0f, 0.10f, 0.02f },
    {  30.0f, 0.05f, 0.20f },
};

/**
 * @brief 喂 freq 正弦，测稳态增益 |H(freq)|。
 *
 * 相位累加器 + **RMS 比值**，理由见 test_notch.c 中同名函数的注释：
 * 峰值在低采样密度下被低估（f=200Hz@1kHz 只有 5 点/周期 → 峰值 0.951），
 * 会凭空造出 5% 的假偏差；RMS 比值对采样密度免疫。
 *
 * @param[in] f      滤波器对象指针。
 * @param[in] freq   测试频率（Hz）。
 * @return           后四分之一样本上的 |H(freq)|。
 */
static float measure(const notch_filter_t *f, float freq)
{
    notch_filter_t local = *f; /* 每个频点从零状态起跑，避免上一段的残余 */
    notch_reset(&local, 0.0f);

    const float two_pi = 2.0f * (float)M_PI;
    const float dphase = two_pi * freq / FS;
    float phase = 0.0f;
    float sum_x2 = 0.0f, sum_y2 = 0.0f; /* 零 double；精度余量见 test_notch.c */

    for (int n = 0; n < STEPS; n++) {
        float x = sinf(phase);
        float y = notch_update(&local, x);
        phase += dphase;
        if (phase >= two_pi) phase -= two_pi;
        if (n >= STEPS - STEPS / 4) {
            sum_x2 += x * x;
            sum_y2 += y * y;
        }
    }
    if (sum_x2 <= 1e-6f) return 1.0f;
    return sqrtf(sum_y2 / sum_x2);
}

int main(int argc, char **argv)
{
    /* 输出路径经 argv[1] 传入，让 ctest 把 CSV 固定到 build 目录
       （Python 对比脚本从那里读取，绝不会读到源码目录里的陈旧副本）。 */
    const char *path = (argc > 1) ? argv[1] : "test_notch_data.csv";
    FILE *out = fopen(path, "w");
    if (!out) { perror(path); return 1; }

    const unsigned ncfg = sizeof(CFG) / sizeof(CFG[0]);
    const unsigned nfrq = sizeof(FREQ_REL) / sizeof(FREQ_REL[0]);

    notch_filter_t nf[sizeof(CFG) / sizeof(CFG[0])];
    for (unsigned c = 0; c < ncfg; c++) {
        nf[c].valid = 0;
        notch_init(&nf[c], CFG[c][0], CFG[c][1], CFG[c][2], FS);
        /* 任何一条设计被闸拒掉都不该吐出看似正常的 CSV——否则 Python 侧
           会拿它去校验一个根本不存在的设计。 */
        if (!nf[c].valid) {
            fprintf(stderr, "notch config %u (f0=%g xi=%g g=%g) failed (valid=0)"
                            " — refusing to emit CSV\n",
                    c, CFG[c][0], CFG[c][1], CFG[c][2]);
            fclose(out);
            return 1;
        }
    }

    fprintf(out, "config,f0,xi,g,freq,measured\n");
    for (unsigned c = 0; c < ncfg; c++) {
        for (unsigned i = 0; i < nfrq; i++) {
            float fr = FREQ_REL[i] * CFG[c][0];
            fprintf(out, "%u,%.6f,%.6f,%.6f,%.8f,%.8f\n",
                    c, CFG[c][0], CFG[c][1], CFG[c][2], fr,
                    measure(&nf[c], fr));
        }
    }

    fclose(out);
    return 0;
}
