#!/usr/bin/env python3
"""峰值滤波器的 f64 参考实现与对比。

独立复现整条设计链，与被测的 C/f32 实现逐频点对拍：

    模拟原型 → 预畸 → scipy.signal.bilinear → 频响求值

被测原型（与 C 侧同一套系数）：

    H(s) = (s² + 2ξω₀ s + ω₀²) / (s² + 2ξ(1/g) ω₀ s + ω₀²)

scipy 的 bilinear 是独立实现（已确认与本库 biquad_c2d_bilinear 逐位吻合，
输入按**降幂** [s^2, s^1, s^0] 给出）。C 侧的时域递推则是 f32，
两边能对上说明设计管线与部署都对。

numpy/scipy 缺失时以 77 退出（由 ctest 的 SKIP_RETURN_CODE 接住）。

用法：test_peak_use_py.py <数据目录>   （目录里须有 test_peak_data.csv）
"""

import csv
import math
import os
import sys

# 相对容差：与陷波器不同，峰值高到 50，绝对容差没有可比性，故按
# max(|H_ref|) 归一。实测四个配置的最坏相对偏差 2.19e-3（出现在
# f0=100、xi=0.02、g=20 的峰顶——等效陷波深度 d=0.05，是四个配置里
# 数值最难的；深陷波的病态在这里以峰高的形式出现，与 notch 镜像），
# 故 0.01 留了约 4.6 倍余量，同时仍是有牙齿的上限：把 num/den 传反、
# 或用 g 而不是 1/g 都会立刻顶穿它。
TOL_REL = 0.01
# 峰高容差：f/f0 == 1.0 那一点必须落回 g 的 ±2% 内
TOL_PEAK = 0.02


def reference_response(f0, xi, g, fs, freqs):
    """返回每个 freq 处的 |H(f)|。

    与 C 侧同一套参数语义：ω₀ 预畸到数字频率 f0，峰高 g，带宽因子 xi。
    """
    import numpy as np
    from scipy.signal import bilinear

    # 预畸：ω0 = 2*fs*tan(pi*f0/fs)，与 filter_utils.c 的 prewarp 一致
    w0 = 2.0 * fs * math.tan(math.pi * f0 / fs)

    # bilinear 取降幂系数；等效陷波深度 d = 1/g
    d = 1.0 / g
    num_s = [1.0, 2.0 * xi * w0, w0 ** 2]
    den_s = [1.0, 2.0 * xi * d * w0, w0 ** 2]
    b, a = bilinear(num_s, den_s, fs)

    z = np.exp(-2j * np.pi * np.asarray(freqs, dtype=float) / fs)
    num = b[0] + b[1] * z + b[2] * z * z
    den = a[0] + a[1] * z + a[2] * z * z
    return np.abs(num / den)


def main():
    data_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    csv_path = os.path.join(data_dir, "test_peak_data.csv")

    try:
        import numpy as np
    except ImportError:
        print("numpy 不可用，跳过 peak 对比")
        return 77
    try:
        from scipy.signal import bilinear  # noqa: F401
    except ImportError:
        print("scipy 不可用，跳过 peak 对比")
        return 77

    if not os.path.exists(csv_path):
        print(f"缺少 {csv_path}——peak_csv_gen 没跑过？")
        return 1

    rows = {}
    with open(csv_path, newline="") as fh:
        for r in csv.DictReader(fh):
            key = (float(r["f0"]), float(r["xi"]), float(r["g"]))
            rows.setdefault(key, []).append(
                (float(r["freq"]), float(r["measured"]))
            )

    if not rows:
        print("CSV 为空")
        return 1

    fs = 1000.0  # 与 test_peak_with_py.c 的 FS 一致
    failures = 0

    print(f"{'f0':>7} {'xi':>6} {'g':>6} {'最大相对偏差':>13} {'(在 f/f0)':>10} "
          f"{'峰顶 |H(f0)|':>13} {'目标 g':>8}")
    for (f0, xi, g), pts in sorted(rows.items()):
        pts.sort()
        freqs = [p[0] for p in pts]
        measured = np.array([p[1] for p in pts])
        ref = reference_response(f0, xi, g, fs, freqs)

        rel = np.abs(measured - ref) / ref.max()
        k = int(np.argmax(rel))
        worst = rel[k]
        worst_rel = freqs[k] / f0

        # 峰高：网格里有 f/f0 == 1.0 这一点
        at_f0 = measured[[i for i, f in enumerate(freqs) if abs(f / f0 - 1.0) < 1e-9][0]]

        if worst > TOL_REL:
            failures += 1
            print(f"  FAIL: 最大相对偏差 {worst:.5f} > {TOL_REL}")
        if abs(at_f0 - g) > TOL_PEAK * g:
            failures += 1
            print(f"  FAIL: 峰高 {at_f0:.5f} 偏离目标 g={g} 超过 {TOL_PEAK:.0%}")

        print(f"{f0:>7g} {xi:>6g} {g:>6g} {worst:>13.5f} {worst_rel:>10.4f} "
              f"{at_f0:>13.6f} {g:>8g}")

    if failures:
        print(f"{failures} 项 peak 对比失败")
        return 1
    print("peak 频响与 f64 参考实现一致")
    return 0


if __name__ == "__main__":
    sys.exit(main())
