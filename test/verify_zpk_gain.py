"""
verify_zpk_gain.py — 验证 zpk 变换中增益 k 对 float32 精度的需求。

模拟 scipy 的完整 zpk 设计管线，分别用 float64 和 float32 执行，
对比各阶段增益 k 的相对误差、最终 SOS 系数差异和频率响应差异。

覆盖: butter, cheby1, cheby2 × LP, HP, BP, BS
"""

import sys
import math
import cmath
import numpy as np

# =====================================================================
#  scipy 管线复现（独立实现，不依赖 scipy.signal）
# =====================================================================

def _relative_degree(z, p):
    return len(p) - len(z)


def buttap(N, dtype=np.float64):
    """Butterworth 模拟原型: 返回 (z, p, k)."""
    z = np.array([], dtype=dtype)
    m = np.arange(-N + 1, N, 2, dtype=dtype)
    p = -np.exp(1j * np.pi * m / (2 * N))
    return z, p.astype(np.complex128), dtype(1.0)


def cheb1ap(N, rp, dtype=np.float64):
    """Chebyshev I 模拟原型: 返回 (z, p, k)."""
    if N == 0:
        return (np.array([], dtype=dtype),
                np.array([], dtype=np.complex128),
                dtype(10 ** (-rp / 20)))
    z = np.array([], dtype=dtype)
    eps = math.sqrt(10 ** (0.1 * rp) - 1.0)
    mu = math.asinh(1.0 / eps) / N
    m = np.arange(-N + 1, N, 2, dtype=dtype)
    theta = np.pi * m / (2 * N)
    p = -np.sinh(mu + 1j * theta)
    k = np.real(np.prod(-p))
    if N % 2 == 0:
        k = k / math.sqrt(1 + eps * eps)
    return z, p.astype(np.complex128), dtype(k)


def cheb2ap(N, rs, dtype=np.float64):
    """Chebyshev II 模拟原型: 返回 (z, p, k)."""
    if N == 0:
        return (np.array([], dtype=np.complex128),
                np.array([], dtype=np.complex128),
                dtype(1.0))
    de = 1.0 / math.sqrt(10 ** (0.1 * rs) - 1)
    mu = math.asinh(1.0 / de) / N
    if N % 2:
        m = np.concatenate(
            (np.arange(-N + 1, 0, 2, dtype=dtype),
             np.arange(2, N, 2, dtype=dtype))
        )
    else:
        m = np.arange(-N + 1, N, 2, dtype=dtype)
    z = 1j / np.sin(m * np.pi / (2 * N))
    m1 = np.arange(-N + 1, N, 2, dtype=dtype)
    theta1 = np.pi * m1 / (2 * N)
    p = -1.0 / np.sinh(mu + 1j * theta1)
    k = np.real(np.prod(-p) / np.prod(-z))
    return z.astype(np.complex128), p.astype(np.complex128), dtype(k)


def prewarp(fd, fs):
    """双线性频率预畸变: f_analog = fs/pi * tan(pi * fd / fs)."""
    return fs / np.pi * np.tan(np.pi * fd / fs)


def lp2lp_zpk(z, p, k, wo):
    """LP → LP 频率变换 (zpk 域)."""
    degree = _relative_degree(z, p)
    z_lp = wo * z
    p_lp = wo * p
    k_lp = k * (wo ** degree)
    return z_lp, p_lp, k_lp


def lp2hp_zpk(z, p, k, wo):
    """LP → HP 频率变换 (zpk 域)."""
    degree = _relative_degree(z, p)
    z_hp = wo / z
    p_hp = wo / p
    z_hp = np.concatenate((z_hp, np.zeros(degree, dtype=z_hp.dtype)))
    with np.errstate(invalid='ignore'):
        k_hp = k * np.real(np.prod(-z) / np.prod(-p))
    return z_hp, p_hp, k_hp


def lp2bp_zpk(z, p, k, wo, bw):
    """LP → BP 频率变换 (zpk 域)."""
    degree = _relative_degree(z, p)
    z_lp = z * bw / 2
    p_lp = p * bw / 2
    z_lp = z_lp.astype(np.complex128)
    p_lp = p_lp.astype(np.complex128)
    z_bp = np.concatenate(
        (z_lp + np.sqrt(z_lp ** 2 - wo ** 2),
         z_lp - np.sqrt(z_lp ** 2 - wo ** 2))
    )
    p_bp = np.concatenate(
        (p_lp + np.sqrt(p_lp ** 2 - wo ** 2),
         p_lp - np.sqrt(p_lp ** 2 - wo ** 2))
    )
    z_bp = np.concatenate((z_bp, np.zeros(degree, dtype=z_bp.dtype)))
    k_bp = k * (bw ** degree)
    return z_bp, p_bp, k_bp


def lp2bs_zpk(z, p, k, wo, bw):
    """LP → BS 频率变换 (zpk 域)."""
    degree = _relative_degree(z, p)
    z_hp = (bw / 2) / z
    p_hp = (bw / 2) / p
    z_hp = z_hp.astype(np.complex128)
    p_hp = p_hp.astype(np.complex128)
    z_bs = np.concatenate(
        (z_hp + np.sqrt(z_hp ** 2 - wo ** 2),
         z_hp - np.sqrt(z_hp ** 2 - wo ** 2))
    )
    p_bs = np.concatenate(
        (p_hp + np.sqrt(p_hp ** 2 - wo ** 2),
         p_hp - np.sqrt(p_hp ** 2 - wo ** 2))
    )
    z_bs = np.concatenate(
        (z_bs, np.full(degree, +1j * wo, dtype=np.complex128))
    )
    z_bs = np.concatenate(
        (z_bs, np.full(degree, -1j * wo, dtype=np.complex128))
    )
    with np.errstate(invalid='ignore'):
        k_bs = k * np.real(np.prod(-z) / np.prod(-p))
    return z_bs, p_bs, k_bs


def bilinear_zpk(z, p, k, fs):
    """双线性变换 (zpk 域)."""
    degree = _relative_degree(z, p)
    fs2 = 2.0 * fs
    z_z = (fs2 + z) / (fs2 - z)
    p_z = (fs2 + p) / (fs2 - p)
    z_z = np.concatenate((z_z, -np.ones(degree, dtype=z_z.dtype)))
    k_z = k * np.real(np.prod(fs2 - z) / np.prod(fs2 - p))
    return z_z, p_z, k_z


# =====================================================================
#  zpk2sos — 模拟 scipy 实现
# =====================================================================

def _cplxreal(z, tol=1e-12):
    """将数组拆分为实数和共轭对（返回上半平面元素）。"""
    # 简化版: 返回所有非实数部分的"代表"
    z = np.asarray(z)
    if len(z) == 0:
        return np.array([]), np.array([])
    real_mask = np.abs(z.imag) <= tol * np.abs(z)
    real = z[real_mask].real
    # 取虚部为正的复数
    complex_mask = ~real_mask & (z.imag > 0)
    complex_pairs = z[complex_mask]
    return real, complex_pairs


def _nearest_real_complex_idx(z, p1, kind='any'):
    """找到距离 p1 最近的零点索引。"""
    z = np.asarray(z)
    if kind == 'real':
        mask = np.abs(z.imag) < 1e-12
    elif kind == 'complex':
        mask = np.abs(z.imag) >= 1e-12
    else:
        mask = np.ones(len(z), dtype=bool)
    if not np.any(mask):
        return 0
    candidates = np.where(mask)[0]
    dists = np.abs(z[candidates] - p1)
    return candidates[np.argmin(dists)]


def zpk2tf(z, p, k):
    """zpk → 传递函数系数。"""
    z = np.asarray(z)
    p = np.asarray(p)
    if len(z) == 0:
        b = np.array([k])
    else:
        b = k * np.poly(z)
    if len(p) == 0:
        a = np.array([1.0])
    else:
        a = np.poly(p)
    return b.real, a.real


def zpk2sos(z, p, k):
    """scipy 风格的 zpk → SOS 转换。"""
    z = np.asarray(z, dtype=np.complex128)
    p = np.asarray(p, dtype=np.complex128)

    if len(z) == 0 and len(p) == 0:
        return np.array([[k, 0., 0., 1., 0., 0.]])

    # 补齐到相同数量
    if len(p) < len(z):
        p = np.concatenate((p, np.zeros(len(z) - len(p))))
    elif len(z) < len(p):
        z = np.concatenate((z, np.zeros(len(p) - len(z))))

    n_sections = (max(len(p), len(z)) + 1) // 2
    if len(p) % 2 == 1:
        p = np.concatenate((p, [0. + 0.j]))
        z = np.concatenate((z, [0. + 0.j]))

    z = np.concatenate(_cplxreal(z))
    p = np.concatenate(_cplxreal(p))
    k = float(np.real(k))

    sos = np.zeros((n_sections, 6))

    def idx_worst(p_arr):
        return np.argmin(np.abs(1 - np.abs(p_arr)))

    for si in range(n_sections - 1, -1, -1):
        p1_idx = idx_worst(p)
        p1 = p[p1_idx]
        p = np.delete(p, p1_idx)

        if np.isreal(p1) and np.sum(np.isreal(p)) == 0:
            if len(z) > 0:
                z1_idx = _nearest_real_complex_idx(z, p1, 'real')
                z1 = z[z1_idx]
                z = np.delete(z, z1_idx)
                b, a = zpk2tf([z1, 0], [p1, 0], 1)
            else:
                b, a = zpk2tf([], [p1], 1)
        else:
            if np.isreal(p1):
                prealidx = np.flatnonzero(np.isreal(p))
                p2_idx = prealidx[idx_worst(p[prealidx])]
                p2 = p[p2_idx]
                p = np.delete(p, p2_idx)
            else:
                p2 = np.conj(p1)

            if len(z) > 0:
                z1_idx = _nearest_real_complex_idx(z, p1, 'any')
                z1 = z[z1_idx]
                z = np.delete(z, z1_idx)
                if not np.isreal(z1):
                    b, a = zpk2tf([z1, np.conj(z1)], [p1, p2], 1)
                else:
                    if len(z) > 0:
                        z2_idx = _nearest_real_complex_idx(z, p1, 'real')
                        z2 = z[z2_idx]
                        z = np.delete(z, z2_idx)
                        b, a = zpk2tf([z1, z2], [p1, p2], 1)
                    else:
                        b, a = zpk2tf([z1], [p1, p2], 1)
            else:
                b, a = zpk2tf([], [p1, p2], 1)

        # 填充到 sos 数组
        pad_b = 3 - len(b)
        pad_a = 3 - len(a)
        sos[si, :3] = np.concatenate((np.zeros(pad_b), b))
        sos[si, 3:] = np.concatenate((np.zeros(pad_a), a))

    # 全局增益 k 乘到第一节分子
    sos[0, :3] *= k
    return sos


# =====================================================================
#  完整设计管线
# =====================================================================

def design_filter(ftype, btype, N, fc, fs, ripple=1.0, dtype=np.float64):
    """
    运行完整 zpk 设计管线，返回各阶段状态。

    返回 dict:
        proto_k, after_analog_k, after_bilinear_k, final_sos
    """
    # 前置: 数字频率预畸变
    if btype in ('lp', 'hp'):
        wc = 2 * np.pi * prewarp(fc, fs)
        w0_norm = 0.0 if btype == 'lp' else np.pi
    else:
        fc1, fc2 = fc
        wc1 = 2 * np.pi * prewarp(fc1, fs)
        wc2 = 2 * np.pi * prewarp(fc2, fs)
        w0_analog = np.sqrt(wc1 * wc2)
        xi = wc2 - wc1
        if btype == 'bp':
            w0_norm = 2 * np.pi * np.sqrt(fc1 * fc2) / fs
        else:
            w0_norm = 0.0

    # 1. 原型
    if ftype == 'butter':
        z, p, k = buttap(N, dtype=dtype)
    elif ftype == 'cheby1':
        z, p, k = cheb1ap(N, ripple, dtype=dtype)
    elif ftype == 'cheby2':
        z, p, k = cheb2ap(N, ripple, dtype=dtype)
    else:
        raise ValueError(f"Unknown ftype: {ftype}")

    proto_k = dtype(k)

    # 2. 模拟频率变换 (需要保存原始 z,p 用于 HP/BS 的增益计算)
    z_orig = z.copy() if len(z) > 0 else z
    p_orig = p.copy()

    if btype == 'lp':
        z, p, k = lp2lp_zpk(z, p, k, wc)
    elif btype == 'hp':
        z, p, k = lp2hp_zpk(z, p, k, wc)
    elif btype == 'bp':
        z, p, k = lp2bp_zpk(z, p, k, w0_analog, xi)
    elif btype == 'bs':
        z, p, k = lp2bs_zpk(z, p, k, w0_analog, xi)

    after_analog_k = dtype(k)

    # 3. 双线性变换 (保存 s 域零极点用于增益计算)
    z_analog = z.copy()
    p_analog = p.copy()
    z, p, k = bilinear_zpk(z, p, k, fs)

    after_bilinear_k = dtype(k)

    # 4. zpk → sos
    sos = zpk2sos(z, p, k)

    return {
        'proto_k': proto_k,
        'after_analog_k': after_analog_k,
        'after_bilinear_k': after_bilinear_k,
        'sos': sos,
        'z': z, 'p': p,
        'w0_norm': w0_norm,
    }


# =====================================================================
#  频率响应
# =====================================================================

def sos_freqz(sos, worN=512, whole=False):
    """计算 SOS 滤波器的频率响应。"""
    worN = np.asarray(worN)
    if worN.ndim == 0:
        w = np.linspace(0, np.pi if not whole else 2 * np.pi, worN)
    else:
        w = worN

    h = np.ones(len(w), dtype=np.complex128)
    for section in sos:
        b0, b1, b2, a0, a1, a2 = section
        # a0 应归一化为 1
        if a0 != 1.0:
            b0, b1, b2 = b0 / a0, b1 / a0, b2 / a0
            a1, a2 = a1 / a0, a2 / a0
        z = np.exp(-1j * w)
        h_section = (b0 + b1 * z + b2 * z ** 2) / (1 + a1 * z + a2 * z ** 2)
        h *= h_section

    return w, h


# =====================================================================
#  主验证逻辑
# =====================================================================

def run_verification():
    """对所有滤波器类型运行精度对比。"""
    fs = 1000.0
    test_cases = [
        # (ftype, btype, order, fc, ripple)
        ('butter',  'lp', 4, 100.0,    1.0),
        ('butter',  'hp', 4, 100.0,    1.0),
        ('butter',  'bp', 4, (50.0, 150.0), 1.0),
        ('butter',  'bs', 4, (50.0, 150.0), 1.0),
        ('cheby1',  'lp', 4, 100.0,    1.0),
        ('cheby1',  'hp', 4, 100.0,    1.0),
        ('cheby1',  'bp', 4, (50.0, 150.0), 1.0),
        ('cheby1',  'bs', 4, (50.0, 150.0), 1.0),
        ('cheby2',  'lp', 4, 100.0,    40.0),
        ('cheby2',  'hp', 4, 100.0,    40.0),
        ('cheby2',  'bp', 4, (50.0, 150.0), 40.0),
        ('cheby2',  'bs', 4, (50.0, 150.0), 40.0),
        # 高 order 测试
        ('cheby1',  'lp', 8, 100.0,    1.0),
        ('cheby2',  'hp', 8, 100.0,    40.0),
        ('butter',  'bp', 8, (50.0, 150.0), 1.0),
    ]

    all_pass = True

    for ftype, btype, order, fc, ripple in test_cases:
        print(f"\n{'='*60}")
        fc_str = f"{fc}" if isinstance(fc, (int, float)) else f"[{fc[0]}, {fc[1]}]"
        print(f"{ftype:8s} | {btype:3s} | N={order} | fc={fc_str} | ripple={ripple}")

        # float64 参考
        ref = design_filter(ftype, btype, order, fc, fs, ripple, dtype=np.float64)

        # float32 模拟
        f32 = design_filter(ftype, btype, order, fc, fs, ripple, dtype=np.float32)

        # --- 增益对比 ---
        k64 = ref['after_bilinear_k']
        k32 = float(f32['after_bilinear_k'])
        k_relerr = abs(k64 - k32) / max(abs(k64), 1e-30) if abs(k64) > 1e-30 else abs(k64 - k32)

        print(f"  proto_k:     f64={ref['proto_k']:.6e}   f32={float(f32['proto_k']):.6e}")
        print(f"  analog_k:    f64={ref['after_analog_k']:.6e}   f32={float(f32['after_analog_k']):.6e}")
        print(f"  bilinear_k:  f64={k64:.6e}   f32={k32:.6e}")
        print(f"  k relerr:    {k_relerr:.2e}")

        # --- SOS 系数对比 ---
        sos64 = np.asarray(ref['sos'], dtype=np.float64)
        sos32 = np.asarray(f32['sos'], dtype=np.float64)
        sos_max_diff = np.max(np.abs(sos64 - sos32))
        sos_relerr = sos_max_diff / max(np.max(np.abs(sos64)), 1e-30)

        print(f"  sos max diff:  {sos_max_diff:.2e}")
        print(f"  sos relerr:    {sos_relerr:.2e}")

        # --- 频率响应对比 ---
        w, h64 = sos_freqz(sos64)
        _, h32 = sos_freqz(sos32)
        h_diff = np.abs(h64 - h32)
        h_max_diff = np.max(h_diff)
        h_max_relerr = np.max(h_diff / np.maximum(np.abs(h64), 1e-15))

        db64 = 20 * np.log10(np.maximum(np.abs(h64), 1e-15))
        db32 = 20 * np.log10(np.maximum(np.abs(h32), 1e-15))
        db_max_diff = np.max(np.abs(db64 - db32))

        print(f"  freq resp max diff:     {h_max_diff:.2e}")
        print(f"  freq resp max relerr:   {h_max_relerr:.2e}")
        print(f"  freq resp dB max diff:  {db_max_diff:.4f} dB")

        # --- 判据 ---
        passes = True
        if k_relerr > 1e-5 and abs(k64) > 1e-10:
            print(f"  ⚠ WARNING: k relative error > 1e-5")
            passes = False
        if db_max_diff > 0.1:
            print(f"  ⚠ WARNING: dB response diff > 0.1 dB")
            passes = False
        if sos_relerr > 1e-4:
            print(f"  ⚠ WARNING: SOS coefficient relative error > 1e-4")
            passes = False

        if passes:
            print(f"  ✅ PASS")
        else:
            print(f"  ❌ FAIL")
            all_pass = False

    print(f"\n{'='*60}")
    if all_pass:
        print("✅ ALL TESTS PASS — float32 is adequate for zpk gain tracking.")
    else:
        print("❌ SOME TESTS FAILED — may need double precision for gain.")

    return 0 if all_pass else 1


if __name__ == '__main__':
    sys.exit(run_verification())
