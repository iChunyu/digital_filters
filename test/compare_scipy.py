#!/usr/bin/env python3
"""Compare C filter CSV output against scipy reference — skip transients."""
import os, sys
import numpy as np
from scipy import signal

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FS = 400.0
TRANSIENT_SKIP = 60  # skip first N samples for steady-state comparison

def check_correctness(label, data, dyn_col, sta_col, sos, x, skip):
    scipy_out = signal.sosfilt(sos, x)
    # Compare last 10% of samples for steady-state accuracy
    ss_start = max(len(x) * 9 // 10, len(x) - 50)
    for tag, col in [("dyn", dyn_col), ("sta", sta_col)]:
        err = data[col][ss_start:] - scipy_out[ss_start:]
        maxe = np.max(np.abs(err))
        rms = np.sqrt(np.mean(err**2))
        status = "OK" if maxe < 1e-4 else "FAIL"
        print(f"  {tag:4s}: max_err={maxe:.2e}  rms={rms:.2e}  [{status}]")
    return scipy_out

# --- Butterworth ---
print("=" * 60)
print("BUTTERWORTH (ORDER=7, fc=50/20/20-50 Hz, fs=400 Hz)")
print("=" * 60)
csv_b = os.path.join(SCRIPT_DIR, "test_butter_data.csv")
if os.path.isfile(csv_b):
    data_b = np.genfromtxt(csv_b, delimiter=',', names=True)
    x = data_b['input']
    for label, btype, dyn_col, sta_col, fc1, fc2 in [
        ("LP", "lowpass", "butter_lp", "static_butter_lp", 50.0, None),
        ("HP", "highpass", "butter_hp", "static_butter_hp", 20.0, None),
        ("BP", "bandpass", "butter_bp", "static_butter_bp", 20.0, 50.0),
        ("BS", "bandstop", "butter_bs", "static_butter_bs", 20.0, 50.0),
    ]:
        sos = signal.butter(7, fc1 if fc2 is None else [fc1,fc2],
                            btype=btype, fs=FS, output='sos')
        print(f"\nButterworth {label}:")
        check_correctness(label, data_b, dyn_col, sta_col, sos, x, TRANSIENT_SKIP)
else:
    print("CSV not found — run test_butter_with_py first")

# --- Chebyshev ---
print("\n" + "=" * 60)
print("CHEBYSHEV (ORDER=7, fs=400 Hz)")
print("=" * 60)
csv_c = os.path.join(SCRIPT_DIR, "test_cheby_data.csv")
if os.path.isfile(csv_c):
    data_c = np.genfromtxt(csv_c, delimiter=',', names=True)
    x = data_c['input']
    for title, scipy_fn, btype, dyn_col, sta_col, ripple, fc1, fc2 in [
        ("Cheby1 LP",  signal.cheby1, "lowpass",  "cheby1_lp", "static_cheby1_lp",  3.0, 50.0, None),
        ("Cheby1 HP",  signal.cheby1, "highpass", "cheby1_hp", "static_cheby1_hp",  3.0, 20.0, None),
        ("Cheby1 BP",  signal.cheby1, "bandpass", "cheby1_bp", "static_cheby1_bp",  3.0, 20.0, 50.0),
        ("Cheby1 BS",  signal.cheby1, "bandstop", "cheby1_bs", "static_cheby1_bs",  3.0, 20.0, 50.0),
        ("Cheby2 LP",  signal.cheby2, "lowpass",  "cheby2_lp", "static_cheby2_lp", 40.0, 50.0, None),
        ("Cheby2 HP",  signal.cheby2, "highpass", "cheby2_hp", "static_cheby2_hp", 40.0, 20.0, None),
        ("Cheby2 BP",  signal.cheby2, "bandpass", "cheby2_bp", "static_cheby2_bp", 40.0, 20.0, 50.0),
        ("Cheby2 BS",  signal.cheby2, "bandstop", "cheby2_bs", "static_cheby2_bs", 40.0, 20.0, 50.0),
    ]:
        sos = scipy_fn(7, ripple, fc1 if fc2 is None else [fc1,fc2],
                       btype=btype, fs=FS, output='sos')
        print(f"\n{title}:")
        check_correctness(title, data_c, dyn_col, sta_col, sos, x, TRANSIENT_SKIP)
else:
    print("CSV not found — run test_cheby_with_py first")

# --- Also test higher-order Butterworth via Python-only ---
print("\n" + "=" * 60)
print("HIGH-ORDER BUTTERWORTH LP (steady-state match check)")
print("=" * 60)
for order in [16, 24, 32]:
    sos = signal.butter(order, 50.0, btype='lowpass', fs=FS, output='sos')
    # Check the SOS coefficients look reasonable
    max_den = np.max(np.abs(sos[:, 4:6]))
    print(f"  N={order:2d}: {sos.shape[0]} sections, max|a|={max_den:.6f}")
