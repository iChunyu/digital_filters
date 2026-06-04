#!/usr/bin/env python3
"""Compare C filter CSV output against scipy reference — skip transients."""
import os
import numpy as np
from scipy import signal

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FS = 400.0

def check_correctness(label, data, col, sos, x):
    scipy_out = signal.sosfilt(sos, x)
    # Compare last 10% of samples for steady-state accuracy
    ss_start = max(len(x) * 9 // 10, len(x) - 50)
    err = data[col][ss_start:] - scipy_out[ss_start:]
    maxe = np.max(np.abs(err))
    rms = np.sqrt(np.mean(err**2))
    status = "OK" if maxe < 1e-4 else "FAIL"
    print(f"  {label:8s}: max_err={maxe:.2e}  rms={rms:.2e}  [{status}]")
    return scipy_out

# --- Butterworth ---
print("=" * 60)
print("BUTTERWORTH (ORDER=7, fc=50/20/20-50 Hz, fs=400 Hz)")
print("=" * 60)
csv_b = os.path.join(SCRIPT_DIR, "test_butter_data.csv")
if os.path.isfile(csv_b):
    data_b = np.genfromtxt(csv_b, delimiter=',', names=True)
    x = data_b['input']
    for label, btype, col, fc1, fc2 in [
        ("LP", "lowpass", "butter_lp", 50.0, None),
        ("HP", "highpass", "butter_hp", 20.0, None),
        ("BP", "bandpass", "butter_bp", 20.0, 50.0),
        ("BS", "bandstop", "butter_bs", 20.0, 50.0),
    ]:
        sos = signal.butter(7, fc1 if fc2 is None else [fc1,fc2],
                            btype=btype, fs=FS, output='sos')
        check_correctness(label, data_b, col, sos, x)
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
    for title, scipy_fn, btype, col, ripple, fc1, fc2 in [
        ("Cheby1 LP",  signal.cheby1, "lowpass",  "cheby1_lp",  3.0, 50.0, None),
        ("Cheby1 HP",  signal.cheby1, "highpass", "cheby1_hp",  3.0, 20.0, None),
        ("Cheby1 BP",  signal.cheby1, "bandpass", "cheby1_bp",  3.0, 20.0, 50.0),
        ("Cheby1 BS",  signal.cheby1, "bandstop", "cheby1_bs",  3.0, 20.0, 50.0),
        ("Cheby2 LP",  signal.cheby2, "lowpass",  "cheby2_lp", 40.0, 50.0, None),
        ("Cheby2 HP",  signal.cheby2, "highpass", "cheby2_hp", 40.0, 20.0, None),
        ("Cheby2 BP",  signal.cheby2, "bandpass", "cheby2_bp", 40.0, 20.0, 50.0),
        ("Cheby2 BS",  signal.cheby2, "bandstop", "cheby2_bs", 40.0, 20.0, 50.0),
    ]:
        sos = scipy_fn(7, ripple, fc1 if fc2 is None else [fc1,fc2],
                       btype=btype, fs=FS, output='sos')
        check_correctness(title, data_c, col, sos, x)
else:
    print("CSV not found — run test_cheby_with_py first")
