#!/usr/bin/env python3
"""Compare Butterworth filter C outputs with scipy.

Usage:
  1. Build: cmake -B build && cmake --build build
  2. Generate CSV: ./build/test/test_butter_with_py
  3. Run this script: python3 test/test_butter_use_py.py
"""

import os
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

FS = 400.0
FC_LP = 50.0
FC_HP = 20.0
FC1_BP = 20.0
FC2_BP = 50.0
ORDER = 7
INPUT_FREQ = 30.0

FILTER_CONFIGS = [
    ("LP", "lowpass", "butter_lp", "static_butter_lp", FC_LP, None),
    ("HP", "highpass", "butter_hp", "static_butter_hp", FC_HP, None),
    ("BP", "bandpass", "butter_bp", "static_butter_bp", FC1_BP, FC2_BP),
    ("BS", "bandstop", "butter_bs", "static_butter_bs", FC1_BP, FC2_BP),
]

CSV_FILE = os.path.join(os.path.dirname(__file__), "test_butter_data.csv")


def make_pickable_legend(fig):
    """Wire up pick events so clicking a legend entry toggles the line."""
    for ax in fig.axes:
        leg = ax.get_legend()
        if leg is None:
            continue
        lined = {}
        for legline, origline in zip(leg.get_lines(), ax.get_lines()):
            legline.set_picker(True)
            lined[legline] = origline

        def on_pick(event):
            legline = event.artist
            origline = lined.get(legline)
            if origline is None:
                return
            visible = not origline.get_visible()
            origline.set_visible(visible)
            legline.set_alpha(1.0 if visible else 0.2)
            fig.canvas.draw()

        fig.canvas.mpl_connect("pick_event", on_pick)


def main():
    data = np.genfromtxt(CSV_FILE, delimiter=",", names=True)
    t = data["timestamp"]
    x = data["input"]

    for label, btype, dyn_col, sta_col, fc1, fc2 in FILTER_CONFIGS:
        if fc2 is None:
            sos = signal.butter(ORDER, fc1, btype=btype, fs=FS, output="sos")
        else:
            sos = signal.butter(ORDER, [fc1, fc2], btype=btype, fs=FS, output="sos")

        scipy_out = signal.sosfilt(sos, x)

        c_dyn = data[dyn_col]
        c_sta = data[sta_col]

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
        fig.suptitle(f"Butterworth {label}  (N={ORDER}, fs={FS} Hz, input={INPUT_FREQ} Hz)")

        ax1.plot(t, x, label="input", color="gray", alpha=0.5)
        ax1.plot(t, c_dyn, label="butter_t", linestyle="--")
        ax1.plot(t, c_sta, label="static_butter_7th_t", linestyle=":")
        ax1.plot(t, scipy_out, label="scipy.signal.butter", linestyle="-.")
        ax1.set_ylabel("Amplitude")
        ax1.legend(loc="best")
        ax1.grid(True, alpha=0.3)

        ax2.plot(t, c_dyn - scipy_out, label="butter_t − scipy", linestyle="--")
        ax2.plot(t, c_sta - scipy_out, label="static_butter_7th_t − scipy", linestyle=":")
        ax2.set_xlabel("Time (s)")
        ax2.set_ylabel("Error")
        ax2.legend(loc="best")
        ax2.grid(True, alpha=0.3)

        make_pickable_legend(fig)
        plt.tight_layout()

    if plt.get_backend() != "agg":
        plt.show()
    else:
        print("Non-interactive backend detected — skipping plt.show()")


if __name__ == "__main__":
    main()
