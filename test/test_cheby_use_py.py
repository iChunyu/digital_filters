#!/usr/bin/env python3
"""Compare Chebyshev filter C outputs with scipy.

Usage:
  1. Build: cmake -B build && cmake --build build
  2. Generate CSV: ./build/test/test_cheby_with_py
  3. Run this script: python3 test/test_cheby_use_py.py
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
    # (title, scipy_fn, btype, dyn_col, sta_col, ripple, fc1, fc2)
    ("Chebyshev I LP",  signal.cheby1, "lowpass",  "cheby1_lp", "static_cheby1_lp",  3.0, FC_LP,  None),
    ("Chebyshev I HP",  signal.cheby1, "highpass", "cheby1_hp", "static_cheby1_hp",  3.0, FC_HP,  None),
    ("Chebyshev I BP",  signal.cheby1, "bandpass", "cheby1_bp", "static_cheby1_bp",  3.0, FC1_BP, FC2_BP),
    ("Chebyshev I BS",  signal.cheby1, "bandstop", "cheby1_bs", "static_cheby1_bs",  3.0, FC1_BP, FC2_BP),
    ("Chebyshev II LP", signal.cheby2, "lowpass",  "cheby2_lp", "static_cheby2_lp", 40.0, FC_LP,  None),
    ("Chebyshev II HP", signal.cheby2, "highpass", "cheby2_hp", "static_cheby2_hp", 40.0, FC_HP,  None),
    ("Chebyshev II BP", signal.cheby2, "bandpass", "cheby2_bp", "static_cheby2_bp", 40.0, FC1_BP, FC2_BP),
    ("Chebyshev II BS", signal.cheby2, "bandstop", "cheby2_bs", "static_cheby2_bs", 40.0, FC1_BP, FC2_BP),
]

CSV_FILE = os.path.join(os.path.dirname(__file__), "test_cheby_data.csv")


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

    for title, scipy_fn, btype, dyn_col, sta_col, ripple, fc1, fc2 in FILTER_CONFIGS:
        if fc2 is None:
            sos = scipy_fn(ORDER, ripple, fc1, btype=btype, fs=FS, output="sos")
        else:
            sos = scipy_fn(ORDER, ripple, [fc1, fc2], btype=btype, fs=FS, output="sos")

        scipy_out = signal.sosfilt(sos, x)

        c_dyn = data[dyn_col]
        c_sta = data[sta_col]

        family = scipy_fn.__name__
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
        fig.suptitle(f"{title}  (N={ORDER}, ripple={ripple} dB, fs={FS} Hz, input={INPUT_FREQ} Hz)")

        ax1.plot(t, x, label="input", color="gray", alpha=0.5)
        ax1.plot(t, c_dyn, label=f"{family}_t", linestyle="--")
        ax1.plot(t, c_sta, label=f"static_{family}_7th_t", linestyle=":")
        ax1.plot(t, scipy_out, label=f"scipy.signal.{family}", linestyle="-.")
        ax1.set_ylabel("Amplitude")
        ax1.legend(loc="best")
        ax1.grid(True, alpha=0.3)

        ax2.plot(t, c_dyn - scipy_out, label=f"{family}_t − scipy", linestyle="--")
        ax2.plot(t, c_sta - scipy_out, label=f"static_{family}_7th_t − scipy", linestyle=":")
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
