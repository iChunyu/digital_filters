# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Test

```bash
cmake -B build -DBUILD_TESTS=ON
cmake --build build
cd build && ctest --output-on-failure
```

`BUILD_TESTS` is ON by default. The library is a static archive (`libbiquad_filter.a`).

## Architecture

This is a **digital filter library** in C, implementing **biquad (second-order IIR)** filters and their **cascade (SOS — second-order sections)** to form higher-order filters. All filters use **Direct Form II (canonical form)**.

Three filter families are available, all sharing the same design pipeline (analog prototype → frequency transform → bilinear discretisation → pole-zero pairing → SOS deployment):
- **Butterworth** — maximally flat passband
- **Chebyshev Type I** — equiripple passband, monotonic stopband
- **Chebyshev Type II** — monotonic passband, equiripple stopband

**Biquad transfer function:** `H(z) = (b0 + b1·z⁻¹ + b2·z⁻²) / (1 + a1·z⁻¹ + a2·z⁻²)`

### Key design decisions

**biquad_filter**
- **Coefficient normalisation**: `den_z[0]` is always normalised to 1.0 internally. `biquad_filter_init()` divides all coefficients by `den_z[0]` on entry.
- **Silent fallback to identity**: If the denominator is zero or the poles are unstable, `biquad_filter_init()` silently replaces the filter with a unity pass-through (`H(z) = 1`). No error codes.
- **Stability check** uses all three **Jury conditions** for second-order systems: `|a2| < 1`, `1 + a1 + a2 > 0`, `1 - a1 + a2 > 0`.
- **State vector `w[3]`** holds `w[n], w[n-1], w[n-2]` — the intermediate nodes in the canonical Direct Form II signal flow graph. Only 3 floats of state per biquad.
- **Bilinear transform** (`biquad_c2d_bilnear`) maps s-domain to z-domain using the substitution `s = 2·fs · (1 - z⁻¹)/(1 + z⁻¹)`.
- **Steady-state reset** computes `w_ss = equilibrium / (1 + a1 + a2)`. The denominator is guaranteed non-zero because `biquad_filter_init()` enforces `1 + a1 + a2 > 0` (otherwise it falls back to identity where `a1 = a2 = 0`).
- `biquad_filter_get_input()` reconstructs `x[n]` from state — useful for debugging or cascading filters.

**sos_filter (cascade of biquads)**
- `sos_filter_init` allocates the sections array with `malloc` and initialises all sections to identity. Must pair with `sos_filter_destroy` to free. Does **not** store sampling frequency — that belongs at the design level (butter_filter_t, cheby1_t, etc.).
- `sos_filter_update` processes samples in series: `input → sec[0] → sec[1] → ... → output`.
- `sos_filter_reset` propagates equilibrium through the cascade: each section is reset with the DC output of the previous section (computed via `biquad_filter_get_output`), preserving the overall DC response.
- `sos_filter_get_input` recovers the input to the first section (the overall filter input). `sos_filter_get_output` returns the output of the last section.
- `sos_filter_set_section` silently ignores out-of-bounds indices, consistent with the identity-fallback pattern.

**filter_utils (design-time helpers)**
- `complex_t` — simple `{float re, im}` struct used throughout the design pipeline.
- `prewarp(fd, fs)` — bilinear frequency pre-warping: `f_analog = fs/π · tan(π · fd / fs)`.
- `analog_lp_transform / analog_hp_transform` — s-domain frequency transforms (scale by `wc` for LP, `wc/p` for HP).
- `bilinear_transform` — maps s-domain poles/zeros to z-domain in-place: `z = (2fs + s) / (2fs - s)`.
- `zpk2sos` — converts z-domain pole/zero arrays to SOS coefficient matrix. Uses the "most-unfavorable-pole-first" pairing algorithm (selects the pole closest to the unit circle, pairs with nearest zero, handles real/complex conjugate pairs, reverses section order so fast poles come first). Each section is gain-normalised for unity at the reference frequency (DC for LP, Nyquist for HP).

**butter_filter / cheby_filter (high-level filter objects)**
- Each filter type has its own struct (`butter_t`, `cheby1_t`, `cheby2_t`) wrapping a `sos_filter_t*`.
- `_init` orchestrates the full design pipeline, taking type (LP/HP), order, cutoff(s), fs, and (for Chebyshev) ripple dB. On failure sets `valid = 0` and leaves `sos` NULL.
- `_update` / `_reset` delegate to the underlying SOS cascade; both pass through unchanged if `!valid`.
- `_destroy` frees the SOS cascade and marks the filter invalid.
- Chebyshev structs include a `ripple_db` field (passband ripple for Type I, minimum stopband attenuation for Type II).
- BP and BS transforms are not yet implemented.

### Files

| File | Role |
|---|---|
| `include/biquad_filter.h` | Biquad public API with Doxygen docs |
| `include/sos_filter.h` | SOS cascade public API |
| `include/filter_utils.h` | Design-time utilities (complex_t, prewarp, transforms, zpk2sos) |
| `include/butter_filter.h` | Butterworth filter type and API |
| `include/cheby_filter.h` | Chebyshev I & II filter types and API |
| `src/biquad_filter.c` | Biquad implementation |
| `src/sos_filter.c` | SOS cascade implementation |
| `src/filter_utils.c` | Design utility implementations |
| `src/butter_filter.c` | Butterworth prototype + init/destroy/update/reset |
| `src/cheby_filter.c` | Chebyshev I & II prototypes + init/destroy/update/reset |
| `test/test_biquad.c` | Biquad tests (CHECK/CLOSE harness) |
| `test/test_sos.c` | SOS cascade tests |
| `test/test_butter.c` | Butterworth LP/HP tests |
| `test/test_cheby.c` | Chebyshev I/II LP/HP tests |
| `CMakeLists.txt` | Top-level CMake (static lib, tests, install rules) |
| `test/CMakeLists.txt` | Test executables + CTest registration |
