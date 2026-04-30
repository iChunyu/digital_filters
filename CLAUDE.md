# CLAUDE.md

本文件为 Claude Code (claude.ai/code) 在此仓库中工作时提供指引。

## 构建与测试

```bash
cmake -B build -DBUILD_TESTS=ON
cmake --build build
cd build && ctest --output-on-failure
```

`BUILD_TESTS` 默认 ON。库为静态归档 (`libdigital_filters.a`)。

## 架构

这是一个 **数字滤波器库**，用 C 实现，基于 **双二阶 (biquad, 二阶 IIR)** 滤波器及其 **级联 (SOS — 二阶节)** 构成高阶滤波器。所有滤波器均采用 **Direct Form II（规范型）**。

三种滤波器族共享同一条设计流水线（模拟原型 → 频率变换 → 双线性离散化 → 零极点配对 → SOS 部署）：
- **Butterworth** — 通带最大平坦
- **Chebyshev Type I** — 通带等波纹，阻带单调
- **Chebyshev Type II** — 通带单调，阻带等波纹

支持四种滤波器类型：`FILTER_LOWPASS`, `FILTER_HIGHPASS`, `FILTER_BANDPASS`, `FILTER_BANDSTOP`。

**双二阶传递函数:** `H(z) = (b0 + b1·z⁻¹ + b2·z⁻²) / (1 + a1·z⁻¹ + a2·z⁻²)`

### 两套 API

项目提供两套 API 共存：

| API | 内存模型 | 适用场景 |
|---|---|---|
| 动态 API (`butter_t`, `cheby1_t`, `cheby2_t`) | `malloc` / `free` | 桌面端，阶数运行时决定 |
| 静态 API (`butter_lp_2nd_t`, `cheby1_hp_3rd_t` 等) | 结构体内嵌数组，零 `malloc` | MCU，阶数编译时确定 |

静态 API 的关键设计：
- 所有 Butterworth 静态结构体共享统一前缀 (`STATIC_BUTTER_FIELDS`)，因此 **update/reset 只需一套函数**：`static_butter_update` / `static_butter_reset`，通过 `(static_butter_t *)` 强制转换即可通用于全部阶数和类型。Chebyshev 同理 (`static_cheby_t`、`static_cheby_update` / `static_cheby_reset`)。
- 只有 `_init` 按阶数分别定义，因为不同阶数需要不同大小的栈上临时数组。
- Butterworth 原型极点仅依赖阶数 N，预计算为 `static const complex_t butter_proto[8][8]` 存入 ROM（~512 字节），init 时无需调用 `cosf`/`sinf`。
- Chebyshev 原型依赖 ripple，在 init 时运行时计算。
- 节数规则：LP/HP 需 `ceil(N/2)` 节，BP/BS 需 `N` 节（频率变换使阶数翻倍）。8 阶 BP/BS 最多 8 节。

### 关键设计决策

**biquad_filter**
- **系数归一化**: `den_z[0]` 内部始终归一化为 1.0。`biquad_filter_init()` 会在入口处将所有系数除以 `den_z[0]`。
- **静默降级为直通**: 若分母为零或极点不稳定，`biquad_filter_init()` 静默地将滤波器替换为单位直通 (`H(z) = 1`)。无错误码。
- **稳定性检测** 使用二阶系统的全部三个 **Jury 条件**: `|a2| < 1`, `1 + a1 + a2 > 0`, `1 - a1 + a2 > 0`。
- **状态向量 `w[3]`** 保存 `w[n], w[n-1], w[n-2]` — 规范型 Direct Form II 信号流图中的中间节点。每个 biquad 仅需 3 个 float 状态。
- **双线性变换** (`biquad_c2d_bilnear`) 使用代换 `s = 2·fs · (1 - z⁻¹)/(1 + z⁻¹)` 将 s 域映射到 z 域。
- **稳态复位** 计算 `w_ss = equilibrium / (1 + a1 + a2)`。分母保证非零，因为 `biquad_filter_init()` 强制 `1 + a1 + a2 > 0`（否则降级为直通，此时 `a1 = a2 = 0`）。
- `biquad_filter_get_input()` 从状态重建 `x[n]` — 用于调试或滤波器级联。

**sos_filter（biquad 级联，动态 API）**
- `sos_filter_init` 通过 `malloc` 分配节数组并初始化所有节为直通。必须与 `sos_filter_destroy` 配对以释放。**不** 存储采样频率 — 这属于设计层面（butter_t, cheby1_t 等）。
- `sos_filter_update` 串行处理采样：`input → sec[0] → sec[1] → ... → output`。
- `sos_filter_reset` 将稳态值传播通过级联：每节以前一节的 DC 输出（通过 `biquad_filter_get_output` 计算）复位，保持整体 DC 响应。
- `sos_filter_get_input` 恢复第一节的输入（整体滤波器输入）。`sos_filter_get_output` 返回最后一节的输出。
- `sos_filter_set_section` 静默忽略越界索引，与直通降级模式一致。

**filter_utils（设计时工具）**
- `complex_t` — 简单的 `{float re, im}` 结构体，贯穿整个设计流水线。
- `prewarp(fd, fs)` — 双线性频率预畸变：`f_analog = fs/π · tan(π · fd / fs)`。
- `analog_lp_transform / analog_hp_transform` — s 域频率变换（LP 缩放 `wc`，HP 为 `wc/p`）。
- `analog_bp_transform / analog_bs_transform` — s 域带通/带阻变换（代换 `s → (s²+ω₀²)/(ξ·s)` 或反向，阶数翻倍）。
- `bilinear_transform` — 原地将 s 域极点/零点映射到 z 域：`z = (2fs + s) / (2fs - s)`。
- `zpk2sos` — 将 z 域零极点数组转为 SOS 系数矩阵。使用"最不利极点优先"配对算法（选择最靠近单位圆的极点，与最近的零点配对，处理实/复共轭对，反转节序使快极点在前）。每节在参考频率处增益归一化（LP/BS 为 DC，HP 为 Nyquist，BP 为中心频率）。

**动态 API 滤波器对象 (butter_filter / cheby_filter)**
- 每种滤波器类型有各自的结构体 (`butter_t`, `cheby1_t`, `cheby2_t`)，包装一个 `sos_filter_t*`。
- `_init` 编排完整设计流水线，接收 type (LP/HP/BP/BS)、order、截止频率、fs，以及（对于 Chebyshev）ripple dB。失败时设 `valid = 0` 并将 `sos` 置 NULL。
- `_update` / `_reset` 委托给底层 SOS 级联；均会在 `!valid` 时直通返回。
- `_destroy` 释放 SOS 级联并将滤波器标记为无效。
- Chebyshev 结构体包含 `ripple_db` 字段（Type I 为通带纹波，Type II 为最小阻带衰减）。
- Chenyshev II 的 ε 参数公式：`ε = 1 / sqrt(10^(dB/10) - 1)`（与 Type I 的 `ε = sqrt(10^(dB/10) - 1)` 互为倒数）。

**静态 API 滤波器对象 (static_butter_filter / static_cheby_filter)**
- 每个滤波器族有统一基类型 (`static_butter_t`, `static_cheby_t`)，所有按阶数/类型定义的结构体共享相同前缀布局，可安全强制转换。
- 结构体通过 X-macro 生成：`FOR_EACH_STATIC_BUTTER_LP_ORDER` (order, sections, ordinal) 和 `FOR_EACH_STATIC_BUTTER_BP_ORDER`。每个滤波器族 32 个结构体（4 种类型 × 8 阶），Chebyshev 为 64 个（I 和 II 各 32 个）。
- `_init` 函数接收截止频率和（对于 Chebyshev）ripple，在栈上完成设计流水线，将系数部署到内嵌的 `biquad_filter_t sections[]` 中。
- 无需 `_destroy` — 结构体离开作用域即自动回收。

### 文件

| 文件 | 用途 |
|---|---|
| `include/biquad_filter.h` | Biquad 公开 API 含 Doxygen 文档 |
| `include/sos_filter.h` | SOS 级联公开 API（动态，malloc） |
| `include/filter_utils.h` | 设计时工具 (complex_t, prewarp, 变换, zpk2sos) |
| `include/butter_filter.h` | Butterworth 动态 API |
| `include/cheby_filter.h` | Chebyshev I & II 动态 API |
| `include/static_butter_filter.h` | Butterworth 静态 API（32 个按阶结构体，统一 update/reset） |
| `include/static_cheby_filter.h` | Chebyshev I & II 静态 API（64 个按阶结构体，统一 update/reset） |
| `src/biquad_filter.c` | Biquad 实现 |
| `src/sos_filter.c` | SOS 级联实现（含 malloc/free） |
| `src/filter_utils.c` | 设计工具实现 |
| `src/butter_filter.c` | Butterworth 原型 + 动态 init/destroy/update/reset |
| `src/cheby_filter.c` | Chebyshev I & II 原型 + 动态 init/destroy/update/reset |
| `src/static_butter_filter.c` | 预计算极点表、共享流水线、按阶 init、统一 update/reset |
| `src/static_cheby_filter.c` | 运行时原型计算、共享流水线、按阶 init、统一 update/reset |
| `test/test_biquad.c` | Biquad 测试 (CHECK/CLOSE 框架) |
| `test/test_sos.c` | SOS 级联测试 |
| `test/test_butter.c` | Butterworth 动态 API 测试 (LP/HP/BP/BS) |
| `test/test_cheby.c` | Chebyshev I/II 动态 API 测试 |
| `test/test_static_butter.c` | Butterworth 静态 API 测试（全类型、多阶数） |
| `test/test_static_cheby.c` | Chebyshev 静态 API 测试 |
| `CMakeLists.txt` | 顶层 CMake（静态库、测试、安装规则） |
| `test/CMakeLists.txt` | 测试可执行文件 + CTest 注册 |
