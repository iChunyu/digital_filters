# CLAUDE.md

本文件为 Claude Code (claude.ai/code) 在此仓库中工作时提供指引。

## 构建与测试

```bash
cmake -B build -DBUILD_TESTS=ON
cmake --build build
cd build && ctest --output-on-failure
```

`BUILD_TESTS` 默认 ON。库为静态归档 (`libdigital_filters.a`)。

### MCU / 嵌入式 使用注意事项

**FPU 要求**: 所有滤波器 `_update` 路径执行密集 float 运算。建议使用带硬件 FPU 的 MCU（Cortex-M4/M7 及以上）。

**栈需求**: init 期间峰值约 800 字节（`butter_design`/`cheby_design` → `zpk2sos` 调用链）。运行时 `_update` 为全 inline，无额外栈开销。建议 MCU 主栈 ≥ 2 KB。

**中断安全**: `_update` 和 `_reset` 不可重入。同一个滤波器结构体如果被 ISR 和主循环共享，需在调用 `_update` 前关中断或使用双缓冲。


## 架构

这是一个面向 **MCU / 嵌入式** 的 **IIR 数字滤波器库**，纯 C，基于 **双二阶 (biquad)** 滤波器及其级联（SOS，二阶节）构成高阶滤波器。所有滤波器采用 **Direct Form II（规范型）**。**零 `malloc`，零 `double`，全部 `float`。**

三种滤波器族共享同一条设计流水线（模拟原型 → 频率变换 → 双线性离散化 → 零极点配对 → SOS 部署）：
- **Butterworth** — 通带最大平坦
- **Chebyshev Type I** — 通带等波纹，阻带单调
- **Chebyshev Type II** — 通带单调，阻带等波纹

支持四种滤波器类型：`FILTER_LOWPASS`, `FILTER_HIGHPASS`, `FILTER_BANDPASS`, `FILTER_BANDSTOP`。
最大原型阶数 **8 阶**（BP/BS 有效 16 阶）。

**双二阶传递函数:** `H(z) = (b0 + b1·z⁻¹ + b2·z⁻²) / (1 + a1·z⁻¹ + a2·z⁻²)`

### API 设计

所有结构体由 X-macro 生成，阶数编译时确定，内嵌 `biquad_filter_t sections[]`，
零堆分配。离开作用域自动回收，无需 `_destroy`。

```
butter_lp_2nd_t   cheby1_hp_3rd_t   cheby2_bs_5th_t  ...
```

- **Butterworth**：32 个结构体（4 类型 × 8 阶），共享 `BUTTER_FIELDS` 前缀
- **Chebyshev**：64 个结构体（Type I/II × 4 类型 × 8 阶），共享 `CHEBY_FIELDS` 前缀

内部辅助函数直接接收 `biquad_filter_t *sections` + 参数，宏生成的 init/update/reset
在调用侧传递对应字段，无需强转。内部函数命名遵循 `{butter,cheby{1,2}}_{type}_{func}` 规则
（如 `butter_lp_init`、`butter_design`、`cheby2_bp_init`）。

对外 exposed 的 update/reset 按阶数/类型分别定义（如 `butter_lp_2nd_update`、
`cheby1_bp_5th_reset`），调用时无需强转。

- 只有 `_init` 按阶数分别定义，因为不同阶数需要不同大小的栈上临时数组
- Butterworth 原型极点预计算为 `static const complex_t butter_proto[8][8]` 存入 ROM（~512 字节），
  init 时无需调用 `cosf`/`sinf`
- Chebyshev 原型依赖 ripple，在 init 时运行时计算
- 节数规则：LP/HP 需 `ceil(N/2)` 节，BP/BS 需 `N` 节（频率变换使阶数翻倍）

### 关键设计决策

**biquad_filter**
- **系数归一化**: `den_z[0]` 内部始终归一化为 1.0
- **静默降级为直通**: 分母为零或极点不稳定 → 替换为单位直通 (`H(z)=1`)
- **稳定性检测**：全部三个 Jury 条件 — `|a2| < 1`, `1 + a1 + a2 > 0`, `1 - a1 + a2 > 0`
- **状态向量 `w[3]`**：Direct Form II，每 biquad 仅需 3 个 `float` 状态
- **`biquad_filter_update` 为 `static inline`**（header-only）。状态更新与输出计算融合——
  先快照 `w[0]`、`w[1]` 到寄存器，一次性缓存全部 5 个系数，计算完 `w0` 后批量写回状态 + 计算输出
- **稳态复位**：`w_ss = equilibrium / (1 + a1 + a2)`，分母保证非零

**filter_utils（设计时工具）**
- `complex_t` — `{float re, im}`
- `prewarp(fd, fs)` — `f_analog = fs/π · tan(π · fd / fs)`
- `analog_lp/hp/bp/bs_transform` — s 域频率变换
- `bilinear_transform` — 原地映射 `z = (2fs + s) / (2fs - s)`
- 增益追踪全部使用 `float`，交替乘除避免中间溢出
- `zpk2sos` — "最不利极点优先"配对算法。工作数组 `ZPK2SOS_MAX_N = 48`

### 文件

| 文件 | 用途 |
|---|---|
| `include/biquad_filter.h` | Biquad 公开 API（含 inline update） |
| `include/filter_utils.h` | 设计工具 (complex_t, prewarp, 变换, zpk2sos) |
| `include/butter_filter.h` | Butterworth API（X-macro 生成 32 结构体） |
| `include/cheby_filter.h` | Chebyshev I & II API（X-macro 生成 64 结构体） |
| `src/biquad_filter.c` | Biquad 实现（init, reset, get_output, get_input, c2d_bilinear） |
| `src/filter_utils.c` | 设计工具实现 |
| `src/butter_filter.c` | 预计算极点表、设计流水线、按阶 init/update/reset |
| `src/cheby_filter.c` | 运行时原型计算、设计流水线、按阶 init/update/reset |
| `test/test_biquad.c` | Biquad 测试 |
| `test/test_butter.c` | Butterworth 测试（全类型、多阶数） |
| `test/test_cheby.c` | Chebyshev I/II 测试（全类型、多阶数） |
| `test/test_butter_with_py.c` | 生成 CSV 与 scipy 参考对比 |
| `test/test_cheby_with_py.c` | 生成 CSV 与 scipy 参考对比 |
| `test/test_butter_use_py.py` | Python 参考滤波器（Butterworth，含绘图） |
| `test/test_cheby_use_py.py` | Python 参考滤波器（Chebyshev，含绘图） |
| `test/compare_scipy.py` | 稳态精度对比脚本（跳过瞬态取最后 10%） |
| `test/verify_zpk_gain.py` | 验证 float32 精度足够 zpk 增益追踪（独立复现 scipy 管线，f64 vs f32） |
| `CMakeLists.txt` | 顶层 CMake |
| `test/CMakeLists.txt` | 测试可执行文件 + CTest 注册 |
