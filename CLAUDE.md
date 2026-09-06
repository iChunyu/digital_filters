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

**栈需求**: init 期间峰值约 800 字节（`design_filter` → `zpk2sos` 调用链：poles 128 B
+ zeros 128 B + sos 192 B ≈ 448 B，加上调用侧原型数组 ~128 B；该数值在 -O2 下实测成立，
-O0 下约 1.1 KB）。运行时 `_update` 全 inline（biquad → 级联 → 逐阶函数三层全部
头文件内联），无额外栈开销。建议 MCU 主栈 ≥ 2 KB。

**libm 依赖**: `biquad_filter_init`/`_update`/`_reset` 路径零 libm 调用（裕量为纯乘加
代数式），裸机不链 libm 也可用 biquad 层。Butterworth init 无需 libm（ROM 查表）。
Chebyshev init 需要 `logf/sqrtf/sinhf/coshf/cosf/sinf`（仅 init 时）。

**Flash 粒度**: 库以 `-ffunction-sections/-fdata-sections` 编译；链接时加
`--gc-sections` 只拉入用到的滤波器族（实测 butter_lp_2nd 单独使用 29.4 → 20.9 KB text）。
update/reset 为头文件 inline，不占库 text。

**中断安全**: `_update` 和 `_reset` 不可重入。同一个滤波器结构体如果被 ISR 和主循环共享，需在调用 `_update` 前关中断或使用双缓冲。

**非有限输入**: 热路径无 NaN/Inf 防护（每样本分支开销），非有限输入毒化状态直至 reset；调用方在源头清洗。


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

内部辅助函数直接接收 `biquad_filter_t *sections` + 参数，宏生成的 init 在调用侧
传递对应字段，无需强转。内部函数命名遵循 `{butter,cheby{1,2}}_{type}_{func}` 规则
（如 `butter_lp_init`、`cheby2_bp_init`）。

设计管线只有一份：`design_filter()`（`filter_utils.c`）——模拟原型 → 频率变换 →
增益折叠 → 双线性 → `zpk2sos` → 部署。Butterworth 传 ROM 极点表 + `k=1, nz=0`；
Chebyshev 运行时算原型传自己的 `k` 和有限零点。级联增益校验 `check_cascade_gains()`
同样只有一份。fail-closed 咽喉点（k 有限非零、节数上限、逐节 init、增益窗口）
全部单点维护。

对外 exposed 的 update/reset 按阶数/类型分别定义（如 `butter_lp_2nd_update`、
`cheby1_bp_5th_reset`），调用时无需强转。它们是**头文件 X-macro 生成的
`static inline`**：调 `biquad_cascade_update(f->sections, ns, x)`，`ns` 是字面量
（循环边界编译期折叠）。每样本路径无函数调用、无运行时节数装载、无链接符号
（96 个字节级相同的 out-of-line 副本曾占 ~12.8 KB flash + 每样本一次 BL）。

- 只有 `_init` 按阶数分别定义，因为不同阶数需要不同大小的栈上临时数组
- Butterworth 原型极点预计算为 `static const complex_t butter_proto[8][8]` 存入 ROM（~512 字节），
  init 时无需调用 `cosf`/`sinf`
- Chebyshev 原型依赖 ripple，在 init 时运行时计算
- 节数规则：LP/HP 需 `ceil(N/2)` 节，BP/BS 需 `N` 节（频率变换使阶数翻倍）

### 关键设计决策

**biquad_filter**
- **系数归一化**: `den_z[0]` 内部始终归一化为 1.0
- **静默降级为直通**: 分母为零或非有限、任一系数非有限、极点不稳定或裕量不足
  → 替换为单位直通 (`H(z)=1`)，`biquad_filter_init` 返回 0。宁可直通也不要发散/静音
- **稳定性检测**：全部三个 Jury 条件 — `|a2| < 1`, `1 + a1 + a2 > 0`, `1 - a1 + a2 > 0`。
  Jury 和用**补偿求和**（TwoSum 式 `sum3f`）——裸 f32 求和在窄带设计上会恰好
  归零误拒合法滤波器；补偿和**恰好为 0** 则说明 f32 系数把极点放到了单位圆上
  （量化塌缩到 z=±1）→ 拒绝。
  外加**极点半径裕量检查**：半径 > 0.99995（距单位圆 < 5e-5，振铃 ≥ 10⁴ 采样）拒绝。
  按极点半径本身判定，不用 `a2`（`a2 = r1·r2` 乘积对非等实根对的主导极点有盲区，
  单侧 `a2 > 0.9999` 漏掉异号实根对）：共轭对 `a2 = r²`，实根
  `r_max = (|a1| + √disc)/2` 化为纯乘加代数式（`√disc > 1.9999 − |a1|` 平方化），
  **biquad init 路径零 libm 依赖**（裸机不链 libm 也能用）；`a1²−4a2` 用
  Dekker 分裂补偿计算（超宽带近实对的裸 f32 判别式噪声 ~5e-7 与真值同量级，
  会让共轭/实根分支翻转）。实测接受的合法设计最小裕量 ~5.2e-5
  （cheby1_bs_8th [100,200]@48k）、~6.5e-5（lp_1st fc=0.5 Hz@48k），
  阈值 5e-5 恰好卡在其下方
- **状态向量 `w[3]`**：Direct Form II，每 biquad 仅需 3 个 `float` 状态
- **`biquad_filter_update` 为 `static inline`**（header-only）。状态更新与输出计算融合——
  先快照 `w[0]`、`w[1]` 到寄存器，一次性缓存全部 5 个系数，计算完 `w0` 后批量写回状态 + 计算输出。
  热路径**刻意无 NaN/Inf 输入防护**（每样本分支开销）——非有限输入毒化状态，
  reset 恢复；调用方在源头清洗（文档已写明）
- **稳态复位**：`w_ss = equilibrium / (1 + a1 + a2)`（补偿求和）；分母恰好为 0
  （纯积分器）或非有限、或 `w_ss` 溢出 → 状态清零（对齐 @note，永不 inf/NaN）。
  reset 是公开结构体上的公开 API，系数无需先过 init

**filter_utils（设计时工具）**
- `complex_t` — `{float re, im}`
- `prewarp(fd, fs)` — `f_analog = fs/π · tan(π · fd / fs)`
- `analog_lp/hp/bp/bs_transform` — s 域频率变换。BP/BS 用**稳定二次公式**
  （小根 = `C / root1`）避免宽频带实极点的大数消减噪声（噪声曾破坏共轭配对）
- `c_sqrt` — **缩放幅值**计算（`m·sqrt((re/m)²+(im/m)²)`）：裸平方在近 Nyquist
  BP/BS 判别式（|re| ~ 1e20 > √FLT_MAX）上溢出 → inf/NaN 级联
- `bilinear_transform` — 原地映射 `z = (2fs + s) / (2fs - s)`
- 增益追踪全部使用 `float`，交替乘除避免中间溢出；LP/BP 的 `s^degree` 因子
  经 `bilinear_zpk_gain_scaled` 与 `(K−p)` 除法**交错折叠**——近 Nyquist 设计
  （wc^8 ≈ FLT_MAX）不再先溢出后抵消
- `zpk2sos` — "最不利极点优先"配对算法。**两个容差分开**：实/复分类 `eps_class
  = 1e-5`（真·实极点 im 恒为 0，真·共轭对 |im| ≥ ~1e-5；1e-3 会把宽频带近实
  共轭对误判为实、跨对错配出极点恰在 z=1 的节）；共轭认领盒 `eps_claim = 1e-3`
  （须容纳双线性 f32 舍入 ~2e-4）。`claim_conjugate` 在盒内取**最近**未用极点
  （原"第一个盒内"在高 Q 簇、对间距 8e-4 < 盒 1e-3 时偷别对的伴侣，部署出完全
  重复的两节 + 丢一对，终态不变量检测不到——实测 cheby1_bs_8th [100,200]@48k
  通带内 +8.7 dB 尖峰）。循环结束强制**不变量**：全部零极点被认领 + 每节根与
  输入零极点**多重集匹配**，任一失败返回 0 → 调用方直通（曾拦截：误配导致
  350× 谐振的静默错误滤波器）
- 工作数组 `ZPK2SOS_MAX_N = 16`
- `design_filter` / `check_cascade_gains` — 两族共享的设计管线与增益校验，
  fail-closed 咽喉点单点维护（见上文 API 设计）

**设计管线级联校验（design_filter 部署后）**
- `k` 必须有限且非零（k=0 会部署出输出恒 0 的"静音"滤波器）
- 解析级联 DC/Nyquist 增益必须匹配族/类型/阶数奇偶的期望值：
  精确结构增益（0/1，由 z=±1 零点保证）窗口 ±0.1；纹波边缘增益
  （cheby1/2 偶数阶 `10^(−rp/rs/20)`）窗口 ±0.25。实测接受的合法设计
  实现误差 ≤ ~3.4e-2（窄带 fc≈6 Hz@48k 的 f32 系数量化真实影响），
  配错对缺陷偏差 ≥ 0.45——窗口两侧都有 ≥ 3× 分离。失败 → 整链直通

**频率包络（fail-closed 拒绝 = 直通，属预期）**
- 极窄带/近 DC 设计极点进入单位圆 5e-5 内被拒（例：lp_2nd fc=1 Hz@48k 被拒、
  fc=6 Hz 通过）；极近 Nyquist 同理（bs_1st [100,23999.8]@48k 被拒）
- 超宽带 BP/BS（fc2/fc1 ≳ 1000）节内状态巨大（w ≈ 1/(1+a1+a2)），f32 状态
  更新噪声把阻带衰减地板抬到 ~−12 dB 量级（DF-II f32 固有，非缺陷）
- cheby2 奇数阶零点跳过阈值 1e-4（合法零点 |sinθ| ≥ 0.38，余量 3 个数量级；
  对软浮点 libm 的 sinf(π) 误差也留足 3 个数量级）

### 文件

| 文件 | 用途 |
|---|---|
| `include/biquad_filter.h` | Biquad 公开 API（inline update + 级联核心 cascade_update/reset） |
| `include/filter_utils.h` | 设计工具 (complex_t, prewarp, 变换, zpk2sos, design_filter, check_cascade_gains) |
| `include/butter_filter.h` | Butterworth API（X-macro 生成 32 结构体 + inline update/reset） |
| `include/cheby_filter.h` | Chebyshev I & II API（X-macro 生成 64 结构体 + inline update/reset） |
| `src/biquad_filter.c` | Biquad 实现（init, reset, get_output, get_input, c2d_bilinear） |
| `src/filter_utils.c` | 设计工具实现 + 共享设计管线 design_filter + 增益校验 |
| `src/butter_filter.c` | 预计算极点表、按阶 init（管线走 design_filter） |
| `src/cheby_filter.c` | 运行时原型计算、按阶 init（管线走 design_filter） |
| `test/test_biquad.c` | Biquad 测试（含裕量对称性、积分器 reset、NaN 语义） |
| `test/test_butter.c` | Butterworth 测试（全类型、多阶数 + 近实配对/窄带回归） |
| `test/test_cheby.c` | Chebyshev I/II 测试（全类型、多阶数 + 高 Q 配对回归） |
| `test/test_butter_with_py.c` | 生成 CSV（argv 指定路径 + valid 检查）与 scipy 对比 |
| `test/test_cheby_with_py.c` | 生成 CSV（argv 指定路径 + valid 检查）与 scipy 对比 |
| `test/test_butter_use_py.py` | Python 参考滤波器（Butterworth，含绘图，argv 指定数据目录） |
| `test/test_cheby_use_py.py` | Python 参考滤波器（Chebyshev，含绘图，argv 指定数据目录） |
| `test/compare_scipy.py` | 稳态精度对比（已注册 ctest，无 numpy/scipy 时 exit 77 SKIP） |
| `test/verify_zpk_gain.py` | 验证 float32 精度足够 zpk 增益追踪（独立复现 scipy 管线，f64 vs f32） |
| `CMakeLists.txt` | 顶层 CMake（库带 -ffunction-sections/-fdata-sections） |
| `test/CMakeLists.txt` | 测试可执行文件 + CTest 注册（6 项，含 CSV 生成 + scipy 对比） |
