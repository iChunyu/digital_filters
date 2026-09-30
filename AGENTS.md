# AGENTS.md

给在本仓库改代码的 agent 的工作手册：本文件是**跨文件约定与"不要踩的坑"的唯一总账**，
每次会话自动加载。用户向 API 用法见 `README.md`；什么内容该写在哪一层见「知识分层约定」。

## 硬约束

- 纯 C，面向 MCU / 嵌入式；**只有 `float`**，零 `double`，零 `malloc`。
- 所有结构体阶数编译期确定，内嵌 `biquad_filter_t sections[]`，离开作用域即回收。
- `_update` / `_reset` 是头文件 `static inline`；热路径无函数调用，也**刻意无 NaN/Inf
  防护**（每样本分支开销），非有限输入毒化状态直至 reset，调用方在源头清洗。
- `_update` / `_reset` 不可重入：同一结构体被 ISR 与主循环共享时需关中断或双缓冲。

## 构建与测试

```bash
cmake -B build -DBUILD_TESTS=ON
cmake --build build
cd build && ctest --output-on-failure
```

预期 **13/13 通过**：5 项 C 单元/回归（`test_biquad`、`test_butter`、`test_cheby`、
`test_notch`、`test_peak`）+ 4 项 CSV 生成 + scipy/陷波/峰值三项黄金参考对比 +
`verify_zpk_gain`。后 4 项需要 Python3 + numpy/scipy（`verify_zpk_gain` 只需 numpy），
缺库时以退出码 77 自动 SKIP；配置期找不到 Python3 时 CSV 生成器不注册，并在配置输出
打印 WARNING——那时 ctest 只有 5 项，**不要把它当成 13/13**。

另有 ASan+UBSan 变体：`cmake -B build-asan -DBUILD_SANITIZE=ON`；`design_filter` 的
BP/BS 越界写、`proto_zeros` 未初始化读这类缺陷只有它抓得到（返回值恰好是 0 会把 UB
掩盖成"测试通过"）。改动系数、配对或阈值后，黄金参考对比会量化响应偏差。

`test_biquad` 不再有不链 `m` 的特殊地位——库已把 `m` 作为 `PUBLIC` 依赖声明，各测试
无需自己写。新增源文件或测试要同步 `CMakeLists.txt` / `tests/CMakeLists.txt`。

## 知识分层约定

知识按来源分流；改动时把内容放进对应的层：

| 层 | 放什么 | 判据 |
|---|---|---|
| `include/`、`src/` | 只留**"删了会静默坏"的陷阱一行** | 见下面三条准入判据 |
| `docs/<源文件名>.md` | 推导、设计理由、被否掉的备选方案、公式（数学用 LaTeX） | 需要回答"为什么" |
| `tests/*.c` | 实测标定、回归档位、实测统计表、缺陷复盘 | 来自运行结果，不是推理 |
| `AGENTS.md` | 构建/测试、硬约束、跨文件约定、"别去修"清单、知识地图 | AI 每次会话都要知道 |
| `README.md` | 用户向用法、API 一览、精度对比、MCU 部署、文件树 | 用户需要知道 |

**新增源码注释的准入判据**（三条同时成立才留在源码，否则进 `docs/`）：① **非显然**——
读者看代码本身推不出来；② 写错了是**静默失败**——编译器与测试都不会红；③ **≤2 行**。
跨文件约定进本文件，公式进 `docs/`，实测数字进 `tests/`。

注释一律中文（`/* */`、构建脚本 `#`、Python `#` 与 docstring）；代码标识符、数学记号、
单位、专有名词（`scipy`、`Butterworth`、`Dekker`、函数名、类型名）保持原样，不要硬译。
`include/` + `src/` **禁 doxygen 标记**（`@brief` / `@param` / `@return` / `@note` /
`@ref` / `@f$` / `@file` / `@c` / `@p`）：符号语义靠命名与 `docs/` 对应篇目承载。
（`tests/*.c` 的注释是标定存档、**一字不动**，不受此约束。）

## 架构地图

五族，两条设计流水线：

- **Butterworth / Chebyshev I / Chebyshev II** — 高阶，走共享管线 `design_filter()`：
  模拟原型 → 频率变换 → 增益折叠 → 双线性 → `zpk2sos()` → 部署。Butterworth 传 ROM
  极点表（`butter_proto[8][8]`，~512 B）、`k=1`、`nz=0`；Chebyshev 运行时算原型
  （依赖 ripple），传自己的 `k` 与有限零点。结构体与 init/update/reset 由 X-macro
  生成（`BUTTER_FIELDS` / `CHEBY_FIELDS`、`FOR_EACH_*_ORDER`），共 32 + 64 个结构体。
- **Notch / Peak** — 固定二阶、单节、无类型变体，因此**不用 X-macro**：预畸 →
  `biquad_c2d_bilinear` → `biquad_filter_init` 直写一个 `biquad_filter_t`，不走零极点
  配对。两族严格互为倒数（`H_peak(ξ,g) ≡ 1/H_notch(ξ,1/g)`），`peak_filter.c` 把分母
  写成 `2·xi·d·w0` 正是为了让这条对偶不变量在 s 域同形——改公式时别弄丢它。

| 文件 | 职责 | 知识篇目 |
|---|---|---|
| `include/biquad_filter.h`、`src/biquad_filter.c` | 单节 DF-II 核、`biquad_c2d_bilinear`；系数归一化、三道 fail-closed 闸、补偿求和 | `docs/biquad_filter.md` |
| `include/filter_utils.h`、`src/filter_utils.c` | 复数工具、四种频率变换、`zpk2sos`、共享管线 `design_filter` / `check_cascade_gains` | `docs/filter_utils.md` |
| `include/butter_filter.h`、`src/butter_filter.c` | Butterworth 原型表、各类型 init 辅助、X-macro 族 API | `docs/butter_filter.md` |
| `include/cheby_filter.h`、`src/cheby_filter.c` | Chebyshev 原型计算、边缘增益、X-macro 族 API | `docs/cheby_filter.md` |
| `include/notch_filter.h`、`src/notch_filter.c` | 二阶陷波：参数闸 + f32 范围闸 + 逐级闸 | `docs/notch_filter.md` |
| `include/peak_filter.h`、`src/peak_filter.c` | 二阶峰值：同三道闸；与 notch 严格互为倒数 | `docs/peak_filter.md` |
| `tests/` | C 单元测试、CSV 生成器、Python 黄金参考与精度验证；**实测标定与回归依据的存档处** | 见「知识地图」 |

命名：公开 API `{族}_{类型}_{阶数序数}_{init,update,reset}`；内部辅助函数
`{butter,cheby{1,2}}_{类型}_{func}`，直接收 `biquad_filter_t *sections`。
节数规则：LP/HP 为 `ceil(N/2)`，BP/BS 为 `N`（频率变换使有效阶数翻倍）。

## 不变量清单

**fail-closed 与 `num_sections`**：设计失败一律 `valid = 0` **且** `num_sections = 0`，
`_update` 直通返回输入。init 的每条早退路径都要显式清零 `num_sections`（调用方可能漏检
`valid`）——新增分支时照抄该模式。

**NaN 拦截写法**：入口闸写 `!(x > 0)`，不要写 `x <= 0`（NaN 比较恒假，后者放行 NaN）。
notch/peak 没有下游 `design_filter` 的有限性闸兜底，入口必须拦；butter/cheby 有
`design_filter` 的 `k` 有限性闸兜底。

**biquad 的三层闸**（`src/biquad_filter.c`）：系数有限性 → Jury 三条件（补偿求和 `sum3f`）
→ 极点半径裕量（半径 > 0.99995 拒绝；按半径本身判、不用 `a2`；判别式用 Dekker 分裂补偿，
两边平方避免 `sqrtf`）。任一失败 → 换单位直通并返回 0。这三条 Jury 判据在当前标定下被
半径裕量闸蕴含（Jury 变红需要一个极点落在单位圆外，而那种半径必被裕量闸拒），所以它们的
价值是**教科书写法与独立契约**，不是第二道有效防线：别为了"让 Jury 也能被单独测到"去改实现。

**`biquad_filter_reset` 是公开结构体上的公开 API**：系数不必先过 init。分母为 0 / 非有限、
或 `w_ss` 溢出时必须清零状态，永不产生 inf/NaN。

**`design_filter` 是导出符号**：栈上工作数组容量固定 `ZPK2SOS_MAX_N = 16`，
`degree = np − nz` 是 `uint8_t`。入口必须拦 `np == 0 || np > 16 || nz > np`，否则数组越界 /
无符号下溢。**BP/BS 另有更紧的上限 `np > ZPK2SOS_MAX_N/2`**：这两种类型会把每个原型零极点
一分为二原地写回同一个 16 元素数组，原型 9~16 阶会直接越界写栈（ASan 已复现）。
`proto_zeros == NULL` 只在 `nz == 0` 时合法。

**`zpk2sos` 的两个容差不可合并**：`eps_class = 1e-5`（实/复分类）、`eps_claim = 1e-3`
（共轭认领盒）。认领取盒内**最近**未用点，不是第一个（否则高 Q 簇里会偷别对的伴侣）。
循环末尾两条不变量——全部零极点被认领 + 每节根与输入多重集匹配——失败即返回 0，
不要"修好继续"。

**热路径保持 inline**：`biquad_filter_update`、`biquad_cascade_*` 与各阶 `_update` /
`_reset` 全在头文件内 `static inline`，节数以编译期字面量传入。不要挪成 out-of-line 函数。

**libm 只在非热路径用**：`_update` 路径零 libm（含全族头文件 inline 与共享
`biquad_cascade_*`）；`init` / `reset` 需要 libm（`tanf` / `powf` / `fabsf` / `isfinite`
等）。新增热路径代码不要引入 libm。

**禁止 `-ffast-math`**：`-ffast-math` 隐含 `-ffinite-math-only`，会把 `isfinite` 折叠成
恒真，让全部有限性闸静默失效（实测：NaN 系数被当作 `valid=1` 部署出 nan 输出）。
`include/biquad_filter.h` 有编译期 `#error` 守卫，是六个头的唯一公共顶点。

**改阈值先读标定**：闸门常数（`NOTCH_GATE` / `PEAK_GATE`、5e-5 半径裕量、增益窗口
±0.1 / ±0.25）的标定依据与边界档位在 `tests/test_notch.c` / `tests/test_peak.c` /
`tests/test_cheby.c`，改常数就要同步它们。`NOTCH_GATE`（`src/notch_filter.c`）与
`PEAK_GATE`（`src/peak_filter.c`）**各自独立定义、取值恰好相同（5e-8）**——改一个必须
同步另一个，否则两条对偶路径的拒绝域会分叉。

**测试约定**：频响测量用相位累加器 + RMS 比值（理由见 `tests/test_notch.c` 注释），不要改回
`sinf(2πfn/fs)` 或改用峰值；频响测量函数里无效滤波器返回 `1.0f` 让检查失败，而不是
`0.0f` 空过。

## 已知边界（被拒是预期行为，别去"修"）

- **极窄带 / 近 DC / 近 Nyquist** 会被拒（`valid = 0`，直通），主因是 f32 系数量化让
  DC/纹波边缘增益塌缩，不是极点半径。边界非单调，以 `valid` 为准。
- **Notch / Peak 另有 f32 设计参数范围闸**（`ξ·g·(ω0/K)² ≥ 5e-8`，Peak 用 `d = 1/g`）：
  低于阈值的参数会被**拒绝**（`valid = 0`），不是静默失真放行。这是**预测式**闸门：在被拒的
  参数上求 |H(jω0)| 本身就是那个病态运算（谷底/峰顶失真而 DC/Nyquist 增益仍然精确是 1），
  所以只能按参数判、且必须在设计之前判。Peak 通常由极点半径闸先触发。
- **`g` 浅于 0.001 dB**（notch `g ≥ 0.9999` / peak `g ≤ 1.0001`）按直通部署——数学上就是
  恒等，不要为了"完整性"烧一个高 Q 节。
- **超宽带 BP/BS**（`fc2/fc1 ≳ 1000`）的阻带地板 ~−12 dB、**reset 后的瞬态**（notch/peak 的
  `w_ss` 达 1e5 量级）都是结构固有，别当 bug 修。

## 知识地图

设计机理、推导、被否方案 → `docs/` 六篇（与源文件同名对应）：`biquad_filter.md`（三层闸与
`sum3f`）、`filter_utils.md`（频率变换、`zpk2sos` 配对、`design_filter` 边界）、
`butter_filter.md` / `cheby_filter.md`（原型、增益折叠、边缘增益）、
`notch_filter.md` / `peak_filter.md`（预测式参数闸、对偶不变量）。

实测标定与回归依据 → `tests/`：

- `tests/test_notch.c`「数值闸边界」段——Q 标定表、标定必须在邻域取最坏的告诫。
- `tests/test_peak.c`「数值闸边界」段——两条闸的先后关系、判据用 `d = 1/g` 的理由、
  "不要用 f64 解析求值做标定"的口径告诫。
- `tests/test_butter.c` / `tests/test_cheby.c`——极点半径裕量、超宽带与近 Nyquist
  回归、配对缺陷（曾产出 DC 增益 97.9 与 350× 谐振）、增益窗口标定。
- `tests/test_biquad.c`——Jury 与极点半径裕量的边界档位、NaN 毒化语义。

用户向用法、与 scipy 的精度对比、MCU 部署、文件树 → `README.md`。
