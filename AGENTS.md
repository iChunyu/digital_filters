# AGENTS.md

给在本仓库改代码的 agent 的指引。用户向的 API 用法与精度表见 `README.md`；
**数值推导与实测标定的权威位置在各源文件头部的注释**（见文末"深挖入口"）。
本文件只写导航、约定和动代码时必须守住的不变量，不重复那些推导。

## 硬约束

- 纯 C，面向 MCU / 嵌入式；**只有 `float`**，零 `double`，零 `malloc`
- 所有结构体阶数编译期确定，内嵌 `biquad_filter_t sections[]`，离开作用域即回收
- `_update` / `_reset` 是头文件 `static inline`；热路径无函数调用，也**刻意无
  NaN/Inf 防护**（每样本分支开销），非有限输入毒化状态直至 reset，调用方在源头清洗
- `_update` / `_reset` 不可重入：同一结构体被 ISR 与主循环共享时需关中断或双缓冲

## 代码约定

注释一律用中文（`/* */` 与 doxygen、构建脚本 `#`、Python `#` 与 docstring）。
代码标识符、数学记号、单位、专有名词（`scipy`、`Butterworth`、`Dekker`、函数名、
类型名）保持原样，不要硬译。

## 构建与测试

```bash
cmake -B build -DBUILD_TESTS=ON
cmake --build build
cd build && ctest --output-on-failure
```

预期 **13/13 通过**：5 项 C 单元/回归（`test_biquad`、`test_butter`、`test_cheby`、
`test_notch`、`test_peak`）+ 4 项 CSV 生成 + scipy/陷波/峰值三项黄金参考对比 +
`verify_zpk_gain`。后 4 项需要 Python3 + numpy/scipy（`verify_zpk_gain` 只需 numpy），
缺库时以退出码 77 自动 SKIP；配置期找不到 Python3 时 CSV 生成器不注册。
改动系数、配对或阈值后，这些对比会量化响应偏差。

`test_biquad` **故意不链接 `m`**——它是"biquad 层零 libm"这条承诺的回归测试，
别给它加 `-lm`。新增源文件或测试要同步 `CMakeLists.txt` / `test/CMakeLists.txt`。

## 架构地图

五族，两条设计流水线：

- **Butterworth / Chebyshev I / Chebyshev II** — 高阶，走共享管线
  `design_filter()`：模拟原型 → 频率变换 → 增益折叠 → 双线性 → `zpk2sos()` → 部署。
  Butterworth 传 ROM 极点表（`butter_proto[8][8]`，~512 B）、`k=1`、`nz=0`；
  Chebyshev 运行时算原型（依赖 ripple），传自己的 `k` 与有限零点。
  结构体与 init/update/reset 由 X-macro 生成（`BUTTER_FIELDS` / `CHEBY_FIELDS`、
  `FOR_EACH_*_ORDER`），共 32 + 64 个结构体。
- **Notch / Peak** — 固定二阶、单节、无类型变体，因此**不用 X-macro**：
  预畸 → `biquad_c2d_bilinear` → `biquad_filter_init` 直写一个 `biquad_filter_t`，
  不走零极点配对。两族严格互为倒数（`H_peak(ξ,g) ≡ 1/H_notch(ξ,1/g)`），
  `peak_filter.c` 把分母写成 `2·xi·d·w0` 正是为了让这条对偶不变量在 s 域同形——
  改公式时别弄丢它。

| 文件 | 职责 |
|---|---|
| `include/biquad_filter.h`、`src/biquad_filter.c` | 单节 DF-II 核、`biquad_c2d_bilinear`；系数归一化、三道 fail-closed 闸、补偿求和 |
| `include/filter_utils.h`、`src/filter_utils.c` | 复数工具、四种频率变换、`zpk2sos`、共享管线 `design_filter` / `check_cascade_gains` |
| `include/butter_filter.h`、`src/butter_filter.c` | Butterworth 原型表、各类型 init 辅助、X-macro 族 API |
| `include/cheby_filter.h`、`src/cheby_filter.c` | Chebyshev 原型计算、边缘增益、X-macro 族 API |
| `include/{notch,peak}_filter.h`、`src/*.c` | 二阶族：参数闸 + f32 范围闸 + 逐级闸 |
| `test/` | C 单元测试、CSV 生成器、Python 黄金参考与精度验证 |

命名：公开 API `{族}_{类型}_{阶数序数}_{init,update,reset}`；内部辅助函数
`{butter,cheby{1,2}}_{类型}_{func}`，直接收 `biquad_filter_t *sections`。
节数规则：LP/HP 为 `ceil(N/2)`，BP/BS 为 `N`（频率变换使有效阶数翻倍）。

## 动代码时必须守住的

**fail-closed 与 `num_sections`**：设计失败一律 `valid = 0` **且**
`num_sections = 0`，`_update` 直通返回输入。init 的每条早退路径都要显式清零
`num_sections`（调用方可能漏检 `valid`）——新增分支时照抄该模式。

**NaN 拦截写法**：入口闸写 `!(x > 0)`，不要写 `x <= 0`（NaN 比较恒假，后者放行
NaN）。notch/peak 没有下游 `design_filter` 的有限性闸兜底，入口必须拦；
butter/cheby 里的 `fc <= 0.0f` 就是漏 NaN 的反例。

**biquad 的三层闸**（`src/biquad_filter.c`）：系数有限性 → Jury 三条件（补偿求和
`sum3f`）→ 极点半径裕量（半径 > 0.99995 拒绝；按半径本身判、不用 `a2`；判别式用
Dekker 分裂补偿；**init 路径禁止出现 `sqrtf`/libm**，绝对值用本地 `abs_f` 而非
`fabsf`）。任一失败 → 换单位直通并返回 0。

**`biquad_filter_reset` 是公开结构体上的公开 API**：系数不必先过 init。分母为 0 /
非有限、或 `w_ss` 溢出时必须清零状态，永不产生 inf/NaN。

**`design_filter` 是导出符号**：栈上工作数组容量固定 `ZPK2SOS_MAX_N = 16`，
`degree = np − nz` 是 `uint8_t`。入口必须拦 `np == 0 || np > 16 || nz > np`，
否则数组越界 / 无符号下溢。

**`zpk2sos` 的两个容差不可合并**：`eps_class = 1e-5`（实/复分类）与
`eps_claim = 1e-3`（共轭认领盒，须容纳双线性 f32 舍入 ~2e-4）。认领取盒内**最近**
未用点，不是第一个（否则高 Q 簇里会偷别对的伴侣）。循环末尾两条不变量——全部零极点
被认领 + 每节根与输入多重集匹配——失败即返回 0，不要"修好继续"。

**热路径保持 inline**：`biquad_filter_update`、`biquad_cascade_*` 与各阶
`_update` / `_reset` 全在头文件内 `static inline`，节数以编译期字面量传入。
不要挪成 out-of-line 函数（96 个字节级相同的副本曾占 ~12.8 KB flash，且每样本
多一次 BL）。

**libm 只在 init 用**：`_update` 路径零 libm；`biquad_filter.o` 整体零 libm；
其余族只有 init 需要 `tanf` / `powf` 等。新增热路径代码不要引入 libm。

**改阈值先读标定**：闸门常数（`NOTCH_GATE` / `PEAK_GATE` = 5e-8、5e-5 半径裕量、
增益窗口 ±0.1 / ±0.25）都有源码注释里的实测标定依据；`test_notch.c` / `test_peak.c`
有专门档位钉住边界（含让某一道闸"露出来"的构造），改常数就要同步它们。

**测试约定**：频响测量用相位累加器 + RMS 比值（理由见 `test_notch.c` 注释），
不要改回 `sinf(2πfn/fs)` 或改用峰值；频响测量函数里无效滤波器返回 `1.0f` 让检查
失败，而不是 `0.0f` 空过。

## 已知边界（被拒是预期行为，别去"修"）

- **极窄带 / 近 DC / 近 Nyquist** 会被拒，主因是 f32 系数量化让 DC/纹波边缘增益
  塌缩，不是极点半径。可用边界非单调（`lp_2nd`@48k 在 2.5–6 Hz 内以 0.001 Hz
  步长扫会翻转 11 次），以 `valid` 为准，不要用解析公式去预测。
- **Notch / Peak 另有 f32 设计参数范围闸**：`ξ·g·(ω0/K)² ≥ 5e-8`，Peak 用
  `d = 1/g`。被拒的参数上会**静默失真**（谷底变浅、峰高变矮、中心偏移而
  `valid` 仍为 1），这是预测式闸门，不是保守过度。Peak 通常由极点半径闸先触发，
  实用上限 `g ≤ 4e4·ξ·tan(πf0/fs)`。
- **`g` 浅于 0.001 dB**（notch `g ≥ 0.9999` / peak `g ≤ 1.0001`）按直通部署
  （`valid = 0`），因为数学上就是恒等——不要为了"完整性"烧一个高 Q 节去实现单位增益。
- **超宽带 BP/BS**（`fc2/fc1 ≳ 1000`）：DF-II 的 f32 状态更新噪声把阻带地板抬到
  ~−12 dB，属该结构固有，换实现也压不掉。
- **reset 后的瞬态**：notch/peak 的 `w_ss` 可达 1e5 量级（DC 增益为 1 的必然结果），
  f32 下复位后头几个输出带绝对误差，随后收敛。

## 深挖入口

- `src/notch_filter.c` / `src/peak_filter.c` 文件头——闸门推导、Q 标定表、
  邻域取最坏的告诫、已被否掉的方案（旁通形式 / 耦合型）。
- `src/biquad_filter.c`——Jury 与极点半径裕量的数值论证。
- `src/filter_utils.c`——`zpk2sos` 容差与配对缺陷史（曾产出 350× 谐振）。
- `README.md`——用户向用法、与 scipy 的精度对比、MCU 部署（libm 依赖、
  flash 分项、--gc-sections 实测）、文件树。
