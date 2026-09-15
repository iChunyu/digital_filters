# CLAUDE.md

本文件为 Claude Code (claude.ai/code) 在此仓库中工作时提供指引。

## 代码约定

**注释一律用中文**——包括 `.c` / `.h` 的 `/* */` 与 doxygen 块、构建脚本的 `#`、
Python 的 `#` 与 docstring。代码标识符、数学记号、单位、外部专有名词（`scipy`、
`Butterworth`、`Dekker`、函数名、类型名）保持原样，不要硬译。

正文用中文说明，标识符/公式原样嵌入，例如：
`/* 极点倒序遍历变换，避免覆盖尚未读取的元素。 */`

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
+ zeros 128 B + sos 192 B ≈ 448 B，加上调用侧原型数组 ~128 B；该数值在 x86-64 -O2 下
实测成立，-O0 下约 1.1 KB）。ARM Cortex-M4 `-Os -fstack-usage` 实测最坏一条链为
`cheby2_lp_init` 208 B + `design_filter` 544 B + `zpk2sos` 192 B ≈ 944 B。
`notch_init` / `peak_init` 走的是直写 biquad 的短路径（无零极点工作数组），
x86-64 -O2 实测均为 **112 B**，是全场最小的一条 init 链。
运行时 `_update` 全 inline（biquad → 级联 → 逐阶函数三层全部
头文件内联），无额外栈开销。建议 MCU 主栈 ≥ 2 KB（两平台实测峰值均有 2× 以上余量）。

**libm 依赖**: `biquad_filter_init`/`_update`/`_reset` 路径零 libm 调用（裕量为纯乘加
代数式），裸机不链 libm 也可用 biquad 层。Butterworth init 无需 libm（ROM 查表）。
Chebyshev init 需要 `logf/sqrtf/sinhf/coshf/cosf/sinf`（仅 init 时）。
Notch / Peak init 需要 `prewarp` 的 `tanf`（仅 init 时，用于预畸）。

**Flash 粒度**: 库以 `-ffunction-sections/-fdata-sections` 编译；链接时加
`--gc-sections` 只拉入用到的滤波器族（实测 arm-none-eabi-gcc / Cortex-M4 / -O2：
库四个目标文件 .text 合计 19.8 KB，而只用 `butter_lp_2nd` 的程序链接后 text 仅
1.4 KB）。缺这些标志时链接器按目标文件粒度拉取，`butter_filter.o` 会整体
（4.3 KB，含全部 32 个 init）被拉入。
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

**Notch（陷波）**是独立的第四族，**不走这条流水线**：它是二阶的，只有一对
一阶分子/分母，没有零极点配对问题，init 直接 `biquad_c2d_bilinear` → 
`biquad_filter_init` 写进单个 `biquad_filter_t`。详见下文"Notch"一节。

**Peak（峰值）**是第五族，与 Notch **严格互为倒数**（`H_peak(ξ,g) ≡ 1/H_notch(ξ,1/g)`，
逐点恒等）。结构、三道闸、测试骨架全部与 Notch 同构，区别只有两处：`g > 1`，
以及全部数值判据用**等效陷波深度** `d = 1/g`。详见下文"Peak"一节。

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
同样只有一份。fail-closed 咽喉点（输入边界 np/nz、k 有限非零、节数上限、
逐节 init、增益窗口）全部单点维护。输入边界指 `np > ZPK2SOS_MAX_N` 或
`nz > np` 直接返回 0——`design_filter` 是导出符号，其栈上工作数组容量固定
（16）且 `degree = np - nz` 是 uint8_t，越界/下溢必须在入口拦下。

对外 exposed 的 update/reset 按阶数/类型分别定义（如 `butter_lp_2nd_update`、
`cheby1_bp_5th_reset`），调用时无需强转。它们是**头文件 X-macro 生成的
`static inline`**：调 `biquad_cascade_update(f->sections, ns, x)`，`ns` 是字面量
（循环边界编译期折叠）。每样本路径无函数调用、无运行时节数装载、无链接符号
（96 个字节级相同的 out-of-line 副本曾占 ~12.8 KB flash + 每样本一次 BL）。

- 只有 `_init` 按阶数分别定义，因为不同阶数需要不同大小的栈上临时数组
- `valid == 0` 时 `num_sections` 恒为 0——init 的每条早退路径都显式清零，
  调用方即使漏检 `valid` 也不会读到未初始化的节数
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

**Notch（`notch_filter.h/c`）**

唯一单阶数、无类型变体的滤波器，因此**不用 X-macro**：一个 `notch_filter_t`、
一个 `notch_init(f, f0, xi, g, fs)`。`sections[1]` 仍是标准 `biquad_filter_t`，
update/reset 复用 `biquad_cascade_*`。

```
H(s) = (s² + 2ξgω₀s + ω₀²) / (s² + 2ξω₀s + ω₀²)     ω₀ = 2πf₀
```

- **预畸只保中心频率与深度**。`ω₀ = 2π·prewarp(f0, fs)` 使数字域谷底精确落在
  `f0`、深度精确等于 `g`（全频段，含近 Nyquist）。但双线性对整个频率轴非线性，
  **宽度被扭曲**：数字域宽度 = `2ξg·f₀` 乘以随 `f0/fs` 变化的因子——`f0/fs ≲ 0.02`
  时 < 1%，`f0/fs = 0.2` 时窄 25%，再高完全走样。`xi` 始终按模拟原型带宽因子解释。
- **`g` 上限 `0.9999`**：`g = 1` 时分子分母逐项相等、`H ≡ 1`，数学上就是直通，
  按 `valid=0` 直通部署而不是烧一个节去实现恒等——后者要把高 Q 极点结构的状态
  推到 1.5e5 量级再精确相消回单位增益，纯属自找数值麻烦。
- **三道闸**（全部 fail-closed → `valid=0`、`num_sections=0`）：参数闸
  （`!(x > 0)` 式写法，同时拦 NaN）、f32 设计参数范围闸、`biquad_filter_init`
  + `check_cascade_gains` 的逐级闸。注意 butter 用的 `fc <= 0.0f` 拦不住 NaN
  （NaN 比较恒假），notch 没有下游 `design_filter` 兜底，必须在入口用 `!(x>0)`。

**f32 设计参数范围闸的来历**（这一段是本族最贵的知识，别丢）：

> **范围界定**：本库**不负责**修复/补偿/消除 f32 系数量化误差，那是数值类型的
> 固有代价、属于调用方的选择范围。这道闸不试图修好什么，只拒绝已知会**静默
> 失真**的参数范围，性质与 `biquad_filter_init` 的极点半径裕量闸相同。下面整段
> 只服务于"阈值该取多少"，不是"库要不要管量化误差"。

DF-II 下的陷波器有一个窄而可精确刻画的失效区。分子分母的 s⁰、s² 系数相同，
只有 s¹ 不同，故 `b1 == a1` **逐位相等**，于是

```
H(z) = 1 + c(z⁻² − 1) / (1 + a1z⁻¹ + a2z⁻²),   c = 2ξω₀K(1−g)/a0
```

**陷波深度只活在 c 里**，而 `c ~ 1e-4` 量级；b0/b2/a2 各带 ±1 ulp（≈1.2e-7）
绝对误差，对 c 的扰动是 0.1%。这点扰动被陷波频率处的极小分母放大：
`|D(e^{jθ₀})| ≈ 8π²ξ(f0/fs)²`（f0=20Hz/48kHz 算出 6.96e-7，数值求值实测
6.955e-7）。相对深度误差 `≈ δ/(g·|D|) ∝ 1/[ξg(f0/fs)²]`，放大 ~1e7 倍，
正好吃掉 f32 的全部 7.2 位十进制。

时域上等价：DF-II 状态放大 `1/(1+a1+a2) ≈ (fs/πf0)²` 倍（f0=20Hz/48kHz 时
1.5e5），输出要把 ~1.45e6 量级的三项相消到 `g`，而系数误差 `δb·|w| ≈ 0.17`
与真实输出同量级。**窄带低通没这个问题**——分子三项同号相加、输出与状态同量级；
陷波与峰值（互为倒数的同一病态）是全库仅有的两族在 DF-II 下病态的滤波器。

判据 `Q = ξ·g·(ω₀/K)² ≥ 5e-8`（纯乘除，无新 libm）。**标定必须在邻域上取最坏，
不能采孤立点**：失效在 0.002 Hz 尺度上剧烈震荡（f0 = 19.95~20.05 Hz @48kHz,
ξ=0.05, g=0.1 区间内比值在 1.04 和 3.04 之间随机跳），那些"好"的值只是舍入运气。
按孤立点标定会得出乐观 25~50 倍的阈值。实测邻域最坏值：Q=2e-9 → 12.1×，
5e-8 → 1.17×，1e-7 → 1.04×。实用刻度：ξ=0.05、g=0.1 时需 `f0/fs ≥ 1.007e-3`
（48kHz 下 f0 ≥ 48 Hz，8kHz 下 f0 ≥ 8 Hz）；工频 50/60Hz @48kHz 在内。

**已被否掉的方案**（别再试）：旁通形式（深度改由 c 单点承载）与耦合型
（coupled form）在**同一个** f0/fs 处一起失效；耦合型把状态从 1.4e6 压到
3.7e3，包络分毫未动——病根是系数量化，不是状态溢出。耦合型若用 `acos` 从
a1 反解角度会更加病态（θ 误差 0.9%），必须直接用模拟原型算 θ、r，即使用对了
可用下限也只比 DF-II 低约 25 倍（f0/fs ≳ 4e-5，代价是库中第一个非 DF-II
结构、每样本 7 次乘法 vs 5 次）。真需要那么低的 f0/fs 时再回头捡这段。

**Peak（`peak_filter.h/c`）**

Notch 的**严格对偶**，全库唯一能与另一族做系数级对照的滤波器：

```
H(s) = (s² + 2ξω₀s + ω₀²) / (s² + 2ξ/g·ω₀s + ω₀²)     ω₀ = 2πf₀，g > 1
```

- **为什么分母带 g（H1）而不是分子带 g（H2）**：`H1(ξ,g) ≡ 1/H2(ξ,1/g)` 逐点恒等
  （数值验证 2e-15），所以 g 的 boost 就是已验证的 (1/g) 陷波的镜像——同一个 ξ 下
  boost/cut 的相对带宽一致。H2 做 boost 会宽 g 倍：g=10、ξ=0.05、f₀=20 Hz 时
  √g 电平全宽 H1 是 0.632 Hz、H2 是 6.32 Hz，且 H2 在 f₀/2 处已 +1.58 dB
  （H1 只有 +0.02 dB）。若哪天要"Q 与增益无关"的经典 constant-Q peaking EQ，
  两个都不是，得换 RBJ 对称型（√g 平分给分子分母，√g 电平处全宽恒为
  `2ξ·f₀`，与增益无关）。
- **全部数值判据用统一参数 `d = 1/g`（等效陷波深度）**："g 的 boost 有多难做"与
  "深度 d = 1/g 的陷波有多难做"是同一个问题。闸门判据
  `Q = ξ·d·(ω₀/K)² ≥ 5e-8`，**不能用 g**：用 g 会被放大 g² 倍，g=10 时 48 kHz 下
  20 Hz 会"过闸"，而它的镜像陷波在同一点是被拒的。
- **半径闸几乎处处比数值闸紧**——与 Notch 恰好相反。Notch 的极点阻尼是 ξ（与 g
  无关），Peak 的是 ξ/g。两条边界 `2ξd·r ≥ 5e-5`（半径）与 `ξd·r² ≥ 5e-8`（数值）
  在 `ξd = 0.0125` 处交叉：ξd < 0.0125（绝大多数实用参数）半径闸先触发，
  ξd > 0.0125 数值闸先触发。所以 boost 的实际上限是
  `g ≤ 4e4·ξ·tan(πf0/fs)`（ξ=0.05：f0/fs=0.05 → g≈317；f0/fs=2e-3 → g≈12.6）。
  测试里两条闸各有专用档位（ξ=0.5、f0≈0.28~0.42 Hz 才让数值闸露出来）。
- **`g ≤ 1.0001` 按直通**（镜像 notch 的 `g ≥ 0.9999`）。判据 `!(d < 0.9999f)`；
  `1.0f/1.0001f` 在 f32 下恰好舍入到 `0.9999f`，等号即拒，与 notch 严格镜像。
- **f32 实测标定**（f32 系数 + f32 DF-II 状态递推，扫 f0 的 ±1% 邻域取最坏；
  f0=20 Hz、ξ=0.05、g=10）：1 kHz 误差 < 3e-5、8 kHz 2.8e-4、48 kHz 峰高
  7.11~11.06（理想 10，被半径闸拒）。镜像核对：同点 g=0.1 的陷波深度
  0.126~0.166。**标定必须走时域**——"f32 系数 + f64 解析求值"给 0.4%，低估两个
  数量级：它漏掉 f32 **状态**舍入，而这里 `|w| ≈ 1/|D(e^{jθ0})|` 到 1e7 量级，
  其 ulp 本身就和输出同量级。系数误差只是病态的一半。
- **对偶不变量测试**：`peak_init(f0,ξ,g,fs)` 与 `notch_init(f0,ξ,1/g,fs)` 的 z 域
  系数向量互为交换（只差公共归一化因子 `1/notch.num_z[0]`），在 50~400 Hz 扫掠上
  逐点校验；`src/peak_filter.c` 把分母写成 `2·xi·d·w0`（而非 `2·xi·w0/g`）正是
  为了让它在 f32 下逐位成立。CSV 层面同样对得上：两族配置表互为倒数，
  f32 实测 `max|peak·notch − 1| = 2.2e-3`。
- init 走与 notch 相同的直写 biquad 短路径，栈同为 112 B。

### 文件

| 文件 | 用途 |
|---|---|
| `include/biquad_filter.h` | Biquad 公开 API（inline update + 级联核心 cascade_update/reset） |
| `include/filter_utils.h` | 设计工具 (complex_t, prewarp, 变换, zpk2sos, design_filter, check_cascade_gains) |
| `include/butter_filter.h` | Butterworth API（X-macro 生成 32 结构体 + inline update/reset） |
| `include/cheby_filter.h` | Chebyshev I & II API（X-macro 生成 64 结构体 + inline update/reset） |
| `include/notch_filter.h` | Notch API（单结构体 + inline update/reset；预畸语义与三道闸的文档） |
| `include/peak_filter.h` | Peak API（与 Notch 同构；H1 形式的来历、d = 1/g 判据、boost 上限的文档） |
| `src/biquad_filter.c` | Biquad 实现（init, reset, get_output, get_input, c2d_bilinear） |
| `src/filter_utils.c` | 设计工具实现 + 共享设计管线 design_filter + 增益校验 |
| `src/butter_filter.c` | 预计算极点表、按阶 init（管线走 design_filter） |
| `src/cheby_filter.c` | 运行时原型计算、按阶 init（管线走 design_filter） |
| `src/notch_filter.c` | 参数闸 + f32 数值闸 + 直写 biquad（不走 design_filter） |
| `src/peak_filter.c` | 同上，判据用 d = 1/g；含 H1 vs H2 的选型理由与 f32 时域标定 |
| `test/test_biquad.c` | Biquad 测试（含裕量对称性、积分器 reset、NaN 语义） |
| `test/test_butter.c` | Butterworth 测试（全类型、多阶数 + 全阶扫掠 + 近实配对/窄带回归） |
| `test/test_cheby.c` | Chebyshev I/II 测试（全类型、多阶数 + 全阶扫掠 + 高 Q 配对回归） |
| `test/test_notch.c` | Notch 测试（参数闸、数值闸边界、深度/DC/Nyquist、g 上限、reset、f 与 xi 扫掠） |
| `test/test_peak.c` | 同上镜像 + **对偶不变量**（与 notch_init(1/g) 的系数互为交换）+ 半径闸/数值闸各一档 |
| `test/test_butter_with_py.c` | 生成 CSV（argv 指定路径 + valid 检查）与 scipy 对比 |
| `test/test_cheby_with_py.c` | 生成 CSV（argv 指定路径 + valid 检查）与 scipy 对比 |
| `test/test_notch_with_py.c` | 生成**频响** CSV（25 个相对 f0 的频点 × 4 配置；valid 失败拒绝出表） |
| `test/test_peak_with_py.c` | 同上，配置表是 notch 表的逐项倒数（g → 1/g），频率网格相同 |
| `test/test_butter_use_py.py` | Python 参考滤波器（Butterworth，含绘图，argv 指定数据目录） |
| `test/test_cheby_use_py.py` | Python 参考滤波器（Chebyshev，含绘图，argv 指定数据目录） |
| `test/test_notch_use_py.py` | Notch f64 参考频响（scipy.signal.bilinear 独立实现 c2d）对拍；无 numpy/scipy 时 exit 77 |
| `test/test_peak_use_py.py` | Peak f64 参考频响同上；容差按 max\|H_ref\| 归一（峰高到 50，绝对容差无意义） |
| `test/compare_scipy.py` | 稳态精度对比（已注册 ctest，无 numpy/scipy 时 exit 77 SKIP） |
| `test/verify_zpk_gain.py` | 验证 float32 精度足够 zpk 增益追踪（独立复现 scipy 管线，f64 vs f32） |
| `CMakeLists.txt` | 顶层 CMake（库带 -ffunction-sections/-fdata-sections） |
| `test/CMakeLists.txt` | 测试可执行文件 + CTest 注册（13 项，含 CSV 生成 + scipy 对比 + notch/peak 频响对拍 + zpk 增益精度验证） |
