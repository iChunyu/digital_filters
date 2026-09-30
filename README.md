# digital_filters — 数字滤波器库（C 语言）

> MCU 用的 IIR 滤波器库。零 `malloc`，零 `double`，全部 `float` 一把梭。
> 和 scipy 对过答案了，稳态误差 1e-6 ~ 1.5e-5 🤏

## 项目简介

基于 **双二阶 (biquad)** 滤波器及其级联（SOS，二阶节）构建高阶 IIR 滤波器。
全部采用 **Direct Form II（规范型）**——每个 biquad 只要 3 个 `float` 状态变量，
是你能写出来的最省内存的 IIR 实现。

所有结构体由 X-macro 生成，阶数编译时确定，biquad 节内嵌在结构体里。
不需要 `malloc`，不需要 `free`，离开作用域自动回收。`<stdlib.h>` 都不用 include。
（Notch / Peak 是固定二阶，各一个结构体，不走 X-macro。）

### 支持的滤波器

| 族 | 通带 | 阻带 | 一句话 |
|---|---|---|---|
| **Butterworth** | 最大平坦 | 单调衰减 | 老实人，不搞花活 |
| **Chebyshev Type I** | 等波纹（指定纹波 dB） | 单调衰减 | 通带里蹦迪，阻带装死 |
| **Chebyshev Type II** | 单调 | 等波纹（指定最小衰减 dB） | 反过来，通带佛系阻带蹦迪 |
| **Notch（陷波）** | 二阶，f0 处衰减到 g | — | 只掐一个频点，别处不动 |
| **Peak（峰值）** | 二阶，f0 处提升到 g | — | 陷波器的倒数，专治"这一段没劲" |

四种类型：**低通 (LP)**、**高通 (HP)**、**带通 (BP)**、**带阻 (BS)**。
原型阶数 1~8（BP/BS 有效阶数翻倍，最高 16 阶）。

Notch / Peak 是独立的二阶族，参数是 `(f0, xi, g, fs)`：`g < 1` 用
`notch_init`（线性衰减），`g > 1` 用 `peak_init`（线性增益）；两者满足
`H_peak(ξ,g) ≡ 1 / H_notch(ξ,1/g)`，boost 与 cut 的带宽语义对称。

## 快速开始

### 构建

```bash
cmake -B build -DBUILD_TESTS=ON
cmake --build build
```

### 运行测试

```bash
cd build && ctest --output-on-failure
```

预期输出：**13 项测试全部通过**（5 项 C 单元/回归测试 + 4 项 CSV 生成 +
3 项黄金参考对比 + 1 项 zpk 增益精度验证）。4 项 CSV 生成只要求配置期能找到
Python3；4 项 Python 对比在缺 numpy/scipy 时以退出码 77 自动 SKIP。

### 基本用法

```c
#include "butter_filter.h"

// 二阶低通，截止 2 Hz，采样 20 Hz
butter_lp_2nd_t filter;
butter_lp_2nd_init(&filter, 2.0f, 20.0f);

// 必须检查 valid：设计失败时滤波器是直通（H(z)=1），
// update 会原样返回输入。宁可直接过，也不要发散的输出。
if (!filter.valid) {
    // 按产品策略处理：换参数重试 / 报错 / 继续跑直通
}

// 跳过烦人的起振瞬态
butter_lp_2nd_reset(&filter, 1.0f);

for (int i = 0; i < 1000; i++) {
    float y = butter_lp_2nd_update(&filter, input[i]);
}

// 不用 destroy，离开作用域自动释放 😌
```

带通：

```c
butter_bp_2nd_t bp;
butter_bp_2nd_init(&bp, 2.0f, 5.0f, 40.0f);  // fc1, fc2, fs
butter_bp_2nd_reset(&bp, 0.0f);
float y = butter_bp_2nd_update(&bp, x);
```

Chebyshev Type I（带纹波）：

```c
cheby1_lp_3rd_t c1;
cheby1_lp_3rd_init(&c1, 3.0f, 20.0f, 0.5f);  // fc, fs, ripple_dB
cheby1_lp_3rd_reset(&c1, 1.0f);
float y = cheby1_lp_3rd_update(&c1, x);
```

Chebyshev Type II（阻带衰减）：

```c
cheby2_hp_2nd_t c2;
cheby2_hp_2nd_init(&c2, 5.0f, 40.0f, 40.0f);  // fc, fs, stopband_dB
cheby2_hp_2nd_reset(&c2, 0.0f);
float y = cheby2_hp_2nd_update(&c2, x);
```

陷波（工频 50 Hz @ 1 kHz，深 −20 dB ⇒ `g = 0.1`）：

```c
#include "notch_filter.h"

notch_filter_t n;
notch_init(&n, 50.0f, 0.05f, 0.1f, 1000.0f);  // f0, xi, g, fs
if (!n.valid) { /* 参数被闸拒（如 f0/fs 太小）→ 直通 */ }
notch_reset(&n, 0.0f);
float y = notch_update(&n, x);
```

峰值（同频点 +20 dB ⇒ `g = 10`）。与陷波**逐点互为倒数**，语义完全对称：

```c
#include "peak_filter.h"

peak_filter_t p;
peak_init(&p, 50.0f, 0.05f, 10.0f, 1000.0f);  // f0, xi, g, fs
if (!p.valid) { /* 参数被闸拒 → 直通 */ }
peak_reset(&p, 0.0f);
float y = peak_update(&p, x);
```

### 命名规则

```
{族}_{类型}_{阶数序数}_t

族:   butter, cheby1, cheby2
类型: lp, hp, bp, bs
阶数: 1st ~ 8th
```

示例：`butter_lp_2nd_t`, `cheby1_bp_5th_t`, `cheby2_bs_3rd_t`

对应的函数：
- `{族}_{类型}_{阶数序数}_init(f, ...)`
- `{族}_{类型}_{阶数序数}_update(f, input)`
- `{族}_{类型}_{阶数序数}_reset(f, equilibrium)`

## 阶数与节数

| 类型 | 原型阶数 N | 有效阶数 | biquad 节数 |
|---|---|---|---|
| LP, HP | N | N | ceil(N/2) |
| BP, BS | N | 2N | N |

8 阶 BP → 16 阶有效滤波器，8 个 biquad 节。

## 设计流水线

```
1. 模拟原型极点（Butterworth：ROM 查表；Chebyshev：运行时算）
2. 模拟频率变换（LP: 缩放; HP: 倒数; BP/BS: 阶数翻倍）
3. 双线性变换 (s → z)
4. 零点补齐
5. 零极点配对 → biquad 系数
6. 部署到内嵌 biquad 节
```

预畸变：`f_analog = fs/π · tan(π · f_digital / fs)`

## 与 scipy 的对比验证

C 代码生成 CSV → scipy 做黄金参考 → 对比稳态精度（跳过瞬态取最后 10%，
已注册进 ctest，装好 numpy/scipy 即自动运行）：

| 滤波器 | LP | HP | BP | BS |
|---|---|---|---|---|
| Butterworth | 0.9e-6 | 2.1e-6 | 2.7e-6 | 1.1e-6 |
| Chebyshev I | 1.2e-6 | 3.6e-6 | 1.5e-5 | 1.1e-5 |
| Chebyshev II | 1.7e-6 | 2.1e-6 | 2.5e-6 | 1.1e-6 |

> C 的 `reset(equilibrium)` 和 scipy 的零初态起点不同，前面几百个采样对不上是
> 正常的（瞬态响应差异）。上表取的是稳态数据。

## 注意事项

### MCU 使用

- **零 `malloc`**：`<stdlib.h>` 不需要，堆管理器关掉照样跑
- **只有 biquad 层不依赖 libm**：`biquad_filter_init`/`update` 是纯乘加代数
  形式，绝对值用本地 `abs_f` 而非 `fabsf`，不依赖编译器内建——裸机固件不链 libm
  也能用它（`-fno-builtin -ffreestanding` 下 `nm -u biquad_filter.o` 实测为空）
- **所有族的 init 都要 libm**：设计管线经 `prewarp()` → `tanf()` 预畸，
  `filter_utils.o` 另需 `sqrtf`/`fmaxf`/`memcpy`/`memset`。Butterworth 的 ROM
  极点表只免掉原型计算的 `cosf`/`sinf`，**不减免 libm 链接**
- Chebyshev init 另需 `powf`/`logf`/`cosf`/`sinf`/`sinhf`/`coshf`（仅 init 一次，
  非逐采样）；`powf` 会连带拉入 errno 版数学内核（`__errno` 等），不含 reent
  支持的裸机工程需自行提供这些符号
- `_update`/`_reset` 全部为**头文件 `static inline`**：每样本路径是带字面量
  节数的 biquad 级联循环，没有函数调用、没有运行时节数装载（需 −O1 及以上；
  −O0 下函数体仍会生成，每样本 1~3 次 BL）
- 全部 `float`，零 `double`
- **flash 粒度**：库以 `-ffunction-sections/-fdata-sections` 编译，链接时
  加 `--gc-sections`（或对应链接器选项）后按函数粒度拉取。实测
  arm-none-eabi-gcc 16.2 / Cortex-M4 / `.text`+`.rodata`：六个目标文件
  −Os 合计 **14.7 KB**、−O2 合计 **18.2 KB**；分项（−Os / −O2）——
  `biquad_filter.o` 0.9/1.0 KB、`filter_utils.o` 4.9/5.8 KB、
  `butter_filter.o` 2.6/4.5 KB、`cheby_filter.o` 5.7/6.2 KB、
  `notch_filter.o` 与 `peak_filter.o` 各 0.3 KB
- **`--gc-sections` 裁不掉设计管线**：`design_filter` 的运行时 `switch(type)`
  引用全部四种频率变换，只用 `butter_lp_2nd` 的最小程序也要保留
  `biquad`+`filter_utils`+`butter` 三个目标文件（−Os 全量 8.4 KB，链接并 GC 后
  实测 ≈6.8 KB；链 newlib libm 后 ≈11.3 KB）。缺 `--gc-sections` 时链接器按
  目标文件粒度拉取，`butter_filter.o` 整体（2.6 KB@−Os / 4.5 KB@−O2，含全部
  32 个 init）被拉入

### 参数校验与 fail-closed 语义

- 原型阶数 1~8
- 截止频率 `0 < fc < fs/2`
- BP/BS 需 `fc1 < fc2` 且 `fc2 < fs/2`
- Chebyshev 的 `ripple_db` > 0
- Notch：`0 < f0 < fs/2`、`xi > 0`、`0 < g < 0.9999`；Peak：同前但 `g > 1.0001`
  （`g` 贴近 1 时峰/谷浅于 0.001 dB，按直通部署）。两者另有 f32 数值包络闸：
  可用区间随 `f0/fs` 与 `g` 收窄，被拒时 `valid = 0`（见头文件里的实测刻度）
- 任一校验失败 → `valid = 0`，update 直通返回输入。**调用方必须检查
  `valid`**——直通是"宁可不过滤也不要错误输出"的兜底，不是静默成功的保证

### 稳定性（三层 fail-closed 防线）

1. 每节 biquad：三个 Jury 条件（`sum3f` 补偿求和只回收两次加法的舍入，
   窄带系数下与裸和逐位相同；和恰好为 0 = 量化后的极点确实落在 z=±1 → 拒绝）
2. 极点半径裕量：任何极点半径 > 0.99995（距单位圆 < 5e-5，会振铃
   ≥ 10⁴ 采样）→ 拒绝。按极点半径本身判定（对称覆盖共轭对、异号实根
   对、非等实根对的主导极点——`a2 = r1·r2` 乘积检查有盲区）
3. 级联 DC/Nyquist 增益校验：部署后解析 H(0)/H(π) 必须落在期望窗口内
   （结构增益 ±0.1、纹波边缘 ±0.25）——拦截"稳定但配错对"的静默错误
   （实测曾拦下 DC 增益 97.9、350× 谐振的错误滤波器）

### 频率包络（设计被拒 = 直通，属预期行为）

- **极窄带 / 近 DC**：被拒的主因是 f32 系数量化让 DC/纹波边缘增益塌缩，
   **不是**极点半径——`butter_lp_2nd`@48k 的 2 万点扫掠里极点半径闸一次都没
  触发。实测（fs=48 kHz）：LP2 可用下界 ≈2.5 Hz、被拒设计全部落在 5.6 Hz 以下，
  LP1 在 0.5 Hz 已可用。边界不是单调整数——LP2 的 2.5~6 Hz 区间以 0.001 Hz
  步长扫会出现 valid 翻转（3500 点中 1469 个被拒）。不同阶数/族/类型边界各不
  相同，一律以 `valid` 为准
- **极近 Nyquist**：例如 BS1 [100, 23990]@48k 可用，[100, 23999.8]@48k
  被拒（极点半径 0.99997）
- **超宽带 BP/BS**（fc2/fc1 ≳ 1000）：极点贴近 z=1 使内部状态巨大
  （w ≈ 1/(1+a1+a2)），f32 状态更新的固有噪声把阻带衰减地板抬高到
  ~−12 dB 量级（本设计类固有，非缺陷）
- **非有限输入**：NaN/Inf 会毒化状态、持续输出 NaN 直到 reset。热路径
  刻意不做防护（每样本分支开销）；传感器/不可信数据请在源头清洗

### 数值精度

- 全 `float`。8 阶以内完全够用，加 double 只增加 ROM/RAM 负担
- "最不利极点优先"配对策略 + 最近共轭认领，减少有限精度舍入噪声
- 窄带高 Q 场景建议实测评估
- scipy 对比（ctest 自动跑）：稳态误差 9e-7 ~ 1.5e-5

## 已知局限

- 原型阶数 8 阶上限（ROM 表 + 结构体枚举的工程约束）
- f32 系数在极窄带/极近 Nyquist 设计上无法忠实表达（见"频率包络"）——
  这些设计被 fail-closed 拒绝而不是硬部署
- 不支持椭圆滤波器（Jacobi 椭圆函数写起来太抽象了，下次一定）

## 文件结构

```
digital_filters/
├── include/
│   ├── biquad_filter.h          # 单节 biquad（含 inline update）
│   ├── filter_utils.h           # 设计工具 + 类型定义
│   ├── butter_filter.h          # Butterworth API
│   ├── cheby_filter.h           # Chebyshev I & II API
│   ├── notch_filter.h           # 陷波 API（f0, xi, g<1, fs）
│   └── peak_filter.h            # 峰值 API（f0, xi, g>1, fs；陷波的倒数）
├── src/
│   ├── biquad_filter.c
│   ├── filter_utils.c
│   ├── butter_filter.c
│   ├── cheby_filter.c
│   ├── notch_filter.c
│   └── peak_filter.c
├── test/
│   ├── CMakeLists.txt
│   ├── test_biquad.c
│   ├── test_butter.c
│   ├── test_cheby.c
│   ├── test_notch.c
│   ├── test_peak.c
│   ├── test_butter_with_py.c     # 生成 CSV 与 scipy 对比
│   ├── test_cheby_with_py.c
│   ├── test_notch_with_py.c
│   ├── test_peak_with_py.c
│   ├── test_butter_use_py.py     # Python 参考滤波器（含绘图）
│   ├── test_cheby_use_py.py
│   ├── test_notch_use_py.py
│   ├── test_peak_use_py.py
│   ├── compare_scipy.py          # 稳态精度对比
│   └── verify_zpk_gain.py        # float32 精度验证（独立复现 scipy 管线）
├── CMakeLists.txt
├── AGENTS.md                    # agent 指引
└── README.md
```
