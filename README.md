# digital_filters — 数字滤波器库（C 语言）

> MCU 用的 IIR 滤波器库。零 `malloc`，零 `double`，全部 `float` 一把梭。
> 和 scipy 对过答案了，稳态误差在 1e-6 量级 🤏

## 项目简介

基于 **双二阶 (biquad)** 滤波器及其级联（SOS，二阶节）构建高阶 IIR 滤波器。
全部采用 **Direct Form II（规范型）**——每个 biquad 只要 3 个 `float` 状态变量，
是你能写出来的最省内存的 IIR 实现。

所有结构体由 X-macro 生成，阶数编译时确定，biquad 节内嵌在结构体里。
不需要 `malloc`，不需要 `free`，离开作用域自动回收。`<stdlib.h>` 都不用 include。

### 支持的滤波器

| 族 | 通带 | 阻带 | 一句话 |
|---|---|---|---|
| **Butterworth** | 最大平坦 | 单调衰减 | 老实人，不搞花活 |
| **Chebyshev Type I** | 等波纹（指定纹波 dB） | 单调衰减 | 通带里蹦迪，阻带装死 |
| **Chebyshev Type II** | 单调 | 等波纹（指定最小衰减 dB） | 反过来，通带佛系阻带蹦迪 |

四种类型：**低通 (LP)**、**高通 (HP)**、**带通 (BP)**、**带阻 (BS)**。
原型阶数 1~8（BP/BS 有效阶数翻倍，最高 16 阶）。

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

预期输出：**7 项测试全部通过**（前 3 项为 C 单元/回归测试；后 4 项为 CSV
生成 + scipy 黄金参考对比 + zpk 增益精度验证，未安装 numpy/scipy 时自动 SKIP）。

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
- Butterworth init **不调 `cosf`/`sinf`**（ROM 查表）；`biquad_filter_init`
  / `update` 路径**完全不依赖 libm**（裕量检查为纯乘加代数形式，绝对值用本地
  `abs_f` 而非 `fabsf`，不依赖编译器内建）——裸机固件不链 libm 也能用
  biquad 层（`-fno-builtin -ffreestanding` 下 `nm -u biquad_filter.o` 实测为空）
- Chebyshev init 需 `logf`/`sqrtf`/`sinhf`/`coshf`（仅 init 一次，非逐采样）
- `_update`/`_reset` 全部为**头文件 `static inline`**：每样本路径是带字面量
  节数的 biquad 级联循环，没有函数调用、没有运行时节数装载
- 全部 `float`，零 `double`
- **flash 粒度**：库以 `-ffunction-sections/-fdata-sections` 编译，链接时
  加 `--gc-sections`（或对应链接器选项）后，只拉入实际用到的滤波器族
  （实测 arm-none-eabi-gcc / Cortex-M4 / -O2：库四个目标文件 .text 合计
  19.8 KB，而只用 `butter_lp_2nd` 的程序链接后 text 仅 1.4 KB）。缺这些标志
  时链接器按目标文件粒度拉取，`butter_filter.o` 会整体（4.3 KB，含全部
  32 个 init）被拉入

### 参数校验与 fail-closed 语义

- 原型阶数 1~8
- 截止频率 `0 < fc < fs/2`
- BP/BS 需 `fc1 < fc2` 且 `fc2 < fs/2`
- Chebyshev 的 `ripple_db` > 0
- 任一校验失败 → `valid = 0`，update 直通返回输入。**调用方必须检查
  `valid`**——直通是"宁可不过滤也不要错误输出"的兜底，不是静默成功的保证

### 稳定性（三层 fail-closed 防线）

1. 每节 biquad：三个 Jury 条件（f32 补偿求和消除窄带设计的误拒；
   和恰好为 0 = 极点在单位圆上 → 拒绝）
2. 极点半径裕量：任何极点半径 > 0.99995（距单位圆 < 5e-5，会振铃
   ≥ 10⁴ 采样）→ 拒绝。按极点半径本身判定（对称覆盖共轭对、异号实根
   对、非等实根对的主导极点——`a2 = r1·r2` 乘积检查有盲区）
3. 级联 DC/Nyquist 增益校验：部署后解析 H(0)/H(π) 必须落在期望窗口内
   （结构增益 ±0.1、纹波边缘 ±0.25）——拦截"稳定但配错对"的静默错误
   （实测曾拦下 DC 增益 97.9、350× 谐振的错误滤波器）

### 频率包络（设计被拒 = 直通，属预期行为）

- **极窄带 / 近 DC**：极点进入单位圆 5e-5 内会被拒。示例（fs=48 kHz）：
  LP2 fc ≥ ~1.2 Hz 可用、fc = 1 Hz 被拒；LP1 fc ≥ ~0.45 Hz 可用。
  不同阶数/族/类型边界不同，以 `valid` 为准
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
│   └── cheby_filter.h           # Chebyshev I & II API
├── src/
│   ├── biquad_filter.c
│   ├── filter_utils.c
│   ├── butter_filter.c
│   └── cheby_filter.c
├── test/
│   ├── CMakeLists.txt
│   ├── test_biquad.c
│   ├── test_butter.c
│   ├── test_cheby.c
│   ├── test_butter_with_py.c     # 生成 CSV 与 scipy 对比
│   ├── test_cheby_with_py.c
│   ├── test_butter_use_py.py     # Python 参考滤波器（含绘图）
│   ├── test_cheby_use_py.py
│   ├── compare_scipy.py          # 稳态精度对比
│   └── verify_zpk_gain.py        # float32 精度验证（独立复现 scipy 管线）
├── CMakeLists.txt
├── CLAUDE.md
└── README.md
```
