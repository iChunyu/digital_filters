# digital_filters — 数字滤波器库（C 语言）

面向 MCU 的 IIR 滤波器库，纯 C 实现，零 `malloc`、零 `double`，全部 `float`。库基于
**双二阶 (biquad)** 滤波器及其级联（SOS，二阶节）构建高阶 IIR 滤波器，全部采用
**Direct Form II（规范型）**——每个 biquad 只要 3 个 `float` 状态变量。结构体阶数编译期
确定、biquad 节内嵌在结构体里，不需要 `malloc`/`free`，离开作用域自动回收。

每种类型都有 **低通 (LP)**、**高通 (HP)**、**带通 (BP)**、**带阻 (BS)**，原型阶数 1~8 阶
（BP/BS 有效阶数翻倍，最高 16 阶）；LP/HP 的 biquad 节数为 `ceil(N/2)`，BP/BS 为 N
（8 阶 BP → 16 阶有效、8 个 biquad 节）。

| 族 | 通带 | 阻带 |
|---|---|---|
| **Butterworth** | 最大平坦 | 单调衰减 |
| **Chebyshev Type I** | 等波纹（指定纹波 dB） | 单调衰减 |
| **Chebyshev Type II** | 单调 | 等波纹（指定最小衰减 dB） |
| **Notch（陷波）** | 二阶，f0 处衰减到 g | — |
| **Peak（峰值）** | 二阶，f0 处提升到 g | — |

Notch / Peak 是独立的固定二阶族，参数是 `(f0, xi, g, fs)`：`g < 1` 用 `notch_init`
（线性衰减），`g > 1` 用 `peak_init`（线性增益），两者满足 `H_peak(ξ,g) ≡ 1/H_notch(ξ,1/g)`。

## 快速开始

### 构建与测试

```bash
cmake -B build -DBUILD_TESTS=ON && cmake --build build
cd build && ctest --output-on-failure
```

预期 **13 项测试全部通过**（5 项 C 单元/回归 + 4 项 CSV 生成 + 3 项黄金参考对比 +
1 项 zpk 增益精度验证）。4 项 Python 对比在缺 numpy/scipy 时以退出码 77 自动 SKIP；
配置期找不到 Python3 时后面 8 项不注册、ctest 只有 5 项，**不要当成 13/13**。
另有 ASan+UBSan 变体：`cmake -B build-asan -DBUILD_SANITIZE=ON`。

### 基本用法

```c
#include "butter_filter.h"

// 二阶低通，截止 2 Hz，采样 20 Hz
butter_lp_2nd_t filter;
butter_lp_2nd_init(&filter, 2.0f, 20.0f);

// 必须检查 valid：设计失败时滤波器按直通（H(z)=1）部署，update 会原样返回
// 输入。宁可直接过，也不要发散的输出。
if (!filter.valid) {
    // 按产品策略处理：换参数重试 / 报错 / 继续跑直通
}

butter_lp_2nd_reset(&filter, 1.0f);   // 跳过起振瞬态

for (int i = 0; i < 1000; i++) {
    float y = butter_lp_2nd_update(&filter, input[i]);
}
// 不需要 destroy，离开作用域自动回收
```

其余族：

```c
butter_bp_2nd_t bp;
butter_bp_2nd_init(&bp, 2.0f, 5.0f, 40.0f);    // fc1, fc2, fs
cheby1_lp_3rd_t c1;
cheby1_lp_3rd_init(&c1, 3.0f, 20.0f, 0.5f);    // fc, fs, ripple_dB

// 工频陷波 50 Hz @ 1 kHz，深 −20 dB ⇒ g = 0.1
notch_filter_t n;
notch_init(&n, 50.0f, 0.05f, 0.1f, 1000.0f);   // f0, xi, g, fs
if (!n.valid) { /* 参数被闸拒（如 f0/fs 太小）→ 直通 */ }
notch_reset(&n, 0.0f);

// 同频点 +20 dB ⇒ g = 10；与陷波逐点互为倒数
peak_filter_t p;
peak_init(&p, 50.0f, 0.05f, 10.0f, 1000.0f);
if (!p.valid) { /* 参数被闸拒 → 直通 */ }
peak_reset(&p, 0.0f);
```

## API 一览

命名规则：`{族}_{类型}_{阶数序数}_t`，配套 `_init` / `_update` / `_reset` 三个函数；
阶数序数为 `1st` ~ `8th`（如 `butter_lp_2nd_t`、`cheby1_bp_5th_t`）。

| 族 | 头文件 | 结构体 | `_init` 参数（除首参指针外） |
|---|---|---|---|
| Butterworth | `butter_filter.h` | `butter_{lp,hp,bp,bs}_{1st..8th}_t` | LP/HP: `(fc, fs)`；BP/BS: `(fc1, fc2, fs)` |
| Chebyshev I | `cheby_filter.h` | `cheby1_{lp,hp,bp,bs}_{1st..8th}_t` | LP/HP: `(fc, fs, ripple_db)`；BP/BS: `(fc1, fc2, fs, ripple_db)` |
| Chebyshev II | `cheby_filter.h` | `cheby2_{lp,hp,bp,bs}_{1st..8th}_t` | LP/HP: `(fc, fs, stopband_db)`；BP/BS: `(fc1, fc2, fs, stopband_db)` |
| Notch | `notch_filter.h` | `notch_filter_t` | `(f0, xi, g, fs)`，`0 < g < 0.9999` |
| Peak | `peak_filter.h` | `peak_filter_t` | `(f0, xi, g, fs)`，`g > 1.0001` |

各族三个函数的形状一致（以二阶低通为例）：`_init(&f, ...)` 设计滤波器，失败时
`f.valid == 0` 按直通部署；`_update(&f, x)` 推一个样本返回一个样本（头文件
`static inline`，每样本零函数调用）；`_reset(&f, equilibrium)` 把状态置到直流稳态，
跳过起振瞬态。

## 与 scipy 的对比验证

C 代码生成 CSV → scipy 做黄金参考 → 对比稳态精度（跳过瞬态取最后 10%），已注册进
ctest，装好 numpy/scipy 即自动运行，实测稳态误差 **1e-6 ~ 1.5e-5**。C 的
`reset(equilibrium)` 与 scipy 的零初态起点不同，前面几百个采样不一致属正常。

## MCU 部署

- **零 `malloc`**：不需要 `<stdlib.h>`，关闭堆管理器也能运行
- **零 `double`**：全部 `float`
- **不可重入**：`_update`/`_reset` 直接读写滤波器结构体内的状态，不做任何同步。
  同一个滤波器对象被 ISR 与主循环共享时，必须自行关中断、双缓冲，或让 ISR 只置
  标志、由主循环统一推状态
- **全库需要 libm**：设计管线经 `prewarp()` → `tanf()` 预畸，`filter_utils.o` 另需
  `sqrtf`/`fabsf`/`fmaxf`/`memcpy`/`memset`；biquad 的 `init`/`reset` 需 `fabsf`
  （`biquad_filter_init` 另需 `isfinite`）。库以 `PUBLIC` 声明 `m`，使用者无需自行
  链接。Butterworth 的 ROM 极点表只免掉原型计算的 `cosf`/`sinf`，不减免 libm 链接；
  Chebyshev init 另需 `powf`/`logf`/`sinhf`/`coshf`，`powf` 会连带拉入 errno 版
  数学内核，不含 reent 支持的裸机工程需自行提供这些符号
- **`_update` 路径不需要 libm**：`biquad_filter_update` 与全族 `_update` 均为纯
  乘加代数形式且是头文件 `static inline`，节数是编译期字面量——每样本零函数调用、
  零运行时节数装载（需 −O1 及以上），代价只落在 `init`（及 `reset` 的 `fabsf`）上
- **flash 实测**（arm-none-eabi-gcc 16.2 / Cortex-M4，`.text + .rodata`）：库以
  `-ffunction-sections -fdata-sections` 编译，**必须配合链接器 `--gc-sections`
  才生效**。六个目标文件 `-Os` 合计 14686 B、`-O2` 合计 18204 B；只用
  `butter_lp_2nd` 的最小程序（`init` 一次 + `update` 八次，`--gc-sections -e main`）
  链接后 text + rodata 实测 11308 B（链 newlib libm）——设计管线无法裁除，因为
  `design_filter` 的运行时 `switch(type)` 引用全部四种频率变换。缺这些标志时链接器
  按目标文件粒度工作，会把 `butter_filter.o` 整体拉入（`-Os` 2590 B / `-O2`
  4496 B，含全部 32 个 init）

## 参数校验与 fail-closed 语义

- 原型阶数 1~8；截止频率 `0 < fc < fs/2`；BP/BS 需 `fc1 < fc2` 且 `fc2 < fs/2`；
  Chebyshev 的 `ripple_db` > 0
- Notch：`0 < f0 < fs/2`、`xi > 0`、`0 < g < 0.9999`；Peak 同前但 `g > 1.0001`
  （`g` 贴近 1 时峰/谷浅于 0.001 dB，按直通部署）。两者另有 f32 数值包络闸，
  可用区间随 `f0/fs` 与 `g` 收窄
- 任一校验失败 → `valid = 0`，update 直通返回输入。**调用方必须检查 `valid`**
  ——直通是"宁可不过滤也不要错误输出"的兜底，不是静默成功的保证

可用频率区间是三道闸（每节 Jury 条件 → 极点半径裕量 → 级联 DC/Nyquist 增益窗口）的
交集，**一律以 `valid` 为准**，不要用解析公式预测。极窄带 / 近 DC / 近 Nyquist /
超宽带 BP/BS（`fc2/fc1` ≳ 1000）被拒都是预期行为，机理与标定见 `docs/` 各知识篇目。

## 已知局限

- 原型阶数上限 8 阶（ROM 表 + 结构体枚举的工程约束）；不支持椭圆滤波器
- f32 系数在极窄带 / 极近 Nyquist 上无法忠实表达，这些设计被 fail-closed 拒绝而非强行部署
- 超宽带 BP/BS 的 DF-II 状态更新噪声把阻带衰减地板抬到 ~−12 dB 量级，为结构固有、
  无法消除
- **非有限输入**：NaN/Inf 会毒化状态、持续输出 NaN 直到 reset（热路径刻意不做防护，
  每样本分支开销）；传感器 / 不可信数据请在源头清洗

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
├── docs/                        # 设计机理、推导、被否方案（与源文件同名）
│   ├── biquad_filter.md         # 三层闸与补偿求和
│   ├── filter_utils.md          # 频率变换、zpk2sos 配对、design_filter 边界
│   ├── butter_filter.md         # 原型 ROM 表、增益折叠
│   ├── cheby_filter.md          # 原型计算、边缘增益
│   ├── notch_filter.md          # 预测式参数闸
│   └── peak_filter.md           # 对偶不变量
├── tests/                       # C 单元测试、CSV 生成器、Python 黄金参考
│   ├── CMakeLists.txt
│   ├── test_biquad.c
│   ├── test_butter.c
│   ├── test_cheby.c
│   ├── test_notch.c
│   ├── test_peak.c
│   ├── test_butter_with_py.c    # 生成 CSV 与 scipy 对比（butter/cheby/notch/peak）
│   ├── test_cheby_with_py.c
│   ├── test_notch_with_py.c
│   ├── test_peak_with_py.c
│   ├── test_butter_use_py.py    # Python 参考滤波器（butter/cheby/notch/peak）
│   ├── test_cheby_use_py.py
│   ├── test_notch_use_py.py
│   ├── test_peak_use_py.py
│   ├── compare_scipy.py         # 稳态精度对比
│   └── verify_zpk_gain.py       # float32 精度验证（独立复现 scipy 管线）
├── CMakeLists.txt
├── AGENTS.md                    # agent 工作手册（不变量清单 + 知识地图）
└── README.md
```
