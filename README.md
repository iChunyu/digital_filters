# digital_filters — 数字滤波器库（C 语言）

> **声明：本仓库所有内容（代码、文档、测试、构建脚本）均由 DeepSeek V4 Pro 生成。**

## 项目简介

这是一个用 C 语言实现的 **IIR 数字滤波器库**，基于 **双二阶 (biquad)** 滤波器及其级联（SOS，二阶节）构建高阶滤波器。所有滤波器均采用 **Direct Form II（规范型）** 结构，仅需 3 个状态变量，是内存效率最高的 IIR 实现方式。

### 支持的滤波器

| 族 | 通带特性 | 阻带特性 |
|---|---|---|
| **Butterworth** | 最大平坦 | 单调衰减 |
| **Chebyshev Type I** | 等波纹（指定纹波 dB） | 单调衰减 |
| **Chebyshev Type II** | 单调 | 等波纹（指定最小衰减 dB） |

每种族支持四种类型：**低通 (LP)**、**高通 (HP)**、**带通 (BP)**、**带阻 (BS)**，最大支持 **8 阶**。

### 两套 API

项目提供两套 API 共存，适用于不同场景：

| | 动态 API | 静态 API |
|---|---|---|
| **结构体** | `butter_t`, `cheby1_t`, `cheby2_t` | `butter_lp_2nd_t`, `cheby1_hp_3rd_t` 等 |
| **内存** | `malloc` / `free` | 内嵌数组，零堆分配 |
| **阶数** | 运行时指定 | 编译时固定 |
| **释放** | 必须调用 `_destroy` | 离开作用域自动回收 |
| **适用** | 桌面端、阶数可变 | MCU、裸机、堆不可用 |

静态 API 的关键设计：
- **按阶 update/reset**：每个阶数/类型组合有各自的 `_update` / `_reset` 函数（如 `butter_lp_2nd_update`），调用时无需强转。内部通过共享前缀 `STATIC_BUTTER_FIELDS` 复用实现。Chebyshev 同理。
- **预计算极点表**：Butterworth 原型极点仅与阶数有关，以 `static const` 存入 ROM（~512 字节），init 时无需调用三角函数。
- **栈上设计流水线**：init 时整个模拟原型→频率变换→双线性变换→零极点配对过程在栈上完成，最大约 448 字节临时空间（8 阶 BP/BS）。

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

预期输出：6 项测试全部通过。

### 基本用法（静态 API — 推荐 MCU 使用）

```c
#include "static_butter_filter.h"

// 二阶低通 Butterworth，截止频率 2 Hz，采样率 20 Hz
butter_lp_2nd_t filter;
butter_lp_2nd_init(&filter, 2.0f, 20.0f);

// 设置为稳态（DC 输入为 1.0）
butter_lp_2nd_reset(&filter, 1.0f);

// 处理采样
for (int i = 0; i < 1000; i++) {
    float output = butter_lp_2nd_update(&filter, input[i]);
    // 使用 output ...
}

// 无需 destroy，filter 离开作用域即释放
```

带通示例：

```c
butter_bp_2nd_t bp;
butter_bp_2nd_init(&bp, 2.0f, 5.0f, 40.0f);  // fc1, fc2, fs

butter_bp_2nd_reset(&bp, 0.0f);
float y = butter_bp_2nd_update(&bp, x);
```

Chebyshev Type I 示例（引入纹波参数）：

```c
cheby1_lp_3rd_t c1;
cheby1_lp_3rd_init(&c1, 3.0f, 20.0f, 0.5f);  // fc, fs, ripple_db

cheby1_lp_3rd_reset(&c1, 1.0f);
float y = cheby1_lp_3rd_update(&c1, x);
```

Chebyshev Type II 示例（阻带衰减）：

```c
cheby2_hp_2nd_t c2;
cheby2_hp_2nd_init(&c2, 5.0f, 40.0f, 40.0f);  // fc, fs, stopband_dB

cheby2_hp_2nd_reset(&c2, 0.0f);
float y = cheby2_hp_2nd_update(&c2, x);
```

### 基本用法（动态 API — 桌面端使用）

```c
#include "butter_filter.h"

butter_t b;
butter_init(&b, FILTER_LOWPASS, 2, 2.0f, 0.0f, 20.0f);

if (b.valid) {
    butter_reset(&b, 1.0f);
    float y = butter_update(&b, 0.5f);
}

butter_destroy(&b);  // 必须调用！
```

### 可用的静态结构体命名规则

```
{族}_{类型}_{阶数序数}_t

族:   butter, cheby1, cheby2
类型: lp, hp, bp, bs
阶数: 1st, 2nd, 3rd, 4th, 5th, 6th, 7th, 8th
```

示例：`butter_lp_2nd_t`, `cheby1_bp_5th_t`, `cheby2_bs_3rd_t`

对应的 init / update / reset 函数命名规则一致：
- `{族}_{类型}_{阶数序数}_init(f, ...)`
- `{族}_{类型}_{阶数序数}_update(f, input)`
- `{族}_{类型}_{阶数序数}_reset(f, equilibrium)`

## 阶数与节数

不同配置下 biquad 节的数量：

| 类型 | 阶数 N | 节数 |
|---|---|---|
| LP, HP | N | ceil(N/2) |
| BP, BS | N | N（频率变换使阶数翻倍）|

以 8 阶为例：LP 需 4 节，BP 需 8 节。

## 设计流水线

所有滤波器共享同一条设计流水线：

```
1. 模拟原型极点（Butterworth：查表；Chebyshev：运行时计算）
2. 模拟频率变换（LP: 缩放; HP: 倒数; BP/BS: 阶数翻倍变换）
3. 双线性变换 (s → z)
4. 零点补齐（s=∞ 的零点映射到 z=-1；BS 除外）
5. 零极点配对 → biquad 系数 (zpk2sos)
6. 部署到 biquad 节
```

### 预畸变

所有截止频率在进入模拟域设计之前均经过双线性预畸变，确保数字域截止位置准确：

```
f_analog = fs/π · tan(π · f_digital / fs)
```

## 注意事项

### MCU 使用

- 静态 API 文件（`static_butter_filter.c` / `static_cheby_filter.c`）**不包含任何 `malloc`/`free` 调用**，无需链接堆管理器。
- 最小栈需求：约 448 字节（8 阶 BP/BS init 时）。更低阶的滤波器占用更少。
- Butterworth init **不调用 `cosf`/`sinf`**（预计算查表），适合 `libm` 不可用或代价高的场景。
- Chebyshev init 需 `asinhf`/`sinhf`/`coshf`（仅 init 时调用一次，非逐采样）。
- `#include <stdlib.h>` 不出现在任何静态 API 文件中。

### 参数校验

- 阶数必须 ≥ 1 且 ≤ 8（静态 API 编译时保证；动态 API 运行时检查）。
- 截止频率必须满足 `0 < fc < fs/2`。
- BP/BS 的 `fc1` 必须 `< fc2` 且 `fc2 < fs/2`。
- Chebyshev 的 `ripple_db` 必须 > 0。
- 非法参数导致 `valid = 0`，此时 update 直通返回输入值（无滤波）。

### 稳定性

- 每个 biquad 节在初始化时经过全部三个 Jury 稳定性条件检验。
- 不稳定的节被静默替换为直通 (`H(z) = 1`)。
- 由于 Butterworth 和 Chebyshev 原型本身稳定，正常参数下不应触发此降级。

### 数值精度

- 全部使用 `float`（32 位单精度）。对于极高阶窄带滤波器，双精度可能更优，但超出本库范围。
- 零极点配对算法使用"最不利极点优先"策略，最大程度降低有限精度下的舍入噪声。
- 8 阶 BP/BS 的极点密度较高，窄带场景下建议评估实际精度是否满足需求。

## 文件结构

```
digital_filters/
├── include/
│   ├── biquad_filter.h          # 单节 biquad
│   ├── sos_filter.h             # SOS 级联（动态）
│   ├── filter_utils.h           # 设计工具 + 类型定义
│   ├── butter_filter.h          # Butterworth 动态 API
│   ├── cheby_filter.h           # Chebyshev 动态 API
│   ├── static_butter_filter.h   # Butterworth 静态 API
│   └── static_cheby_filter.h    # Chebyshev 静态 API
├── src/
│   ├── biquad_filter.c
│   ├── sos_filter.c
│   ├── filter_utils.c
│   ├── butter_filter.c
│   ├── cheby_filter.c
│   ├── static_butter_filter.c
│   └── static_cheby_filter.c
├── test/
│   ├── CMakeLists.txt
│   ├── test_biquad.c
│   ├── test_sos.c
│   ├── test_butter.c
│   ├── test_cheby.c
│   ├── test_static_butter.c
│   ├── test_static_cheby.c
│   ├── test_butter_with_py.c     # 生成 CSV 与 Python 参考对比
│   └── test_cheby_with_py.c      # 生成 CSV 与 Python 参考对比
├── CMakeLists.txt
├── CLAUDE.md
└── README.md
```

## 许可

本项目仅作示例用途。
