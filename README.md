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

预期输出：**3 项测试全部通过**。

### 基本用法

```c
#include "butter_filter.h"

// 二阶低通，截止 2 Hz，采样 20 Hz
butter_lp_2nd_t filter;
butter_lp_2nd_init(&filter, 2.0f, 20.0f);

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

C 代码生成 CSV → scipy 做黄金参考 → 对比稳态精度（跳过瞬态取最后 10%）：

| 滤波器 | LP | HP | BP | BS |
|---|---|---|---|---|
| Butterworth | 0.9e-6 | 2.1e-6 | 3.2e-6 | 1.1e-6 |
| Chebyshev I | 1.2e-6 | 3.6e-6 | 1.5e-5 | 3.7e-5 |
| Chebyshev II | 1.4e-6 | 3.2e-6 | 4.5e-6 | 1.0e-6 |

> C 的 `reset(equilibrium)` 和 scipy 的零初态起点不同，前面几百个采样对不上是
> 正常的（瞬态响应差异）。上表取的是稳态数据。

## 注意事项

### MCU 使用

- **零 `malloc`**：`<stdlib.h>` 不需要，堆管理器关掉照样跑
- Butterworth init **不调 `cosf`/`sinf`**（ROM 查表），`libm` 不是必需品
- Chebyshev init 需 `asinhf`/`sinhf`/`coshf`（仅 init 一次，非逐采样）
- `biquad_filter_update` 是 `static inline`，编译器直接内联——没有函数调用开销
- 全部 `float`，零 `double`

### 参数校验

- 原型阶数 1~8
- 截止频率 `0 < fc < fs/2`
- BP/BS 需 `fc1 < fc2` 且 `fc2 < fs/2`
- Chebyshev 的 `ripple_db` > 0
- 非法参数 → `valid = 0`，update 直通返回输入

### 稳定性

- 每节 biquad 过全部三个 Jury 条件
- 不稳定 → 静默换直通 (`H(z)=1`)，不报错

### 数值精度

- 全 `float`。8 阶以内完全够用，加 double 只增加 ROM/RAM 负担
- "最不利极点优先"配对策略，减少有限精度舍入噪声
- 窄带高 Q 场景建议实测评估

## 已知局限

- 原型阶数 8 阶上限（ROM 表 + 结构体枚举的工程约束）
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
