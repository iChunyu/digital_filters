# digital_filters — 数字滤波器库（C 语言）

> 一个给 MCU 用的 IIR 滤波器库。没有 `double`，没有 `malloc`（静态 API），
> 全部 `float` 一把梭。和 scipy 对过答案了，稳态误差在 1e-6 量级 🤏

## 项目简介

基于 **双二阶 (biquad)** 滤波器及其级联（SOS，二阶节）构建高阶 IIR 滤波器。
全部采用 **Direct Form II（规范型）**——每个 biquad 只要 3 个 `float` 状态变量，
是你能写出来的最省内存的 IIR 实现。再多省一个字节就得去改物理定律了。

### 支持的滤波器

| 族 | 通带特性 | 阻带特性 | 一句话点评 |
|---|---|---|---|
| **Butterworth** | 最大平坦 | 单调衰减 | 老实人，不搞花活 |
| **Chebyshev Type I** | 等波纹（指定纹波 dB） | 单调衰减 | 通带里蹦迪，阻带装死 |
| **Chebyshev Type II** | 单调 | 等波纹（指定最小衰减 dB） | 反过来，通带佛系阻带蹦迪 |

每种支持四种类型：**低通 (LP)**、**高通 (HP)**、**带通 (BP)**、**带阻 (BS)**。

### 两套 API

| | 动态 API | 静态 API |
|---|---|---|
| **结构体** | `butter_t`, `cheby1_t`, `cheby2_t` | `butter_lp_2nd_t`, `cheby1_hp_3rd_t` 等 |
| **内存** | `malloc` / `free` | 内嵌数组，零堆分配 |
| **最大阶数（原型）** | **12 阶** | **8 阶** |
| **BP/BS 有效阶数** | 原型×2（最高 24 阶） | 原型×2（最高 16 阶） |
| **释放** | 必须调用 `_destroy` | 离开作用域自动回收 |
| **适用** | 桌面端、阶数可变 | MCU、裸机、堆不可用 |
| **增益精度** | `float` | `float` |

> **关于阶数**：API 里的 `order` 参数是模拟原型阶数。对 LP/HP，这就是最终阶数；
> 对 BP/BS，频率变换会翻倍（12 阶原型 → 24 阶有效滤波器，12 个 biquad 节）。
>
> 为什么动态 API 只到 12 阶（原型）？因为增益追踪全程 `float`，`wc^N` 在原型 N=13 时
> 可能溢出 float 上限（~3.4e38）。12 阶对所有常见采样率和截止频率都是安全的。
> 真需要更高阶？级联两个 12 阶的，自己拼。

静态 API 亮点：
- **零 `malloc`**：`<stdlib.h>` 都不需要 `#include`，堆管理器关掉照样跑
- **预计算极点表**：Butterworth 原型极点存 ROM（~512 字节），init 不用调 `cosf`/`sinf`
- **栈上设计流水线**：init 全程在栈上完成，不用堆

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

预期输出：**7 项测试全部通过**。

### 基本用法（静态 API — 推荐 MCU 使用）

```c
#include "static_butter_filter.h"

// 二阶低通 Butterworth，截止频率 2 Hz，采样率 20 Hz
butter_lp_2nd_t filter;
butter_lp_2nd_init(&filter, 2.0f, 20.0f);

// 设置为稳态（跳过烦人的起振瞬态）
butter_lp_2nd_reset(&filter, 1.0f);

// 处理采样
for (int i = 0; i < 1000; i++) {
    float output = butter_lp_2nd_update(&filter, input[i]);
    // 干点啥...
}

// 不用 destroy，离开作用域自动释放。就是这么省心 😌
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

### 基本用法（动态 API — 桌面端）

```c
#include "butter_filter.h"

butter_t b;
butter_init(&b, FILTER_LOWPASS, 12, 50.0f, 0.0f, 400.0f);  // 12 阶 LP

if (b.valid) {
    butter_reset(&b, 1.0f);
    float y = butter_update(&b, 0.5f);
}

butter_destroy(&b);  // 别忘了！malloc 出来的得还回去
```

### 静态结构体命名规则

```
{族}_{类型}_{阶数序数}_t

族:   butter, cheby1, cheby2
类型: lp, hp, bp, bs
阶数: 1st ~ 8th
```

示例：`butter_lp_2nd_t`, `cheby1_bp_5th_t`, `cheby2_bs_3rd_t`

## 阶数与节数

| 类型 | 阶数 N | 节数 |
|---|---|---|
| LP, HP | N | ceil(N/2) |
| BP, BS | N | N（频率变换使阶数翻倍）|

动态 API 12 阶 BP → 12 节 biquad；静态 API 8 阶 BP → 8 节。

## 设计流水线

```
1. 模拟原型极点（Butterworth：查表；Chebyshev：运行时算）
2. 模拟频率变换（LP: 缩放; HP: 倒数; BP/BS: 阶数翻倍）
3. 双线性变换 (s → z)
4. 零点补齐（s=∞ 的零点映射到 z=-1；BS 除外）
5. 零极点配对 → biquad 系数 (zpk2sos)
6. 部署到 biquad 节
```

预畸变公式：`f_analog = fs/π · tan(π · f_digital / fs)`

## 与 scipy 的对比验证

仓库自带 Python 对比脚本。C 代码生成 CSV → scipy 做黄金参考 → 对比稳态精度：

| 滤波器 | LP | HP | BP | BS |
|---|---|---|---|---|
| Butterworth | 1.2e-6 | 2.1e-6 | 3.5e-6 | 1.1e-6 |
| Chebyshev I | 1.2e-6 | 3.6e-6 | 1.5e-5 | 3.7e-5 |
| Chebyshev II | 1.4e-6 | 3.2e-6 | 4.5e-6 | 9.8e-7 |

> 注意：C 的 `reset(equilibrium)` 和 scipy 的零初态起点不同，前面几百个采样点对不上
> 是正常的（瞬态响应差异）。上表取的是最后 10% 的稳态数据。

## 注意事项

### MCU 使用

- 静态 API **不含 `malloc`/`free`**，无需堆管理器
- Butterworth init **不调 `cosf`/`sinf`**（查表），`libm` 不是必需品
- Chebyshev init 需 `asinhf`/`sinhf`/`coshf`（仅 init 一次）
- `biquad_filter_update` 是 `static inline`，编译器会直接内联到调用点——没有函数调用开销

### 参数校验

- 动态 API 阶数 1~12，静态 API 阶数 1~8
- 截止频率 `0 < fc < fs/2`
- BP/BS 需 `fc1 < fc2` 且 `fc2 < fs/2`
- Chebyshev 的 `ripple_db` > 0
- 非法参数 → `valid = 0`，update 直通返回输入

### 稳定性

- 每节 biquad 过全部三个 Jury 条件：`|a2|<1`, `1+a1+a2>0`, `1-a1+a2>0`
- 不稳定 → 静默换直通 (`H(z)=1`)，不报错（嵌入式环境你报给谁看？）

### 数值精度

- **全 `float`，零 `double`**。这是故意的——不是不能加，是不想加。
  12 阶以内 float 完全够用，加 double 只增加 ROM/RAM 负担。
- "最不利极点优先"配对策略，减少有限精度舍入噪声
- 窄带 BP/BS 高 Q 值场景建议实测评估精度

## 已知局限

- 动态 API 原型阶数上限 12（BP/BS 有效阶数 24）。需要更高就级联。
- 静态 API 8 阶上限（ROM 极点表 + 结构体枚举的工程约束）。
- 不支持椭圆滤波器（那货的 Jacobi 椭圆函数写起来太抽象了，下次一定）。

## 文件结构

```
digital_filters/
├── include/
│   ├── biquad_filter.h          # 单节 biquad（含 inline update）
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
│   ├── test_high_order.c         # 高阶测试（8/10/12 阶）
│   ├── test_butter_with_py.c     # 生成 CSV 与 scipy 对比
│   ├── test_cheby_with_py.c
│   ├── test_butter_use_py.py     # Python 参考滤波器（含绘图）
│   ├── test_cheby_use_py.py
│   └── compare_scipy.py          # 稳态精度对比（跳过瞬态）
├── CMakeLists.txt
├── CLAUDE.md
└── README.md
```
