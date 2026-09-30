# Butterworth 族（butter_filter）

**源文件**：`include/butter_filter.h`、`src/butter_filter.c`
**测试**：`tests/test_butter.c`（标定表、回归档位与实测统计只在此维护，本文只给结论与指针）

最大平坦幅频响应的四类型滤波器（LP / HP / BP / BS）× 八个原型阶数，全部走
`docs/filter_utils.md` 里的共享管线 `design_filter()`。族特有的知识只有两件：
ROM 原型表与 X-macro 族生成。

## 原型极点与 ROM 表

$$
\theta_k = \frac{\pi(2k + N + 1)}{2N}, \qquad
p_k = \cos\theta_k - j\,\sin\theta_k, \qquad k = 0, 1, \dots, N - 1
$$

极点全在左半平面（$\mathrm{Re}\,p_k < 0$），原型按 $H(0) = 1$ 归一，所以本族传
给管线的增益常数恒为 $k = 1$、有限零点恒为 0 个（`proto_zeros == NULL` 在
`nz == 0` 时合法）。

八个阶数的极点预计算成 `static const` 表 `butter_proto[order − 1][pole_index]`：
行宽固定 8 列、未用到的尾部零填充，`complex_t` 是两个 `float`，于是
$8 \times 8 \times 8\ \mathrm{B} \approx 512\ \mathrm{B}$——放 flash，不占 RAM。

### 为什么 Butterworth 查 ROM 表，Chebyshev 却运行时算

Butterworth 的原型极点**只依赖阶数** $N$：可能的原型只有 8 个，是有限集，直接列进
`static const` 表最省——init 里省掉 `cosf` / `sinf`，也不需要在栈上算。

Chebyshev 的原型还依赖 `ripple_db`（连续量），无法枚举成有限表，只能每次 init
现算（见 `docs/cheby_filter.md`）。这是两族在原型层的唯一差异，不是随意选择。

注意 ROM 表只免掉原型计算的三角函数，**不减免 libm 链接**：`prewarp()` 的 `tanf`、
双线性与配对里的 `sqrtf` / `fabsf` 仍在 init 路径上；只有 `_update` 路径才是零
libm。

## X-macro 族生成

32 个结构体（4 类型 × 8 阶）与它们的 init / update / reset 全部由宏展开，没有手抄：

- `BUTTER_FIELDS` 是四种类型结构体的公共前缀（`valid` / `type` / `order` /
  `num_sections` / `fc1` / `fc2` / `fs`），后面接内嵌的
  `biquad_filter_t sections[ns]`。`ns` 是编译期字面量，结构体离开作用域即回收
  （无 `malloc`）。
- 两张阶数表：`FOR_EACH_BUTTER_LP_ORDER` 与 `FOR_EACH_BUTTER_BP_ORDER`，各列
  `X(order, sections, ordinal_label)`。**两张表不能合并**——LP/HP 与 BP/BS 的节数
  规则不同（见下节），共用一张表会让结构体的 `sections[]` 容量错配。
- 结构体 typedef、init 声明与定义、`static inline` 的 update / reset 由同一张表
  展开，所以结构与容量不会互相漂移。

命名：`butter_{lp,hp,bp,bs}_{1st..8th}_{init,update,reset}`，类型
`butter_{lp,hp,bp,bs}_{1st..8th}_t`。

## 节数规则

| 类型 | 原型阶数 $N$ | 有效阶数 | biquad 节数 |
|---|---|---|---|
| LP, HP | $N$ | $N$ | $\lceil N/2 \rceil$ |
| BP, BS | $N$ | $2N$ | $N$ |

LP/HP 一节装一对共轭极点，奇数阶剩一个实极点，故 $\lceil N/2 \rceil$。BP/BS 的
频率变换把每个原型极点一分为二，有效阶数翻倍，节数因此等于 $N$（8 阶 BP → 16 阶
有效滤波器 → 8 节）。宏表里的 ns 是容量，运行时的 `num_sections` 由
`design_filter()` 返回并写入，失败时为 0。

## 期望 DC/Nyquist 增益：只用 ±0.1 窗口

四种类型在 $z = \pm 1$ 处的零点/极点位置是结构确定的，DC 与 Nyquist 增益因此是
**精确**的 0 或 1：

| 类型 | DC 期望 | Nyquist 期望 |
|---|---|---|
| LP | 1 | 0 |
| HP | 0 | 1 |
| BP | 0 | 0 |
| BS | 1 | 1 |

`check_cascade_gains` 的窗口按期望值分工：期望 ∈ {0, 1}（精确结构增益）用 ±0.1，
只有非 0/1 的纹波边缘才用 ±0.25。Butterworth 四种类型没有任何一个边缘期望落在
纹波档，**本族只用 ±0.1、从不使用 ±0.25**——若哪天 butter 走了纹波档，说明期望值
传错了。

这道闸拦的是"稳定但配错对 / 增益缩放错"的静默错误（缺陷签名与实测见
`tests/test_cheby.c` 的「回归：过去会静默部署错误滤波器的宽带 cheby2 配置」段）。
机制与窗口分工见 `docs/filter_utils.md`，窗口两侧边界档位见 `tests/test_butter.c`
的「回归：design_filter 输入边界闸（导出的管线）」段。

## 设计流程与 fail-closed 契约

每个类型的 `static` 辅助函数（`butter_{lp,hp,bp,bs}_init`）三步：

1. 入口闸：`order ∈ 1..8`、$0 < f_c < f_s/2$（BP/BS 另需 $f_{c1} < f_{c2} < f_s/2$）。
   一律写成 `!(fc > 0.0f)` 而不是 `fc <= 0.0f`——后者与 NaN 比较恒假，会把 NaN
   放行。
2. 查 ROM 表 → `design_filter(..., k = 1.0f, butter_proto[order − 1], order,
   NULL, 0)`。
3. `check_cascade_gains` 按上表四种类型的期望校验 DC/Nyquist，不通过即返回 0。

X-macro 宏体先写 `valid = 0` **且** `num_sections = 0`，再按辅助函数的返回值决定
是否置 `valid = 1` 与写入节数。`num_sections` 必须显式清零：调用方可能漏检
`valid`，垃圾节数会让后续按节数遍历的代码读到越界节（`_update` 另有 `!valid` 判据
先返回输入，节数按编译期字面量传入、不读 `num_sections`）；清成 0 后 0 节级联天然
直通。失败后 `sections[]` 的内容不保证，不得直接使用。

`_update` / `_reset` 是头文件 `static inline`（节数编译期字面量、每样本零函数
调用），且**不做 NaN/Inf 防护、不可重入**；`!valid` 时 update 原样返回输入、reset
为空操作。这两条语义与理由见 `docs/biquad_filter.md`。

init 成功后内部状态为零而不是稳态：带直流的信号从零相位开始时要先调
`_reset(equilibrium)`。

## 已知边界：被拒是预期行为

近 DC / 近 Nyquist / 极窄带 / 超宽带设计返回 `valid = 0`（直通）是**预期**行为，
不是 bug：

- 主因是 **f32 系数量化**让 DC / 边缘增益塌缩（有时是极点半径裕量），不是算法
  错误；
- 边界**非单调**：相邻的一小段频率里 `valid` 会反复翻转，且随阶数 / 类型不同。
  不要用解析公式预测可用区间，**一律以 `valid` 为准**；
- 超宽带 BP/BS（$f_{c2}/f_{c1} \gtrsim 1000$）的阻带地板 ~−12 dB 是 DF-II 结构
  固有，换实现也压不掉。

已标定的两侧档位见 `tests/test_butter.c`：「回归：可表示的窄带 LP 接受，超窄带
（极点落在距 z=1 的 5e-5 裕量之内）fail-closed 拒绝」段，与「回归：近 Nyquist
带边——在极点裕量之内设计正常部署；超出则确定性地拒绝（没有 NaN 垃圾，相邻规格
之间也没有接受悬崖）」段。

## 测试指针

标定表、回归档位与实测统计只维护在 `tests/test_butter.c`（本文只给结论与指针）：

- 「参数非法 → valid = 0，直通」、「回归：init 被拒后 num_sections == 0，而不是
  垃圾值」、「入口闸的 NaN 行为：必须 fail-closed」——fail-closed 契约的三面；
- 「回归：超宽带 BP8（增益链 f32 溢出）」、「回归：近 Nyquist LP8（此前会产生
  NaN）」、「回归：含近实极点对的宽带 BP」、「c_sqrt 的缩放计算（防近 Nyquist
  BP/BS 的平方幅值溢出）」——增益折叠、双线性与求根的 f32 极限；
- 「回归：design_filter 输入边界闸（导出的管线）」——管线入口闸与增益窗口两侧
  边界；
- 「Sweep: butter 全阶矩阵（LP/HP/BP/BS × 全 8 阶）」——X-macro 全阶覆盖与
  `max|H|` 约束；
- 「结构体尺寸」——内嵌 `sections[]` 随阶数增长。

## 为什么不这样写（被否方案）

- **像 Chebyshev 那样运行时算原型**：Butterworth 原型只依赖阶数、是有限集，查表
  省掉 init 里的 `cosf` / `sinf`，代价只是一张 ~512 B 的 const 表。
- **用一张 X-macro 阶数表覆盖 LP/HP 与 BP/BS**：两种节数规则（$\lceil N/2 \rceil$
  vs $N$）会让 `sections[]` 容量错配。
- **辅助函数只写 `valid`、不返回节数**：`num_sections` 会留下垃圾值，而调用方
  可能漏检 `valid`。保持 `uint8_t` 返回值 + 宏体先双双清零。
