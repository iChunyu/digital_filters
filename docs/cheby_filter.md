# Chebyshev I / II 族（cheby_filter）

**源文件**：`include/cheby_filter.h`、`src/cheby_filter.c`
**测试**：`tests/test_cheby.c`（标定表、回归档位与实测统计只在此维护，本文只给结论与指针）

等纹波的两族：Chebyshev I 是通带等纹波 + 单调阻带，Chebyshev II 是单调通带 +
阻带等纹波。两族 × 四类型（LP / HP / BP / BS）× 八个阶数 = 64 个结构体，全部走
`docs/filter_utils.md` 里的共享管线 `design_filter()`。与 Butterworth 的结构差异
只有两点：原型要运行时算（依赖 `ripple_db`），以及 Chebyshev II 带有限零点。

## Chebyshev I 原型（通带等纹波）

$$
\theta_k = \frac{\pi(2k + N + 1)}{2N}, \qquad k = 0, 1, \dots, N - 1
$$

$$
\varepsilon = \sqrt{10^{r_p/10} - 1}, \qquad
\mu = \frac{\mathrm{asinh}(1/\varepsilon)}{N}, \qquad
p_k = \sinh\mu \cdot \cos\theta_k + j\,\cosh\mu \cdot \sin\theta_k
$$

$r_p$ 是 `ripple_db`（通带纹波 dB）。角分布与 Butterworth 相同，但极点被推到
椭圆上：实部乘 $\sinh\mu$、虚部乘 $\cosh\mu$。纹波越深（$r_p$ 越大）$\varepsilon$
越大、$\mu$ 越小，椭圆越扁——极点越靠虚轴。`ripple_db <= 0` 被入口闸拒绝（含
NaN），见下文「设计流程与 fail-closed 契约」。

## Chebyshev II 原型（阻带等纹波）

$$
\varepsilon = \frac{1}{\sqrt{10^{r_s/10} - 1}}, \qquad
\mu = \frac{\mathrm{asinh}(1/\varepsilon)}{N}, \qquad
z_k = \frac{j}{\sin\theta_k}
$$

$r_s$ 是 `ripple_db`（阻带最小衰减 dB）。极点取 Chebyshev I 极点的**倒数**：记
Chebyshev I 的极点为 $p_k^{(\mathrm{I})}$，则

$$
p_k^{(\mathrm{II})} = \frac{\overline{p_k^{(\mathrm{I})}}}{\left|p_k^{(\mathrm{I})}\right|^2}
= \frac{1}{p_k^{(\mathrm{I})}}
$$

这正是"频率变量取倒数"（$\omega \mapsto 1/\omega$）的后果：极点 $p \mapsto 1/p$，
同时把原型在 $s = \infty$ 的零点映到有限点 $z_k = j/\sin\theta_k$。写共轭除以模方
（而不是直接 `1/p`）是为了在实数域完成这一步：`cheby2_proto` 先算分母

$$
\mathrm{den} = \sinh^2\mu\cos^2\theta + \cosh^2\mu\sin^2\theta
            = \left|p_k^{(\mathrm{I})}\right|^2
$$

再取实部 $\sinh\mu\cos\theta/\mathrm{den}$、虚部
$-\cosh\mu\sin\theta/\mathrm{den}$，全程只有实数运算，不引入复数除法。

零点 $z_k = j/\sin\theta_k$ 落在虚轴上：$\theta_k$ 与 $2\pi - \theta_k$（即 $k$ 与
$N - 1 - k$）的 $\sin$ 异号，给出互为共轭的一对；$\theta_k$ 恰为 $\pi$ 的那个零点
（奇数阶有一个）在 $s = \infty$，不输出——见下节。

## 零跳过阈值 1e-4

$z_k = j/\sin\theta_k$ 在 $\theta = \pi$ 处发散。`cheby2_proto` 用
`fabsf(sinf(theta)) > 1e-4f` 筛掉这些角，返回值 `nz` 就是有限零点个数：偶数阶
$N$、奇数阶 $N - 1$，直接作为 `design_filter` 的 `proto_zeros` / `nz` 传入。

阈值只需吞掉 $\theta = \pi$ 处 `sinf` 的**求值误差**，不需要留"接近 π 的合法角"
的余量：$\theta_k = \pi$ 这个角本身是精确的（$N$ 为奇数时 $k = (N-1)/2$ 取到，
算式退化成 $1 \cdot \pi$），`sinf` 的残留只来自 `(float)M_PI` 与真 π 之差
（$-8.74 \times 10^{-8}$），取倒数就是 $j \cdot 10^7$ 量级的**幻影零点**；而偏离
$\pi$ 最近的合法角相距 $\pi/(2N)$（$N$ 为偶数时取到），于是对 $N \le 8$

$$
|\sin\theta| = \sin\frac{\pi}{2N} \ge \sin\frac{\pi}{16} \approx 0.195
$$

与阈值 $10^{-4}$ 相差约 **1950 倍**。这里用的是 $N = 8$ 的 $\sin(\pi/16)$，
**不是** $\sin(\pi/8) \approx 0.38$——后者对应 $N = 4$，不在最小处取值，拿它估余量
会误以为余量比实际大。

危险的不是发散，是**静默**：几个 $j \cdot 10^7$ 的幻影零点会把响应悄悄重塑成另一个
滤波器（LP/HP 下双线性会把它们映到 $z \approx -1$ 附近），而 `valid = 1`、不报错也
不发散。能抓它的是响应回归——`tests/test_cheby.c` 的「扫掠：cheby2 BP/BS 全矩阵
必须部署出响应正常的结果」段与「Sweep: cheby2 LP/HP 全阶矩阵」段（DC/Nyquist
窗口）。

## 边缘增益期望与 ±0.25 窗口

边缘（纹波）增益随**奇偶阶**不同，所以 `check_cascade_gains` 的期望值要按阶数算：

| 族 | 类型 | DC 期望 | Nyquist 期望 |
|---|---|---|---|
| cheby1 | LP | 奇 $1$ / 偶 $10^{-r_p/20}$ | 0 |
| cheby1 | HP | 0 | 奇 $1$ / 偶 $10^{-r_p/20}$ |
| cheby1 | BP | 0 | 0 |
| cheby1 | BS | 奇 $1$ / 偶 $10^{-r_p/20}$ | 同 DC |
| cheby2 | LP | 1 | 奇 $0$ / 偶 $10^{-r_s/20}$ |
| cheby2 | HP | 奇 $0$ / 偶 $10^{-r_s/20}$ | 1 |
| cheby2 | BP | 奇 $0$ / 偶 $10^{-r_s/20}$ | 同 DC |
| cheby2 | BS | 1 | 1 |

$0$ 与 $1$ 都是**精确**的结构增益（由 $z = \pm 1$ 处的零点/极点位置决定），
其余条目是纹波边缘 $10^{-r_p/20}$ / $10^{-r_s/20}$——表里没有第三种可能。

为什么会有 $10^{-r_p/20}$ 这类期望值：偶数阶 Chebyshev I 原型在 $\omega = 0$ 处
正好落在纹波极值 $|H| = 1/\sqrt{1+\varepsilon^2} = 10^{-r_p/20}$，奇数阶在
$\omega = 0$ 处恰为 1；Chebyshev II 对称地是 $\omega = \infty$ 处的
$10^{-r_s/20}$（偶）与 0（奇）。这是**原型**的归一化事实，传到数字域的哪一端
（DC 还是 Nyquist）由各类型在 $z = \pm 1$ 处的结构零点决定。`design_filter` 的
增益 $k$ 也按同一套期望归一（见下文「设计流程与 fail-closed 契约」），所以本表既是
校验的期望，也是归一化口径。

**这就是纹波档 ±0.25 的唯一来由**：`check_cascade_gains` 的窗口按期望值分工，期望
∈ {0, 1}（精确结构增益）用 ±0.1，其余用 ±0.25。纹波边缘落在带边响应最陡处，f32
量化误差在那里最大，故放宽到 ±0.25。放宽的代价是缺陷侧分离窄（≈1.8×），它主要靠
"期望值本身可能很小"（$r_s = 40$ dB 时期望只有 0.01）而非窗口宽度设防——可接受偏差
与两侧分离倍数见 `tests/test_cheby.c` 的「Sweep: cheby1 HP/BP/BS 全阶矩阵」段。
Butterworth 四种类型的边缘期望全是精确 0/1，因此从不使用纹波档（见
`docs/butter_filter.md`）。

## 为什么 μ 用 `logf` + `sqrtf` 而不是 `asinhf`

两个原型函数都用恒等式

$$
\mathrm{asinh}(x) = \ln\left(x + \sqrt{x^2 + 1}\right), \qquad x = 1/\varepsilon
$$

即 `logf(inv_eps + sqrtf(inv_eps * inv_eps + 1.0f)) / n`，而不是直接调 `asinhf`。
依据在 git 历史（提交 `0f73e3c`「MCU 兼容性优化：asinhf 替换」）：`asinhf` 在
newlib-nano 上可能没有符号，会链接失败；改写后只用到 `logf` / `sqrtf`，二者本库
其它地方已经在用（`prewarp` / `powf` / 双线性都要求 libm），不新增链接依赖。

那条提交信息只记录了"避免 newlib-nano 链接缺失"这一条依据，**没有精度对比的
记录**——不要把它当成有实测依据的精度取舍，也不要替它补理由。

`cheby1_proto` 与 `cheby2_proto` 各自现算这段（`inv_eps` / `mu` / `sinh_mu` /
`cosh_mu` 四行完全相同），没有抽公共函数：两者其余部分（极点倒数、零点筛选）毫无
共同代码，抽出来只会多一层调用。

## 设计流程与 fail-closed 契约

每个类型的 `static` 辅助函数（`cheby{1,2}_{lp,hp,bp,bs}_init`）五步：

1. 入口闸：`order ∈ 1..8`、$0 < f_c < f_s/2$（BP/BS 另需 $f_{c1} < f_{c2} < f_s/2$）、
   `ripple_db > 0`；一律写 `!(x > 0)`（含 `ripple_db`），NaN 一样 fail-closed。
2. 算 $\varepsilon$：cheby1 用 $\sqrt{10^{r_p/10} - 1}$，cheby2 用
   $1/\sqrt{10^{r_s/10} - 1}$。
3. 算原型，写在栈上的 `complex_t poles[8]` / `zeros[8]` 里。
4. 算增益常数 $k = 1/\mathrm{Re}\frac{\prod(-z_i)}{\prod(-p_i)}$
   （`zpk_hp_bs_gain`；`nz == 0` 时分子积取 1）；cheby1 的**偶数阶**再除以
   $\sqrt{1 + \varepsilon^2}$。
5. `design_filter(...)` → `check_cascade_gains(...)` 按上表的奇偶期望校验。

$k$ 为什么要那一步除法：`zpk_hp_bs_gain` 算的是原型零极点在 $s = 0$ 处的比值，取
倒数就是把 $H(0)$ 归一化到 1。cheby1 偶数阶再除 $\sqrt{1 + \varepsilon^2}$：
偶数阶原型在 $\omega = 0$ 处的增益本来就是 $10^{-r_p/20}$
（即 $1/\sqrt{1+\varepsilon^2}$）；不除这一项，DC 会被强行拉回 1，与上表的期望
冲突。cheby2 在 $\omega = 0$ 处恒为 1（$T_N \to \infty$），故不需要这一项。

X-macro 的用法与 Butterworth 相同（`CHEBY_FIELDS` 多一个 `ripple_db` 字段；命名
`cheby{1,2}_{lp,hp,bp,bs}_{1st..8th}_{init,update,reset}`；LP/HP 节数
$\lceil N/2 \rceil$、BP/BS 为 $N$；两族 × 4 类型 × 8 阶共 64 个结构体），机制见
`docs/butter_filter.md`。fail-closed 契约也相同：宏体先 `valid = 0` 且
`num_sections = 0`，之后按辅助函数返回值决定是否置 `valid = 1`；失败后
`sections[]` 内容不保证，`_update` 直通、`_reset` 空操作，**调用方必须检查
`valid`**。init 成功后内部状态为零而不是稳态，带直流信号要先
`_reset(equilibrium)`。热路径的 inline / 无 NaN 防护 / 不可重入语义见
`docs/biquad_filter.md`。

## 已知边界

近 DC / 近 Nyquist / 极窄带 / 超宽带返回 `valid = 0`（直通）是预期行为，主因是
f32 系数量化（有时是极点半径裕量），边界非单调，一律以 `valid` 为准——与
Butterworth 同（见 `docs/butter_filter.md`）。Chebyshev 特有的被拒档位：

- **退化纹波**：`ripple_db` 很大时设计被拒——先是极点半径裕量闸，再往上则算出的
  增益常数 $k$ 变为 0，被 `design_filter` 的 `k` 有限非零闸拒绝（该测试的判据就
  写作 `k=0`）。这是预期行为，不是参数校验缺失；档位见 `tests/test_cheby.c` 的
  「退化纹波 → fail-closed（直通）」段。
- **cheby2 HP 的近 Nyquist 带边**可能被拒（极点进入 $z = -1$ 近旁的 5e-5 裕量内）。
  `tests/test_cheby.c` 的「Sweep: cheby2 LP/HP 全阶矩阵」段里有一个被白名单扣掉的
  点——标定档位 5 阶 × $f_c = 0.96$·Nyquist，预期拒绝，不是回归。
- **cheby1 BP 的近 Nyquist 带边**可能被拒（极点落在 $z = +1$ 近旁）。`tests/test_cheby.c`
  的「Sweep: cheby1 HP/BP/BS 全阶矩阵」段里有一个被白名单扣掉的点——标定档位 7 阶 ×
  $[20, 480]$ Hz（$f_s = 1000$，上沿 480 Hz 逼近 Nyquist 500 Hz）× `rp = 3` dB，
  预期拒绝，不是回归。

## 测试指针

标定表、回归档位与实测统计只维护在 `tests/test_cheby.c`（本文只给结论与指针）：

- 「参数非法 → valid = 0」、「回归：init 被拒后 num_sections == 0，而不是垃圾
  值」、「入口闸的 NaN 行为：必须 fail-closed」、「无效滤波器的直通行为」——
  fail-closed 契约；
- 「回归：cheby2 BP 五阶 [50,120]@1000，rs=0.5——实/复混合零点必须按类型配对
  （zpk2sos 配对缺陷）」、「回归：cheby1 BS 五阶 [50,120]@1000，rp=3——实极点必须
  与实极点配对（zpk2sos 配对缺陷）」、「回归：高 Q BS 极点簇——共轭配对」——
  零极点配对缺陷（配对算法本身见 `docs/filter_utils.md`）；
- 「回归：过去会静默部署错误滤波器的宽带 cheby2 配置」——增益校验闸拦下的静默
  错误；
- 「近 Nyquist LP 的抢救（增益链 f32 溢出已修）」——增益折叠的 f32 极限；
- 「退化纹波 → fail-closed（直通）」——`ripple_db` 过大时两道闸的先后；
- 「扫掠：cheby2 BP/BS 全矩阵必须部署出响应正常的结果」、「扫掠：cheby1 LP 全阶
  矩阵」、「Sweep: cheby1 HP/BP/BS 全阶矩阵」、「Sweep: cheby2 LP/HP 全阶矩阵」
  ——X-macro 全阶覆盖、节数断言、DC/Nyquist 窗口与 `max|H|` 约束；
- 「解析级联增益（不含瞬态影响）」——`cascade_dc_gain` / `cascade_nyq_gain` /
  `max_gain_over` 三个本地解析辅助。

## 为什么不这样写（被否方案）

- **直接 `1/p` 求 Chebyshev II 极点**：需要复数除法。改写为
  $\overline{p}/|p|^2$，分母 $\mathrm{den}$ 一次算好、实虚部各除一次即可。
- **把零跳过阈值抬到最小合法 $|\sin\theta| \approx 0.195$ 之上**（例如按"接近 π
  就算发散"的直觉取 0.5）：偶数阶的那对合法零点会被一并吞掉，静默丢掉一对有限
  零点、把滤波器变成另一个响应。1e-4 已比 `sinf` 的求值误差大三个数量级，不需要
  抬高。
- **纹波边缘也用 ±0.1**：纹波边缘落在带边响应最陡处，f32 量化误差在那里最大，
  故放宽到 ±0.25（代价与现状见上文「边缘增益期望与 ±0.25 窗口」节）。
- **μ 用 `asinhf`**：newlib-nano 上可能缺符号导致链接失败，故改写为
  `logf` + `sqrtf` 的恒等式（见上文「为什么 μ 用 `logf` + `sqrtf` 而不是
  `asinhf`」节）。
