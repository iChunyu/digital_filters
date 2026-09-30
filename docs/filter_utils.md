# filter_utils 共享管线（filter_utils）

**源文件**：`include/filter_utils.h`、`src/filter_utils.c`
**测试**：`tests/test_butter.c`、`tests/test_cheby.c`（标定表、回归档位与实测统计只在此维护，本文只给结论与指针）

Butterworth 与 Chebyshev 两族共用的设计管线与零极点工具：预畸、四种频率变换、
增益折叠、`zpk2sos` 配对、`design_filter` 六步管线、级联增益校验。Notch / Peak 不用
`zpk2sos`（固定单节），但复用 `prewarp`、`biquad_c2d_bilinear` 与
`check_cascade_gains`。

## 预畸

$$
f_{\mathrm{analog}} = \frac{f_s}{\pi}\,\tan\!\left(\frac{\pi\, f_{\mathrm{digital}}}{f_s}\right)
$$

区间 $(0, f_s/2)$ 之外双线性映射无定义（过了 Nyquist，`tan` 变号，根本不存在调用方
能用的模拟频率），`prewarp` 直接返回 NaN。入口写成取反比较
`!(fd > 0) || !(fd < fs*0.5)`，因此 NaN 的 `fd` / `fs` 也 fail-closed 而不是漏进
`tanf`。返回的 NaN 沿 `wc` 传播，最终使 `design_filter` 的增益 `k` 变 NaN，被其
有限性闸（`!isfinite(k)`）兜成直通，而不是悄悄错掉的滤波器——族 API 的入口闸已挡掉
越界，这里只是第二层防线。

双线性变换本身（`bilinear_transform`）的映射式是 $z = (2f_s + s)/(2f_s - s)$；
实现按 $K = 2 f_s$ 写成 $z = (K + s)(K - \bar{s})/|K - s|^2$，直接给出实虚部。

## 四种频率变换

四个函数都原地改写数组，`degree = np - nz` 是 `uint8_t`，容量规则必须由调用方遵守。

| 变换 | 代数 | 极点输出 | 零点输出与追加 |
|---|---|---|---|
| `analog_lp_transform` | $p \mapsto w_c\,p$，$z \mapsto w_c\,z$ | `np` | `nz`，不追加 |
| `analog_hp_transform` | $p \mapsto w_c / p$，$z \mapsto w_c / z$ | `np` | `nz`，调用方追加 `np - nz` 个原点零点 |
| `analog_bp_transform` | $s \mapsto (s^2 + \omega_0^2)/(\xi s)$ | `2·np_old` | `np_old + nz_old`，追加 `np_old - nz_old` 个原点零点 |
| `analog_bs_transform` | $s \mapsto \xi s / (s^2 + \omega_0^2)$ | `2·np_old` | `2·np_old`，追加 `2·(np_old - nz_old)` 个 $\pm j\omega_0$ 零点 |

其中 $\omega_0 = \sqrt{\omega_1 \omega_2}$、$\xi = \omega_2 - \omega_1$。BP/BS 阶数
翻倍：每个极点/零点经二次公式一分为二（BP 解 $s^2 - \xi p s + \omega_0^2 = 0$，BS 解
$s^2 - (\xi/p) s + \omega_0^2 = 0$），极点原地双写要求数组容量 $\ge 2\,np$。

**各自要追加的原型无穷远零点**：

- LP 不追加（$s = \infty$ 的零点留在 $s = \infty$，双线性后由 `design_filter` 补
  $z = -1$）；
- HP 追加 `np - nz` 个原点零点（$s = \infty$ 在 $s \mapsto w_c/s$ 下变 $s = 0$，
  双线性后到 $z = +1$）；
- BP 追加 `np_old - nz_old` 个原点零点（双线性后到 $z = +1$），另一半 $s = \infty$
  的零点仍由 `design_filter` 补 $z = -1$；
- BS 追加 $2\,(np_old - nz_old)$ 个 $\pm j\omega_0$ 零点（有限，落在阻带中心），
  **不做** $z = -1$ 补零。

## 三式增益折叠

zpk 形式 $H(s) = k \prod(s - z_i) / \prod(s - p_i)$ 在各变换下，把新的主项系数折回
$k$：

- `zpk_hp_bs_gain`（LP→HP 与 LP→BS 共用）：在**变换前**的原型零极点上
  $k' = k \cdot \mathrm{Re}\!\left(\frac{\prod(-z_i)}{\prod(-p_i)}\right)$；`nz == 0`
  时 $\prod(-z_i)$ 取 1。
- `bilinear_zpk_gain`（双线性）：在**双线性代换前**的 s 域零极点上
  $k_z = k \cdot \mathrm{Re}\!\left(\frac{\prod(K - z_i)}{\prod(K - p_i)}\right)$，
  $K = 2 f_s$。
- `bilinear_zpk_gain_scaled`（LP/BP 的频率变换 + 双线性合并）：
  $k' = k \cdot s^{degree} \cdot \mathrm{Re}\!\left(\frac{\prod(K - z_i)}{\prod(K - p_i)}\right)$，
  $s = w_c$（LP）或 $s = \xi$（BP）。$s^{degree}$ 就是 LP/BP 变换本身带出的
  $w_c^{degree}$ / $\xi^{degree}$ 增益。

三条都用增量乘除、逐对折叠，共轭对下虚部相消，只取实部。f32 精度足以支撑这条增益
跟踪的量化验证见 `tests/verify_zpk_gain.py`。

### `s^degree` 必须与 `∏(K−p)` 交错

`bilinear_zpk_gain_scaled` 里**不能**先单独算 $s^{degree}$ 再乘 $\prod(K-p)$：近
Nyquist 的 $w_c^8 \approx \mathrm{FLT\_MAX}$ 会先溢出，而中间量一旦溢出成 inf，后面
同量级的 $\prod(K-p)$ 因子再也抵消不回来。实现按"每个 $s$ 乘紧跟一个 $(K-p)$ 除、再
处理 $(K-z)/(K-p)$ 对、最后清尾"的顺序交错，让运行乘积始终有界。回归档位见
`tests/test_butter.c` 的「回归：超宽带 BP8（增益链 f32 溢出）」段。

## c_sqrt 与稳定二次求根

`c_sqrt` 求 $\sqrt{re + j\,im}$ 的主支（实部非负）。幅值按最大分量缩放：
$m = \max(|re|, |im|)$，$\sqrt{re^2 + im^2} = m\sqrt{(re/m)^2 + (im/m)^2}$。
近 Nyquist 的 BP/BS 判别式量级 $\sim 1.5\times10^{10}$，裸平方在任一分量超过
$\sqrt{\mathrm{FLT\_MAX}} \approx 1.8\times10^{19}$ 时溢出；缩放把每个中间量压到
$\le \sqrt{2}$。另有 $re < 0$ 且 $|im| \ll |re|$ 的分支，避免 $(mag + re)/2$ 的灾难性
消减。档位见 `tests/test_butter.c` 的「c_sqrt 的缩放计算」段（该段当前只有间接覆盖）。

`stable_roots` 求 $s^2 + B s + C = 0$（$B, C$ 复数）的稳定根：
$r_1 = (-B - \sqrt{B^2 - 4C})/2$，$r_2 = C / r_1$。直接公式 $(-B + \sqrt{\cdots})/2$
在 $|B|^2 \gg 4|C|$ 时灾难性消减（宽带 BP/BS 的实原型极点），小根改用 $C / r_1$ 绕开；
$r_1$ 退化（$|r_1|^2 < 10^{-20}$）时退回直接公式。档位见 `tests/test_butter.c` 的
「回归：含近实极点对的宽带 BP」段。

复数除法 `c_div`（`stable_roots` 与频率变换的除数）在分母为零时把结果置 0。这是
防御路径，正常流程不会发生；即便去掉这层保护，算出的非有限量也会沿增益 `k` 被
`design_filter` 的 `!isfinite(k)` 兜成直通，而不是静默部署。

## zpk2sos

把 $n$ 个 z 域零点与 $n$ 个极点配成 $\lceil n/2 \rceil$ 个二阶节（$n \le 16$，对应
8 阶原型经 BP/BS 翻倍）。算法：

1. **最不利极点优先**：先挑 $|p|$ 最大（离单位圆最近）的未用极点，再找最近的零点；
   最后整体反转节顺序，慢节排后、快节排前。
2. 实数与实数配对，复数与共轭配对；只剩一个实极点时用实零点组成一阶节。
3. **两个容差不可合并**：`eps_class = 1e-5`（实/复分类）与 `eps_claim = 1e-3`
   （共轭认领盒）是不同量。真·实根虚部严格为 `0.0f`，真·共轭对 $|im| \gtrsim 10^{-5}$；
   用 `1e-3` 做分类会把宽带 BP/BS 的合法近实共轭对压成"实"，与另一对成员跨对配对，
   造出极点恰在 $z = 1$ 的节，整个设计随之失败。认领盒则必须吞下双线性变换在单位圆
   附近的 f32 舍入（可把共轭对劈开 $\sim 2\times10^{-4}$）。
4. **认领取盒内最近、不是第一个**：否则高 Q 簇里会偷走别对的伴侣，部署出两节完全
   重复、丢失一对（旧实现曾静默 `valid=1` 部署出重复节；实测缺陷签名见
   `tests/test_cheby.c` 的「回归：高 Q BS 极点簇——共轭配对」段）。
5. **两条收尾不变量，失败即返回 0**：(a) 全部零极点都被认领（有残留说明记账与数组
   脱节、有元素被丢出级联）；(b) 每节的根能复现输入零极点**多重集**（`poly_roots` +
   `match_roots`，防跨对误认领——`used[]` 全置位了但配错了元素）。
6. 总增益 $k$ **只施加到第一节分子**（`sos[0]` 的 $b_0,b_1,b_2$ 乘 $k$），与 scipy
   `zpk2sos` 约定一致；每节分子内部 $b_0 = 1$ 由 `make_biquad` 保证。

配对缺陷与认领盒回归档位见 `tests/test_cheby.c`（「回归：cheby2 BP 五阶…」、「回归：
cheby1 BS 五阶…」、「回归：高 Q BS 极点簇——共轭配对」）。

## design_filter 六步管线

`design_filter` 是导出符号，栈上工作数组固定 `ZPK2SOS_MAX_N = 16`。步骤：

1. 模拟频率变换（`switch` 四分支，含 HP/BS 的 `zpk_hp_bs_gain` 与 BP/BS 的
   `w0 = sqrt(wc1*wc2)`、`xi = wc2 - wc1`）；
2. 双线性增益（在 s 域零极点上算，须在双线性覆写它们之前）；
3. 双线性变换 $s \to z$；
4. 补 $z = -1$ 零点（$s = \infty$ 的零点；BS 不做）；
5. `zpk2sos` 配对成 SOS，节数 $ns = \lceil np/2 \rceil$ 超过 `max_sections` 即拒；
6. 逐节 `biquad_filter_init` 部署，任一节失败即拒。

每个 fail-closed 咽喉点在整条管线里恰好存在一次：

- **入口边界**：`np == 0`、`np > 16`、`nz > np` 直接返回 0——数组越界 / `uint8_t`
  相对阶数下溢。**BP/BS 另有更紧上限 `np > 8`**：这两种类型把每个原型零极点一分为二
  **原地双写**同一个 16 元素数组，原型 9~16 阶会直接越界写栈（ASan 已复现写穿
  `poles[16]`，发生在 BP/BS 分支，而唯一能拦它的 `ns > max_sections` 在越界之后）。
- **`proto_zeros == NULL` 仅当 `nz == 0` 合法**：`nz > 0` 却传 NULL 直接返回 0，不静默
  跳过拷贝、不在**未初始化**的栈数组上继续频率变换/增益折叠/配对（读未初始化值是 UB，
  很可能部署出看似有效的垃圾设计）。
- **`k` 有限非零**：`!isfinite(k) || k == 0` 即拒——`k == 0` 会部署出分子全零的静音
  滤波器，是退化结果不是有效设计。
- **节数上限**：`ns > max_sections` 即拒。
- **逐节 `init`**：任一节过不了 biquad 三道闸即拒。

各咽喉点的标定与回归档位见 `tests/test_butter.c`（尤其「回归：design_filter 输入
边界闸」段）与 `tests/test_cheby.c`。

## check_cascade_gains

错但稳定的零极点配对与增益缩放错误能通过所有逐节检查，却毁掉响应形状。
`check_cascade_gains` 在部署后的 f32 系数上**精确有理求值** H(0) 与 H(π)（无采样、
无三角，`init` 时 $O(\text{节数})$）：

$$
H(0) = \prod_i \frac{b_{0,i} + b_{1,i} + b_{2,i}}{1 + a_{1,i} + a_{2,i}}, \qquad H(\pi) = \prod_i \frac{b_{0,i} - b_{1,i} + b_{2,i}}{1 - a_{1,i} + a_{2,i}}
$$

窗口按期望值分工：期望为精确结构增益（0 或 1，由 $z = \pm 1$ 处的结构零点保证）用
$\pm 0.1$；纹波边缘（cheby1/2 偶数阶的 $10^{-rp/20}$ / $10^{-rs/20}$，落在带边响应
最陡处）放宽到 $\pm 0.25$。$\pm 0.25$ 档的缺陷侧分离只有 $\sim 1.8\times$（$\pm 0.1$
档两侧分离更宽，具体倍数见 `tests/test_cheby.c`），它主要靠"期望值本身可能很小"而非
窗口宽度设防。

`num_sections == 0` 返回 0：空级联两端乘积都是 1.0，带阻类的 $(1, 1)$ 期望会被空洞
通过——0 节不是设计，拒绝而不是盖章放行。

已确认的缺陷实测：配错极点导致 DC 增益 97.9、带阻上 350× 谐振（`tests/test_cheby.c`
的「回归：过去会静默部署错误滤波器的宽带 cheby2 配置」段）。窗口两侧边界见
`tests/test_butter.c` 的「回归：design_filter 输入边界闸」段，标定与分离倍数见
`tests/test_cheby.c` 的「Sweep: cheby1 HP/BP/BS 全阶矩阵」段。

## 测试指针

- 预畸 / 频率变换 / 增益折叠 / `design_filter` 边界：`tests/test_butter.c`（超宽带
  BP8、近 Nyquist LP8、c_sqrt 缩放、近实极点对宽带 BP、design_filter 输入边界闸、
  增益窗口两侧边界）。
- 配对与增益窗口：`tests/test_cheby.c`（实/复混合零点配对、实极点配对、高 Q BS 极点
  簇共轭配对、宽带 cheby2 缺陷、sweep 与窗口标定）。
- `tests/verify_zpk_gain.py`：zpk 增益折叠的 Python 黄金参考（f32 精度充分性）。
