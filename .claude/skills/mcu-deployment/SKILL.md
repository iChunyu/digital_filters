---
name: mcu-deployment
description: MCU / 嵌入式部署本库时的约束与实测数据：硬件 FPU 要求、init 栈峰值、libm 依赖、Flash 粒度（--gc-sections）。当用户问及在 Cortex-M / 裸机上使用、栈或 flash 占用、是否依赖 libm、逐阶链接体积时使用。
---

# MCU / 嵌入式 使用注意事项

**FPU 要求**: 所有滤波器 `_update` 路径执行密集 float 运算。建议使用带硬件 FPU 的 MCU（Cortex-M4/M7 及以上）。

**栈需求**: init 峰值出现在 Chebyshev II 的 BP/BS 调用链上。`-fstack-usage` 实测
（arm-none-eabi-gcc 16.2 / Cortex-M4 / -Os）：`cheby2_bp_init` 200 B +
`design_filter` 528 B + `zpk2sos` 216 B + `match_roots` 48 B = **992 B**；
同一条链在 Cortex-M4 -O0 为 1080 B、x86-64 -O0 为 1136 B。`design_filter` 的帧是
结构性的：poles 128 B + zeros 128 B + sos 192 B，再加逐级调用的帧。
`notch_init` / `peak_init` 走的是直写 biquad 的短路径（无零极点工作数组），
Cortex-M4 -Os 实测各 **80 B**、x86-64 -O2 各 **112 B**，是全场最小的一条 init 链。
运行时 `_update` 全 inline（biquad → 级联 → 逐阶函数三层全部
头文件内联），无额外栈开销。建议 MCU 主栈 ≥ 2 KB。

**libm 依赖**: 只有 biquad 层免 libm——`biquad_filter_init`/`_update`/`_reset` 是纯
乘加代数式，`nm -u biquad_filter.o` 实测为空。**所有族的 init 都要 libm**：设计管线
经 `prewarp()` → `tanf()` 预畸；`nm -u` 实测 `filter_utils.o` 需要
`tanf sqrtf fmaxf memcpy memset`，`cheby_filter.o` 另有
`powf logf cosf sinf sinhf coshf`。Butterworth 的 ROM 极点表只免掉原型计算的
`cosf`/`sinf`，不减免链接——只用 `butter_lp_2nd` 的最小程序不链 libm 时链接失败。

**Flash 粒度**: 库以 `-ffunction-sections/-fdata-sections` 编译；链接时加
`--gc-sections` 按函数粒度拉取。实测（arm-none-eabi-gcc 16.2 / Cortex-M4，
`.text` + `.rodata`）：六个目标文件 -Os 合计 **14686 B**、-O2 合计 **18204 B**；
分项 -Os/-O2 为 biquad_filter 928/1028、filter_utils 4874/5772、butter_filter
2590/4496、cheby_filter 5670/6236、notch_filter 312/332、peak_filter 312/340。
只用 `butter_lp_2nd` 的最小程序（init 一次 + update 八次，`--gc-sections -e main`）
链接后 text + .rodata 实测 **11308 B**（链 newlib libm）——**`--gc-sections` 裁不掉
设计管线**，因为 `design_filter` 的运行时 `switch(type)` 引用全部四种频率变换。
缺这些标志时链接器按目标文件粒度拉取，`butter_filter.o` 会整体
（-Os 2590 B / -O2 4496 B，含全部 32 个 init）被拉入。
update/reset 为头文件 inline，库侧不占 text；每个被调用的 update 会在调用点生成
代码——实测 98 个 update 全部被调用（不调 init）的编译单元 text + .rodata 为
2120 B，只调 1 个时为 128 B（-Os）。

**中断安全**: `_update` 和 `_reset` 不可重入。同一个滤波器结构体如果被 ISR 和主循环共享，需在调用 `_update` 前关中断或使用双缓冲。

**非有限输入**: 热路径无 NaN/Inf 防护（每样本分支开销），非有限输入毒化状态直至 reset；调用方在源头清洗。
