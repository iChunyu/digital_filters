cmake_minimum_required(VERSION 3.14)
project(digital_filters VERSION 1.0.0 LANGUAGES C)

# 可选：把库与全部测试一起编成 ASan+UBSan 版本
# （cmake -B build-asan -DBUILD_SANITIZE=ON && cmake --build build-asan &&
#   cd build-asan && ctest）。必须在任何 add_library/add_executable 之前调用
# add_compile_options 才能生效。design_filter 的 BP/BS 栈越界、proto_zeros
# 未初始化读这类缺陷只有开 sanitizer 才能在常规测试里被抓到——把“返回值恰好
# 是 0”当成测试通过会掩盖 UB。
option(BUILD_SANITIZE "Build library and tests with ASan+UBSan" OFF)
if(BUILD_SANITIZE AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer -g)
    add_link_options(-fsanitize=address,undefined)
endif()

# 严格警告：目录级，覆盖库与全部测试（不会随安装头文件传递给使用者）。
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra)
endif()

# ── Library ────────────────────────────────────────────────────────────────

add_library(digital_filters STATIC
    src/biquad_filter.c
    src/filter_utils.c
    src/butter_filter.c
    src/cheby_filter.c
    src/notch_filter.c
    src/peak_filter.c
)

target_include_directories(digital_filters
    PUBLIC include
)

# 库需要 libm：prewarp → tanf，Chebyshev 原型 → powf/logf/sinhf/coshf，
# filter_utils 的 sqrtf/fabsf，以及 init/reset 路径的其余数学函数。
# 以 PUBLIC 声明，使用者无需自行知道要链 m。
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_link_libraries(digital_filters PUBLIC m)
endif()

# 把每个函数/数据项放进各自的段，使带 --gc-sections 的 MCU 链接器按函数粒度
# 拉取。实测（arm-none-eabi-gcc 16.2 / Cortex-M4，.text + .rodata）：
# 六个目标文件 -Os 合计 14686 B、-O2 合计 18204 B；只用 butter_lp_2nd 的最小
# 程序（init 一次 + update 八次，--gc-sections -e main）链接后 text + .rodata
# 实测 11308 B（链 newlib libm）——设计管线裁不掉，因为 design_filter 的运行时
# switch(type) 引用全部四种频率变换。缺这些标志时链接器按目标文件粒度工作，
# 会把 butter_filter.o 整体拉入（-Os 2590 B / -O2 4496 B，含全部 32 个 init）。
# 使用方必须加 --gc-sections（或对应链接器选项）才能受益；
# 否则该标志不起作用。
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(digital_filters PRIVATE -ffunction-sections -fdata-sections)
endif()

# ── Options ────────────────────────────────────────────────────────────────

option(BUILD_TESTS "Build test programs" ON)

# ── Tests ──────────────────────────────────────────────────────────────────

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# ── Installation ───────────────────────────────────────────────────────────

include(GNUInstallDirs)

install(TARGETS digital_filters
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
