cmake_minimum_required(VERSION 3.14)
project(digital_filters VERSION 1.0.0 LANGUAGES C)

# ── Library ────────────────────────────────────────────────────────────────

add_library(digital_filters STATIC
    src/biquad_filter.c
    src/filter_utils.c
    src/butter_filter.c
    src/cheby_filter.c
)

target_include_directories(digital_filters
    PUBLIC include
)

# 把每个函数/数据项放进各自的段，使带 --gc-sections 的 MCU 链接器
# 只拉入实际用到的滤波器族。实测（arm-none-eabi-gcc / Cortex-M4 / -O2）：
# 库的四个目标文件 .text 合计约 19.8 KB，而只用 butter_lp_2nd 的程序
# 链接后 text 仅约 1.4 KB。缺这些标志时链接器按目标文件粒度工作，
# 会把 butter_filter.o 整体拉入（4.3 KB，含全部 32 个 init）。
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
    add_subdirectory(test)
endif()

# ── Installation ───────────────────────────────────────────────────────────

include(GNUInstallDirs)

install(TARGETS digital_filters
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
