cmake_minimum_required(VERSION 3.14)
project(digital_filters VERSION 1.0.0 LANGUAGES C)

# 可选：把库与全部测试一起编成 ASan+UBSan 版本
# （cmake -B build-asan -DBUILD_SANITIZE=ON）。必须在任何 add_library/add_executable
# 之前调用 add_compile_options 才生效。未开 sanitizer 时，栈越界 / 未初始化读会被
# “返回值恰好是 0”掩盖（缺陷实录见 AGENTS.md）。
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

# 库需要 libm（设计管线 tanf/powf/logf/sinhf/coshf，filter_utils 的
# sqrtf/fabsf，init/reset 的 fabsf/isfinite）。以 PUBLIC 声明，使用者无需自行链 m。
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_link_libraries(digital_filters PUBLIC m)
endif()

# 把每个函数/数据项放进各自的段，MCU 链接器才能按函数粒度拉取；这依赖使用方
# 加 --gc-sections（否则该标志不起作用）。flash 实测数字与用法见 README「MCU 部署」。
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
