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

# Emit each function/data item in its own section so MCU linkers with
# --gc-sections only pull in the design pipeline of the families actually
# used (measured: butter_lp_2nd alone 36.7 KB → 14.4 KB of flash).
# Consumers must link with --gc-sections (or the equivalent linker flag)
# to benefit; the flag is inert otherwise.
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
