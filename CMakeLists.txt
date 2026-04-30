cmake_minimum_required(VERSION 3.14)
project(digital_filters VERSION 1.0.0 LANGUAGES C)

# ── Library ────────────────────────────────────────────────────────────────

add_library(digital_filters STATIC
    src/biquad_filter.c
    src/sos_filter.c
    src/filter_utils.c
    src/butter_filter.c
    src/cheby_filter.c
    src/static_butter_filter.c
    src/static_cheby_filter.c
)

target_include_directories(digital_filters
    PUBLIC include
)

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
