cmake_minimum_required(VERSION 3.10)
project(MemoriaLCA_LZ VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -msse4.2 -mbmi -mbmi2")
set(CMAKE_CXX_FLAGS_DEBUG   "-g -Wall -Wextra -pedantic")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# ── Paths de dependencias ─────────────────────────────────────────────────────
set(SDSL_INCLUDE     $ENV{HOME}/.local/include)
set(SR_INDEX_DIR     ${CMAKE_SOURCE_DIR}/external/sr-index)
set(SR_BUILD_DIR     ${SR_INDEX_DIR}/build)

include_directories(${CMAKE_SOURCE_DIR}/src ${SDSL_INCLUDE}
    ${SR_INDEX_DIR}/include
    ${SR_BUILD_DIR}/_deps/json-src/single_include)

# ── Disponibilidad del sr-index ───────────────────────────────────────────────
if(EXISTS ${SR_INDEX_DIR}/include AND EXISTS ${SDSL_INCLUDE}/sdsl)
    set(SR_AVAILABLE TRUE)
    message(STATUS "sr-index: disponible")
else()
    set(SR_AVAILABLE FALSE)
    message(STATUS "sr-index: no encontrado — targets sr_* deshabilitados")
endif()

# ── Google Test (vía FetchContent) ────────────────────────────────────────────
include(FetchContent)
FetchContent_Declare(googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.14.0
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
include(GoogleTest)
enable_testing()

include(cmake/helpers.cmake)

# ── PIMPL del sr-index ────────────────────────────────────────────────────────
# Compilado una sola vez y referenciado con $<TARGET_OBJECTS:sr_index_locator_obj>
# para evitar múltiples definiciones de símbolos del sr-index al linkar.
if(SR_AVAILABLE)
    add_library(sr_index_locator_obj OBJECT src/baseline/sr_index_locator.cpp)
    target_include_directories(sr_index_locator_obj PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${SR_INDEX_DIR}/include
        ${SDSL_INCLUDE}
        ${SR_BUILD_DIR}/_deps/json-src/single_include
    )
    target_compile_options(sr_index_locator_obj PRIVATE -O3)
endif()

# ── Tests LZ77 ───────────────────────────────────────────────────────────────
lz_test(test_lz77            SOURCES src/lz77/parser.cpp)
lz_test(test_wm_rmq_min      SOURCES src/wavelet/wm_rmq_min.cpp)
lz_test(test_grid            SOURCES src/lz77/grid.cpp
                                      src/wavelet/wm_rmq_min.cpp)
lz_test(test_index           SOURCES src/index.cpp src/lz77/parser.cpp src/lz77/grid.cpp
                                      src/wavelet/wm_rmq_min.cpp)
lz_test(test_classify        SOURCES src/index.cpp src/mem/extractor.cpp
                                      src/lz77/parser.cpp src/lz77/grid.cpp
                                      src/wavelet/wm_rmq_min.cpp)
lz_test(test_mem             SOURCES src/index.cpp src/mem/extractor.cpp
                                      src/lz77/parser.cpp src/lz77/grid.cpp
                                      src/wavelet/wm_rmq_min.cpp)
lz_test(test_index_roundtrip SOURCES src/index.cpp src/lz77/parser.cpp src/lz77/grid.cpp
                                      src/wavelet/wm_rmq_min.cpp)

# ── Tools LZ77 ───────────────────────────────────────────────────────────────
lz_tool(measure_index)
lz_tool(build_lz_index)
lz_tool(inspect_lz_index)
lz_tool(verify_lca_equiv)

# ── sr_locate: caso especial (usa gflags, sin sr_index_locator_obj) ───────────
if(SR_AVAILABLE)
    add_executable(sr_locate tools/sr_locate.cpp)
    target_include_directories(sr_locate PRIVATE
        ${SR_INDEX_DIR}/include ${SDSL_INCLUDE}
        ${SR_BUILD_DIR}/_deps/gflags-build/include
        ${SR_BUILD_DIR}/_deps/json-src/single_include
    )
    target_link_libraries(sr_locate PRIVATE
        ${SR_BUILD_DIR}/_deps/gflags-build/libgflags_nothreads.a pthread
    )
    message(STATUS "sr_locate: habilitado")
else()
    message(STATUS "sr_locate: deshabilitado (sr-index no encontrado)")
endif()

# ── Tools sr-index ────────────────────────────────────────────────────────────
sr_tool(build_sr_index)
sr_tool(inspect_sr_index)
sr_tool(compare_classifiers     LZ MULTIDEF SOURCES src/mem/extractor.cpp)
sr_tool(bench_classify_compare  LZ MULTIDEF SOURCES src/mem/extractor.cpp)
sr_tool(bench_locate_compare    LZ MULTIDEF)

# ── Tests sr-index ────────────────────────────────────────────────────────────
sr_test(test_baseline_sr     SOURCES src/mem/extractor.cpp)
sr_test(test_compare_indexes LZ MULTIDEF SOURCES src/mem/extractor.cpp)
