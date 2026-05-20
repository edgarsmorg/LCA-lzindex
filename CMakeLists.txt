cmake_minimum_required(VERSION 3.10)

project(MemoriaLCA_LZ VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Flags de compilación obligatorios
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -msse4.2 -mbmi -mbmi2")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra -pedantic")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Headers del proyecto + sdsl (header-only, duscob fork)
set(SDSL_INCLUDE $ENV{HOME}/.local/include)

# sr-index headers (header-only) + nlohmann json ya descargado por su build/
set(SR_INDEX_INCLUDE ${CMAKE_SOURCE_DIR}/external/sr-index/include)
set(SR_INDEX_JSON    ${CMAKE_SOURCE_DIR}/external/sr-index/build/_deps/json-src/single_include)

include_directories(${CMAKE_SOURCE_DIR}/src ${SDSL_INCLUDE} ${SR_INDEX_INCLUDE} ${SR_INDEX_JSON})

# Archivos fuente del proyecto
set(PROJECT_SOURCES
    src/main.cpp
    src/index.cpp
    src/mem/extractor.cpp
    src/lz77/parser.cpp
    src/lz77/grid.cpp
    src/wavelet/wt_rmq_min.cpp
    src/wavelet/wm_rmq_min.cpp
)

# Ejecutable principal
add_executable(memoria_lca ${PROJECT_SOURCES})

# ── Tests con Google Test ─────────────────────────────────────────────────────
include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.14.0
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(test_lz77 tests/test_lz77.cpp src/lz77/parser.cpp)
target_include_directories(test_lz77 PRIVATE ${CMAKE_SOURCE_DIR}/src ${SDSL_INCLUDE})
target_link_libraries(test_lz77 PRIVATE GTest::gtest_main)
target_compile_options(test_lz77 PRIVATE -O2)

include(GoogleTest)
gtest_discover_tests(test_lz77)

add_executable(test_wt_rmq_min
    tests/test_wt_rmq_min.cpp
    src/wavelet/wt_rmq_min.cpp
)
target_include_directories(test_wt_rmq_min PRIVATE ${CMAKE_SOURCE_DIR}/src ${SDSL_INCLUDE})
target_link_libraries(test_wt_rmq_min PRIVATE GTest::gtest_main)
target_compile_options(test_wt_rmq_min PRIVATE -O2)
gtest_discover_tests(test_wt_rmq_min)

add_executable(test_wm_rmq_min
    tests/test_wm_rmq_min.cpp
    src/wavelet/wt_rmq_min.cpp
    src/wavelet/wm_rmq_min.cpp
)
target_include_directories(test_wm_rmq_min PRIVATE ${CMAKE_SOURCE_DIR}/src ${SDSL_INCLUDE})
target_link_libraries(test_wm_rmq_min PRIVATE GTest::gtest_main)
target_compile_options(test_wm_rmq_min PRIVATE -O2)
gtest_discover_tests(test_wm_rmq_min)

add_executable(test_grid
    tests/test_grid.cpp
    src/lz77/parser.cpp
    src/lz77/grid.cpp
    src/wavelet/wt_rmq_min.cpp
    src/wavelet/wm_rmq_min.cpp
)
target_include_directories(test_grid PRIVATE ${CMAKE_SOURCE_DIR}/src ${SDSL_INCLUDE})
target_link_libraries(test_grid PRIVATE GTest::gtest_main)
target_compile_options(test_grid PRIVATE -O2)
gtest_discover_tests(test_grid)

add_executable(test_index
    tests/test_index.cpp
    src/index.cpp
    src/lz77/parser.cpp
    src/lz77/grid.cpp
    src/wavelet/wt_rmq_min.cpp
    src/wavelet/wm_rmq_min.cpp
)
target_include_directories(test_index PRIVATE ${CMAKE_SOURCE_DIR}/src ${SDSL_INCLUDE})
target_link_libraries(test_index PRIVATE GTest::gtest_main)
target_compile_options(test_index PRIVATE -O2)
gtest_discover_tests(test_index)

add_executable(test_classify
    tests/test_classify.cpp
    src/index.cpp
    src/mem/extractor.cpp
    src/lz77/parser.cpp
    src/lz77/grid.cpp
    src/wavelet/wt_rmq_min.cpp
    src/wavelet/wm_rmq_min.cpp
)
target_include_directories(test_classify PRIVATE ${CMAKE_SOURCE_DIR}/src ${SDSL_INCLUDE})
target_link_libraries(test_classify PRIVATE GTest::gtest_main)
target_compile_options(test_classify PRIVATE -O2)
gtest_discover_tests(test_classify)

add_executable(test_mem
    tests/test_mem.cpp
    src/index.cpp
    src/mem/extractor.cpp
    src/lz77/parser.cpp
    src/lz77/grid.cpp
    src/wavelet/wt_rmq_min.cpp
    src/wavelet/wm_rmq_min.cpp
)
target_include_directories(test_mem PRIVATE ${CMAKE_SOURCE_DIR}/src ${SDSL_INCLUDE})
target_link_libraries(test_mem PRIVATE GTest::gtest_main)
target_compile_options(test_mem PRIVATE -O2)
gtest_discover_tests(test_mem)

# Habilitamos testing
enable_testing()

# ── sr_locate: tool de locate sobre el sr-index ───────────────────────────
# duscob/sdsl-lite es header-only (divsufsort incluido inline), no hay .a
set(SR_INDEX_DIR ${CMAKE_SOURCE_DIR}/external/sr-index)
set(SDSL_INCLUDE  $ENV{HOME}/.local/include)
set(SR_BUILD_DIR  ${SR_INDEX_DIR}/build)

if(EXISTS ${SR_INDEX_DIR}/include AND EXISTS ${SDSL_INCLUDE}/sdsl)
    add_executable(sr_locate tools/sr_locate.cpp)
    target_include_directories(sr_locate PRIVATE
        ${SR_INDEX_DIR}/include
        ${SDSL_INCLUDE}
        ${SR_BUILD_DIR}/_deps/gflags-build/include
        ${SR_BUILD_DIR}/_deps/json-src/single_include
    )
    target_link_libraries(sr_locate PRIVATE
        ${SR_BUILD_DIR}/_deps/gflags-build/libgflags_nothreads.a
        pthread
    )
    message(STATUS "sr_locate: habilitado")
else()
    message(STATUS "sr_locate: deshabilitado (sdsl o sr-index no encontrado)")
endif()

# ── measure_index: mide tamaño del LZ77-index por componente ─────────────────
add_executable(measure_index tools/measure_index.cpp
    src/index.cpp
    src/lz77/parser.cpp
    src/lz77/grid.cpp
    src/wavelet/wt_rmq_min.cpp
    src/wavelet/wm_rmq_min.cpp
)
target_include_directories(measure_index PRIVATE ${CMAKE_SOURCE_DIR}/src ${SDSL_INCLUDE})
target_compile_options(measure_index PRIVATE -O3)
message(STATUS "measure_index: habilitado")

# ── sr_index_locator: objeto PIMPL que aísla los headers del sr-index ─────────
# Se compila UNA sola vez y se reutiliza en baseline_lca_sr y test_baseline_sr,
# evitando múltiples definiciones de símbolos de sr-index al linkar con index.cpp.
if(EXISTS ${SR_INDEX_DIR}/include AND EXISTS ${SDSL_INCLUDE}/sdsl)
    add_library(sr_index_locator_obj OBJECT
        src/baseline/sr_index_locator.cpp
    )
    target_include_directories(sr_index_locator_obj PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${SR_INDEX_DIR}/include
        ${SDSL_INCLUDE}
        ${SR_BUILD_DIR}/_deps/json-src/single_include
    )
    target_compile_options(sr_index_locator_obj PRIVATE -O3)
endif()

# ── baseline_lca_sr: benchmark del pipeline LCA usando sr-index ───────────────
if(EXISTS ${SR_INDEX_DIR}/include AND EXISTS ${SDSL_INCLUDE}/sdsl)
    add_executable(baseline_lca_sr
        tools/baseline_lca_sr.cpp
        $<TARGET_OBJECTS:sr_index_locator_obj>
    )
    target_include_directories(baseline_lca_sr PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${SR_INDEX_DIR}/include
        ${SDSL_INCLUDE}
        ${SR_BUILD_DIR}/_deps/json-src/single_include
    )
    target_compile_options(baseline_lca_sr PRIVATE -O3)
    message(STATUS "baseline_lca_sr: habilitado")
else()
    message(STATUS "baseline_lca_sr: deshabilitado (sdsl o sr-index no encontrado)")
endif()

# ── test_baseline_sr: tests del pipeline LCA con SrIndexLocator ──────────────
# No linkea con index.cpp/lz77/ para evitar colisiones de símbolos sr-index.
# El cross-check LZ77 vs sr-index se realiza en baseline_lca_sr (herramienta).
if(EXISTS ${SR_INDEX_DIR}/include AND EXISTS ${SDSL_INCLUDE}/sdsl)
    add_executable(test_baseline_sr
        tests/test_baseline_sr.cpp
        src/mem/extractor.cpp
        $<TARGET_OBJECTS:sr_index_locator_obj>
    )
    target_include_directories(test_baseline_sr PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${SR_INDEX_DIR}/include
        ${SDSL_INCLUDE}
        ${SR_BUILD_DIR}/_deps/json-src/single_include
    )
    target_link_libraries(test_baseline_sr PRIVATE GTest::gtest_main)
    target_compile_options(test_baseline_sr PRIVATE -O2)
    gtest_discover_tests(test_baseline_sr)
    message(STATUS "test_baseline_sr: habilitado")
else()
    message(STATUS "test_baseline_sr: deshabilitado (sdsl o sr-index no encontrado)")
endif()

# Directorios de tests y benchmarks (a agregar después)
# add_subdirectory(tests)
# add_subdirectory(bench)
