# cmake/helpers.cmake — funciones helper para definir targets LZ77 y sr-index.
#
# Requiere que las variables SDSL_INCLUDE, SR_INDEX_DIR, SR_BUILD_DIR,
# SR_AVAILABLE estén definidas en el CMakeLists.txt padre antes de incluir
# este archivo.

# ── Source lists compartidas ──────────────────────────────────────────────────

set(_LZ_CORE
    ${CMAKE_SOURCE_DIR}/src/index.cpp
    ${CMAKE_SOURCE_DIR}/src/lz77/parser.cpp
    ${CMAKE_SOURCE_DIR}/src/lz77/grid.cpp
    ${CMAKE_SOURCE_DIR}/src/wavelet/wt_rmq_min.cpp
    ${CMAKE_SOURCE_DIR}/src/wavelet/wm_rmq_min.cpp
)

set(_LZ_INCS
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/tools
    ${SDSL_INCLUDE}
)

set(_SR_INCS
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/tools
    ${SR_INDEX_DIR}/include
    ${SDSL_INCLUDE}
    ${SR_BUILD_DIR}/_deps/json-src/single_include
)

# ── lz_tool(<name> [SOURCES src...]) ─────────────────────────────────────────
# Tool en tools/<name>.cpp que usa _LZ_CORE + -O3.
function(lz_tool NAME)
    cmake_parse_arguments(A "" "" "SOURCES" ${ARGN})
    add_executable(${NAME} tools/${NAME}.cpp ${_LZ_CORE} ${A_SOURCES})
    target_include_directories(${NAME} PRIVATE ${_LZ_INCS})
    target_compile_options(${NAME} PRIVATE -O3)
    message(STATUS "${NAME}: habilitado")
endfunction()

# ── lz_test(<name> [SOURCES src...]) ─────────────────────────────────────────
# Test en tests/<name>.cpp con GTest + -O2. SOURCES son archivos extra (no usa
# _LZ_CORE automáticamente: cada test declara solo lo que necesita).
function(lz_test NAME)
    cmake_parse_arguments(A "" "" "SOURCES" ${ARGN})
    add_executable(${NAME} tests/${NAME}.cpp ${A_SOURCES})
    target_include_directories(${NAME} PRIVATE ${CMAKE_SOURCE_DIR}/src ${SDSL_INCLUDE})
    target_link_libraries(${NAME} PRIVATE GTest::gtest_main)
    target_compile_options(${NAME} PRIVATE -O2)
    gtest_discover_tests(${NAME})
endfunction()

# ── sr_tool(<name> [LZ] [MULTIDEF] [SOURCES src...]) ─────────────────────────
# Tool en tools/<name>.cpp que depende del sr-index.
#   LZ       — incluye también _LZ_CORE (para tools que mezclan ambos índices)
#   MULTIDEF — añade -Wl,--allow-multiple-definition (necesario cuando LZ y sr-index
#              expanden los mismos templates inline)
#   SOURCES  — fuentes extra (p.ej. src/mem/extractor.cpp)
function(sr_tool NAME)
    if(NOT SR_AVAILABLE)
        message(STATUS "${NAME}: deshabilitado (sr-index no encontrado)")
        return()
    endif()
    cmake_parse_arguments(A "LZ;MULTIDEF" "" "SOURCES" ${ARGN})
    set(_srcs tools/${NAME}.cpp $<TARGET_OBJECTS:sr_index_locator_obj> ${A_SOURCES})
    if(A_LZ)
        list(APPEND _srcs ${_LZ_CORE})
    endif()
    add_executable(${NAME} ${_srcs})
    target_include_directories(${NAME} PRIVATE ${_SR_INCS})
    if(A_MULTIDEF)
        target_link_libraries(${NAME} PRIVATE -Wl,--allow-multiple-definition)
    endif()
    target_compile_options(${NAME} PRIVATE -O3)
    message(STATUS "${NAME}: habilitado")
endfunction()

# ── sr_test(<name> [LZ] [MULTIDEF] [SOURCES src...]) ─────────────────────────
# Test en tests/<name>.cpp que depende del sr-index + GTest + -O2.
function(sr_test NAME)
    if(NOT SR_AVAILABLE)
        message(STATUS "${NAME}: deshabilitado (sr-index no encontrado)")
        return()
    endif()
    cmake_parse_arguments(A "LZ;MULTIDEF" "" "SOURCES" ${ARGN})
    set(_srcs tests/${NAME}.cpp $<TARGET_OBJECTS:sr_index_locator_obj> ${A_SOURCES})
    if(A_LZ)
        list(APPEND _srcs ${_LZ_CORE})
    endif()
    add_executable(${NAME} ${_srcs})
    target_include_directories(${NAME} PRIVATE ${_SR_INCS})
    target_link_libraries(${NAME} PRIVATE GTest::gtest_main)
    if(A_MULTIDEF)
        target_link_libraries(${NAME} PRIVATE -Wl,--allow-multiple-definition)
    endif()
    target_compile_options(${NAME} PRIVATE -O2)
    gtest_discover_tests(${NAME})
endfunction()
