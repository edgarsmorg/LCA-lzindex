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

# Headers del proyecto
include_directories(${CMAKE_SOURCE_DIR}/src)

# Archivos fuente del proyecto
set(PROJECT_SOURCES
    src/main.cpp
    src/index.cpp
    src/lz77/parser.cpp
    src/lz77/grid.cpp
)

# Ejecutable principal
add_executable(memoria_lca ${PROJECT_SOURCES})

# Habilitamos testing
enable_testing()

# Directorios de tests y benchmarks (a agregar después)
# add_subdirectory(tests)
# add_subdirectory(bench)
