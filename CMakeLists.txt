cmake_minimum_required(VERSION 3.10)

project(MemoriaLCA_LZ VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(memoria_lca src/main.cpp)

# agregamos headers
target_include_directories(memoria_lca PUBLIC include)
