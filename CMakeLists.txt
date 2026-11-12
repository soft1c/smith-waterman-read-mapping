cmake_minimum_required(VERSION 3.16)
project(sequence_alignment LANGUAGES CXX)

find_package(Threads REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(sw_mapper src/main.cpp)
target_link_libraries(sw_mapper PRIVATE Threads::Threads)
target_compile_options(sw_mapper PRIVATE
    "$<$<CXX_COMPILER_ID:GNU,Clang>:-Wall;-Wextra;-Wpedantic>"
    "$<$<CXX_COMPILER_ID:MSVC>:/W4>"
)
