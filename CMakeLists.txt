cmake_minimum_required(VERSION 3.15)
project(StressApp CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(mohr main.cpp)

target_include_directories(mohr PRIVATE "C:/msys64/mingw64/include/eigen3")

target_compile_options(mohr PRIVATE -O3 -march=native)

target_link_libraries(mohr PRIVATE
    user32
    gdi32
    opengl32
    gdiplus
    Shlwapi
    dwmapi
    stdc++fs
)