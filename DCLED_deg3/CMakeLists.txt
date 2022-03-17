cmake_minimum_required(VERSION 3.17)
project(DCLED_deg3 C)

set(CMAKE_C_STANDARD 99)

add_executable(DCLED_deg3 main.c tool.h Eval.h)

target_link_libraries(DCLED_deg3 libgmp.so libflint.so m -lpthread)