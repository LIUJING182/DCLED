cmake_minimum_required(VERSION 3.17)
project(DCLED_deg5 C)

set(CMAKE_C_STANDARD 99)

add_executable(DCLED_deg5 main.c tool.h Eval.h)

target_link_libraries(DCLED_deg5 libgmp.so /usr/local/lib/libflint.so m -lpthread)