cmake_minimum_required(VERSION 3.15)
project(pitch C)

set(CMAKE_C_STANDARD 99)

add_executable(pitch main.c tinywav.c kiss_fft.c kiss_fftr.c)