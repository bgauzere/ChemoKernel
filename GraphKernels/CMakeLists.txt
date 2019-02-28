cmake_minimum_required(VERSION 2.6)

project(GraphKernels)

add_definitions("-O3 -std=c++11 -ansi -pedantic -Dcimg_use_vt100 -Dcimg_display=1 -I/usr/X11R6/include  -lm -L/usr/X11R6/lib -lpthread -lX11")
add_subdirectory(src/pandore)
add_subdirectory(src/gtb)
add_subdirectory(src/libsvm)
add_subdirectory(src/GraphKernels)

