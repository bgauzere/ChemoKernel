cmake_minimum_required(VERSION 2.6)

# Project : Classification/regression on molecules

project(Molecules)

SET(CMAKE_CXX_FLAGS_DEBUG "-g")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -DVISU")

SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_module/")
MESSAGE("Path to Cmake Modules is " ${CMAKE_MODULE_PATH})

find_package(X11 REQUIRED)
find_package(Threads REQUIRED)
#find_package(GSL REQUIRED)
#enable_language(Fortran)
find_package(LAPACK REQUIRED)
#find_package(gomp REQUIRED)

SET(LINK_CORE_LIBRARIES
	${X11_LIBRARIES}
	${CMAKE_THREAD_LIBS_INIT}
	${LAPACK_LIBRARIES}
	gomp
	#${GSL_LIBRARIES}

)

SET(CMAKE_EXE_LINKER_FLAGS "-lX11 -lpthread -L/usr/lib -lm -llapack -lstdc++ -lgomp ")

add_definitions("-ansi -Wall -Wextra -pedantic -Dcimg_verbosity=3 -fopenmp -Dcimg_use_lapack")


add_subdirectory(src)
