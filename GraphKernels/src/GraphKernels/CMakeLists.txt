cmake_minimum_required(VERSION 2.6)

# Headers inclusion

include_directories(../../include/GraphKernels)
include_directories(../../include/pandore ../../include/gtb ../../include/cimg ../../include/libsvm)

# GraphKernels library creation

file(GLOB_RECURSE sources
    *.cpp
)

file(GLOB_RECURSE headers
    ../../include/GraphKernels/*.h
    ../../include/GraphKernels/*.hpp
)

add_library(GraphKernels STATIC
  ${sources}
  ${headers}
  )

target_link_libraries(GraphKernels gtb svm)

INSTALL(FILES ${headers} DESTINATION include/GraphKernels)
INSTALL(TARGETS GraphKernels DESTINATION lib/)
