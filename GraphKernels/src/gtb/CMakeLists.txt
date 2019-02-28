# General stuff
include_directories(../../include/pandore ../pandore/operators ../../include/cimg)

# Library for Graph Tool Box
include_directories(../../include/gtb)

# CImg stuff
add_definitions(${GTB_CIMG_FLAGS})
include_directories(${GTB_CIMG_INCPATH})

file(GLOB_RECURSE sources
    *.cpp
)

file(GLOB_RECURSE headers
    ../../include/gtb/*.h
)

add_library(gtb STATIC
  ${sources}
  ${headers}
)

target_link_libraries(gtb pandore)

INSTALL(FILES ${headers} DESTINATION include/gtb)
INSTALL(FILES ../../include/cimg/CImg.h DESTINATION include/cimg)
INSTALL(TARGETS gtb DESTINATION lib/)
