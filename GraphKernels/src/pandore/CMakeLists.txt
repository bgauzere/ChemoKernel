include_directories(../../include/pandore)

file(GLOB_RECURSE headers
    ../../include/pandore/*.h
)

add_library(pandore STATIC
  objects/bundled.cpp
  objects/collection.cpp
  objects/errc.cpp
  objects/graph.cpp
  objects/main.cpp
  objects/pobject.cpp
  objects/region.cpp)

INSTALL(FILES ${headers} DESTINATION include/pandore)
INSTALL(TARGETS pandore DESTINATION lib/)
