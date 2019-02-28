cmake_minimum_required(VERSION 2.6)

# Headers inclusion

include_directories(../../include/libsvm)

# libsvm library creation

add_library(svm STATIC
  libsvm.cpp
  ../../include/libsvm/libsvm.h
)

target_link_libraries(svm)

INSTALL(FILES ../../include/libsvm/libsvm.h DESTINATION include/libsvm)
INSTALL(TARGETS svm DESTINATION lib/)
