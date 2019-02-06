# - Basic CMake setup
# Check version meets your requirements
# Declare project, which will configure compiler for us
cmake_minimum_required(VERSION 2.8)
project(SNDrift)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# enable threading
list(APPEND CMAKE_CXX_FLAGS "-pthread -std=c++11 ${CMAKE_CXX_FLAGS}")

#Require ROOT, initially try finding previously installed root
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR} ${CMAKE_MODULE_PATH})

find_package(ROOT REQUIRED Tree Hist Spectrum Geom Physics MathCore)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")

# START of simulation part
# Extra libraries from ROOT
set(ROOT_LIBRARIES ${ROOT_LIBRARIES} -lSpectrum -lGeom -lMathCore)
include_directories(include ${ROOT_INCLUDES})

# Build a dynamic library from our sources
add_library(transportlib SHARED 
  include/thread_pool.hpp
  include/fields.hh 
  include/geomodel.hh
  include/utils.hh
  include/electrode.hh 
  include/getopt_pp.h
  include/ctransport.hh
  src/collection.cpp
  src/getopt_pp.cpp
  src/thread_pool.cpp
  src/fields.cpp 
  src/geomodel.cpp 
  src/utils.cpp 
  src/electrode.cpp )
target_link_libraries(transportlib ${ROOT_LIBRARIES})

#Executables
add_executable(scan.exe examples/scan.cpp)
target_link_libraries(scan.exe ${ROOT_LIBRARIES} transportlib)

add_executable(mcdrift.exe examples/mcdrift.cpp)
target_link_libraries(mcdrift.exe ${ROOT_LIBRARIES} transportlib)

# Build the testing code, tell CTest about it
enable_testing()
set(CMAKE_CXX_STANDARD 11)

# main catch library compile only once
add_library(Catch testing/catch-main.cpp)
target_include_directories(Catch INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/testing)

# add test executable
include_directories(. testing)
add_executable(trial testing/trial.cpp)
target_link_libraries(trial PUBLIC Catch transportlib)

# CTest comand
add_test(NAME trial
  COMMAND trial -s
)
