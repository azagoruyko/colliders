cmake_minimum_required(VERSION 3.10)

project (colliders) 
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

find_package( Maya REQUIRED )

set(sources
  sources/main.cpp 
  sources/bellCollider.cpp
  sources/bellCollider.h
  sources/bellColliderSolver.cpp
  sources/bellColliderSolver.h
  sources/planeCollider.cpp
  sources/planeCollider.h
  sources/skirtBellCollider.cpp
  sources/skirtBellCollider.h
  sources/utils.hpp)

add_library(colliders SHARED ${sources})

MAYA_PLUGIN( colliders )