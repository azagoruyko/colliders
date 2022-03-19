cmake_minimum_required(VERSION 2.6)

if(CMAKE_CONFIGURATION_TYPES)
   set(CMAKE_CONFIGURATION_TYPES Debug Release RelWithDebInfo)
   set(CMAKE_CONFIGURATION_TYPES "${CMAKE_CONFIGURATION_TYPES}" CACHE STRING
     "Reset the configurations to Debug and Release"
     FORCE)
endif()

if (NOT CMAKE_BUILD_TYPE) 
    set (CMAKE_BUILD_TYPE "Release") 
endif ()

project (colliders) 
set(CMAKE_MODULE_PATH $ENV{MYTONA_PATH}/resources/cpp/cmake)
SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)

# Packages
if (MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
endif ()

find_package( Maya REQUIRED )

ADD_DEFINITIONS(-D__TBB_MAKE_EXCEPTION_PTR_PRESENT -D__TBB_NO_IMPLICIT_LINKAGE)

set(sources
  sources/main.cpp 
  sources/bellCollider.cpp
  sources/bellCollider.h
  sources/utils.hpp)

add_library(colliders SHARED ${sources})

target_include_directories(colliders PUBLIC 
  ${typeIdDir} 
  ${MAYA_INCLUDE_DIR})

target_link_libraries( 
  colliders 
  ${MAYA_Foundation_LIBRARY} 
  ${MAYA_OpenMaya_LIBRARY} 
  ${MAYA_OpenMayaUI_LIBRARY} 
  ${MAYA_OpenMayaRender_LIBRARY} 
  ${MAYA_tbb_LIBRARY})

set_target_properties( colliders PROPERTIES OUTPUT_NAME colliders CLEAN_DIRECT_OUTPUT 1 SUFFIX ".mll" )

MAYA_SET_PLUGIN_PROPERTIES( colliders )

install(TARGETS colliders RUNTIME DESTINATION .)

