# Copyright 2017-2019 University of Huddersfield.
# Licensed under the BSD-3 License.
# See license.md file in the project root for full license information.
# This project has received funding from the European Research Council (ERC)
# under the European Union’s Horizon 2020 research and innovation programme
# (grant agreement No 725899).

cmake_minimum_required(VERSION 3.11)

################################################################################
# Paths
set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/install" CACHE PATH "")

set(SC_LIBRARY_OUTPUT_PREFIX "release-packaging" CACHE STRING "Where in the hierarchy to write sc plugins")

set(SC_PATH "" CACHE PATH "Path to the top of SuperCollider source tree")
if(NOT SC_PATH)
  message(FATAL_ERROR "SuperCollider source path is not set")
endif()

set(FLUID_PATH "" CACHE PATH "Optional path to the Fluid Decomposition repo")

if (APPLE)
  set(CMAKE_XCODE_GENERATE_SCHEME ON)
	set(CMAKE_OSX_ARCHITECTURES "x86_64" CACHE STRING "")
  set(CMAKE_OSX_DEPLOYMENT_TARGET "10.8" CACHE STRING "")
  #A consequence of targetting 10.8. Needs to be set globally from 10.15 onwards in order for the test program to compile successfully during configure
  string(APPEND CMAKE_CXX_FLAGS " -stdlib=libc++")
endif()

################################################################################
# Main project
project (graph-loop-grain LANGUAGES CXX)

if(NOT MSVC)
  add_compile_options(-fdiagnostics-color=always)
endif()

#set correct std lib linking for Windows (in CMake 3.15 this has a native function)
if(MSVC) #needs to be after project()
  foreach(flag_var
      CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
      CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
    if(${flag_var} MATCHES "/MD")
      string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
    endif()
  endforeach()
endif()

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_XCODE_GENERATE_TOP_LEVEL_PROJECT_ONLY ON)

MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()

include(FetchContent)

set(FETCHCONTENT_QUIET FALSE)

FetchContent_Declare(
  flucoma-core
  GIT_REPOSITORY https://github.com/flucoma/flucoma-core.git
  GIT_PROGRESS TRUE
  GIT_TAG origin/main
)

if(FLUID_PATH)
  get_filename_component(
    FETCHCONTENT_SOURCE_DIR_FLUCOMA-CORE ${FLUID_PATH} ABSOLUTE
  )
endif()

FetchContent_GetProperties(flucoma-core)
if(NOT flucoma-core_POPULATED)
  FetchContent_Populate(flucoma-core)
  add_subdirectory(${flucoma-core_SOURCE_DIR} ${flucoma-core_BINARY_DIR})
  include(flucoma_version)
  include(flucoma-buildtools)
  include(flucoma-buildtype)
endif()


FetchContent_Declare(
  flucoma-sc
  GIT_REPOSITORY https://github.com/flucoma/flucoma-sc.git
  GIT_PROGRESS TRUE
  GIT_TAG origin/main
)

FetchContent_GetProperties(flucoma-sc)

if(NOT flucoma-sc_POPULATED)
  FetchContent_Populate(flucoma-sc)
endif()


set(CMAKE_SHARED_MODULE_PREFIX "")
if(APPLE OR WIN32)
  set(CMAKE_SHARED_MODULE_SUFFIX ".scx")
endif()

#needed for complaint-free static linking with GCC
if(CMAKE_COMPILER_IS_GNUCXX)
  target_compile_options( HISSTools_FFT PUBLIC -fPIC )
endif()

#sandbox regrettable dependency on SC internals for SendReply()
add_library(FLUID_SC_COPYREPLYADDR STATIC
  "${flucoma-sc_SOURCE_DIR}/include/wrapper/CopyReplyAddress.cpp"
  "${SC_PATH}/common/SC_Reply.cpp"
  "${SC_PATH}/external_libraries/boost/libs/system/src/error_code.cpp"
)

target_include_directories(FLUID_SC_COPYREPLYADDR SYSTEM PRIVATE
  "${SC_PATH}/include/plugin_interface"
  "${SC_PATH}/include/common"
  "${SC_PATH}/common"
  "${SC_PATH}/external_libraries/boost"
)
set_target_properties(FLUID_SC_COPYREPLYADDR PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)

if(CMAKE_COMPILER_IS_GNUCXX)
  target_compile_options(FLUID_SC_COPYREPLYADDR PUBLIC -fPIC )
endif()

if(APPLE)
  target_compile_options(FLUID_SC_COPYREPLYADDR PRIVATE -stdlib=libc++)
endif()

target_compile_definitions(FLUID_SC_COPYREPLYADDR PRIVATE BOOST_ALL_NO_LIB  BOOST_CONFIG_SUPPRESS_OUTDATED_MESSAGE)

add_library(FLUID_SC_WRAPPER INTERFACE)
target_include_directories(FLUID_SC_WRAPPER
  INTERFACE
  "${flucoma-sc_SOURCE_DIR}/include/"
)

target_sources(FLUID_SC_WRAPPER
  INTERFACE
  "${flucoma-sc_SOURCE_DIR}/include/FluidSCWrapper.hpp"
  "${flucoma-sc_SOURCE_DIR}/include/SCBufferAdaptor.hpp"
)

target_link_libraries(FLUID_SC_WRAPPER INTERFACE FLUID_SC_COPYREPLYADDR)

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/FluidGraphGrain")
source_group("" FILES "FluidGraphGrain.cpp")

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/FluidGraphLoop")
source_group("" FILES "FluidGraphLoop.cpp")

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/FluidGraphPlay")
source_group("" FILES "FluidGraphPlay.cpp")


install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Classes
        DESTINATION ${CMAKE_BINARY_DIR}/dist
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/HelpSource
        DESTINATION ${CMAKE_BINARY_DIR}/dist
)

install(TARGETS FluidGraphGrain FluidGraphLoop FluidGraphPlay
        DESTINATION ${CMAKE_BINARY_DIR}/dist/plugins
)


#[[
#install bits.


install(TARGETS fluid.graphgrain_tilde fluid.graphloop_tilde fluid.graphplay_tilde
        DESTINATION ${CMAKE_BINARY_DIR}/dist
)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/demo
        DESTINATION ${CMAKE_BINARY_DIR}/dist
)



set(SC_INSTALL_PREFIX "." CACHE PATH "Prefix for assembling SC packages")
set(FLUID_PACKAGE_NAME FluidCorpusManipulation CACHE STRING "Name for published package")
set(SC_PACKAGE_ROOT ${SC_INSTALL_PREFIX}/${FLUID_PACKAGE_NAME})

foreach(PACKAGE_DIRECTORY Classes;HelpSource;Examples)
  install(DIRECTORY "release-packaging/${PACKAGE_DIRECTORY}"
          DESTINATION ${SC_PACKAGE_ROOT})
endforeach()

install(DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/
        DESTINATION ${SC_PACKAGE_ROOT}/plugins
        PATTERN "*.ilk" EXCLUDE
        PATTERN "*.PDB" EXCLUDE)
install(DIRECTORY "${flucoma-core_SOURCE_DIR}/AudioFiles"
        DESTINATION ${SC_PACKAGE_ROOT})
install(FILES QuickStart.md
        DESTINATION ${SC_PACKAGE_ROOT})
install(FILES ${flucoma-core_SOURCE_DIR}/distribution.lic
        DESTINATION ${SC_PACKAGE_ROOT}
        RENAME LICENSE.md)
]]
