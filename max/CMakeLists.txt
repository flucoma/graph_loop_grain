# Part of the Fluid Corpus Manipulation Project (http://www.flucoma.org/)
# Copyright 2017-2019 University of Huddersfield.
# Licensed under the BSD-3 License.
# See license.md file in the project root for full license information.
# This project has received funding from the European Research Council (ERC)
# under the European Union’s Horizon 2020 research and innovation programme
# (grant agreement No 725899).

cmake_minimum_required(VERSION 3.11)
################################################################################
# Paths
set(MAX_TARGET_NAMESPACE "" CACHE STRING "Optional namesapce to prepend to targets to prevent collision in superbuilds")
set(MAX_SDK_PATH "" CACHE PATH "Path to the Max SDK")
if(NOT MAX_SDK_PATH)
  message(FATAL_ERROR "Max SDK Path is not set")
endif()

set(FLUID_PATH "" CACHE PATH "Optional path to the flucoma-core repo")# TODO: set to FLUCOMA_CORE_PATh

if (APPLE)
  set(CMAKE_XCODE_GENERATE_SCHEME ON)
  set(CMAKE_XCODE_SCHEME_EXECUTABLE "/Applications/Max.app")
  set(CMAKE_OSX_DEPLOYMENT_TARGET "10.8" CACHE STRING "")
  #A consequence of targetting 10.8. Needs to be set globally from 10.15 onwards in order for the test program to compile successfully during configure
  string(APPEND CMAKE_CXX_FLAGS " -stdlib=libc++")
endif()

################################################################################
project (graph-loop-grain LANGUAGES CXX)

set(CMAKE_XCODE_GENERATE_TOP_LEVEL_PROJECT_ONLY ON)
set(AUTHOR_DOMAIN "org.flucoma")
string(TIMESTAMP "%Y" THISYEAR)
set(COPYRIGHT_STRING "Copyright (c) 2017-${THISYEAR} University of Huddersfield")
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
################################################################################
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
  flucoma-max
  GIT_REPOSITORY https://github.com/flucoma/flucoma-max.git
  GIT_PROGRESS TRUE
  GIT_TAG origin/main
)

FetchContent_GetProperties(flucoma-max)

if(NOT flucoma-max_POPULATED)
  FetchContent_Populate(flucoma-max)
endif()

# Misc setup and subroutines
include("${flucoma-max_SOURCE_DIR}/source/script/max-package.cmake")

# Max Wrapper
add_library( FLUID_MAX INTERFACE )
target_sources(
	FLUID_MAX INTERFACE
	"${flucoma-max_SOURCE_DIR}/source/include/FluidMaxWrapper.hpp"
	"${flucoma-max_SOURCE_DIR}/source/include/MaxBufferAdaptor.hpp"
)

################################################################################
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/fluid.graphgrain_tilde")
source_group("" FILES "fluid.graphgrain_tilde.cpp")

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/fluid.graphloop_tilde")
source_group("" FILES "fluid.graphloop_tilde.cpp")

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/fluid.graphplay_tilde")
source_group("" FILES "fluid.graphplay_tilde.cpp")

install(TARGETS fluid.graphgrain_tilde fluid.graphloop_tilde fluid.graphplay_tilde
        DESTINATION ${CMAKE_BINARY_DIR}/dist/externals
)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/help
        DESTINATION ${CMAKE_BINARY_DIR}/dist/
)
