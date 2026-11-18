cmake_minimum_required(VERSION 3.24)

cmake_policy(SET CMP0135 NEW)

########################################################################
####                        Project metadata                        ####
########################################################################

project(LNOT
	VERSION 1.0.0
	DESCRIPTION "Lightweight Numerical Optimization Tools"
	LANGUAGES CXX)

########################################################################
####                           Variables                            ####
########################################################################

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build.")
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif()

set(LNOT_SAN_FLAGS -fsanitize=address -fsanitize=undefined -fno-omit-frame-pointer)

set(LNOT_COMMON_COMPILE_OPTIONS
	-Wall
	-Wcast-align
	-Wconversion
	-Wdouble-promotion
	-Werror
	-Wextra
	-Wformat=2
	-Winit-self
	-Wmissing-declarations
	-Wno-c99-extensions
	-Wno-error=array-bounds
	-Wnon-virtual-dtor
	-Wnull-dereference
	-Wold-style-cast
	-Woverloaded-virtual
	-Wpedantic
	-Wrange-loop-construct
	-Wshadow
	-Wsign-conversion
	-Wswitch-enum
	-Wundef
	-Wuninitialized
	-Wvarargs
	-Wvexing-parse
	-flax-vector-conversions
	-fopenmp-simd
	-pedantic
	-pedantic-errors)

set(LNOT_GNU_CXX_COMPILE_OPTIONS ${LNOT_COMMON_COMPILE_OPTIONS})

set(LNOT_Clang_CXX_COMPILE_OPTIONS ${LNOT_COMMON_COMPILE_OPTIONS})

########################################################################
####                            Options                             ####
########################################################################

option(LNOT_ENABLE_SANITIZERS "Enable sanitizers in Debug builds"  OFF)
option(LNOT_BUILD_DOC         "Build Doxygen documentation"        OFF)
option(LNOT_USE_NATIVE_OPT    "Enable -march=native -mtune=native" OFF)
option(LNOT_BUILD_TESTS       "Build tests for LNOT"               OFF)

########################################################################
####                             Tests                              ####
########################################################################

if(${LNOT_BUILD_TESTS})
	enable_testing()
endif()

########################################################################
####                          Dependencies                          ####
########################################################################

########################################################################
####                         SubDirectories                         ####
########################################################################

if(${LNOT_BUILD_TESTS})
	add_subdirectory(tests)
endif()

########################################################################
####                            Targets                             ####
########################################################################

add_library(LNOT_include INTERFACE)

add_library(LNOT_options INTERFACE)

add_library(LNOT)

########################################################################
####                            Sources                             ####
########################################################################

target_sources(LNOT PRIVATE
	src/LNOT/Containers/ArrayHandle.cpp
	src/LNOT/Containers/CircularBuffer.cpp
	src/LNOT/Containers/HostAllocator.cpp)

target_include_directories(LNOT_include INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

########################################################################
####                              Link                              ####
########################################################################

target_link_libraries(LNOT
	PUBLIC LNOT_include
	PUBLIC LNOT_options)

########################################################################
####                        Compile options                         ####
########################################################################

target_compile_features(LNOT_options INTERFACE cxx_std_23)

target_compile_options(LNOT_options INTERFACE
	$<$<CXX_COMPILER_ID:GNU>:${LNOT_GNU_CXX_COMPILE_OPTIONS}>
	$<$<CXX_COMPILER_ID:Clang>:${LNOT_Clang_CXX_COMPILE_OPTIONS}>)

if(LNOT_ENABLE_SANITIZERS)
	target_compile_options(LNOT_options INTERFACE $<$<CONFIG:Debug>:${SAN_FLAGS}>)
endif()

if(LNOT_USE_NATIVE_OPT)
	target_compile_options(LNOT_options INTERFACE -march=native -mtune=native)
endif()
