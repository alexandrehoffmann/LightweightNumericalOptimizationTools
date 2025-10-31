cmake_minimum_required(VERSION 3.14)
project(LightweightNumericalOptimizationTools VERSION 1.0.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# List of compiler flags
set(LNOT_GENERAL_COMPILE_COMPILE_OPTIONS
    -fopenmp-simd
    -Wundef
    -Wvarargs
    -Wall
    -Wextra
    -Winit-self
    -Wpedantic
    -Werror
    -Wconversion
    -Wuninitialized
    -Wmissing-declarations
    -Wsign-conversion
    -Wshadow
    -Wcast-align
    -Wnull-dereference
    -Wformat=2
    -flax-vector-conversions
    -pedantic
    -pedantic-errors
    -Wno-error=array-bounds
    -Wno-c99-extensions
)

set(LNOT_COMPILE_C_COMPILE_OPTIONS
    -Wstrict-prototypes
    -Wmissing-prototypes
    -Wnested-externs
    -Wbad-function-cast
    -Wmissing-parameter-type
    -Wimplicit-function-declaration
    -Wreturn-type
)

set(LNOT_COMPILE_CXX_COMPILE_OPTIONS
	-Wrange-loop-construct
	-Wnon-virtual-dtor
	-Wold-style-cast
	-Woverloaded-virtual
)

# Set default build type if not specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)

    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif()

# === Options ===
option(LNOT_WITH_BLAS  "Enable BLAS support via pkg-config" OFF)
option(LNOT_NATIVE_OPT "Enable CPU-specific optimizations (-march=native -mtune=native)" OFF)
option(LNOT_BUILD_DEMO "Build demo executable" OFF)
option(LNOT_BUILD_DOC  "Build Doxygen documentation" OFF)

# === Dependencies ===
find_package(fmt REQUIRED)

if(LNOT_WITH_BLAS)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(OpenBLAS REQUIRED openblas)

    message(STATUS "OpenBLAS found: ${OpenBLAS_INCLUDE_DIRS}")
endif()

# === Doxygen ===
if(LNOT_BUILD_DOC)
    find_package(Doxygen REQUIRED dot)

    # Set default values for Doxyfile
    set(DOXYGEN_PROJECT_NAME "Lightweight Numerical Optimization Tools (LNOT)")
    set(DOXYGEN_INPUT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
    set(DOXYGEN_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/doc")

    set(DOXYGEN_FILE ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in ${DOXYGEN_FILE} @ONLY)

    add_custom_target(doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_FILE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM
    )
endif()

# === Demo subdir ===
if(LNOT_BUILD_DEMO)
    add_subdirectory(demo)
endif()

# === Library Target ===
set(LNOT_SRC
	src/LNOT/Containers/CircularBuffer.cpp
	#
	src/LNOT/BasicLinalg/BasicLinalg.c
	src/LNOT/BasicLinalg/SymmetricDenseMatrixOp.cpp
	src/LNOT/BasicLinalg/DiagonalPreconditionerOp.cpp
	#
	src/LNOT/LinearSolvers/ConjugateGradient.cpp
	src/LNOT/LinearSolvers/LanczosSolver.cpp
	#
	src/LNOT/TRSSolvers/TruncatedConjugateGradient.cpp
	src/LNOT/TRSSolvers/LanczosTRSSolver.cpp
	#
	src/LNOT/LineSearch/BacktrackingLineSearch.cpp
	src/LNOT/LineSearch/BisectionLineSearch.cpp
	#
	src/LNOT/GeneralSolvers/NewtonSolver.cpp
	src/LNOT/GeneralSolvers/BFGS.cpp
	src/LNOT/GeneralSolvers/LBFGS.cpp
	src/LNOT/GeneralSolvers/NonLinearConjugateGradient.cpp
	src/LNOT/GeneralSolvers/NewtonTrustRegionSolver.cpp
	src/LNOT/GeneralSolvers/SR1TrustRegionSolver.cpp
	src/LNOT/GeneralSolvers/LSR1TrustRegionSolver.cpp
)

add_library(LNOT STATIC ${LNOT_SRC})

target_include_directories(LNOT
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(LNOT PUBLIC fmt::fmt)

if(LNOT_WITH_BLAS)
    target_include_directories(LNOT PUBLIC ${OpenBLAS_INCLUDE_DIRS})
    target_link_libraries(LNOT PUBLIC ${OpenBLAS_LIBRARIES})
    target_compile_definitions(LNOT PUBLIC LNOT_WITH_BLAS)
endif()

target_compile_options(LNOT PUBLIC
	$<$<COMPILE_LANGUAGE:C>:${LNOT_GENERAL_COMPILE_COMPILE_OPTIONS} ${LNOT_COMPILE_C_COMPILE_OPTIONS}>
    $<$<COMPILE_LANGUAGE:CXX>:${LNOT_GENERAL_COMPILE_COMPILE_OPTIONS} ${LNOT_COMPILE_CXX_COMPILE_OPTIONS}>
)

if(LNOT_NATIVE_OPT)
    message(STATUS "Enabling native CPU optimizations (-march=native -mtune=native)")
    target_compile_options(LNOT PRIVATE -march=native -mtune=native)
endif()

set_source_files_properties(
    src/LNOT/BasicLinalg/BasicLinalg.c
    PROPERTIES
        COMPILE_OPTIONS "-ffast-math"
)

# === Installation ===
include(GNUInstallDirs)

install(DIRECTORY include/LNOT DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(TARGETS LNOT
    EXPORT LNOTTargets
)

install(EXPORT LNOTTargets
    FILE LNOTTargets.cmake
    NAMESPACE LNOT::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/LNOT
)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/LNOTConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/LNOTConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/LNOTConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/LNOT
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/LNOTConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/LNOTConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/LNOT
)

# === pkg-config ===
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/lnot.pc.in
    ${CMAKE_CURRENT_BINARY_DIR}/lnot.pc
    @ONLY
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/lnot.pc
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
)
