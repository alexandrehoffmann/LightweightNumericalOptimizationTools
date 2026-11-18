########################################################################
####                       tests subdirectory                       ####
########################################################################

########################################################################
####                          Dependencies                          ####
########################################################################

find_package(GTest REQUIRED)

include(FetchContent)


########################################################################
####                            Targets                             ####
########################################################################

add_library(tests_include INTERFACE)

add_executable(LNOTContainersTests)

########################################################################
####                            Sources                             ####
########################################################################

target_sources(LNOTContainersTests PRIVATE
	src/Containers/ArrayHandleTests.cpp
	src/Containers/CircularBufferTests.cpp)

target_include_directories(tests_include INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

########################################################################
####                      Compile Definitions                       ####
########################################################################

########################################################################
####                              Link                              ####
########################################################################

target_link_libraries(LNOTContainersTests
	PRIVATE tests_include
	PRIVATE LNOT
	PRIVATE GTest::gtest_main)

########################################################################
####                             Tests                              ####
########################################################################

add_test(NAME TestLNOTContainers COMMAND LNOTContainersTests)

