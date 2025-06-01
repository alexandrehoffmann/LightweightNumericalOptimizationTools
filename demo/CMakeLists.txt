add_executable(demo_LinearSolver                demo_LinearSolver.cpp)
add_executable(demo_TrustRegionSubproblemSolver demo_TrustRegionSubproblemSolver.cpp)

target_include_directories(demo_LinearSolver PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(demo_TrustRegionSubproblemSolver PRIVATE ${PROJECT_SOURCE_DIR}/include)

target_link_libraries(demo_LinearSolver                PRIVATE LNOT)
target_link_libraries(demo_TrustRegionSubproblemSolver PRIVATE LNOT)

