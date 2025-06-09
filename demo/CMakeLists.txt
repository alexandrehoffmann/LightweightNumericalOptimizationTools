find_package(PkgConfig REQUIRED)

add_executable(demo_LinearSolver                demo_LinearSolver.cpp)
add_executable(demo_PrecondLinearSolver         demo_PrecondLinearSolver.cpp)
add_executable(demo_TrustRegionSubproblemSolver demo_TrustRegionSubproblemSolver.cpp)
add_executable(demo_Newton                      demo_Newton.cpp)
add_executable(demo_QuasiNewton                 demo_QuasiNewton.cpp)

target_include_directories(demo_LinearSolver                PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(demo_PrecondLinearSolver         PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(demo_TrustRegionSubproblemSolver PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(demo_Newton                      PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(demo_QuasiNewton                 PRIVATE ${PROJECT_SOURCE_DIR}/include)

target_link_libraries(demo_LinearSolver                PRIVATE LNOT)
target_link_libraries(demo_PrecondLinearSolver         PRIVATE LNOT)
target_link_libraries(demo_TrustRegionSubproblemSolver PRIVATE LNOT)
target_link_libraries(demo_Newton                      PRIVATE LNOT)
target_link_libraries(demo_QuasiNewton                 PRIVATE LNOT)

