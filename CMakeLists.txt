cmake_minimum_required(VERSION 3.10)
project(MolecularDynamics)

# Require C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -Wall")

# Find Required Libraries
find_package(Boost REQUIRED COMPONENTS program_options)
include_directories(${Boost_INCLUDE_DIRS})  
find_package(BLAS REQUIRED)
find_package(OpenMP REQUIRED)

# Serial Version
add_executable(md mdmS.cpp compute_serial.cpp init_testing.cpp)
target_link_libraries(md PUBLIC ${Boost_LIBRARIES} ${BLAS_LIBRARIES})

# Parallel Version
add_executable(mdpar mdmPsingNP.cpp compute_parallel.cpp init_testing.cpp)
target_link_libraries(mdpar PUBLIC OpenMP::OpenMP_CXX ${Boost_LIBRARIES} ${BLAS_LIBRARIES})

# Default Build Target
add_custom_target(default ALL DEPENDS md mdpar)

# Clean Target
add_custom_target(clean_all 
    COMMAND rm -rf CMakeFiles CMakeCache.txt cmake_install.cmake Makefile md mdpar build
    COMMENT "Cleaning build files"
)

# Serial Test Cases
add_custom_target(run_tests
    COMMAND ./md --ic-one
    COMMAND ./md --ic-one-vel
    COMMAND ./md --ic-two
    COMMAND ./md --ic-two-pass1
    COMMAND ./md --ic-two-pass2
    COMMAND ./md --ic-two-pass3
    DEPENDS md
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running all test cases for Molecular Dynamics Simulation"
)

# Parallel Test Cases
add_custom_target(run_tests_parallel
    COMMAND OMP_NUM_THREADS=1 ./mdpar --ic-one
    COMMAND OMP_NUM_THREADS=1 ./mdpar --ic-one-vel
    COMMAND OMP_NUM_THREADS=1 ./mdpar --ic-two
    COMMAND OMP_NUM_THREADS=1 ./mdpar --ic-two-pass1
    COMMAND OMP_NUM_THREADS=1 ./mdpar --ic-two-pass2
    COMMAND OMP_NUM_THREADS=1 ./mdpar --ic-two-pass3
    DEPENDS mdpar
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running all test cases for Parallel Molecular Dynamics Simulation"
)

# Doxygen Documentation
find_package(Doxygen REQUIRED)
set(DOXYGEN_IN ${CMAKE_SOURCE_DIR}/Doxyfile)
set(DOXYGEN_OUT ${CMAKE_BINARY_DIR}/docs)

add_custom_target(doc
    COMMAND ${CMAKE_COMMAND} -E make_directory ${DOXYGEN_OUT} 
    COMMAND doxygen ${DOXYGEN_IN}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM
)
