cmake_minimum_required(VERSION 3.10)
project(MolecularDynamics)

# Require C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -Wall")

# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}) - This moves md to the source file so the text outputs can be plotted

find_package(Boost REQUIRED COMPONENTS program_options)
include_directories(${Boost_INCLUDE_DIRS})  
find_package(BLAS REQUIRED)

add_executable(md maincmdtests.cpp)
target_link_libraries(md ${Boost_LIBRARIES} ${BLAS_LIBRARIES})

add_custom_target(default ALL DEPENDS md)

add_custom_target(clean_all 
    COMMAND rm -rf CMakeFiles CMakeCache.txt cmake_install.cmake Makefile md build
    COMMENT "Cleaning build files"
)

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
# add_custom_target(run
#     COMMAND ./md ${ARGS}
#     DEPENDS md
#     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#     COMMENT "Running Molecular Dynamics Simulation"
# )
# 🔹 Doxygen Documentation

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
