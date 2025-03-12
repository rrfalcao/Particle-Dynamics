cmake_minimum_required(VERSION 3.10)
project(MolecularDynamics)

# Require C++17
set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -Wall")

find_package(Boost REQUIRED COMPONENTS program_options)
include_directories(${Boost_INCLUDE_DIRS})  
find_package(BLAS REQUIRED)

add_executable(md maincmd.cpp)
target_link_libraries(md ${Boost_LIBRARIES} ${BLAS_LIBRARIES})

add_custom_target(default ALL DEPENDS md)

add_custom_target(clean_all 
    COMMAND rm -rf CMakeFiles CMakeCache.txt md
    COMMENT "Cleaning build files"
)


# add_custom_target(run
#     COMMAND ./md ${ARGS}
#     DEPENDS md
#     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#     COMMENT "Running Molecular Dynamics Simulation"
# )