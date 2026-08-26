cmake_minimum_required(VERSION 3.10)
project(seg_trajopt_cpp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Ipopt QUIET)
if(NOT Ipopt_FOUND)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(IPOPT ipopt)
    if(IPOPT_FOUND)
        set(IPOPT_INC_DIR ${IPOPT_INCLUDE_DIRS})
        set(IPOPT_LIBRARIES ${IPOPT_LIBRARIES})
    endif()
endif()

if(NOT Ipopt_FOUND AND NOT IPOPT_FOUND)
    set(IPOPT_POSSIBLE_INC_DIRS /usr/include/coin-or /usr/local/include/coin-or)
    set(IPOPT_POSSIBLE_LIB_DIRS /usr/lib /usr/local/lib /usr/lib/x86_64-linux-gnu)
    find_path(IPOPT_INC_DIR NAMES IpTNLP.hpp PATHS ${IPOPT_POSSIBLE_INC_DIRS})
    find_library(IPOPT_LIBRARIES NAMES ipopt PATHS ${IPOPT_POSSIBLE_LIB_DIRS})
    if(IPOPT_INC_DIR AND IPOPT_LIBRARIES)
        set(IPOPT_FOUND TRUE)
    endif()
endif()

if(NOT IPOPT_FOUND AND NOT Ipopt_FOUND)
    message(FATAL_ERROR "Ipopt not found! Install: sudo apt install coinor-libipopt-dev")
endif()

get_filename_component(TRAJOPT_CPP_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../trajopt_cpp" ABSOLUTE)
if(NOT EXISTS ${TRAJOPT_CPP_DIR}/CMakeLists.txt)
    set(TRAJOPT_CPP_DIR ${CMAKE_SOURCE_DIR}/trajopt_cpp)
endif()
include_directories(${TRAJOPT_CPP_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
if(IPOPT_INC_DIR)
    include_directories(${IPOPT_INC_DIR})
elseif(Ipopt_FOUND)
    include_directories(${IPOPT_INC_DIR})
endif()

# trajopt_cpp: 基函数、多项式、RBF、JSON；偏航-路径对齐目标（与 8 变量布局一致）
set(SHARED_SOURCES
    ${TRAJOPT_CPP_DIR}/src/chebyshev_basis.cpp
    ${TRAJOPT_CPP_DIR}/src/polynomial_basis.cpp
    ${TRAJOPT_CPP_DIR}/src/polynomial.cpp
    ${TRAJOPT_CPP_DIR}/src/rbf_terrain.cpp
    ${TRAJOPT_CPP_DIR}/src/json_reader.cpp
    ${TRAJOPT_CPP_DIR}/src/objectives/yaw_path_alignment_objective.cpp
)

set(SEG_SOURCES
    src/seg_trajopt_nlp.cpp
    src/segmentation_reader.cpp
    src/constraints/boundary_constraint_8var.cpp
    src/constraints/wheel_base_distance_terrain_constraint.cpp
    src/constraints/wheels_distance_terrain_constraint.cpp
    src/constraints/trajectory_bounds_constraint_8var.cpp
    src/constraints/rollable_region_constraint.cpp
    src/constraints/coefficient_norm_constraint.cpp
    src/objectives/base_height_stability_objective.cpp
    src/objectives/initial_guess_proximity_objective.cpp
    src/waypoint_reader.cpp
    main.cpp
)

add_executable(seg_trajopt_cpp ${SHARED_SOURCES} ${SEG_SOURCES})

if(Ipopt_FOUND)
    target_link_libraries(seg_trajopt_cpp ${IPOPT_LIBRARIES})
else()
    target_link_libraries(seg_trajopt_cpp ${IPOPT_LIBRARIES})
    target_include_directories(seg_trajopt_cpp PRIVATE ${IPOPT_INCLUDE_DIRS})
    target_link_directories(seg_trajopt_cpp PRIVATE ${IPOPT_LIBRARY_DIRS})
endif()
target_link_libraries(seg_trajopt_cpp m)

target_compile_options(seg_trajopt_cpp PRIVATE -Wall -Wextra -O3)
target_compile_definitions(seg_trajopt_cpp PRIVATE HAVE_STDDEF_H)
