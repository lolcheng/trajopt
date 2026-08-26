cmake_minimum_required(VERSION 3.10)
project(trajopt_cpp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 查找Ipopt
# 方法1: 尝试使用CMake的find_package
find_package(Ipopt QUIET)

# 方法2: 如果找不到，尝试使用pkg-config
if(NOT Ipopt_FOUND)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(IPOPT ipopt)
    if(IPOPT_FOUND)
        set(IPOPT_INC_DIR ${IPOPT_INCLUDE_DIRS})
        set(IPOPT_LIBRARIES ${IPOPT_LIBRARIES})
    endif()
endif()

# 方法3: 如果还是找不到，尝试手动查找
if(NOT Ipopt_FOUND AND NOT IPOPT_FOUND)
    # 常见的安装路径
    set(IPOPT_POSSIBLE_INC_DIRS
        /usr/include/coin-or
        /usr/local/include/coin-or
        /usr/include
        /usr/local/include
    )
    set(IPOPT_POSSIBLE_LIB_DIRS
        /usr/lib
        /usr/local/lib
        /usr/lib/x86_64-linux-gnu
    )
    
    find_path(IPOPT_INC_DIR
        NAMES IpTNLP.hpp
        PATHS ${IPOPT_POSSIBLE_INC_DIRS}
    )
    
    find_library(IPOPT_LIBRARIES
        NAMES ipopt
        PATHS ${IPOPT_POSSIBLE_LIB_DIRS}
    )
    
    if(IPOPT_INC_DIR AND IPOPT_LIBRARIES)
        set(IPOPT_FOUND TRUE)
        message(STATUS "Found Ipopt manually: ${IPOPT_INC_DIR}, ${IPOPT_LIBRARIES}")
    endif()
endif()

# 检查是否找到Ipopt
if(NOT IPOPT_FOUND AND NOT Ipopt_FOUND)
    message(FATAL_ERROR "Ipopt not found! Please install it using: sudo apt install coinor-libipopt-dev")
endif()

# 包含目录
include_directories(${CMAKE_SOURCE_DIR}/include)
if(IPOPT_INC_DIR)
    include_directories(${IPOPT_INC_DIR})
elseif(Ipopt_FOUND)
    include_directories(${IPOPT_INC_DIR})
endif()

# 源文件
set(SOURCES
    src/rbf_terrain.cpp
    src/polynomial.cpp
    src/polynomial_basis.cpp
    src/chebyshev_basis.cpp
    src/constraints/boundary_constraint.cpp
    src/constraints/wheel_base_distance_constraint.cpp
    src/constraints/wheels_distance_constraint.cpp
    src/constraints/anti_crossing_constraint.cpp
    src/constraints/trajectory_bounds_constraint.cpp
    src/constraints/normal_force_upper_bound_constraint.cpp
    src/constraints/wheel_airborne_constraint.cpp
    src/trajopt_nlp.cpp
    src/json_reader.cpp
    src/terrain_visualizer.cpp
    src/trajectory_visualizer.cpp
    src/force_data_saver.cpp
    src/objectives/wheel_terrain_distance_objective.cpp
    src/objectives/wheel_airborne_penalty_objective.cpp
    src/objectives/path_length_objective.cpp
    src/objectives/forward_progress_objective.cpp
    src/objectives/regularization_objective.cpp
    src/objectives/smoothing_objective.cpp
    src/objectives/wheel_base_distance_nominal_objective.cpp
    src/objectives/wheels_distance_nominal_objective.cpp
    src/objectives/yaw_path_alignment_objective.cpp
    src/objectives/contact_stiffness_objective.cpp
    src/objectives/force_balance_objective.cpp
    src/objectives/torque_balance_objective.cpp
    src/objectives/friction_cone_objective.cpp
)

# 头文件
set(HEADERS
    include/rbf_terrain.hpp
    include/polynomial.hpp
    include/basis_function.hpp
    include/polynomial_basis.hpp
    include/chebyshev_basis.hpp
    include/constraint_base.hpp
    include/objective_base.hpp
    include/constraints/boundary_constraint.hpp
    include/constraints/wheel_base_distance_constraint.hpp
    include/constraints/wheels_distance_constraint.hpp
    include/constraints/anti_crossing_constraint.hpp
    include/constraints/trajectory_bounds_constraint.hpp
    include/constraints/normal_force_upper_bound_constraint.hpp
    include/constraints/wheel_airborne_constraint.hpp
    include/objectives/wheel_terrain_distance_objective.hpp
    include/objectives/forward_progress_objective.hpp
    include/objectives/regularization_objective.hpp
    include/objectives/smoothing_objective.hpp
    include/objectives/wheel_base_distance_nominal_objective.hpp
    include/objectives/wheels_distance_nominal_objective.hpp
    include/objectives/yaw_path_alignment_objective.hpp
    include/objectives/contact_stiffness_objective.hpp
    include/objectives/force_balance_objective.hpp
    include/objectives/torque_balance_objective.hpp
    include/objectives/friction_cone_objective.hpp
    include/trajopt_nlp.hpp
    include/json_reader.hpp
    include/terrain_visualizer.hpp
    include/trajectory_visualizer.hpp
)

# 可执行文件
add_executable(trajopt_cpp
    main.cpp
    ${SOURCES}
)

# 链接库
if(Ipopt_FOUND)
    target_link_libraries(trajopt_cpp ${IPOPT_LIBRARIES})
elseif(IPOPT_FOUND)
    target_link_libraries(trajopt_cpp ${IPOPT_LIBRARIES})
    target_include_directories(trajopt_cpp PRIVATE ${IPOPT_INCLUDE_DIRS})
    target_link_directories(trajopt_cpp PRIVATE ${IPOPT_LIBRARY_DIRS})
endif()

# 链接数学库（Ipopt可能需要）
target_link_libraries(trajopt_cpp m)

# 编译选项
target_compile_options(trajopt_cpp PRIVATE
    -Wall
    -Wextra
    -O3
)

# 定义Ipopt需要的宏
target_compile_definitions(trajopt_cpp PRIVATE
    HAVE_STDDEF_H
)

