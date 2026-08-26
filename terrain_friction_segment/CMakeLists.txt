cmake_minimum_required(VERSION 3.18)
project(terrain_friction_segment LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)

find_package(CUDAToolkit REQUIRED)

# CUDA 架构（可按机器修改）
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES "70;80;86" CACHE STRING "CUDA architectures")
endif()

# 包含目录
include_directories(${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/external)

# 主机源：JSON 读取 + CUDA 调用封装 + main
set(HOST_SOURCES
  src/main.cpp
  src/json_reader.cpp
  src/rbf_friction_cuda.cpp
)

# CUDA 源
set(CUDA_SOURCES
  cuda/rbf_friction_kernel.cu
)

# 可执行文件（CUDA 与 C++ 混合）
add_executable(terrain_friction_segment
  ${HOST_SOURCES}
  ${CUDA_SOURCES}
)

# 链接 CUDA 运行时
target_link_libraries(terrain_friction_segment
  CUDA::cudart_static
)

target_compile_options(terrain_friction_segment PRIVATE
  $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -O3>
  $<$<COMPILE_LANGUAGE:CUDA>:-O3 --use_fast_math>
)

# --- RRT* with region cost (C++ only, no CUDA) ---
add_executable(rrt_star_region
  src/rrt_star_main.cpp
  src/rrt_star_region.cpp
)
target_include_directories(rrt_star_region PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/external)
target_compile_options(rrt_star_region PRIVATE -Wall -Wextra -O3)
