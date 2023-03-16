cmake_minimum_required(VERSION 3.24)
project(FLIPonCUDA CUDA)

include(CheckLanguage)
check_language(CUDA)

option(ENABLE_GPU "Enable GPU acceleration" ON)
if(ENABLE_GPU)
    add_definitions(-DENABLE_GPU)
endif()

set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CUDA_STANDARD 11)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
set(CMAKE_CUDA_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CUDA_ARCHITECTURES 35 50 70 72)
set(CMAKE_CUDA_FLAGS_DEBUG "-G")
set(CMAKE_CUDA_FLAG_RELEASE "-Wno-deprecated-gpu-targets")
set(CMAKE_CUDA_ARCHITECTURES 70 72)

find_package(CUDAToolkit REQUIRED)

set( FLIPSIM_HEADERS
        src/vecmath.h
        )

add_library( flipfluidsim_lib STATIC
        ${FLIPSIM_HEADERS}
        )

add_executable( flipfluidsim_exe src/main.cu)

target_link_libraries ( flipfluidsim_exe flipfluidsim_lib )

set_target_properties( flipfluidsim_exe PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)

set_target_properties( flipfluidsim_exe PROPERTIES
        OUTPUT_NAME flipsim)

install ( TARGETS
        flipfluidsim_exe
        DESTINATION
        bin
        )

install ( TARGETS
        flipfluidsim_lib
        DESTINATION
        lib
        )
