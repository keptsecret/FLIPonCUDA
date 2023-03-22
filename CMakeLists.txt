cmake_minimum_required(VERSION 3.17)
project(FLIPonCUDA LANGUAGES CXX CUDA)

include(CheckLanguage)
check_language(CUDA)

if (MSVC)
    list (APPEND FOC_DEFINITIONS "FOC_IS_MSVC" "_CRT_SECURE_NO_WARNINGS")
    list (APPEND FOC_DEFINITIONS "FOC_IS_MSVC" "_ENABLE_EXTENDED_ALIGNED_STORAGE")
endif ()

option(CUDA_ENABLED "Enable GPU acceleration" OFF)
if(CUDA_ENABLED)
    add_definitions(-DFOC_BUILD_GPU)
endif()

set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

###########################################################################
# CUDA

find_package(CUDAToolkit REQUIRED)

include_directories (${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

add_library (foc_warnings INTERFACE)
target_compile_options (
        foc_warnings
        INTERFACE
        #"$<$<COMPILE_LANGUAGE:CUDA>:SHELL:-Xptxas --warn-on-double-precision-use>"
        "$<$<COMPILE_LANGUAGE:CUDA>:SHELL:-Xcudafe --diag_suppress=partial_override>"
        "$<$<COMPILE_LANGUAGE:CUDA>:SHELL:-Xcudafe --diag_suppress=virtual_function_decl_hidden>"
        "$<$<COMPILE_LANGUAGE:CUDA>:SHELL:-Xcudafe --diag_suppress=integer_sign_change>"
        "$<$<COMPILE_LANGUAGE:CUDA>:SHELL:-Xcudafe --diag_suppress=declared_but_not_referenced>"
        # WAR invalid warnings about this with "if constexpr"
        "$<$<COMPILE_LANGUAGE:CUDA>:SHELL:-Xcudafe --diag_suppress=implicit_return_from_non_void_function>"
)

string (APPEND CMAKE_CUDA_FLAGS " -Xnvlink -suppress-stack-size-warning")

add_library (cuda_build_configuration INTERFACE)
target_compile_options (
        cuda_build_configuration
        INTERFACE
        "$<$<COMPILE_LANGUAGE:CUDA>:--std=c++17;--use_fast_math;--expt-relaxed-constexpr;--extended-lambda;--forward-unknown-to-host-compiler>"
        "$<$<COMPILE_LANGUAGE:CUDA>:$<IF:$<AND:$<CONFIG:Debug>,$<NOT:$<BOOL:$<TARGET_PROPERTY:CUDA_PTX_COMPILATION>>>>,-G;-g,-lineinfo;-maxrregcount;128>>"
)

###########################################################################

set( FOC_HEADERS
        src/foc.h
        src/vecmath.h
        src/array3d.h
        src/markerparticle.h
        src/macgrid.h
        src/cellmaterialgrid.h
        src/scalarfield.h
        src/fluidsim.h
        )

set( FOC_SOURCE
        src/macgrid.cpp
        src/cellmaterialgrid.cpp
        src/scalarfield.cpp
        src/fluidsim.cpp
        )

if (CUDA_ENABLED)
    set_source_files_properties(
            ${FOC_HEADERS}
            ${FOC_SOURCE}

            PROPERTIES LANGUAGE CUDA
    )
endif()

add_library( flipfluidsim_lib STATIC
        ${FOC_HEADERS}
        ${FOC_SOURCE}
        )

target_compile_definitions ( flipfluidsim_lib PRIVATE ${FOC_DEFINITIONS} )

set_property ( TARGET flipfluidsim_lib PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON )
set_property ( TARGET flipfluidsim_lib PROPERTY CUDA_ARCHITECTURES all-major )

target_include_directories( flipfluidsim_lib PUBLIC
        src
        )

target_compile_options ( flipfluidsim_lib PUBLIC ${FOC_CXX_FLAGS} )
target_link_libraries( flipfluidsim_lib PRIVATE foc_warnings $<$<BOOL:CUDA_ENABLED>:cuda_build_configuration>)

if (WIN32)
    # Avoid a name clash when building on Visual Studio
    set_target_properties ( flipfluidsim_lib PROPERTIES OUTPUT_NAME libfoc )
endif()

set ( ALL_FOC_LIBS
        flipfluidsim_lib
        )

if (CUDA_ENABLED)
    set_target_properties( flipfluidsim_lib PROPERTIES
            CUDA_SEPARABLE_COMPILATION ON)
endif()

add_executable( flipfluidsim_exe src/foc.cu )

target_compile_definitions ( flipfluidsim_exe PRIVATE ${FOC_DEFINITIONS} )
target_compile_options ( flipfluidsim_exe PRIVATE ${FOC_CXX_FLAGS} )
target_include_directories ( flipfluidsim_exe PRIVATE src )
target_link_libraries ( flipfluidsim_exe ${ALL_FOC_LIBS} foc_warnings )

set_target_properties( flipfluidsim_exe PROPERTIES
        OUTPUT_NAME foc)

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
