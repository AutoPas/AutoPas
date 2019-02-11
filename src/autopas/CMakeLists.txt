file(GLOB_RECURSE MY_SRC
        "*.cpp"
        "*.h"
        )

file(GLOB_RECURSE CU_SRC
        "*.cu"
        "*.cuh"
        )
        
ADD_LIBRARY(autopas STATIC
        ${MY_SRC}
        $<$<BOOL:${ENABLE_CUDA}>:${CU_SRC}>
        )

target_link_libraries(autopas
        PUBLIC
        rt # required for Time.h
        ${CMAKE_THREAD_LIBS_INIT} # required for Logger and ExceptionHandler
        # Do not link against openmp when using archer
        $<$<AND:$<BOOL:${OPENMP}>,$<NOT:$<BOOL:${ARCHER}>>>:OpenMP::OpenMP_CXX>
        )

# Ompstuff needs to be here because OpenMP.cmake needs to run before this file to create to OpenMP target.
# this can be resolved by upgrading to CMake 3.13 and enforcing CMP0079.
target_compile_definitions(autopas
        PUBLIC
        $<$<BOOL:${OPENMP}>:AUTOPAS_OPENMP>
        $<$<BOOL:${ENABLE_CUDA}>:AUTOPAS_CUDA>
        )
        
if(ENABLE_CUDA)
target_compile_options(autopas 
        PUBLIC 
        $<$<COMPILE_LANGUAGE:CUDA>:-gencode arch=compute_61,code=sm_61>
        )
        
TARGET_INCLUDE_DIRECTORIES(autopas SYSTEM PUBLIC ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
endif()

TARGET_INCLUDE_DIRECTORIES(autopas PUBLIC ../.)
