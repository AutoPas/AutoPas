
if(ENABLE_CUDA)
set(CUDA_SEPARABLE_COMPILATION ON)

#Get CUDA compute capability
set(CUDAPROGRAM ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/detectCudaComputeCapabilty)
set(CUDAFILE ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/detectCudaComputeCapabilty.cu)
execute_process(COMMAND nvcc -lcuda ${CUDAFILE} -o ${CUDAPROGRAM})
execute_process(COMMAND ${CUDAPROGRAM}
                RESULT_VARIABLE CUDA_RETURN_CODE
                OUTPUT_VARIABLE DETECTED_CCC)
                
if (NOT (${CUDA_RETURN_CODE} EQUAL 0))
    message(WARNING "Could not detect cuda compute capabilty!")
endif()
set(CUDA_COMPUTE_CAPABILITY ${DETECTED_CCC} CACHE STRING "Cuda compute capability of yout GPU")

file(REMOVE ${CUDAPROGRAM})
target_compile_options(autopas
        PUBLIC
        #architecture flags and -Xcompiler to prepend to -fopenmp
        $<$<COMPILE_LANGUAGE:CUDA>:$<$<STREQUAL:${CMAKE_BUILD_TYPE},Debug>:-lineinfo -pg> -gencode arch=compute_${CUDA_COMPUTE_CAPABILITY},code=sm_${CUDA_COMPUTE_CAPABILITY} $<$<BOOL:${OPENMP}>:-Xcompiler>>
        )
endif()
