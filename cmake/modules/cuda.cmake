if (AUTOPAS_ENABLE_CUDA)
    set(CUDA_SEPARABLE_COMPILATION ON)

    # set auto-detect as the default CUDA_COMPUTE_CAPABILITY if it is not yet set.
    if (NOT CUDA_COMPUTE_CAPABILITY)
        set(
            CUDA_COMPUTE_CAPABILITY
            "auto-detect"
            CACHE
                STRING "Choose the cuda compute capability, options are: auto-detect 60 61 70 75."
                FORCE
        )
    endif ()

    set_property(CACHE CUDA_COMPUTE_CAPABILITY PROPERTY STRINGS "auto-detect;60;61;70;75")

    if (CUDA_COMPUTE_CAPABILITY STREQUAL "auto-detect")
        # Get CUDA compute capability
        set(CUDAPROGRAM ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/detectCudaComputeCapabilty)
        set(CUDAFILE ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/detectCudaComputeCapabilty.cu)
        execute_process(
            COMMAND
                nvcc
                -lcuda
                ${CUDAFILE}
                -o
                ${CUDAPROGRAM}
        )
        execute_process(
            COMMAND ${CUDAPROGRAM}
            RESULT_VARIABLE CUDA_RETURN_CODE
            OUTPUT_VARIABLE DETECTED_CCC
        )

        if (NOT (${CUDA_RETURN_CODE} EQUAL 0))
            message(WARNING "Could not detect cuda compute capabilty! Please Set Manually.")
        endif ()
        set(CUDA_COMPUTE_CAPABILITY ${DETECTED_CCC})
        file(REMOVE ${CUDAPROGRAM})
    endif ()

    target_compile_options(
        autopas
        PUBLIC
            # architecture flags and -Xcompiler to prepend to -fopenmp
            $<$<COMPILE_LANGUAGE:CUDA>:$<$<STREQUAL:${CMAKE_BUILD_TYPE},Debug>:-lineinfo
            -pg>
            -gencode
            arch=compute_${CUDA_COMPUTE_CAPABILITY},code=sm_${CUDA_COMPUTE_CAPABILITY}
            $<$<BOOL:${OPENMP}>:-Xcompiler>>
    )
endif ()
