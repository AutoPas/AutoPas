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
                ${CMAKE_CUDA_COMPILER}
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

    # cuda as of 10.2 does not support c++17, so set this to 14 here. If cmake < 3.17 is used this
    # is set globally in the top level CMakeLists.txt
    if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.17)
        target_compile_features(autopas PRIVATE cuda_std_14)
    endif ()

    target_compile_options(
        autopas
        PRIVATE
            # specify host compiler
            $<$<COMPILE_LANGUAGE:CUDA>:-ccbin=${CMAKE_CXX_COMPILER}>
            # add debug flags
            $<$<COMPILE_LANGUAGE:CUDA>:$<$<STREQUAL:${CMAKE_BUILD_TYPE},Debug>:-lineinfo -pg>>
            # notify nvcc of compute capability
            $<$<COMPILE_LANGUAGE:CUDA>:-gencode
            arch=compute_${CUDA_COMPUTE_CAPABILITY},code=sm_${CUDA_COMPUTE_CAPABILITY}>
    )

endif ()
