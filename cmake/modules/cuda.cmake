
if(ENABLE_CUDA)
target_compile_options(autopas
        PUBLIC
        #architecture flags and -Xcompiler to prepend to -fopenmp
        $<$<COMPILE_LANGUAGE:CUDA>:-gencode arch=compute_61,code=sm_61 $<$<BOOL:${OPENMP}>:-Xcompiler>>
        )
endif()
