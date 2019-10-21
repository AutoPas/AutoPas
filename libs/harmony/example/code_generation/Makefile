TO_BASE=../..

MPICC=mpicc
export OMPI_CC:=$(CC)
export MPICH_CC:=$(CC)

ifeq (0, $(shell which $(firstword $(MPICC)) > /dev/null 2>&1; echo $$?))
    SRCS+=gemm.c
    NO_INST_TGTS+=gemm

    gemm: override CC=$(MPICC)
    gemm: REQ_CPPFLAGS+=-I$(TO_BASE)/src
    gemm: REQ_CFLAGS+=$(MPICFLAGS)
    gemm: REQ_LDLIBS+=-ldl
    gemm: $(TO_BASE)/src/libharmony.a
else
    $(info Skipping code_generation example: $(firstword $(MPICC)) not found.)
endif

# Active Harmony makefiles should always include this file last.
include $(TO_BASE)/make/common.mk
