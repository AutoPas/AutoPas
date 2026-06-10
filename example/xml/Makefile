TO_BASE=../..

SRCS=example.c
NO_INST_TGTS=example

EXTRA_CLEANUP=*.xml

example: REQ_CFLAGS+=-I$(TO_BASE)/src
example: $(TO_BASE)/src/libharmony.a

# Active Harmony makefiles should always include this file last.
include $(TO_BASE)/make/common.mk
