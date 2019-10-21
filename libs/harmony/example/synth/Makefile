TO_BASE=../..

SRCS=synth.c \
     testfunc.c
NO_INST_TGTS=synth

synth: REQ_CPPFLAGS+=-I$(TO_BASE)/src
synth: REQ_LDLIBS+=-lm
synth: $(SRCS:.c=.o) $(TO_BASE)/src/libharmony.a

# Active Harmony makefiles should always include this file last.
include $(TO_BASE)/make/common.mk
