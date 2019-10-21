TO_BASE=../../..

REQ_CFLAGS+=-fPIC
REQ_CPPFLAGS+=-I$(TO_BASE)/src
REQ_LDFLAGS+=$(SHAREDOBJ_FLAG)

SRCS=angel.c \
     exhaustive.c \
     libvertex.c \
     nm.c \
     pro.c \
     random.c

LIBEXEC_TGTS=angel.so \
             exhaustive.so \
             nm.so \
             pro.so \
             random.so

angel.so: REQ_LDLIBS+=-lm
angel.so: libvertex.o

exhaustive.so: REQ_LDLIBS+=-lm

nm.so: REQ_LDLIBS+=-lm
nm.so: libvertex.o

pro.so: REQ_LDLIBS+=-lm
pro.so: libvertex.o

# Active Harmony makefiles should always include this file last.
include $(TO_BASE)/make/common.mk
