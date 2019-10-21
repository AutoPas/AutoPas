TO_BASE=../..

SRCS=join.c \
     loadfile.c \
     minimal.c \
     multi.c \
     users_guide.c
NO_INST_TGTS=join \
             loadfile \
             minimal \
             multi \
             users_guide

join: REQ_CPPFLAGS+=-I$(TO_BASE)/src
join: $(TO_BASE)/src/libharmony.a

loadfile: REQ_CPPFLAGS+=-I$(TO_BASE)/src
loadfile: $(TO_BASE)/src/libharmony.a

minimal: REQ_CPPFLAGS+=-I$(TO_BASE)/src
minimal: $(TO_BASE)/src/libharmony.a

multi: REQ_CPPFLAGS+=-I$(TO_BASE)/src
multi: $(TO_BASE)/src/libharmony.a

users_guide: REQ_CPPFLAGS+=-I$(TO_BASE)/src
users_guide: $(TO_BASE)/src/libharmony.a

# ifeq (0, $(shell which $(FC) > /dev/null 2>&1; echo $$?))
#     SRCS+=fortran.f
#     NO_INST_TGTS+=fortran
#
#     fortran: REQ_FFLAGS+=-fno-underscoring
#     fortran: fortran.o $(TO_BASE)/src/libharmony.a
# else
#     $(info Skipping fortran client_api example: $(FC) not found.)
# endif

# CHPLC=chpl
# ifeq (0, $(shell which $(CHPLC) > /dev/null 2>&1; echo $$?))
#     SRCS+=chapel.chpl
#     NO_INST_TGTS+=chapel
#
#     chapel: chapel.chpl $(TO_BASE)/src/libharmony.a
# 	$(CHPLC) $(CHPLFLAGS) $^ $(TO_BASE)/src/hclient.h -o $@
# else
#     $(info Skipping chapel client_api example: $(CHPLC) not found.)
# endif

# Active Harmony makefiles should always include this file last.
include $(TO_BASE)/make/common.mk
