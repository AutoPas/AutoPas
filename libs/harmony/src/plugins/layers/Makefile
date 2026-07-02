TO_BASE=../../..

SRCS=agg.c \
     cache.c \
     codegen.c \
     codegen-helper.c \
     constraint.c \
     group.c \
     log.c \
     TAUdb.c \
     xmlWriter.c

LIBEXEC_TGTS=agg.so \
             cache.so \
             codegen.so \
             codegen-helper \
             constraint.so \
             group.so \
             log.so
SHARED_OBJS=$(filter %.so, $(LIBEXEC_TGTS))

REQ_CPPFLAGS+=-I$(TO_BASE)/src

LIBTAUDB=
LIBTAUDB_INCDIR=$(LIBTAUDB)/include
LIBTAUDB_LIBDIR=$(LIBTAUDB)/lib
ifneq (, $(wildcard $(LIBTAUDB_LIBDIR)/libtaudb.a))
    LIBEXEC_TGTS += TAUdb.so

    TAUdb.so: REQ_CPPFLAGS+=-I$(LIBTAUDB_INCDIR)
    TAUdb.so: REQ_LDFLAGS+=-L$(LIBTAUDB_LIBDIR)
    TAUdb.so: REQ_LDLIBS+=-ltaudb
endif

LIBXML2=
LIBXML2_INCDIR=$(LIBXML2)/include/libxml2
LIBXML2_LIBDIR=$(LIBXML2)/lib
ifneq (, $(wildcard $(LIBXML2_INCDIR)/libxml/xmlwriter.h))
    LIBEXEC_TGTS += xmlWriter.so

    xmlWriter.so: REQ_CPPFLAGS+=-I${LIBXML2_INCDIR}
    xmlWriter.so: REQ_LDFLAGS+=-L${LIBXML2_LIBDIR}
    xmlWriter.so: REQ_LDLIBS+=-lxml2
endif

$(SHARED_OBJS:%.so=%.o): REQ_CFLAGS+=-fPIC
$(SHARED_OBJS): REQ_LDFLAGS+=$(SHAREDOBJ_FLAG)

codegen-helper: $(TO_BASE)/src/libharmony.a

# Active Harmony makefiles should always include this file last.
include $(TO_BASE)/make/common.mk
