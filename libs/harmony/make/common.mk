ifndef TO_BASE
    $(error Makefile does not define TO_BASE variables.)
endif

#
# Variables that control the build system.
#
.DEFAULT_GOAL=all

prefix=$(realpath .)/$(TO_BASE)
exec_prefix=$(prefix)
bindir=$(exec_prefix)/bin
libdir=$(exec_prefix)/lib
libexecdir=$(exec_prefix)/libexec
includedir=$(prefix)/include
TGTS=$(BIN_TGTS) $(LIB_TGTS) $(INCLUDE_TGTS) $(LIBEXEC_TGTS) $(NO_INST_TGTS)

CC_DEFS=$(shell eval "echo | $(CC) -dM -E - 2>&1; $(CC) -dM 2>&1")
ifeq (__ICC, $(findstring __ICC, $(CC_DEFS)))
    C_COMPILER_MAKE=intel
else ifeq (__PGI, $(findstring __PGI, $(CC_DEFS)))
    C_COMPILER_MAKE=pgi
else ifeq (__clang__, $(findstring __clang__, $(CC_DEFS)))
    C_COMPILER_MAKE=clang
else
    C_COMPILER_MAKE=gcc
endif

ifeq (__linux__, $(findstring __linux__, $(CC_DEFS)))
    OS_MAKE=linux
else ifeq (__APPLE__, $(findstring __APPLE__, $(CC_DEFS)))
    OS_MAKE=darwin
endif

CXX_DEFS=$(shell eval "echo | $(CXX) -dM -E - 2>&1; $(CXX) -dM 2>&1")
ifeq (__ICC, $(findstring __ICC, $(CXX_DEFS)))
    CXX_COMPILER_MAKE+=intel
else ifeq (__PGI, $(findstring __PGI, $(CXX_DEFS)))
    CXX_COMPILER_MAKE+=pgi
else ifeq (__clang__, $(findstring __clang__, $(CXX_DEFS)))
    CXX_COMPILER_MAKE+=clang
else
    CXX_COMPILER_MAKE+=gcc
endif

COMPILER_INCLUDES=$(sort $(C_COMPILER_MAKE) $(CXX_COMPILER_MAKE))
include $(COMPILER_INCLUDES:%=$(TO_BASE)/make/%.mk)
include $(TO_BASE)/make/$(OS_MAKE).mk

#
# Standard rules to make available in all subsystems.
#
.PHONY: all \
        clean \
        distclean \
        install \
        subdirs \
        $(SUBDIRS) \
        $(TO_BASE)/src/libharmony.a

all: $(TGTS) subdirs

_BIN_DIRS=$(addprefix $(DESTDIR)$(bindir)/, $(dir $(BIN_TGTS))) \
          $(addprefix $(DESTDIR)$(bindir)/, $(dir $(BIN_COPY)))
_LIB_DIRS=$(addprefix $(DESTDIR)$(libdir)/, $(dir $(LIB_TGTS))) \
          $(addprefix $(DESTDIR)$(libdir)/, $(dir $(LIB_COPY)))
_LIBEXEC_DIRS=$(addprefix $(DESTDIR)$(libexecdir)/, $(dir $(LIBEXEC_TGTS))) \
              $(addprefix $(DESTDIR)$(libexecdir)/, $(dir $(LIBEXEC_COPY)))
_INCLUDE_DIRS=$(addprefix $(DESTDIR)$(includedir)/, $(dir $(INCLUDE_TGTS))) \
              $(addprefix $(DESTDIR)$(includedir)/, $(dir $(INCLUDE_COPY)))
install: all
	@if [ -n "$(BIN_TGTS)" -o -n "$(BIN_COPY)" ]; then \
	    echo mkdir -p $(sort $(_BIN_DIRS)) && \
	         mkdir -p $(sort $(_BIN_DIRS)) && \
	    for i in $(BIN_TGTS); do \
	        echo cp $$i $(DESTDIR)$(bindir)/$$i && \
	             cp $$i $(DESTDIR)$(bindir)/$$i; \
	    done && \
	    for i in $(BIN_COPY); do \
	        echo cp $$i $(DESTDIR)$(bindir)/$$i && \
	             cp $$i $(DESTDIR)$(bindir)/$$i; \
	    done; \
	fi

	@if [ -n "$(LIB_TGTS)" -o -n "$(LIB_COPY)" ]; then \
	    echo mkdir -p $(sort $(_LIB_DIRS)) && \
	         mkdir -p $(sort $(_LIB_DIRS)) && \
	    for i in $(LIB_TGTS); do \
	        echo cp $$i $(DESTDIR)$(libdir)/$$i && \
	             cp $$i $(DESTDIR)$(libdir)/$$i; \
	    done && \
	    for i in $(LIB_COPY); do \
	        echo cp $$i $(DESTDIR)$(libdir)/$$i && \
	             cp $$i $(DESTDIR)$(libdir)/$$i; \
	    done; \
	fi

	@if [ -n "$(LIBEXEC_TGTS)" -o -n "$(LIBEXEC_COPY)" ]; then \
	echo mkdir -p $(sort $(_LIBEXEC_DIRS)) && \
	     mkdir -p $(sort $(_LIBEXEC_DIRS)) && \
	    for i in $(LIBEXEC_TGTS); do \
	        echo cp $$i $(DESTDIR)$(libexecdir)/$$i && \
	             cp $$i $(DESTDIR)$(libexecdir)/$$i; \
	    done && \
	    for i in $(LIBEXEC_COPY); do \
	        echo cp $$i $(DESTDIR)$(libexecdir)/$$i && \
	             cp $$i $(DESTDIR)$(libexecdir)/$$i; \
	    done; \
	fi

	@if [ -n "$(INCLUDE_TGTS)" -o -n "$(INCLUDE_COPY)" ]; then \
	    echo mkdir -p $(sort $(_INCLUDE_DIRS)) && \
	         mkdir -p $(sort $(_INCLUDE_DIRS)) && \
	    for i in $(INCLUDE_TGTS); do \
	        echo cp $$i $(DESTDIR)$(includedir)/$$i && \
	             cp $$i $(DESTDIR)$(includedir)/$$i; \
	    done && \
	    for i in $(INCLUDE_COPY); do \
	        echo cp $$i $(DESTDIR)$(includedir)/$$i && \
	             cp $$i $(DESTDIR)$(includedir)/$$i; \
	    done; \
	fi

	@if [ -n "$(NO_INST_TGTS)" ]; then \
	    for i in $(NO_INST_TGTS); do \
	        echo "$$i must be executed from source directory (for now)."; \
	    done; \
	fi

clean: subdirs
	$(RM) core a.out *.o $(TGTS) $(EXTRA_CLEANUP)

distclean: clean
	$(RM) *~ *.d

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

$(TO_BASE)/src/libharmony.a:
	$(MAKE) -C $(TO_BASE)/src libharmony.a

#
# Standard pattern rules to make available in all subsystems.
#
_FLAGS_CC=$(REQ_CFLAGS) $(CFLAGS)
_FLAGS_CPP=$(REQ_CPPFLAGS) $(CPPFLAGS)
_FLAGS_CXX=$(REQ_CXXFLAGS) $(CXXFLAGS)
_FLAGS_FC=$(REQ_FFLAGS) $(FFLAGS)
_FLAGS_LD=$(REQ_LDFLAGS) $(LDFLAGS)
_FLAGS_LIBS=$(REQ_LDLIBS) $(LOADLIBES) $(LDLIBS)
_FLAGS_AR=$(REQ_ARFLAGS) $(ARFLAGS)

%.o: %.c
	$(CC) -c $(_FLAGS_CPP) $(_FLAGS_CC) $< -o $@

%.o: %.cxx
	$(CXX) -c $(_FLAGS_CPP) $(_FLAGS_CXX) $< -o $@

%.o: %.f
	$(FC) -c $(_FLAGS_FC) $< -o $@

%.o: %.F
	$(FC) -c $(_FLAGS_CPP) $(_FLAGS_FC) $< -o $@

%.so: %.o
	$(CC) -shared $(_FLAGS_CC) $(_FLAGS_LD) $^ $(_FLAGS_LIBS) -o $@

%.a:
	$(AR) $(_FLAGS_AR) $@ $^

%: %.o
	$(CC) $(_FLAGS_LD) $^ $(_FLAGS_LIBS) -o $@

%: %.c
	$(CC) $(_FLAGS_CPP) $(_FLAGS_CC) $(_FLAGS_LD) $^ $(_FLAGS_LIBS) -o $@

%: %.f
	$(FC) $(_FLAGS_FC) $(_FLAGS_LD) $^ $(_FLAGS_LDLIBS) -o $@

%: %.F
	$(FC) $(_FLAGS_CPP) $(_FLAGS_FC) $(_FLAGS_LD) $^ $(_FLAGS_LIBS) -o $@

#
# Auto dependency creation
#
%.d: %.c
	@$(RM) $@; \
	$(CC) -MM $(_FLAGS_CPP) $< > $@.$$$$ 2>/dev/null; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	$(RM) $@.$$$$
-include $(patsubst %.c,%.d,$(filter %.c,$(SRCS)))

%.d: %.cxx
	@$(RM) $@; \
	$(CXX) -MM $(_FLAGS_CPP) $< > $@.$$$$ 2>/dev/null; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	$(RM) $@.$$$$
-include $(patsubst %.cxx,%.d,$(filter %.cxx,$(SRCS)))
