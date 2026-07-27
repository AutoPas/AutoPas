STRICT_FLAGS=$(strip -Weverything -Werror \
                     -Wno-cast-qual \
                     -Wno-disabled-macro-expansion \
                     -Wno-documentation \
                     -Wno-documentation-unknown-command \
                     -Wno-exit-time-destructors \
                     -Wno-float-equal \
                     -Wno-format-nonliteral \
                     -Wno-global-constructors \
                     -Wno-padded \
                     -Wno-missing-field-initializers \
                     -Wno-missing-prototypes \
                     -Wno-missing-noreturn \
                     -Wno-missing-variable-declarations \
                     -Wno-reserved-id-macro \
                     -Wno-shorten-64-to-32 \
                     -Wno-sign-compare \
                     -Wno-sign-conversion \
                     -Wno-switch-enum \
                     -Wno-unknown-attributes \
                     -Wno-unknown-warning-option \
                     -Wno-unused-parameter \
              )

ifeq ($(C_COMPILER_MAKE), clang)
    REQ_CFLAGS+=-std=c99

    ifeq ($(DEBUG), 1)
        REQ_CFLAGS+=-g
    else
        REQ_CFLAGS+=-O2
    endif

    ifeq ($(STRICT), 1)
        REQ_CFLAGS+=$(STRICT_FLAGS)
    else
        REQ_CFLAGS+=-Wno-unknown-attributes
    endif
endif

ifeq ($(CXX_COMPILER_MAKE), clang)
    ifeq ($(DEBUG), 1)
        REQ_CXXFLAGS+=-g
    else
        REQ_CXXFLAGS+=-O2
    endif

    ifeq ($(STRICT), 1)
        REQ_CXXFLAGS+=$(STRICT_FLAGS)
    else
        REQ_CXXFLAGS+=-Wno-unknown-attributes
    endif
endif
