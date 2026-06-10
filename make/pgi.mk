ifeq ($(C_COMPILER_MAKE), pgi)
    REQ_CFLAGS+=-c99

    ifeq ($(DEBUG), 1)
        REQ_CFLAGS+=-g -Minform=severe
    else
        REQ_CFLAGS+=-fast -Minform=fatal
    endif

    ifeq ($(STRICT), 1)
        REQ_CFLAGS+=-Minform=warn
    endif
endif

ifeq ($(CXX_COMPILER_MAKE), pgi)
    ifeq ($(DEBUG), 1)
        REQ_CXXFLAGS+=-g -Minform=severe
    else
        REQ_CXXFLAGS+=-fast -Minform=fatal
    endif

    ifeq ($(STRICT), 1)
        REQ_CXXFLAGS+=-Minform=warn
    endif
endif
