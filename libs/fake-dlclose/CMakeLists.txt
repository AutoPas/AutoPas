add_library(fakedlclose SHARED
        fake-dlclose.c)

SET_TARGET_PROPERTIES(fakedlclose PROPERTIES LINKER_LANGUAGE C)
