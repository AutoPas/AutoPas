TO_BASE=.

SUBDIRS=src \
	code-server \
	example

.PHONY: doc

code-server: src

example: src

doc:
	$(MAKE) -C doc

# Active Harmony makefiles should always include this file last.
include $(TO_BASE)/make/common.mk
