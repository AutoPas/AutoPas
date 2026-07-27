TO_BASE=.

SUBDIRS=src \
	code-server

.PHONY: doc

code-server: src

doc:
	$(MAKE) -C doc

# Active Harmony makefiles should always include this file last.
include $(TO_BASE)/make/common.mk
