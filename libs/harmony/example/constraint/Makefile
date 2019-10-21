TO_BASE=../..

TEST_DIR=$(TO_BASE)/example/client_api
TEST_BIN=minimal

.PHONY: test client

all: client

test: client
	$(TEST_DIR)/$(TEST_BIN) LAYERS=constraint.so OC_FILE=constraint.in

client:
	$(MAKE) -C $(TEST_DIR) $(TEST_BIN)

# Active Harmony makefiles should always include this file last.
include $(TO_BASE)/make/common.mk
