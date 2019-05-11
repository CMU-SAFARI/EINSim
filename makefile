#!/usr/bin/env make
BIN = einsim
BIN_DBG = $(BIN).d
BUILD_DIR = build
SRC_DIR = src
SOURCES := $(wildcard $(SRC_DIR)/*.cpp $(SRC_DIR)/*/*.cpp)
SUBDIRS := $(sort $(dir ${SOURCES}))
OBJS = $(SOURCES:%=$(BUILD_DIR)/%.o)
OBJS_DBG = $(SOURCES:%=$(BUILD_DIR)/%.do)
DEPS = $(OBJS:%=%.d)
DEPS_DBG = $(OBJS_DBG:%=%.d)
CC = g++
LD = g++
INCLUDE_DIRS = $(SRC_DIR) lib/eigen_3.3.7 lib/libtp lib/cxxopts
CFLAGS_OPT = -g -Wfatal-errors -Werror -Wall -O3 $(INCLUDE_DIRS:%=-I%) -std=c++11 -pthread
CFLAGS_DBG = -g -Wfatal-errors -Werror -Wall -O0 $(INCLUDE_DIRS:%=-I%) -std=c++11 -pthread
LDFLAGS = -pthread
LDFLAGS_DBG = -pthread

# $(info SOURCES : $(SOURCES))
# $(info OBJS    : $(OBJS))
# $(info DEPS    : $(DEPS))
# $(info SOURCES : $(SOURCES))

.SUFFIXES:

.PHONY: default all release debug clean
default: release

release: $(BIN)

debug: $(BIN_DBG)

all: release debug

-include $(DEPS) 
-include $(DEPS_DBG)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)
	mkdir -p $(SUBDIRS:%=$(BUILD_DIR)/%)

$(BUILD_DIR)/%.cpp.do: %.cpp | $(BUILD_DIR)
	$(CC) $(CFLAGS_DBG) -c -MMD -MF $@.d -o $@ $<

$(BUILD_DIR)/%.cpp.o: %.cpp | $(BUILD_DIR)
	$(CC) $(CFLAGS_OPT) -c -MMD -MF $@.d -o $@ $<

$(BIN): $(OBJS)
	$(LD) $(LDFLAGS) $^ -o $@

$(BIN_DBG): $(OBJS_DBG)
	$(LD) $(LDFLAGS_DBG) $^ -o $@

clean:
	rm -f $(BIN)
	rm -f $(BIN_DBG)
	rm -rf $(BUILD_DIR)
