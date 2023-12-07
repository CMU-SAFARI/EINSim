#!/usr/bin/env make
BIN = einsim
BIN_DBG = $(BIN).d
BIN_LIB = lib$(BIN).so
BUILD_DIR = build
DOXY_DIR = doxygen
SRC_DIR = src
SOURCES := $(wildcard $(SRC_DIR)/*.cpp $(SRC_DIR)/*/*.cpp)
SUBDIRS := $(sort $(dir ${SOURCES}))
OBJS = $(SOURCES:%=$(BUILD_DIR)/%.o)
OBJS_DBG = $(SOURCES:%=$(BUILD_DIR)/%.do)
OBJS_SO = $(SOURCES:%=$(BUILD_DIR)/%.so)
DEPS = $(OBJS:%=%.d)
DEPS_DBG = $(OBJS_DBG:%=%.d)
DEPS_SO = $(OBJS_SO:%=%.d)
CC = g++
LD = g++
INCLUDE_DIRS = $(SRC_DIR) lib lib/eigen_3.3.7 lib/libtp lib/cxxopts lib/rapidjson
CFLAGS_OPT = -g -Wfatal-errors -Werror -Wall -O3 $(INCLUDE_DIRS:%=-I%) -std=c++11 -pthread
CFLAGS_DBG = -ggdb -Wfatal-errors -Werror -Wall -O0 $(INCLUDE_DIRS:%=-I%) -std=c++11 -pthread
CFLAGS_LIB = $(CFLAGS_OPT) -fPIC -DDYNAMIC_LIB
LDFLAGS = -pthread
LDFLAGS_DBG = -pthread
LDFLAGS_LIB = -shared

# $(info SOURCES : $(SOURCES))
# $(info OBJS    : $(OBJS))
# $(info DEPS    : $(DEPS))
# $(info SOURCES : $(SOURCES))

.SUFFIXES:

.PHONY: default jall all release debug clean doc lib
default: release

release: $(BIN)

doc:
	doxygen

debug: $(BIN_DBG)

lib: $(BIN_LIB)

all: release debug lib

jall: 
	$(MAKE) -j 8 all

-include $(DEPS) 
-include $(DEPS_DBG)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)
	mkdir -p $(SUBDIRS:%=$(BUILD_DIR)/%)

$(BUILD_DIR)/%.cpp.o: %.cpp | $(BUILD_DIR)
	$(CC) $(CFLAGS_OPT) -c -MMD -MF $@.d -o $@ $<

$(BUILD_DIR)/%.cpp.do: %.cpp | $(BUILD_DIR)
	$(CC) $(CFLAGS_DBG) -c -MMD -MF $@.d -o $@ $<

$(BUILD_DIR)/%.cpp.so: %.cpp | $(BUILD_DIR)
	$(CC) $(CFLAGS_LIB) -c -MMD -MF $@.d -o $@ $<

$(BIN): $(OBJS)
	$(LD) $(LDFLAGS) $^ -o $@

$(BIN_DBG): $(OBJS_DBG)
	$(LD) $(LDFLAGS_DBG) $^ -o $@

$(BIN_LIB): $(OBJS_SO)
	$(LD) $(LDFLAGS_LIB) $^ -o $@

clean:
	rm -f $(BIN)
	rm -f $(BIN_DBG)
	rm -f $(BIN_LIB)
	rm -rf $(BUILD_DIR)
	rm -rf $(DOXY_DIR)
