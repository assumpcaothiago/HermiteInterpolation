CC ?= cc
CXX ?= c++
AR ?= ar

CPPFLAGS ?=
CFLAGS ?= -O2
CXXFLAGS ?= -O2

CPPFLAGS += -Iinclude
CFLAGS += -std=c99 -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror
CXXFLAGS += -std=c++11 -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror

BUILD_DIR := build
LIBRARY := $(BUILD_DIR)/libhermite3d.a
LIBRARY_OBJECT := $(BUILD_DIR)/hermite3d.o
C_TEST := $(BUILD_DIR)/test_hermite3d
CXX_SMOKE_TEST := $(BUILD_DIR)/header_smoke

.PHONY: all check clean

all: $(LIBRARY)

$(BUILD_DIR):
	mkdir -p $@

$(LIBRARY_OBJECT): src/hermite3d.c include/hermite3d.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(LIBRARY): $(LIBRARY_OBJECT)
	$(AR) rcs $@ $^

$(C_TEST): tests/test_hermite3d.c include/hermite3d.h $(LIBRARY) | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LIBRARY) -lm -o $@

$(CXX_SMOKE_TEST): tests/header_smoke.cpp include/hermite3d.h $(LIBRARY) | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIBRARY) -lm -o $@

check: $(C_TEST) $(CXX_SMOKE_TEST)
	$(C_TEST)
	$(CXX_SMOKE_TEST)

clean:
	rm -rf $(BUILD_DIR)
