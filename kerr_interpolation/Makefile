CC ?= cc
PYTHON ?= python3

CPPFLAGS ?=
CFLAGS ?= -O2

CPPFLAGS += -Iinclude -Igenerated -I../hermite3d/include
CFLAGS += -std=c99 -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror

BUILD_DIR := build
HERMITE_LIBRARY := ../hermite3d/build/libhermite3d.a
GENERATED_SOURCE := generated/kerr_exact.c
GENERATOR := generate_kerr_exact.py
EXACT_OBJECT := $(BUILD_DIR)/kerr_exact.o
CONVERGENCE_OBJECT := $(BUILD_DIR)/kerr_convergence.o
CONVERGENCE := $(BUILD_DIR)/kerr_convergence
EXACT_TEST := $(BUILD_DIR)/test_kerr_exact

.PHONY: all run regen check-generated check clean help

all: $(CONVERGENCE)

$(BUILD_DIR):
	mkdir -p $@

$(HERMITE_LIBRARY):
	$(MAKE) -C ../hermite3d

$(EXACT_OBJECT): $(GENERATED_SOURCE) include/kerr_exact.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $(GENERATED_SOURCE) -o $@

$(CONVERGENCE_OBJECT): src/kerr_convergence.c include/kerr_exact.h \
                       ../hermite3d/include/hermite3d.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(CONVERGENCE): $(CONVERGENCE_OBJECT) $(EXACT_OBJECT) $(HERMITE_LIBRARY)
	$(CC) $(CFLAGS) $^ -lm -o $@

$(EXACT_TEST): tests/test_kerr_exact.c generated/kerr_exact_fixtures.h \
               $(EXACT_OBJECT) | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(EXACT_OBJECT) -lm -o $@

run: $(CONVERGENCE)
	$(CONVERGENCE) $(ARGS)

regen:
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) $(GENERATOR)

check-generated:
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) $(GENERATOR) --check

check: check-generated $(EXACT_TEST) $(CONVERGENCE)
	$(MAKE) -C ../hermite3d check
	$(MAKE) -C ../kerr_notes check
	$(EXACT_TEST)
	$(CONVERGENCE) --resolutions 8,12 --points 32 \
		--seed 0x4b4552524845524d > $(BUILD_DIR)/smoke-output.txt
	@grep -q "Kerr metric / quintic Hermite interpolation convergence" \
		$(BUILD_DIR)/smoke-output.txt

clean:
	$(RM) -r $(BUILD_DIR) __pycache__

help:
	@echo "make             build the convergence executable"
	@echo "make run ARGS='--resolutions 32,48,64 --points 10000'"
	@echo "make regen       regenerate the committed exact Kerr C source"
	@echo "make check       run symbolic, unit, library, and smoke checks"
	@echo "make clean       remove local build products and Python caches"
