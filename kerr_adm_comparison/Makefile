CC ?= cc
PYTHON ?= python3

CPPFLAGS ?=
CFLAGS ?= -O2
OPENMP_CFLAGS ?= -fopenmp
OPENMP_LDFLAGS ?= -fopenmp

STRICT_CFLAGS := -std=c99 -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror
VENDOR_CFLAGS := -std=gnu99 -Wall -Wextra -Wconversion -Wshadow -Werror
CPPFLAGS += -Iinclude -Igenerated -I../hermite3d/include

BUILD_DIR := build
PROGRAM := $(BUILD_DIR)/kerr_adm_comparison
EXACT_TEST := $(BUILD_DIR)/test_kerr_adm_exact
LAGRANGE_TEST := $(BUILD_DIR)/test_lagrange_reference
NOISE_TEST := $(BUILD_DIR)/test_kerr_adm_noise
HERMITE_LIBRARY := ../hermite3d/build/libhermite3d.a

GENERATOR := generate_kerr_adm_exact.py
GENERATED_SOURCE := generated/kerr_adm_exact.c
GENERATED_FIXTURES := generated/kerr_adm_exact_fixtures.h
CURRENT_METRIC_SOURCE := ../kerr_interpolation/generated/kerr_exact.c

MAIN_OBJECT := $(BUILD_DIR)/kerr_adm_comparison.o
EXACT_OBJECT := $(BUILD_DIR)/kerr_adm_exact.o
CURRENT_METRIC_OBJECT := $(BUILD_DIR)/kerr_exact.o
ADAPTER_OBJECT := $(BUILD_DIR)/lagrange3d_reference.o
VENDOR_OBJECT := $(BUILD_DIR)/interpolation_3d_general__uniform_src_grid.o

VENDOR_DIR := vendor/bhah_lagrange
BHAH_REFERENCE ?= ../changa_uw_new_build/bhah_lib

.PHONY: all run regen check-generated check-vendor check-upstream check clean help

all: $(PROGRAM)

$(BUILD_DIR):
	mkdir -p $@

$(HERMITE_LIBRARY):
	$(MAKE) -C ../hermite3d

$(EXACT_OBJECT): $(GENERATED_SOURCE) include/kerr_adm_exact.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(STRICT_CFLAGS) -c $< -o $@

$(CURRENT_METRIC_OBJECT): $(CURRENT_METRIC_SOURCE) \
                          ../kerr_interpolation/include/kerr_exact.h | $(BUILD_DIR)
	$(CC) -I../kerr_interpolation/include $(CFLAGS) $(STRICT_CFLAGS) -c $< -o $@

$(ADAPTER_OBJECT): src/lagrange3d_reference.c \
                   include/lagrange3d_reference.h include/kerr_adm_exact.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(STRICT_CFLAGS) -c $< -o $@

$(VENDOR_OBJECT): $(VENDOR_DIR)/interpolation_3d_general__uniform_src_grid.c \
                  $(VENDOR_DIR)/interpolation_lagrange_uniform.h \
                  $(VENDOR_DIR)/intrinsics/simd_intrinsics.h | $(BUILD_DIR)
	$(CC) -I$(VENDOR_DIR) $(CFLAGS) $(VENDOR_CFLAGS) $(OPENMP_CFLAGS) \
		-c $< -o $@

$(MAIN_OBJECT): src/kerr_adm_comparison.c src/kerr_adm_noise.h \
                include/kerr_adm_exact.h \
                include/lagrange3d_reference.h \
                ../hermite3d/include/hermite3d.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(STRICT_CFLAGS) $(OPENMP_CFLAGS) -c $< -o $@

$(PROGRAM): $(MAIN_OBJECT) $(EXACT_OBJECT) $(ADAPTER_OBJECT) $(VENDOR_OBJECT) \
            $(HERMITE_LIBRARY)
	$(CC) $(CFLAGS) $^ $(OPENMP_LDFLAGS) -lm -o $@

$(EXACT_TEST): tests/test_kerr_adm_exact.c $(GENERATED_FIXTURES) \
               $(EXACT_OBJECT) $(CURRENT_METRIC_OBJECT) | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) -I../kerr_interpolation/include $(CFLAGS) $(STRICT_CFLAGS) \
		tests/test_kerr_adm_exact.c $(EXACT_OBJECT) $(CURRENT_METRIC_OBJECT) \
		-lm -o $@

$(LAGRANGE_TEST): tests/test_lagrange_reference.c $(ADAPTER_OBJECT) \
                  $(VENDOR_OBJECT) | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(STRICT_CFLAGS) -c \
		tests/test_lagrange_reference.c -o $(BUILD_DIR)/test_lagrange_reference.o
	$(CC) $(CFLAGS) $(BUILD_DIR)/test_lagrange_reference.o $(ADAPTER_OBJECT) \
		$(VENDOR_OBJECT) $(OPENMP_LDFLAGS) -lm -o $@

$(NOISE_TEST): tests/test_kerr_adm_noise.c src/kerr_adm_noise.h | $(BUILD_DIR)
	$(CC) -Isrc $(CFLAGS) $(STRICT_CFLAGS) $< -lm -o $@

run: $(PROGRAM)
	$(PROGRAM) $(ARGS)

regen:
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) $(GENERATOR)

check-generated:
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) $(GENERATOR) --check

check-vendor:
	cd $(VENDOR_DIR) && sha256sum -c SHA256SUMS

check-upstream:
	test -d '$(BHAH_REFERENCE)'
	cmp $(VENDOR_DIR)/interpolation_3d_general__uniform_src_grid.c \
		$(BHAH_REFERENCE)/interpolation_3d_general__uniform_src_grid.c
	cmp $(VENDOR_DIR)/interpolation_lagrange_uniform.h \
		$(BHAH_REFERENCE)/interpolation_lagrange_uniform.h
	cmp $(VENDOR_DIR)/intrinsics/simd_intrinsics.h \
		$(BHAH_REFERENCE)/intrinsics/simd_intrinsics.h

check: check-generated check-vendor $(EXACT_TEST) $(LAGRANGE_TEST) \
       $(NOISE_TEST) $(PROGRAM)
	$(MAKE) -C ../hermite3d check
	$(MAKE) -C ../kerr_notes check
	$(EXACT_TEST)
	$(LAGRANGE_TEST)
	$(NOISE_TEST)
	OMP_NUM_THREADS=2 OMP_DYNAMIC=FALSE $(PROGRAM) \
		--resolutions 8,12 --points 64 --seed 0x4b4552524845524d \
		> $(BUILD_DIR)/smoke-output.txt
	@grep -q "Kerr ADM Hermite-Lagrange interpolation comparison" \
		$(BUILD_DIR)/smoke-output.txt
	@grep -q "HERMITE / VALUE" $(BUILD_DIR)/smoke-output.txt
	@grep -q "LAGRANGE+FD / D/DX" $(BUILD_DIR)/smoke-output.txt
	@grep -q "ANALYTIC-FD / D/DZ" $(BUILD_DIR)/smoke-output.txt
	@awk 'NF == 8 && ($$1 == 8 || $$1 == 12) { if ($$2 != $$1 + 12 || $$8 != 2) exit 1; found = 1 } \
		END { if (!found) exit 1 }' $(BUILD_DIR)/smoke-output.txt
	@awk 'NF == 8 && $$1 == 8 { if ($$7 < 1) exit 1; found = 1 } \
		END { if (!found) exit 1 }' $(BUILD_DIR)/smoke-output.txt
	@awk 'NF == 5 && ($$1 == 8 || $$1 == 12) { if ($$2 <= 0 || $$3 <= 0 || $$4 <= 0 || $$5 <= 0) exit 1; found = 1 } \
		END { if (!found) exit 1 }' $(BUILD_DIR)/smoke-output.txt
	OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE $(PROGRAM) \
		--resolutions 8,12 --points 32 \
		> $(BUILD_DIR)/clean-default.txt
	OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE $(PROGRAM) \
		--resolutions 8,12 --points 32 --noise-epsilon 0 \
		> $(BUILD_DIR)/clean-explicit-zero.txt
	$(PYTHON) tests/compare_noise_runs.py \
		$(BUILD_DIR)/clean-default.txt $(BUILD_DIR)/clean-explicit-zero.txt \
		--errors equal
	OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE $(PROGRAM) \
		--resolutions 8,12 --points 32 --noise-epsilon 1e-8 \
		--noise-seed 0x4e4f49534541444d > $(BUILD_DIR)/noisy-one-thread.txt
	OMP_NUM_THREADS=2 OMP_DYNAMIC=FALSE $(PROGRAM) \
		--resolutions 8,12 --points 32 --noise-epsilon 1e-8 \
		--noise-seed 0x4e4f49534541444d > $(BUILD_DIR)/noisy-two-threads.txt
	$(PYTHON) tests/compare_noise_runs.py \
		$(BUILD_DIR)/noisy-one-thread.txt $(BUILD_DIR)/noisy-two-threads.txt \
		--errors equal
	OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE $(PROGRAM) \
		--resolutions 8,12 --points 32 --noise-epsilon 1e-8 \
		--noise-seed 0x12345678 > $(BUILD_DIR)/noisy-other-seed.txt
	$(PYTHON) tests/compare_noise_runs.py \
		$(BUILD_DIR)/noisy-one-thread.txt $(BUILD_DIR)/noisy-other-seed.txt \
		--errors different
	@! $(PROGRAM) --resolutions 8 --points 1 --noise-epsilon=-1 >/dev/null 2>&1
	@! $(PROGRAM) --resolutions 8 --points 1 --noise-epsilon=1.01 >/dev/null 2>&1
	@! $(PROGRAM) --resolutions 8 --points 1 --noise-epsilon=nan >/dev/null 2>&1
	@! $(PROGRAM) --resolutions 8 --points 1 --noise-epsilon=inf >/dev/null 2>&1
	@! $(PROGRAM) --resolutions 8 --points 1 --noise-epsilon=garbage >/dev/null 2>&1

clean:
	$(RM) -r $(BUILD_DIR) __pycache__

help:
	@echo "make                 build the ADM comparison executable"
	@echo "make run ARGS='--resolutions 32,48,64 --points 10000'"
	@echo "OMP_NUM_THREADS=8 make run ARGS='--resolutions 64 --points 10000'"
	@echo "make run ARGS='--resolutions 64 --points 10000 --noise-epsilon 1e-8'"
	@echo "make regen           regenerate committed ADM expressions and fixtures"
	@echo "make check-generated verify generated files are current"
	@echo "make check-vendor    verify copied Lagrange file checksums"
	@echo "make check-upstream  compare copied files with BHAH_REFERENCE"
	@echo "make check           run symbolic, unit, library, and smoke checks"
	@echo "make clean           remove local build products and Python caches"
