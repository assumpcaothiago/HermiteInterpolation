CC ?= cc
PYTHON ?= python3

CPPFLAGS ?=
CFLAGS ?= -O2
OPENMP_CFLAGS ?= -fopenmp
OPENMP_LDFLAGS ?= -fopenmp

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
PARALLEL_HELPER_TEST := $(BUILD_DIR)/test_parallel_helpers
PROFILE ?= xx:value
Z_SAMPLES ?= 2000
PROFILE_CSV ?= $(BUILD_DIR)/z_profile.csv
PROFILE_PNG ?= $(BUILD_DIR)/z_profile.png

.PHONY: all run plot regen check-generated check clean help

all: $(CONVERGENCE)

$(BUILD_DIR):
	mkdir -p $@

$(HERMITE_LIBRARY):
	$(MAKE) -C ../hermite3d

$(EXACT_OBJECT): $(GENERATED_SOURCE) include/kerr_exact.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $(GENERATED_SOURCE) -o $@

$(CONVERGENCE_OBJECT): src/kerr_convergence.c src/kerr_parallel_helpers.h \
                       include/kerr_exact.h \
                       ../hermite3d/include/hermite3d.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPENMP_CFLAGS) -c $< -o $@

$(CONVERGENCE): $(CONVERGENCE_OBJECT) $(EXACT_OBJECT) $(HERMITE_LIBRARY)
	$(CC) $(CFLAGS) $^ $(OPENMP_LDFLAGS) -lm -o $@

$(EXACT_TEST): tests/test_kerr_exact.c generated/kerr_exact_fixtures.h \
               $(EXACT_OBJECT) | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(EXACT_OBJECT) -lm -o $@

$(PARALLEL_HELPER_TEST): tests/test_parallel_helpers.c \
                         src/kerr_parallel_helpers.h | $(BUILD_DIR)
	$(CC) -Isrc $(CFLAGS) $< -lm -o $@

run: $(CONVERGENCE)
	$(CONVERGENCE) $(ARGS)

plot: $(CONVERGENCE) plot_z_profile.py
	$(CONVERGENCE) $(ARGS) --z-profile '$(PROFILE)' \
		--z-samples $(Z_SAMPLES) --z-output '$(PROFILE_CSV)'
	MPLCONFIGDIR='$(abspath $(BUILD_DIR)/matplotlib-cache)' \
		$(PYTHON) plot_z_profile.py '$(PROFILE_CSV)' \
		--output '$(PROFILE_PNG)'

regen:
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) $(GENERATOR)

check-generated:
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) $(GENERATOR) --check

check: check-generated $(EXACT_TEST) $(PARALLEL_HELPER_TEST) $(CONVERGENCE)
	$(MAKE) -C ../hermite3d check
	$(MAKE) -C ../kerr_notes check
	$(EXACT_TEST)
	$(PARALLEL_HELPER_TEST)
	OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE $(CONVERGENCE) \
		--resolutions 8,12 --points 32 \
		--seed 0x4b4552524845524d --z-profile tt:value \
		--z-samples 8 --z-output $(BUILD_DIR)/smoke-profile-serial.csv \
		> $(BUILD_DIR)/smoke-output-serial.txt
	OMP_NUM_THREADS=2 OMP_DYNAMIC=FALSE $(CONVERGENCE) \
		--resolutions 8,12 --points 32 \
		--seed 0x4b4552524845524d --z-profile tt:value \
		--z-samples 8 --z-output $(BUILD_DIR)/smoke-profile-parallel.csv \
		> $(BUILD_DIR)/smoke-output-parallel.txt
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) tests/compare_parallel_output.py \
		$(BUILD_DIR)/smoke-output-serial.txt \
		$(BUILD_DIR)/smoke-output-parallel.txt --parallel-threads 2
	cmp $(BUILD_DIR)/smoke-profile-serial.csv \
		$(BUILD_DIR)/smoke-profile-parallel.csv
	@grep -q "Kerr metric / quintic Hermite interpolation convergence" \
		$(BUILD_DIR)/smoke-output-parallel.txt
	@grep -q "component,quantity,resolution,z,analytic,interpolated,error" \
		$(BUILD_DIR)/smoke-profile-parallel.csv
	@test "$$(wc -l < $(BUILD_DIR)/smoke-profile-parallel.csv)" -eq 17
	@awk -F, 'NR > 1 && $$4 > -0.46650635094610965 && \
		$$4 < 0.46650635094610965 { exit 1 }' \
		$(BUILD_DIR)/smoke-profile-parallel.csv
	@awk '$$1 == 8 || $$1 == 12 { if ($$4 < 0.46650635094610965) exit 1 }' \
		$(BUILD_DIR)/smoke-output-parallel.txt
	OMP_NUM_THREADS=2 OMP_DYNAMIC=FALSE $(CONVERGENCE) \
		--resolutions 8 --points 1 > $(BUILD_DIR)/smoke-points-1.txt
	OMP_NUM_THREADS=2 OMP_DYNAMIC=FALSE $(CONVERGENCE) \
		--resolutions 8 --points 2 > $(BUILD_DIR)/smoke-points-2.txt
	@awk '$$1 == 8 { if ($$7 != 2 || $$8 != 2) exit 1; found = 1 } \
		END { if (!found) exit 1 }' $(BUILD_DIR)/smoke-points-1.txt
	@awk '$$1 == 8 { if ($$7 != 2 || $$8 != 2) exit 1; found = 1 } \
		END { if (!found) exit 1 }' $(BUILD_DIR)/smoke-points-2.txt
	@PYTHONDONTWRITEBYTECODE=1 $(PYTHON) plot_z_profile.py --help > /dev/null

clean:
	$(RM) -r $(BUILD_DIR) __pycache__

help:
	@echo "make             build the convergence executable"
	@echo "make run ARGS='--resolutions 32,48,64 --points 10000'"
	@echo "OMP_NUM_THREADS=8 make run ARGS='--resolutions 64 --points 10000'"
	@echo "make plot PROFILE=xx:value ARGS='--resolutions 32,48,64 --points 10000'"
	@echo "OMP_NUM_THREADS=8 make plot PROFILE=xx:dz ARGS='--resolutions 64 --points 10000'"
	@echo "make regen       regenerate the committed exact Kerr C source"
	@echo "make check       run symbolic, unit, library, and smoke checks"
	@echo "make clean       remove local build products and Python caches"
