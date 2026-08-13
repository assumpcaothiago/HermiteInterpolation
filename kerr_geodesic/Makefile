CC ?= cc
PYTHON ?= python3

CPPFLAGS ?=
CFLAGS ?= -O2
OPENMP_CFLAGS ?= -fopenmp
OPENMP_LDFLAGS ?= -fopenmp

STRICT_CFLAGS := -std=c99 -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror
VENDOR_CFLAGS := -std=gnu99 -Wall -Wextra -Wconversion -Wshadow -Werror
CPPFLAGS += -Iinclude -I../hermite3d/include -I../kerr_adm_comparison/include \
	-I../kerr_adm_comparison/src

BUILD_DIR := build
PROGRAM := $(BUILD_DIR)/kerr_geodesic
CORE_TEST := $(BUILD_DIR)/test_geodesic_core
HERMITE_LIBRARY := ../hermite3d/build/libhermite3d.a

MAIN_OBJECT := $(BUILD_DIR)/kerr_geodesic.o
CORE_OBJECT := $(BUILD_DIR)/kerr_geodesic_core.o
EXACT_OBJECT := $(BUILD_DIR)/kerr_adm_exact.o
ADAPTER_OBJECT := $(BUILD_DIR)/lagrange3d_reference.o
VENDOR_OBJECT := $(BUILD_DIR)/interpolation_3d_general__uniform_src_grid.o

ADM_DIR := ../kerr_adm_comparison
VENDOR_DIR := $(ADM_DIR)/vendor/bhah_lagrange
TRAJECTORY_CSV ?= $(BUILD_DIR)/trajectory.csv
PLOT_OUTPUT ?= $(BUILD_DIR)/trajectory.png
PLOT_ARGS ?=

.PHONY: all run plot check sanitize clean help

all: $(PROGRAM)

$(BUILD_DIR):
	mkdir -p $@

$(HERMITE_LIBRARY):
	$(MAKE) -C ../hermite3d

$(MAIN_OBJECT): src/kerr_geodesic.c include/kerr_geodesic_core.h \
                $(ADM_DIR)/include/kerr_adm_exact.h \
                $(ADM_DIR)/include/lagrange3d_reference.h \
                $(ADM_DIR)/src/kerr_adm_noise.h \
                ../hermite3d/include/hermite3d.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(STRICT_CFLAGS) $(OPENMP_CFLAGS) -c $< -o $@

$(CORE_OBJECT): src/kerr_geodesic_core.c include/kerr_geodesic_core.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(STRICT_CFLAGS) -c $< -o $@

$(EXACT_OBJECT): $(ADM_DIR)/generated/kerr_adm_exact.c \
                 $(ADM_DIR)/include/kerr_adm_exact.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(STRICT_CFLAGS) -c $< -o $@

$(ADAPTER_OBJECT): $(ADM_DIR)/src/lagrange3d_reference.c \
                   $(ADM_DIR)/include/lagrange3d_reference.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(STRICT_CFLAGS) -c $< -o $@

$(VENDOR_OBJECT): $(VENDOR_DIR)/interpolation_3d_general__uniform_src_grid.c \
                  $(VENDOR_DIR)/interpolation_lagrange_uniform.h \
                  $(VENDOR_DIR)/intrinsics/simd_intrinsics.h | $(BUILD_DIR)
	$(CC) -I$(VENDOR_DIR) $(CFLAGS) $(VENDOR_CFLAGS) $(OPENMP_CFLAGS) \
		-c $< -o $@

$(PROGRAM): $(MAIN_OBJECT) $(CORE_OBJECT) $(EXACT_OBJECT) $(ADAPTER_OBJECT) \
            $(VENDOR_OBJECT) $(HERMITE_LIBRARY)
	$(CC) $(CFLAGS) $^ $(OPENMP_LDFLAGS) -lm -o $@

$(CORE_TEST): tests/test_geodesic_core.c $(CORE_OBJECT) $(EXACT_OBJECT) | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(STRICT_CFLAGS) $< $(CORE_OBJECT) \
		$(EXACT_OBJECT) -lm -o $@

run: $(PROGRAM)
	$(PROGRAM) $(ARGS)

plot: $(PROGRAM)
	$(PROGRAM) $(ARGS) --trajectory-csv $(TRAJECTORY_CSV)
	$(PYTHON) plot_trajectory.py $(TRAJECTORY_CSV) \
		--output $(PLOT_OUTPUT) $(PLOT_ARGS)

check: $(PROGRAM) $(CORE_TEST)
	$(MAKE) -C ../hermite3d check
	$(MAKE) -C ../kerr_adm_comparison check
	$(MAKE) -C ../geodesic_notes check
	$(CORE_TEST)
	OMP_NUM_THREADS=2 OMP_DYNAMIC=FALSE $(PROGRAM) \
		--resolutions 16,20 --orbits 1 --steps-per-orbit 128 \
		--output-every 8 --trajectory-csv $(BUILD_DIR)/smoke.csv \
		> $(BUILD_DIR)/smoke.txt
	$(PYTHON) tests/validate_trajectory_csv.py $(BUILD_DIR)/smoke.csv \
		--resolutions 16,20
	$(PYTHON) tests/test_plot_data.py $(BUILD_DIR)/smoke.csv
	@grep -q "Local RHS defects" $(BUILD_DIR)/smoke.txt
	@grep -q "Invariant and Hamiltonian drift" $(BUILD_DIR)/smoke.txt
	OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE $(PROGRAM) \
		--resolutions 16 --orbits 1 --steps-per-orbit 32 \
		--output-every 8 --trajectory-csv $(BUILD_DIR)/clean-default.csv \
		> /dev/null
	OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE $(PROGRAM) \
		--resolutions 16 --orbits 1 --steps-per-orbit 32 \
		--noise-epsilon 0 --output-every 8 \
		--trajectory-csv $(BUILD_DIR)/clean-zero.csv > /dev/null
	cmp $(BUILD_DIR)/clean-default.csv $(BUILD_DIR)/clean-zero.csv
	OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE $(PROGRAM) \
		--resolutions 16 --orbits 1 --steps-per-orbit 32 \
		--noise-epsilon 1e-8 --output-every 8 \
		--trajectory-csv $(BUILD_DIR)/noisy-one.csv > /dev/null
	OMP_NUM_THREADS=2 OMP_DYNAMIC=FALSE $(PROGRAM) \
		--resolutions 16 --orbits 1 --steps-per-orbit 32 \
		--noise-epsilon 1e-8 --output-every 8 \
		--trajectory-csv $(BUILD_DIR)/noisy-two.csv > /dev/null
	cmp $(BUILD_DIR)/noisy-one.csv $(BUILD_DIR)/noisy-two.csv
	$(PYTHON) -m py_compile plot_trajectory.py tests/validate_trajectory_csv.py \
		tests/test_plot_data.py
	@! $(PROGRAM) --resolutions 15 --orbits 1 >/dev/null 2>&1
	@! $(PROGRAM) --resolutions 16,16 --orbits 1 >/dev/null 2>&1
	@! $(PROGRAM) --resolutions 16 --noise-epsilon=-1 >/dev/null 2>&1
	@! $(PROGRAM) --resolutions 16 --steps-per-orbit 0 >/dev/null 2>&1

sanitize: $(CORE_TEST)
	$(CC) $(CPPFLAGS) -O1 -g $(STRICT_CFLAGS) \
		-fsanitize=address,undefined -fno-omit-frame-pointer \
		tests/test_geodesic_core.c src/kerr_geodesic_core.c \
		$(ADM_DIR)/generated/kerr_adm_exact.c -lm \
		-o $(BUILD_DIR)/test_geodesic_core_sanitize
	ASAN_OPTIONS=detect_leaks=0 $(BUILD_DIR)/test_geodesic_core_sanitize
	$(CC) $(CPPFLAGS) -O1 -g $(STRICT_CFLAGS) $(OPENMP_CFLAGS) \
		-fsanitize=address,undefined -fno-omit-frame-pointer \
		-c src/kerr_geodesic.c -o $(BUILD_DIR)/sanitize_main.o
	$(CC) $(CPPFLAGS) -O1 -g $(STRICT_CFLAGS) \
		-fsanitize=address,undefined -fno-omit-frame-pointer \
		-c src/kerr_geodesic_core.c -o $(BUILD_DIR)/sanitize_core.o
	$(CC) $(CPPFLAGS) -O1 -g $(STRICT_CFLAGS) \
		-fsanitize=address,undefined -fno-omit-frame-pointer \
		-c $(ADM_DIR)/generated/kerr_adm_exact.c \
		-o $(BUILD_DIR)/sanitize_exact.o
	$(CC) $(CPPFLAGS) -O1 -g $(STRICT_CFLAGS) \
		-fsanitize=address,undefined -fno-omit-frame-pointer \
		-c $(ADM_DIR)/src/lagrange3d_reference.c \
		-o $(BUILD_DIR)/sanitize_adapter.o
	$(CC) -I$(VENDOR_DIR) -O1 -g $(VENDOR_CFLAGS) $(OPENMP_CFLAGS) \
		-fsanitize=address,undefined -fno-omit-frame-pointer \
		-c $(VENDOR_DIR)/interpolation_3d_general__uniform_src_grid.c \
		-o $(BUILD_DIR)/sanitize_vendor.o
	$(CC) -I../hermite3d/include -O1 -g $(STRICT_CFLAGS) \
		-fsanitize=address,undefined -fno-omit-frame-pointer \
		-c ../hermite3d/src/hermite3d.c -o $(BUILD_DIR)/sanitize_hermite.o
	$(CC) -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer \
		$(BUILD_DIR)/sanitize_main.o $(BUILD_DIR)/sanitize_core.o \
		$(BUILD_DIR)/sanitize_exact.o $(BUILD_DIR)/sanitize_adapter.o \
		$(BUILD_DIR)/sanitize_vendor.o $(BUILD_DIR)/sanitize_hermite.o \
		$(OPENMP_LDFLAGS) -lm -o $(BUILD_DIR)/kerr_geodesic_sanitize
	ASAN_OPTIONS=detect_leaks=0 OMP_NUM_THREADS=2 OMP_DYNAMIC=FALSE \
		$(BUILD_DIR)/kerr_geodesic_sanitize --resolutions 16 --orbits 1 \
		--steps-per-orbit 32 --noise-epsilon 1e-8 > /dev/null

clean:
	$(RM) -r $(BUILD_DIR) __pycache__ tests/__pycache__

help:
	@echo "make                 build the geodesic comparison"
	@echo "make run ARGS='--resolutions 64,96 --orbits 10 --steps-per-orbit 4096'"
	@echo "make plot ARGS='--resolutions 64 --orbits 10 --steps-per-orbit 4096'"
	@echo "OMP_NUM_THREADS=8 make run ARGS='--resolutions 96 --noise-epsilon 1e-8'"
	@echo "make check           run existing dependencies, unit tests, and smoke tests"
	@echo "make sanitize        run the geodesic core with ASan and UBSan"
	@echo "make clean           remove local build and plot products"
