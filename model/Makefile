# Bubble Model (macaque) — build file. Builds as C23 (-std=gnu2x), required by include/gplot.h.
#
#   make        local build (macOS, Homebrew libomp + mpicc; on Linux: make OMP=-fopenmp)
#   make hpc    cluster build (Intel oneAPI: mpiicx), as used for the paper
#   make debug  local build with -g
#   make sanitize  local build with AddressSanitizer + UBSan (checks only)

EXEC   = bubble-model
C_FILES = $(wildcard src/*.c)
WARN   = -Wall -Wextra -Wpedantic -Werror
# C23 for gplot.h; the GNU variant keeps POSIX popen/open_memstream/mkdtemp visible on glibc
STD    = -std=gnu2x
# compile-time overrides, e.g. make DEFS='-DCONS_MIN_MAINAXIS=5.5 -DMICRONOD_MIN_MAINAXIS=0.5'
DEFS  ?=

# --- local (macOS) ---------------------------------------------------------
CC      ?= mpicc
LIBOMP  ?= $(shell brew --prefix libomp 2>/dev/null)
OMP      = -isystem $(LIBOMP)/include -L$(LIBOMP)/lib -lomp -Xpreprocessor -fopenmp

local: $(C_FILES)
	mpicc $(OMP) -Iinclude $(DEFS) $(STD) $(WARN) -O3 -o $(EXEC) $(C_FILES)

debug: $(C_FILES)
	mpicc $(OMP) -Iinclude $(DEFS) $(STD) $(WARN) -g -o $(EXEC) $(C_FILES)

# checks only (memory errors, undefined behaviour), not for runs
sanitize: $(C_FILES)
	mpicc $(OMP) -Iinclude $(DEFS) $(STD) $(WARN) -g -O1 -fno-omit-frame-pointer \
		-fsanitize=address,undefined -o $(EXEC) $(C_FILES)

# --- cluster (Intel oneAPI) -------------------------------------------------
# I_MPI_CC=icx: some site modules set the Intel MPI wrappers to the classic icc backend, which
# has no C23 support; force the oneAPI icx compiler (what mpiicx is meant to call).
hpc: $(C_FILES)
	I_MPI_CC=icx mpiicx -Iinclude $(DEFS) $(STD) -qopenmp $(C_FILES) -o $(EXEC) -O3

clean:
	rm -f $(EXEC)

.PHONY: local debug sanitize hpc clean
