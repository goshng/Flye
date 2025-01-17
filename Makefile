ROOT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
THREADS := 4

export LIBCUCKOO = -I${ROOT_DIR}/lib/libcuckoo
export INTERVAL_TREE = -I${ROOT_DIR}/lib/interval_tree
export LEMON = -I${ROOT_DIR}/lib/lemon
export BIN_DIR = ${ROOT_DIR}/bin
export BUILD_DIR = ${ROOT_DIR}/build
export MINIMAP2_DIR = ${ROOT_DIR}/lib/minimap2
export SAMTOOLS_DIR = ${ROOT_DIR}/lib/samtools-1.9

export CXXFLAGS += ${LIBCUCKOO} ${INTERVAL_TREE} ${LEMON} -I${MINIMAP2_DIR} -I${BUILD_DIR}
export LDFLAGS += -lz -L${MINIMAP2_DIR} -lminimap2

ifeq ($(shell uname -m),arm64)
	export arm_neon=1
	export aarch64=1
endif

.PHONY: clean all profile debug minimap2 samtools

.DEFAULT_GOAL := all

# Ensure the build directory exists
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

${BIN_DIR}/dflye-minimap2: $(BUILD_DIR)
	make -C ${MINIMAP2_DIR} -j ${THREADS}
	cp ${MINIMAP2_DIR}/minimap2 ${BIN_DIR}/dflye-minimap2

minimap2: ${BIN_DIR}/dflye-minimap2

samtools: ${BIN_DIR}/dflye-samtools

${BIN_DIR}/dflye-samtools: $(BUILD_DIR)
	cd ${SAMTOOLS_DIR} && ./configure --without-curses --disable-bz2 --disable-lzma --enable-plugins
	make samtools -C ${SAMTOOLS_DIR} -j ${THREADS}
	cp ${SAMTOOLS_DIR}/samtools ${BIN_DIR}/dflye-samtools

all: $(BUILD_DIR) minimap2 samtools
	make release -C src -j ${THREADS} BUILD_DIR=$(BUILD_DIR)

profile: $(BUILD_DIR) minimap2 samtools
	make profile -C src -j ${THREADS} BUILD_DIR=$(BUILD_DIR)

debug: $(BUILD_DIR) minimap2 samtools
	make debug -C src -j ${THREADS} BUILD_DIR=$(BUILD_DIR)

clean:
	make clean -C src BUILD_DIR=$(BUILD_DIR)
	make clean -C ${MINIMAP2_DIR}
	make clean-all -C ${SAMTOOLS_DIR}
	rm -rf ${BUILD_DIR}
	rm -f ${BIN_DIR}/dflye-minimap2
	rm -f ${BIN_DIR}/dflye-samtools

