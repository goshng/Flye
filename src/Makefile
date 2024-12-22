.PHONY: all clean debug profile

CXXFLAGS += -Wall -Wextra -pthread -std=c++11 -g
CXXFLAGS += -Wno-missing-field-initializers
LDFLAGS += -pthread -std=c++11 -rdynamic

MODULES_BIN := ${BIN_DIR}/aflye-modules
BUILD_DIR := build

profile: CXXFLAGS += -pg
profile: LDFLAGS += -pg
profile: release

release: CXXFLAGS += -O3 -DNDEBUG
release: aflye-modules

SANITIZE_FLAGS += -D_GLIBCXX_SANITIZE_VECTOR -U_FORTIFY_SOURCE -fsanitize=address -fno-omit-frame-pointer -fsanitize=undefined
#SANITIZE_FLAGS += -fsanitize=thread -fsanitize=leak -fsanitize=undefined
debug: CXXFLAGS += -Og ${SANITIZE_FLAGS}
#debug: CXXFLAGS += -D_GLIBCXX_DEBUG
debug: LDFLAGS += ${SANITIZE_FLAGS}
debug: aflye-modules

# Ensure build directory exists
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# sequence module
sequence_obj := ${patsubst sequence/%.cpp,${BUILD_DIR}/sequence_%.o,${wildcard sequence/*.cpp}}

${BUILD_DIR}/sequence_%.o: sequence/%.cpp sequence/*.h common/*.h | $(BUILD_DIR)
	${CXX} -c ${CXXFLAGS} $< -o $@

# aflye-assemble module
assemble_obj := ${patsubst assemble/%.cpp,${BUILD_DIR}/assemble_%.o,${wildcard assemble/*.cpp}}

${BUILD_DIR}/assemble_%.o: assemble/%.cpp assemble/*.h sequence/*.h common/*.h | $(BUILD_DIR)
	${CXX} -c ${CXXFLAGS} $< -o $@

# aflye-repeat module
repeat_obj := ${patsubst repeat_graph/%.cpp,${BUILD_DIR}/repeat_%.o,${wildcard repeat_graph/*.cpp}}

${BUILD_DIR}/repeat_%.o: repeat_graph/%.cpp repeat_graph/*.h sequence/*.h common/*.h | $(BUILD_DIR)
	${CXX} -c ${CXXFLAGS} $< -o $@

# aflye-contigger module
contigger_obj := ${patsubst contigger/%.cpp,${BUILD_DIR}/contigger_%.o,${wildcard contigger/*.cpp}}

${BUILD_DIR}/contigger_%.o: contigger/%.cpp repeat_graph/*.h sequence/*.h common/*.h | $(BUILD_DIR)
	${CXX} -c ${CXXFLAGS} $< -o $@

# aflye-polish module
polish_obj := ${patsubst polishing/%.cpp,${BUILD_DIR}/polish_%.o,${wildcard polishing/*.cpp}}

${BUILD_DIR}/polish_%.o: polishing/%.cpp polishing/*.h common/*.h | $(BUILD_DIR)
	${CXX} -c ${CXXFLAGS} $< -o $@

# main module
main_obj := ${BUILD_DIR}/main.o

${BUILD_DIR}/main.o: main.cpp | $(BUILD_DIR)
	${CXX} -c ${CXXFLAGS} $< -o $@

aflye-modules: ${assemble_obj} ${sequence_obj} ${repeat_obj} ${contigger_obj} ${polish_obj} ${main_obj}
	${CXX} ${assemble_obj} ${sequence_obj} ${repeat_obj} ${contigger_obj} ${polish_obj} ${main_obj} -o ${MODULES_BIN} ${LDFLAGS}

clean:
	rm -rf ${BUILD_DIR}
	rm -f ${MODULES_BIN}

