CXX=clang++

CXX_FLAGS=-m64 -std=c++11 -pedantic-errors -W -Wall -Wextra -Wshadow \
					-Wpointer-arith -Wcast-qual -Wunused -Wwrite-strings -Werror

CLANG_FLAGS=-Wbool-conversions -Wshift-overflow -Wliteral-conversion \
						-Wstrict-prototypes -Wmissing-prototypes

OPT_FLAGS=-O3 -mmmx -msse -msse2 -msse3 -msse4 -msse4.2 -march=native

CXX_FLAGS += $(OPT_FLAGS) $(CLANG_FLAGS)

# For catch.hpp
CXX_FLAGS += -Iinclude -I.

HTS_PATH=../htslib/
HTS_FLAGS=-I$(HTS_PATH) -L$(HTS_PATH) -lhts # s-Wl,-rpath=$(HTS_PATH)
CXX_FLAGS += $(HTS_FLAGS)

SOURCES=src/io.cpp src/hash.cpp src/extract.cpp
TESTS=src/io-test.cpp src/hash-test.cpp src/extract-test.cpp

default: all

all: extract run-tests

catch.hpp:
	curl -O https://raw.githubusercontent.com/philsquared/Catch/master/single_include/catch.hpp

extract:
	$(CXX) $(CXX_FLAGS) -o $@ $(SOURCES) src/main.cpp

extract-test: catch.hpp
	$(CXX) $(CXX_FLAGS) -o $@ $(SOURCES) $(TESTS)

run-tests: extract-test
	./extract-test

clean:
	rm -f extract extract-test *.o *.dSYM
