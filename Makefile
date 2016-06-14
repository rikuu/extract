CXX=clang++

CXX_FLAGS=-m64 -std=c++11 -pedantic-errors -W -Wall -Wextra -Wshadow \
					-Wpointer-arith -Wcast-qual -Wunused -Wwrite-strings -Werror

CLANG_FLAGS=-Wbool-conversions -Wshift-overflow -Wliteral-conversion \
						-Wstrict-prototypes -Wmissing-prototypes

OPT_FLAGS=-O3 -mmmx -msse -msse2 -msse3 -msse4 -msse4.2 -march=native

CXX_FLAGS += $(OPT_FLAGS)
CXX_FLAGS += $(CLANG_FLAGS)

HTS_PATH=../htslib
HTS_FLAGS=-I$(HTS_PATH)/ -L$(HTS_PATH)/ -lhts

SOURCES=io.cpp hash.cpp extract.cpp
TESTS=io-test.cpp hash-test.cpp extract-test.cpp

default: all

all: extract run-tests

catch.hpp:
	curl -O https://raw.githubusercontent.com/philsquared/Catch/master/single_include/catch.hpp

extract:
	$(CXX) $(CXX_FLAGS) -o $@ $(SOURCES) main.cpp $(HTS_FLAGS)

extract-test: catch.hpp
	$(CXX) $(CXX_FLAGS) -o $@ $(SOURCES) $(TESTS) $(HTS_FLAGS)

run-tests: extract-test
	./extract-test

clean:
	rm -f extract extract-test *.o *.dSYM
