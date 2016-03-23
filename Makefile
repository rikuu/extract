CXX=clang++

CPP_FLAGS=-m64 -std=c++11 -pedantic-errors -W -Wall -Wextra -Wshadow -Wpointer-arith -Wcast-qual \
					-Wunused -Wstrict-prototypes -Wmissing-prototypes -Wwrite-strings \
					-Wbool-conversions -Wshift-overflow -Wliteral-conversion \
					-Werror

OPT_FLAGS=-O3 -mmmx -msse -msse2 -msse3 -msse4 -msse4.2 -march=native
CPP_FLAGS+=$(OPT_FLAGS)

HTS_PATH=../htslib
HTS_FLAGS=-I$(HTS_PATH)/ -L$(HTS_PATH)/ -lhts

SOURCES=io.cpp hash.cpp extract.cpp
TESTS=io-test.cpp hash-test.cpp extract-test.cpp

BINARIES=extract extract-test

default: all

extract:
		$(CXX) $(CPP_FLAGS) -o $@ $(SOURCES) $(HTS_FLAGS)

all: $(BINARIES) run-tests

catch.hpp:
		curl -O https://raw.githubusercontent.com/philsquared/Catch/master/single_include/catch.hpp

extract-test: catch.hpp hash.cpp io.cpp $(TESTS)
		$(CXX) $(CPP_FLAGS) -o $@ $(filter-out %.hpp,$^) $(HTS_FLAGS)

run-tests: extract-test
		./extract-test

clean:
		rm -rf $(BINARIES) *.o *.dSYM
