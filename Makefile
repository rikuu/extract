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

default: all

all: extract run-tests

catch.hpp:
	curl -O https://raw.githubusercontent.com/philsquared/Catch/master/single_include/catch.hpp

extract:
	$(CXX) $(CPP_FLAGS) -o $@ $(SOURCES) main.cpp $(HTS_FLAGS)

extract-test: catch.hpp
	$(CXX) $(CPP_FLAGS) -o $@ $(SOURCES) $(TESTS) $(HTS_FLAGS)

run-tests: extract-test
	./extract-test

clean:
	rm -f catch.hpp extract extract-test *.o *.dSYM
