CXX=clang++

CPP_FLAGS=-m64 -std=c++11 -pedantic-errors -W -Wall -Wextra -Wshadow -Wpointer-arith -Wcast-qual \
					-Wunused -Wstrict-prototypes -Wmissing-prototypes -Wwrite-strings \
					-Wbool-conversions -Wshift-overflow -Wliteral-conversion \
					-Werror

OPT_FLAGS=-O3 -mmmx -msse -msse2 -msse3 -msse4 -msse4.2 -march=native
CPP_FLAGS+=$(OPT_FLAGS)

HTS_PATH=../htslib
HTS_FLAGS=-I$(HTS_PATH)/ -L$(HTS_PATH)/ -lhts

SOURCES=io.cpp hash.cpp bam_bed_extract.cpp

default: all

extract:
		$(CXX) $(CPP_FLAGS) -o $@ $(SOURCES) $(HTS_FLAGS)

all: extract

clean:
		rm -rf $(BINARIES) *.o *.dSYM
