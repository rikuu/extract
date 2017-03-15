# Extract

Standalone tool for filtering reads based on paired-end read alignments.

## Compilation

CMake and a C++11 compiler, such as g++ 4.5 or newer are required.
Depends on HTSlib and GATB-core, both of which are included as git submodules.

Compile HTSlib according to their instructions.
GATB-core is compiled as needed by default.

  mkdir build; cd build; cmake ..; make

## Usage
