# Extract

Standalone tool for filtering reads based on paired-end read alignments.

## Compilation

CMake and a C++11 compiler, such as g++ 4.5 or newer are required.
Depends on HTSlib and GATB-core, both of which are included as git submodules.

Compile HTSlib according to their instructions.
GATB-core is compiled as needed by default.

  mkdir build; cd build; cmake ..; make

## Usage

Basic usage is

  extract -bam reads.bam -region contig1:100-110 -read-length 100 -mean 150 -std-dev 15 -reads filtered.fasta

The regions are defined the same way samtools does, i.e. CONTIG:START-END.
Read length should be given as the maximum read length if the reads are different length.

A more detailed list of possible arguments:

  -region [region]   : Region for filtering
  -std-dev [int]     : Insert size standard deviation
  -mean [int]        : Mean insert size
  -read-length [int] : Read length
  -reads [file]      : FASTA-formatted output file
  -bam [file]        : Aligned BAM file
  -unmapped          : Threshold for using unmapped reads  [default '0']
  -no-overlap        : Don't output overlapping reads
  -unmapped-only     : Only output unmapped reads
  -insertion         : Insertion filtering

Unmapped reads are added when the length of the filtered reads are less than the
threshold, thus the default value of 0 means the unmapped reads are not added.

Insertion filtering refers to the way the region is parsed. Assuming the reads
are aligned against the reference, the regions for filtering reads to genotype
the insertions do not technically exist in the reference. With insertion
filtering only the start of the region is used to define the read pair regions.
