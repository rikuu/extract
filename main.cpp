// Copyright 2016 Riku Walve

#include <string>
#include <iostream>
#include <vector>

#include "htslib/sam.h"

#include "io.hpp"
#include "hash.hpp"
#include "extract.hpp"

int main(int argc, char* argv[]) {
  if (argc == 8) {
    io_t io(argv[1]);
    if (!io.loaded) {
      std::cerr << "Error loading alignments" << std::endl;
      return 1;
    }

    const int read_length = std::stoi(argv[2]);
    const int mean_insert = std::stoi(argv[3]);
    const int std_dev = std::stoi(argv[4]);

    const int tid = bam_name2id(io.header, argv[5]);
    const int start = std::stoi(argv[6]);
    const int end = std::stoi(argv[7]);

    run_extract(io,
      read_length, mean_insert,
      std_dev, tid, start, end);

    io.unload();

    return 0;
  } else if (argc == 9) {
    io_t pe1_io(argv[1]), pe2_io(argv[2]);
    if (!pe1_io.loaded || !pe2_io.loaded) {
      std::cerr << "Error loading alignments" << std::endl;
      return 1;
    }

    const int read_length = std::stoi(argv[3]);
    const int mean_insert = std::stoi(argv[4]);
    const int std_dev = std::stoi(argv[5]);

    const int tid = bam_name2id(pe1_io.header, argv[6]);
    const int start = std::stoi(argv[7]);
    const int end = std::stoi(argv[8]);

    run_extract(pe1_io, pe2_io,
      read_length, mean_insert,
      std_dev, tid, start, end);

    pe1_io.unload();
    pe2_io.unload();

    return 0;
  } else {
    std::cerr << argv[0] << " <in1>.bam [<in2>.bam]" <<
      " <read length> <mu> <sd>" <<
      " <scaffold name> <gap start> <gap end>" << std::endl;

    return 1;
  }

  return 0;
}
