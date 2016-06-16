// Copyright 2016 Riku Walve

#include <string>
#include <iostream>
#include <vector>

#include "htslib/sam.h"

#include "io.hpp"
#include "hash.hpp"
#include "extract.hpp"

int main(int argc, char* argv[]) {
  if (argc != 8) {
    std::cerr << argv[0] << " <in1>.bam [<in2>.bam]" <<
      " <read length> <mu> <sd>" <<
      " <scaffold name> <gap start> <gap end>" << std::endl;

    return 1;
  }
  
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
}
