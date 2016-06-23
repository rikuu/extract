// Copyright 2016 Riku Walve

#include <string>
#include <iostream>
#include <vector>

#include "htslib/sam.h"

#include "io.hpp"
#include "hash.hpp"
#include "extract.hpp"

int main(int argc, char* argv[]) {
  // TODO: Actual argument parsing
  if (argc != 12) {
    std::cerr << "Usage: " << argv[0] << " <alignment>.bam" <<
      " <read length> <mu> <sd>" <<
      " <scaffold name> <gap start> <gap end> <gap length>" <<
      " <exact> <unmapped> <threshold>" << std::endl;

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
  const int length = std::stoi(argv[8]);

  const bool exact = std::stoi(argv[9]) == 1;
  const bool unmapped = std::stoi(argv[10]) == 1;
  const int threshold = std::stoi(argv[11]);

  run_extract(io,
    read_length, mean_insert, std_dev,
    tid, start, end, length,
    exact, unmapped, threshold);

  io.unload();

  return 0;
}
