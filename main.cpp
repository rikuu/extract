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
    std::cerr << argv[0] << " <alignments>.bam <read length> <mu> <sd>"
      << " <scaffold name> <gap start> <gap end>" << std::endl;
    return -1;
  }

  const std::string samFilename = argv[1];
  const int read_length = std::stoi(argv[2]);
  const int mean_insert = std::stoi(argv[3]);
  const int std_dev = std::stoi(argv[4]);

  io_t io = load(samFilename);
  if (!io.loaded) {
    return 1;
  }

  const int tid = bam_name2id(io.header, argv[5]);
  const int start = std::stoi(argv[6]);
  const int end = std::stoi(argv[7]);

  run_extract(io, read_length, mean_insert, std_dev, tid, start, end);

  bam_hdr_destroy(io.header);
  sam_close(io.sam);

  return 0;
}
