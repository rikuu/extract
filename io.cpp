// Copyright 2016 Riku Walve

#include <string>
#include <iostream>
#include <vector>

#include <sstream>
#include <fstream>

#include <cstring>

#include "htslib/sam.h"

#include "io.hpp"

// Loads a bam/sam file into an IO object
io_t load(const std::string &samFilename) {
  io_t io;

  io.sam = sam_open(samFilename.c_str(), "r");
  if (io.sam == NULL) {
    std::cerr << "ERROR: SAM file not found!" << std::endl;
    return io;
  }

  io.header = sam_hdr_read(io.sam);
  if (io.header == NULL) {
    std::cerr << "ERROR: SAM header error!" << std::endl;
    return io;
  }

  io.idx = sam_index_load(io.sam, samFilename.c_str());
  if (io.idx == NULL) {
    std::cerr << "ERROR: SAM index not found!" << std::endl;
    return io;
  }

  io.loaded = true;
  return io;
}

// Converts an alignment to std::string. Handles reverse complements.
std::string convertToString(const uint8_t *query, const int32_t length,
    const bool reverse) {
  char* string = new char[length];

  for (int i = 0; i < length; i++) {
    switch (querySequence(query, length, i, reverse)) {
      case 0x1:
        string[i] = 'A';
        break;
      case 0x2:
        string[i] = 'C';
        break;
      case 0x4:
        string[i] = 'G';
        break;
      case 0x8:
        string[i] = 'T';
        break;
      case 0x15:
      default:
        string[i] = 'N';
        break;
    }
  }

  string[length] = '\0';
  std::string result(string);
  delete string;

  return result;
}
