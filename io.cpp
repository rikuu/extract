// Copyright 2016 Riku Walve

#include <cstring>

#include <string>
#include <iostream>
#include <vector>

#include "htslib/sam.h"

#include "io.hpp"

// Loads a bam/sam file into an IO object
io_t::io_t(const std::string &samFilename) {
  this->sam = sam_open(samFilename.c_str(), "r");
  if (this->sam == NULL) {
    // std::cerr << "ERROR: SAM file not found!" << std::endl;
    return;
  }

  this->header = sam_hdr_read(this->sam);
  if (this->header == NULL) {
    // std::cerr << "ERROR: SAM header error!" << std::endl;
    return;
  }

  this->idx = sam_index_load(this->sam, samFilename.c_str());
  if (this->idx == NULL) {
    // std::cerr << "ERROR: SAM index not found!" << std::endl;
    return;
  }

  this->loaded = true;
}

// Converts an alignment to std::string. Handles reverse complements.
std::string convertToString(const uint8_t *query, const int32_t length,
    const bool reverse, char* buffer) {
  for (int i = 0; i < length; i++) {
    // TODO: Improve cache coherence by reversing blocks
    switch (querySequence(query, length, i, reverse)) {
      case 0x1:
        buffer[i] = 'A';
        break;
      case 0x2:
        buffer[i] = 'C';
        break;
      case 0x4:
        buffer[i] = 'G';
        break;
      case 0x8:
        buffer[i] = 'T';
        break;
      case 0x15:
      default:
        buffer[i] = 'N';
        break;
    }
  }

  buffer[length] = '\0';
  return std::string(buffer);
}
