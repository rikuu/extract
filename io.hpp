#include <string>
#include <iostream>
#include <vector>

#include <sstream>
#include <fstream>

#include <cstring>

#include "htslib/sam.h"

#ifndef EXTRACT_IO_HPP
#define EXTRACT_IO_HPP

typedef struct {
  samFile *sam = NULL;
  bam_hdr_t *header = NULL;
  hts_idx_t *idx = NULL;
  bool loaded = false;
} io_t;

io_t load(const std::string &);

inline uint8_t complement(const uint8_t);

inline uint8_t querySequence(const uint8_t *, const int32_t, const int32_t,
  const bool);

std::string convertToString(const uint8_t *, const int32_t, const bool);

#endif
