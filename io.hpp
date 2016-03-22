// Copyright 2016 Riku Walve

#include <string>

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

std::string convertToString(const uint8_t *, const int32_t, const bool);

// These could be blocks, though
inline uint8_t complement(const uint8_t n) {
  switch (n) {
    case 1: return 8; break;
    case 2: return 4; break;
    case 4: return 2; break;
    case 8: return 1; break;
    case 15:
    default: return 15; break;
  }
}

inline uint8_t querySequence(const uint8_t *query, const int32_t length,
    const int32_t index, const bool reverse) {
  if (!reverse)
    return bam_seqi(query, index);

  return complement(bam_seqi(query, length - 1 - index));
}

#endif
