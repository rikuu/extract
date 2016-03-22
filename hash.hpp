// Copyright 2016 Riku Walve

#include <vector>

#include "htslib/sam.h"

#ifndef EXTRACT_HASH_HPP
#define EXTRACT_HASH_HPP

// Hashes a string
size_t hash_str(const char *);

// Checks if an alignment is in a vector.
// This is basically a slower bloom filter.
inline bool in_alignments(const std::vector<size_t> &alignments,
    const size_t hash) {
  for (size_t i = 0; i < alignments.size(); i++) {
    if (alignments[i] == hash) {
      return true;
    }
  }

  return false;
}

// Hash functions for alignments
inline size_t hash_alignment1(const bam1_t *bam) {
  return hash_str(bam_get_qname(bam)) +
    ((((bam->core.flag & BAM_FREAD1) != 0) ? 1 : 2) << 8);
}

inline size_t hash_alignment2(const bam1_t *bam) {
  return hash_str(bam_get_qname(bam)) +
    ((((bam->core.flag & BAM_FREAD2) != 0) ? 1 : 2) << 8);
}

#endif
