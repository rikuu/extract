#include <vector>

#include "htslib/sam.h"

#ifndef EXTRACT_HASH_HPP
#define EXTRACT_HASH_HPP

size_t hash_str(char *);

inline bool in_alignments(std::vector<size_t> &alignments,
    const size_t hash) {
  for (size_t i = 0; i < alignments.size(); i++) {
    if (alignments[i] == hash) {
      return true;
    }
  }

  return false;
}

inline size_t hash_alignment1(const bam1_t *bam) {
  return hash_str(bam_get_qname(bam)) +
    ((((bam->core.flag & BAM_FREAD1) != 0) ? 1 : 2) << 8);
}

inline size_t hash_alignment2(const bam1_t *bam) {
  return hash_str(bam_get_qname(bam)) +
    ((((bam->core.flag & BAM_FREAD2) != 0) ? 1 : 2) << 8);
}

#endif
