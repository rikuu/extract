// Copyright 2016 Riku Walve

#include <vector>

#include "htslib/sam.h"

#ifndef EXTRACT_HASH_HPP
#define EXTRACT_HASH_HPP

#define BLOOM_SIZE 1024

// Hashes a string
size_t hash_str(const char *);

// Hash functions for alignments
inline size_t hash_alignment1(const bam1_t *bam) {
  return hash_str(bam_get_qname(bam)) +
    ((((bam->core.flag & BAM_FREAD1) != 0) ? 1 : 2) << 8);
}

inline size_t hash_alignment2(const bam1_t *bam) {
  return hash_str(bam_get_qname(bam)) +
    ((((bam->core.flag & BAM_FREAD2) != 0) ? 1 : 2) << 8);
}

// TODO: Use multiple hash functions to reduce collisions
class bloom {
private:
  std::vector<bool> m_vector;

public:
  bloom() {
    m_vector.resize(BLOOM_SIZE);
  }

  // Clears the filter
  inline void clear() {
    m_vector.assign(BLOOM_SIZE, false);
  }

  // Adds an alignment to the filter
  inline void push(const bam1_t *bam) {
    const size_t hash = hash_alignment1(bam);
    m_vector[hash % BLOOM_SIZE] = true;
  }

  // Checks if an alignment is contained
  inline bool contains(const bam1_t *bam) const {
    const size_t hash = hash_alignment1(bam);
    return m_vector[hash % BLOOM_SIZE];
  }

  // Checks if the mate of an alignment is contained
  inline bool contains_mate(const bam1_t *bam) const {
    const size_t hash = hash_alignment2(bam);
    return m_vector[hash % BLOOM_SIZE];
  }
};

#endif
