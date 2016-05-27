#include "htslib/sam.h"

#include "io.hpp"

#ifndef ITERATOR_HPP
#define ITERATOR_HPP

class iterator {
private:
  samFile *m_sam;

public:
  hts_itr_t *iter;
  bam1_t *bam;

  iterator(const io_t io, const int tid, const int start, const int end)
      : m_sam(io.sam), bam(bam_init1()) {
    iter = sam_itr_queryi(io.idx, tid, start, end);
    if (iter == NULL) {
      std::cerr << "ERROR: SAM iterator is NULL!" << std::endl;
    }
  }

  iterator(const io_t io, const char *string)
      : m_sam(io.sam), bam(bam_init1()) {
    iter = sam_itr_querys(io.idx, io.header, string);
    if (iter == NULL) {
      std::cerr << "ERROR: SAM iterator is NULL!" << std::endl;
    }
  }

  ~iterator() {
    hts_itr_destroy(iter);
    bam_destroy1(bam);
  }

  inline bool next() {
    if (iter == NULL) return false;
    return (sam_itr_next(m_sam, iter, bam) >= 0);
  }
};

#endif
