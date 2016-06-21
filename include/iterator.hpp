#include <functional>

#include "htslib/sam.h"

#include "io.hpp"

#ifndef ITERATOR_HPP
#define ITERATOR_HPP

class iterator {
private:
  samFile *m_sam;
  hts_itr_t *m_iter;

public:
  bam1_t *bam;

  iterator(const io_t io, const int tid, const int start, const int end)
      : m_sam(io.sam), bam(bam_init1()) {
    m_iter = sam_itr_queryi(io.idx, tid, start, end);
    if (m_iter == NULL) {
      std::cerr << "WARNING: SAM iterator is NULL!" << std::endl;
    }
  }

  iterator(const io_t io, const char *string)
      : m_sam(io.sam), bam(bam_init1()) {
    m_iter = sam_itr_querys(io.idx, io.header, string);
    if (m_iter == NULL) {
      std::cerr << "WARNING: SAM iterator is NULL!" << std::endl;
    }
  }

  // Copy constructor
  iterator(const iterator& other) {
    m_sam = other.m_sam;
    bam = bam_copy1(bam, other.bam);
    memcpy(m_iter, other.m_iter, sizeof(hts_itr_t));
  }

  // Move constructor
  iterator(iterator&& other) noexcept :
      m_sam(other.m_sam), m_iter(other.m_iter), bam(other.bam) {
    other.m_sam = NULL;
    other.m_iter = NULL;
    other.bam = NULL;
  }

  // Copy assignment operator
  iterator& operator=(const iterator& other) {
    iterator tmp(other);
    *this = std::move(tmp);
    return *this;
  }

  // Move assignment operator
  iterator& operator= (iterator&& other) noexcept {
    hts_itr_destroy(m_iter);
    bam_destroy1(bam);

    m_sam = other.m_sam;
    m_iter = other.m_iter;
    bam = other.bam;

    other.m_sam = NULL;
    other.m_iter = NULL;
    other.bam = NULL;

    return *this;
  }

  // Destructor
  ~iterator() {
    hts_itr_destroy(m_iter);
    bam_destroy1(bam);
  }

  inline bool next() {
    if (m_iter == NULL) return false;
    return (sam_itr_next(m_sam, m_iter, bam) >= 0);
  }

  void map(const std::function<void(bam1_t *)> &lambda) {
    #pragma omp parallel
    {
      bam1_t *tbam = bam_init1();

      int ret;
      while (true) {
        #pragma omp critical
        {
          ret = sam_itr_next(m_sam, m_iter, tbam);
        }

        if (ret < 0) {
          break;
        }

        lambda(tbam);
      }

      bam_destroy1(tbam);
    }
  }
};

#endif
