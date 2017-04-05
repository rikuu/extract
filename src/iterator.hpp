/*****************************************************************************
 *  Extract
 *  Copyright (C) Riku Walve 2017
 *
 *  Contact: riku.walve@cs.helsinki.fi
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#include <functional>

#include "htslib/sam.h"

#include "io.hpp"

#ifndef ITERATOR_HPP
#define ITERATOR_HPP

class sam_iterator {
private:
  samFile *m_sam;
  hts_itr_t *m_iter;

public:
  bam1_t *bam;

  sam_iterator(const io_t io, const int tid, const int start, const int end)
      : m_sam(io.sam), bam(bam_init1()) {
    m_iter = sam_itr_queryi(io.idx, tid, start, end);
    if (m_iter == nullptr) {
      std::cerr << "WARNING: SAM iterator is NULL!" << std::endl;
    }
  }

  sam_iterator(const io_t io, const char *string)
      : m_sam(io.sam), bam(bam_init1()) {
    m_iter = sam_itr_querys(io.idx, io.header, string);
    if (m_iter == nullptr) {
      std::cerr << "WARNING: SAM iterator is NULL!" << std::endl;
    }
  }

  // Copy constructor
  sam_iterator(const sam_iterator& other) {
    m_sam = other.m_sam;
    bam = bam_copy1(bam, other.bam);
    memcpy(m_iter, other.m_iter, sizeof(hts_itr_t));
  }

  // Move constructor
  sam_iterator(sam_iterator&& other) noexcept :
      m_sam(other.m_sam), m_iter(other.m_iter), bam(other.bam) {
    other.m_sam = nullptr;
    other.m_iter = nullptr;
    other.bam = nullptr;
  }

  // Copy assignment operator
  sam_iterator& operator=(const sam_iterator& other) {
    sam_iterator tmp(other);
    *this = std::move(tmp);
    return *this;
  }

  // Move assignment operator
  sam_iterator& operator= (sam_iterator&& other) noexcept {
    hts_itr_destroy(m_iter);
    bam_destroy1(bam);

    m_sam = other.m_sam;
    m_iter = other.m_iter;
    bam = other.bam;

    other.m_sam = nullptr;
    other.m_iter = nullptr;
    other.bam = nullptr;

    return *this;
  }

  // Destructor
  ~sam_iterator() {
    hts_itr_destroy(m_iter);
    bam_destroy1(bam);
  }

  inline bool next() {
    if (m_iter == nullptr) return false;
    return (sam_itr_next(m_sam, m_iter, bam) >= 0);
  }
};

#endif
