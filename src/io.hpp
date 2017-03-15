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

#include <string>

#include "htslib/sam.h"

#ifndef EXTRACT_IO_HPP
#define EXTRACT_IO_HPP

class io_t {
 public:
  samFile *sam = NULL;
  bam_hdr_t *header = NULL;
  hts_idx_t *idx = NULL;
  bool loaded = false;

  io_t(const std::string &);

  void unload() {
    bam_hdr_destroy(this->header);
    sam_close(this->sam);
    loaded = false;
  }
};

std::string convertToString(const uint8_t *, const int32_t, const bool, char*);

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
