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
    return;
  }

  this->header = sam_hdr_read(this->sam);
  if (this->header == NULL) {
    return;
  }

  this->idx = sam_index_load(this->sam, samFilename.c_str());
  if (this->idx == NULL) {
    return;
  }
}

// Converts an alignment to std::string. Handles reverse complements.
std::string convertToString(const uint8_t *query, const int32_t length,
    const bool reverse, char* buffer) {
  for (int32_t i = 0; i < length; i++) {
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
