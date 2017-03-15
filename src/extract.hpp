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

#ifndef EXTRACT_HPP
#define EXTRACT_HPP

#include <sys/types.h>

#include <string>
#include <vector>
#include <functional>

// GATB-core Bloom filter requires hash1 function for items
// TODO: Hash bam object directly
inline u_int64_t hash1(const std::string &key, u_int64_t seed=0) {
  return std::hash<std::string>{}(key);
}

#include <gatb/gatb_core.hpp>

#include "htslib/sam.h"

#include "io.hpp"

uint64_t count_reads(const std::string &);

class Extract : public Tool {
 public:
    Extract();
    void execute();

    // Prints an alignment in fasta format
    void print_fasta(const bam1_t*, char*, BankFasta*);

    // Prints all alignments in a region
    void process_region(const io_t&, const int, const int, const int, char*, IBloom<std::string>*, BankFasta*, int*, int*);
    void process_mates(const io_t&, const int, const int, const int, IBloom<std::string>*);
    void find_mates(const io_t&, char*, IBloom<std::string>*, BankFasta*, int*, int*);
    void process_unmapped(const io_t&, char*, IBloom<std::string>*, BankFasta*, int*);
};

#endif
