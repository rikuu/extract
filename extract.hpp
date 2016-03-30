#include <vector>

#include "htslib/sam.h"

#include "io.hpp"

#ifndef EXTRACT_HPP
#define EXTRACT_HPP

void print_fasta(const bam1_t*);

void process_region(const io_t, const int, const int, const int);
std::vector<size_t> process_mates(const io_t, const int, const int, const int);
void find_mates(const io_t, const std::vector<size_t> &);

void run_extract(const io_t, const int, const int, const int, const int,
  const int, const int);

#endif
