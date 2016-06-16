#include <vector>

#include "htslib/sam.h"

#include "io.hpp"

#ifndef EXTRACT_HPP
#define EXTRACT_HPP

// Prints an alignment in fasta format
void print_fasta(const bam1_t*, char*);

// Prints all alignments in a region
size_t process_region(const io_t, const int, const int, const int, char*,
  const bloom&);

void process_mates(const io_t, const int, const int, const int, bloom*);

size_t find_mates(const io_t, char*, const bloom&);

void process_unmapped(const io_t, char*, const bloom&);

void run_extract(const io_t,
    const int, const int, const int,
    const int, const int, const int);

#endif
