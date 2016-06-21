#include <vector>

#include "htslib/sam.h"

#include "io.hpp"

#ifndef EXTRACT_HPP
#define EXTRACT_HPP

// Prints an alignment in fasta format
void print_fasta(const bam1_t*, char*);

// Prints all alignments in a region
size_t process_region(const io_t, const int, const int, const int, char*,
  bloom_filter*);

void process_mates(const io_t, const int, const int, const int, bloom_filter*);

size_t find_mates(const io_t, char*, bloom_filter*);

void process_unmapped(const io_t, char*, bloom_filter*);

void run_extract(const io_t,
    const int, const int, const int,
    const int, const int, const int,
    const bool, const bool, const int);

#endif
