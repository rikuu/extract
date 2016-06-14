#include <vector>

#include "htslib/sam.h"

#include "io.hpp"

#ifndef EXTRACT_HPP
#define EXTRACT_HPP

void print_fasta(const bam1_t*, char*);

size_t process_region(const io_t, const int, const int, const int, char*);
void process_mates(const io_t, const int, const int, const int, bloom*);
size_t find_mates(const io_t, char*, const bloom&);
void process_unmapped(const io_t, char*);

void run_extract(const io_t,
    const int, const int, const int,
    const int, const int, const int);

#endif
