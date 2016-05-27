#include <cmath>

#include <string>
#include <iostream>
#include <vector>

#include "htslib/sam.h"

#include "io.hpp"
#include "hash.hpp"
#include "extract.hpp"
#include "iterator.hpp"

void print_fasta(const bam1_t *bam, char* buffer) {
  const std::string sequence = convertToString(bam_get_seq(bam),
    bam->core.l_qseq, bam_is_rev(bam), buffer);

  // Skip all reads with Ns. These might add gaps
  /*for (size_t i = 0; i < sequence.size(); i++) {
    if (sequence[i] == 'N' || sequence[i] == 'n') {
      return;
    }
  }*/

  std::cout << '>' << bam_get_qname(bam) << '/'
    << (((bam->core.flag & BAM_FREAD1) != 0) ? '1' : '2') << '\n'
    << sequence << std::endl;
}

void process_region(const io_t io, const int tid, const int start,
    const int end, char* buffer) {
  iterator iter(io, tid, start, end);
  while (iter.next()) {
    print_fasta(iter.bam, buffer);
  }
}

std::vector<size_t> process_mates(const io_t io, const int tid, const int start,
    const int end) {
  std::vector<size_t> reads;

  iterator iter(io, tid, start, end);
  while (iter.next()) {
    if ((iter.bam->core.flag & BAM_FMUNMAP) != 0) {
      reads.push_back(hash_alignment1(iter.bam));
    }
  }

  return reads;
}

void find_mates(const io_t io, const std::vector<size_t> &alignments,
    char *buffer) {
  iterator iter(io, ".");
  while (iter.next()) {
    if (in_alignments(alignments, hash_alignment2(iter.bam))) {
      print_fasta(iter.bam, buffer);
    }
  }
}

void run_extract(const io_t pe1_io, const io_t pe2_io,
    const int read_length, const int mean_insert, const int std_dev,
    const int tid, const int start, const int end) {
  // Allocate memory for string conversions
  char *buffer = new char[read_length+1];

  // Extract reads from the overlap
  process_region(pe1_io, tid, start, end, buffer);
  process_region(pe2_io, tid, start, end, buffer);

  // Extract pairs from the left mappings
  {
    const int left_start = start - (int) (ceilf(mean_insert + (1.96f*std_dev)) + 2*read_length);
    const int left_end = end - (int) (floorf(mean_insert - (1.96f*std_dev)) + read_length);
    std::vector<size_t> left = process_mates(pe1_io, tid, left_start, left_end);
    find_mates(pe2_io, left, buffer);

    // Figure this out
    left = process_mates(pe2_io, tid, left_start, left_end);
    find_mates(pe1_io, left, buffer);
  }

  // Extract pairs from the right mappings
  {
    const int right_start = start + (int) (ceilf(mean_insert + (1.96f*std_dev)) + read_length);
    const int right_end = end + (int) (floorf(mean_insert - (1.96f*std_dev)) + read_length);
    std::vector<size_t> right = process_mates(pe2_io, tid, right_start, right_end);
    find_mates(pe1_io, right, buffer);

    // Figure this out
    right = process_mates(pe1_io, tid, right_start, right_end);
    find_mates(pe2_io, right, buffer);
  }

  delete[] buffer;
}

void run_extract(const io_t io,
    const int read_length, const int mean_insert, const int std_dev,
    const int tid, const int start, const int end) {
  // Allocate memory for string conversions
  char *buffer = new char[read_length+1];

  // Extract reads from the overlap
  process_region(io, tid, start, end, buffer);

  // Extract pairs from the left mappings
  {
    const int left_start = start - (int) (ceilf(mean_insert + (1.96f*std_dev)) + 2*read_length);
    const int left_end = end - (int) (floorf(mean_insert - (1.96f*std_dev)) + read_length);
    std::vector<size_t> left = process_mates(io, tid, left_start, left_end);
    find_mates(io, left, buffer);
  }

  // Extract pairs from the right mappings
  {
    const int right_start = start + (int) (ceilf(mean_insert + (1.96f*std_dev)) + read_length);
    const int right_end = end + (int) (floorf(mean_insert - (1.96f*std_dev)) + read_length);
    std::vector<size_t> right = process_mates(io, tid, right_start, right_end);
    find_mates(io, right, buffer);
  }

  delete[] buffer;
}
