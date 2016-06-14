#include <cmath>
#include <cstring>

#include <string>
#include <iostream>
#include <vector>
#include <functional>

#include "htslib/sam.h"

#include "io.hpp"
#include "hash.hpp"
#include "extract.hpp"
#include "iterator.hpp"

void print_fasta(const bam1_t *bam, char* buffer) {
  const std::string sequence = convertToString(bam_get_seq(bam),
    bam->core.l_qseq, bam_is_rev(bam), buffer);

  // Skip all reads with Ns. These might add gaps
  // for (size_t i = 0; i < sequence.size(); i++) {
  //   if (sequence[i] == 'N' || sequence[i] == 'n') {
  //     return;
  //   }
  // }

  std::cout << '>' << bam_get_qname(bam) << '/'
    << (((bam->core.flag & BAM_FREAD1) != 0) ? '1' : '2') << '\n'
    << sequence << std::endl;
}

size_t process_region(const io_t io, const int tid, const int start,
    const int end, char* buffer) {
  size_t seqlen = 0;

  iterator iter(io, tid, start, end);
  while (iter.next()) {
    print_fasta(iter.bam, buffer);
    seqlen += strlen(buffer);
  }

  // iter.map([&buffer, &seqlen] (bam1_t *bam) -> void {
  //   print_fasta(bam, buffer);
  //   seqlen += strlen(buffer);
  // });

  return seqlen;
}

void process_mates(const io_t io, const int tid, const int start,
    const int end, bloom *bloom) {
  iterator iter(io, tid, start, end);
  while (iter.next()) {
    if ((iter.bam->core.flag & BAM_FMUNMAP) != 0) {
      bloom->push(iter.bam);
    }
  }

  // iter.map([&reads](bam1_t *bam) -> void {
  //   if ((bam->core.flag & BAM_FMUNMAP) != 0) {
  //     bloom->push(iter.bam);
  //   }
  // });
}

size_t find_mates(const io_t io, char *buffer, const bloom &bloom) {
  size_t seqlen = 0;

  iterator iter(io, ".");
  while (iter.next()) {
    if (bloom.in_alignments(iter.bam)) {
      print_fasta(iter.bam, buffer);
      seqlen += strlen(buffer);
    }
  }

  return seqlen;

  // iter.map([&alignments, &buffer](bam1_t *bam) -> void {
  //   if (bloom.in_alignments(bam)) {
  //     print_fasta(bam, buffer);
  //   }
  // });
}

void process_unmapped(const io_t io, char* buffer) {
  iterator iter(io, ".");
  while (iter.next()) {
    if ((iter.bam->core.flag & BAM_FUNMAP) != 0) {
      print_fasta(iter.bam, buffer);
    }
  }

  // iter.map([&buffer](bam1_t *bam) -> void { print_fasta(bam, buffer); });
}

void run_extract(const io_t io,
    const int read_length, const int mean_insert, const int std_dev,
    const int tid, const int start, const int end) {
  // Allocate memory for string conversions
  char *buffer = new char[read_length+1];

  size_t seqlen = 0;

  // Extract reads from the overlap
  seqlen += process_region(io, tid, start, end, buffer);

  bloom bloom;

  // Extract pairs from the left mappings
  {
    const int left_start = start - (mean_insert + (3*std_dev) + 2*read_length);
    const int left_end = end - (mean_insert - (3*std_dev) + read_length);
    process_mates(io, tid, left_start, left_end, &bloom);
    seqlen += find_mates(io, buffer, bloom);
  }

  bloom.clear();

  // Extract pairs from the right mappings
  {
    const int right_start = start + (mean_insert + (3*std_dev) + read_length);
    const int right_end = end + (mean_insert - (3*std_dev) + read_length);
    process_mates(io, tid, right_start, right_end, &bloom);
    seqlen += find_mates(io, buffer, bloom);
  }

  // TODO: Make this not hard-coded
  if (seqlen / (end - start) < 25) {
    process_unmapped(io, buffer);
  }

  delete[] buffer;
}
