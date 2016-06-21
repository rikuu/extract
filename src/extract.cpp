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
    const int end, char* buffer, bloom_filter *bloom) {
  size_t seqlen = 0;

  iterator iter(io, tid, start, end);
  while (iter.next()) {
    if (!bloom->contains(iter.bam)) {
      print_fasta(iter.bam, buffer);
      seqlen += strlen(buffer);
    }
  }

  // iter.map([&buffer, &seqlen, &bloom] (bam1_t *bam) -> void {
  //   if (!bloom.contains(iter.bam)) {
  //     print_fasta(bam, buffer);
  //     seqlen += strlen(buffer);
  //   }
  // });

  return seqlen;
}

void process_mates(const io_t io, const int tid, const int start,
    const int end, bloom_filter *bloom) {
  iterator iter(io, tid, start, end);
  while (iter.next()) {
    if ((iter.bam->core.flag & BAM_FMUNMAP) != 0) {
      bloom->push(iter.bam);
    }
  }

  // iter.map([&bloom](bam1_t *bam) -> void {
  //   if ((bam->core.flag & BAM_FMUNMAP) != 0) {
  //     bloom->push(iter.bam);
  //   }
  // });
}

size_t find_mates(const io_t io, char *buffer, bloom_filter *bloom) {
  size_t seqlen = 0;

  iterator iter(io, ".");
  while (iter.next()) {
    if (bloom->contains_mate(iter.bam)) {
      print_fasta(iter.bam, buffer);
      seqlen += strlen(buffer);
    }
  }

  return seqlen;

  // iter.map([&buffer, &bloom](bam1_t *bam) -> void {
  //   if (bloom.in_alignments(bam)) {
  //     print_fasta(bam, buffer);
  //   }
  // });
}

void process_unmapped(const io_t io, char* buffer, bloom_filter *bloom) {
  iterator iter(io, ".");
  while (iter.next()) {
    if ((iter.bam->core.flag & BAM_FUNMAP) != 0 &&
          !bloom->contains(iter.bam)) {
      print_fasta(iter.bam, buffer);
    }
  }

  // iter.map([&buffer, &bloom](bam1_t *bam) -> void {
  //   if ((bam->core.flag & BAM_FUNMAP) != 0 && !bloom.contains(iter.bam)) {
  //     print_fasta(bam, buffer);
  //   }
  // });
}

void run_extract(const io_t io,
    const int read_length, const int mean_insert, const int std_dev,
    const int tid, const int start, const int end,
    const bool exact, const bool unmapped, const int threshold) {
  bloom_filter *bloom = NULL;
  if (exact) {
    bloom_exact *bloom_ = new bloom_exact();
    bloom = bloom_;
  } else {
    bloom_hash *bloom_ = new bloom_hash();
    bloom = bloom_;
  }

  // Extract pairs from the left mappings
  {
    const int left_start = start - (mean_insert + (3*std_dev) + 2*read_length);
    const int left_end = end - (mean_insert - (3*std_dev) + read_length);
    process_mates(io, tid, left_start, left_end, bloom);
  }

  // Extract pairs from the right mappings
  {
    const int right_start = start + (mean_insert + (3*std_dev) + read_length);
    const int right_end = end + (mean_insert - (3*std_dev) + read_length);
    process_mates(io, tid, right_start, right_end, bloom);
  }

  // Allocate memory for string conversions
  char *buffer = new char[read_length+1];

  const size_t seqlen =
      process_region(io, tid, start, end, buffer, bloom) +
      find_mates(io, buffer, bloom);

  // TODO: Make threshold not hard-coded
  if (unmapped && (seqlen / (end - start) < (unsigned) threshold)) {
    process_unmapped(io, buffer, bloom);
  }

  delete[] buffer;
  delete bloom;
}
