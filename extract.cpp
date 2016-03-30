#include <cmath>

#include <string>
#include <iostream>
#include <vector>

#include "htslib/sam.h"

#include "io.hpp"
#include "hash.hpp"
#include "extract.hpp"

void print_fasta(const bam1_t *bam) {
  const std::string sequence = convertToString(bam_get_seq(bam),
    bam->core.l_qseq, bam_is_rev(bam));

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
    const int end) {
  hts_itr_t *iter = sam_itr_queryi(io.idx, tid, start, end);
  if (iter == NULL) {
    std::cerr << "ERROR: SAM iterator is NULL!" << std::endl;
    return;
  }

  bam1_t *bam = bam_init1();
  while (sam_itr_next(io.sam, iter, bam) >= 0) {
    print_fasta(bam);
  }

  hts_itr_destroy(iter);
  bam_destroy1(bam);
}

std::vector<size_t> process_mates(const io_t io, const int tid, const int start,
    const int end) {
  std::vector<size_t> reads;

  hts_itr_t *iter = sam_itr_queryi(io.idx, tid, start, end);
  if (iter == NULL) {
    std::cerr << "ERROR: SAM iterator is NULL!" << std::endl;
    return reads;
  }

  bam1_t *bam = bam_init1();
  while (sam_itr_next(io.sam, iter, bam) >= 0) {
    if ((bam->core.flag & BAM_FMUNMAP) != 0) {
      reads.push_back(hash_alignment1(bam));
    }
  }

  hts_itr_destroy(iter);
  bam_destroy1(bam);

  return reads;
}

void find_mates(const io_t io, const std::vector<size_t> &alignments) {
  hts_itr_t *iter = sam_itr_querys(io.idx, io.header, ".");
  if (iter == NULL) {
    std::cerr << "ERROR: SAM iterator is NULL!" << std::endl;
    return;
  }

  bam1_t *bam = bam_init1();
  while (sam_itr_next(io.sam, iter, bam) >= 0) {
    if (in_alignments(alignments, hash_alignment2(bam))) {
      print_fasta(bam);
    }
  }

  hts_itr_destroy(iter);
  bam_destroy1(bam);
}

void run_extract(const io_t io, const int read_length, const int mean_insert,
    const int std_dev, const int tid, const int start, const int end) {
  // Extract reads from the overlap
  process_region(io, tid, start, end);

  // Extract pairs from the left mappings
  {
    const int left_start = start - (int) (ceilf(mean_insert + (1.96f*std_dev)) + 2*read_length);
    const int left_end = end - (int) (floorf(mean_insert - (1.96f*std_dev)) + read_length);
    std::vector<size_t> left = process_mates(io, tid, left_start, left_end);
    find_mates(io, left);
  }

  // Extract pairs from the right mappings
  {
    const int right_start = start + (int) (ceilf(mean_insert + (1.96f*std_dev)) + read_length);
    const int right_end = end + (int) (floorf(mean_insert - (1.96f*std_dev)) + read_length);
    std::vector<size_t> right = process_mates(io, tid, right_start, right_end);
    find_mates(io, right);
  }
}
