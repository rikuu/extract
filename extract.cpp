#include <string>
#include <iostream>
#include <vector>

#include <sstream>
#include <fstream>

#include <cstring>

#include "htslib/sam.h"

#include "io.hpp"
#include "hash.hpp"

void process_region(const io_t, const int, const int, const int);
std::vector<size_t> process_mates(const io_t, const int, const int, const int);
void find_mates(const io_t, std::vector<size_t> &);

void process_region(const io_t io, const int tid, const int start,
    const int end) {
  hts_itr_t *iter = sam_itr_queryi(io.idx, tid, start, end);
  if (iter == NULL) {
    std::cerr << "ERROR: SAM iterator is NULL!" << std::endl;
    return;
  }

	bam1_t *bam = bam_init1();
	while (sam_itr_next(io.sam, iter, bam) >= 0) {
    std::cout << '>' << bam_get_qname(bam) << '/' <<
      (((bam->core.flag & BAM_FREAD1) != 0) ? '1' : '2') << '\n' <<
      convertToString(bam_get_seq(bam), bam->core.l_qseq, bam_is_rev(bam)) << std::endl;
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
      /*std::cout << '>' << bam_get_qname(bam) << '/' <<
        (((bam->core.flag & BAM_FREAD1) != 0) ? '1' : '2') << std::endl;*/
      reads.push_back(hash_alignment1(bam));
    }
  }

	hts_itr_destroy(iter);
	bam_destroy1(bam);

  return reads;
}

void find_mates(const io_t io, std::vector<size_t> &alignments) {
  hts_itr_t *iter = sam_itr_querys(io.idx, io.header, "."); //"*");
  if (iter == NULL) {
    std::cerr << "ERROR: SAM iterator is NULL!" << std::endl;
    return;
  }

	bam1_t *bam = bam_init1();
	while (sam_itr_next(io.sam, iter, bam) >= 0) {
    if (in_alignments(alignments, hash_alignment2(bam))) {
      std::cout << '>' << bam_get_qname(bam) << '/'
        << (((bam->core.flag & BAM_FREAD1) != 0) ? '1' : '2') << '\n'
        << convertToString(bam_get_seq(bam), bam->core.l_qseq, bam_is_rev(bam))
        << std::endl;
    }
  }

	hts_itr_destroy(iter);
	bam_destroy1(bam);
}

int main(int argc, char* argv[]) {
  if (argc != 8) {
    std::cerr << argv[0] << " <alignments>.bam <read length> <mu> <sd>"
      << " <scaffold name> <gap start> <gap end>" << std::endl;
    return -1;
  }

  const std::string samFilename = argv[1];
  const int read_length = std::stoi(argv[2]);
  const int mean_insert = std::stoi(argv[3]);
  const int std_dev = std::stoi(argv[4]);

  io_t io = load(samFilename);
  if (!io.loaded) {
    return 1;
  }

  const int tid = bam_name2id(io.header, argv[5]);
  const int start = std::stoi(argv[6]);
  const int end = std::stoi(argv[7]);

  process_region(io, tid, start, end);

  {
    const int left_start = start - ((mean_insert + 2*std_dev) + 2*read_length);
    const int left_end = end - ((mean_insert - 2*std_dev) + read_length);
    std::vector<size_t> left = process_mates(io, tid, left_start, left_end);
    find_mates(io, left);
  }

  {
    const int right_start = start + ((mean_insert + 2*std_dev) + read_length);
    const int right_end = end + ((mean_insert - 2*std_dev) + read_length);
    std::vector<size_t> right = process_mates(io, tid, right_start, right_end);
    find_mates(io, right);
  }

  bam_hdr_destroy(io.header);
	sam_close(io.sam);

  return 0;
}
