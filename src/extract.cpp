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

#include <cmath>
#include <cstring>

#include <string>
#include <iostream>
#include <vector>

#include "htslib/sam.h"

#include "extract.hpp"
#include "io.hpp"
#include "iterator.hpp"

#include <gatb/gatb_core.hpp>

/*****************************************************************************/

static const char* STR_ALIGNMENT = "-bam";
static const char* STR_OUTPUT = "-reads";

static const char* STR_READ_LENGTH = "-read-length";
static const char* STR_MEAN = "-mean";
static const char* STR_STD_DEV = "-std-dev";

static const char* STR_REGION = "-region";
static const char* STR_INSERTION = "-insertion";

static const char* STR_THRESHOLD = "-unmapped";
static const char* STR_ONLY_UNMAPPED = "-unmapped-only";
static const char* STR_NO_OVERLAP = "-no-overlap";

/*****************************************************************************/

Extract::Extract() : Tool("Extract") {
  // Input / output
  getParser()->push_front(new OptionOneParam(STR_ALIGNMENT, "Aligned BAM file", true));
  getParser()->push_front(new OptionOneParam(STR_OUTPUT, "FASTA-formatted output file", true));

  // Read library parameters
  getParser()->push_front(new OptionOneParam(STR_READ_LENGTH, "Read length", true));
  getParser()->push_front(new OptionOneParam(STR_MEAN, "Mean insert size", true));
  getParser()->push_front(new OptionOneParam(STR_STD_DEV, "Insert size standard deviation", true));

  // Gap parameters
  getParser()->push_front(new OptionOneParam(STR_REGION, "Region for filtering", true));
  getParser()->push_front(new OptionNoParam(STR_INSERTION, "Insertion filtering"));

  // Unmapped read parameters
  getParser()->push_front(new OptionOneParam(STR_THRESHOLD, "Threshold for using unmapped reads", false, "0"));
  getParser()->push_front(new OptionNoParam(STR_ONLY_UNMAPPED, "Only output unmapped reads"));
  getParser()->push_front(new OptionNoParam(STR_NO_OVERLAP, "Don't output overlapping reads"));
}

/*****************************************************************************/

static inline const std::string bam2string(const bam1_t *bam) {
  return std::string(bam_get_qname(bam)) + (((bam->core.flag & BAM_FREAD1) != 0) ? "/1" : "/2");
}

static inline const std::string bam2string_mate(const bam1_t *bam) {
  return std::string(bam_get_qname(bam)) + (((bam->core.flag & BAM_FREAD1) != 0) ? "/2" : "/1");
}

/*****************************************************************************/

// Output a bam object into GATB fasta bank
void Extract::print_fasta(const bam1_t *bam, char* buffer, BankFasta *bank) {
  const std::string sequence = convertToString(bam_get_seq(bam),
    bam->core.l_qseq, bam_is_rev(bam), buffer);

  Sequence seq(buffer);
  seq._comment = bam2string(bam);
  bank->insert(seq);
}

// Output reads that map to a region in a scaffold and not in the Bloom filter.
void Extract::process_region(const io_t &io, const int tid, const int start,
    const int end, char* buffer, IBloom<std::string> *bloom, BankFasta *bank,
    int *seqlen, int *num_of_reads) {
  sam_iterator iter(io, tid, start, end);
  while (iter.next()) {
    if (!bloom->contains(bam2string(iter.bam))) {
      print_fasta(iter.bam, buffer, bank);
      *seqlen += strlen(buffer);
      (*num_of_reads)++;
    }
  }
}

// Add all reads in a region to a Bloom filter
void Extract::process_mates(const io_t &io, const int tid, const int start,
    const int end, IBloom<std::string> *bloom) {
  sam_iterator iter(io, tid, start, end);
  while (iter.next()) {
    bloom->insert(bam2string(iter.bam));
  }
}

// Output all reads in a Bloom filter
void Extract::find_mates(const io_t &io, char *buffer, IBloom<std::string> *bloom,
    BankFasta *bank, int *seqlen, int *num_of_reads) {
  sam_iterator iter(io, ".");
  while (iter.next()) {
    if (bloom->contains(bam2string_mate(iter.bam))) {
      print_fasta(iter.bam, buffer, bank);
      *seqlen += strlen(buffer);
      (*num_of_reads)++;
    }
  }
}

// Output all unmapped reads
void Extract::process_unmapped(const io_t &io, char* buffer,
    IBloom<std::string> *bloom, BankFasta *bank, int* num_of_reads) {
  // Go through entire BAM file rather than the unmapped reads only.
  // BWA gives unmapped reads with a mapped mate a position, essentially making
  // it a mapped read. Outputting these reads is debateble.
  sam_iterator iter(io, ".");
  while (iter.next()) {
    if ((iter.bam->core.flag & BAM_FUNMAP) != 0 &&
          !bloom->contains(bam2string(iter.bam))) {
      print_fasta(iter.bam, buffer, bank);
      (*num_of_reads)++;
    }
  }
}

// Counts the number of reads by iterating through the alignment file
uint64_t count_reads(const std::string &filename) {
  io_t io(filename);
  sam_iterator iter(io, ".");

  uint64_t count = 0;
  while (iter.next()) {
    count++;
  }

  io.unload();

  return count;
}

// Execute read extraction
void Extract::execute() {
  const std::string alignment = getInput()->getStr(STR_ALIGNMENT);
  const std::string output = getInput()->getStr(STR_OUTPUT);

  // TODO: Infer read length (assume all reads are same length or nah?)
  const int read_length = static_cast<int>(getInput()->getInt(STR_READ_LENGTH));
  const int mean_insert = static_cast<int>(getInput()->getInt(STR_MEAN));
  const int std_dev = static_cast<int>(getInput()->getInt(STR_STD_DEV));

  const std::string region = getInput()->getStr(STR_REGION);
  const bool insertion = getParser()->saw(STR_INSERTION);

  const int threshold = static_cast<int>(getInput()->getInt(STR_THRESHOLD));
  const bool unmapped_only = getParser()->saw(STR_ONLY_UNMAPPED);
  const bool no_overlap = getParser()->saw(STR_NO_OVERLAP);

  // Parse samtools style region
  const size_t colon = region.find(":");
  const size_t dash = region.find("-");
  if (colon == std::string::npos || dash == std::string::npos) {
    std::cerr << "Malformed region, should be CONTIG:START-END" << std::endl;
    exit(1);
  }

  const std::string scaffold = region.substr(0, colon);
  const int gap_start = std::stoi(region.substr(colon+1, dash-colon));
  const int gap_end = std::stoi(region.substr(dash+1));
  const int length = gap_end - gap_start;

  // TODO: Move start by flank length in insertion?
  const int start = gap_start;
  const int end = insertion ? gap_start : gap_end;

  // Load alignment file
  io_t io(alignment);
  if (!io.loaded()) {
    std::cerr << "Error loading alignments" << std::endl;
    exit(1);
  }

  // Allocate memory for string conversions
  char *buffer = new char[read_length+1];

  // Use basic Bloom filter from GATB
  const uint64_t num_of_reads = count_reads(alignment);
  IBloom<std::string> *bloom = new BloomSynchronized<std::string>(5 * num_of_reads);

  // Open output file
  BankFasta reads(output);
  int seqlen = 0, reads_extracted = 0, unmapped_extracted = 0;

  if (!unmapped_only) {
    // Compute scaffold id from scaffold name
    const int tid = bam_name2id(io.header, scaffold.c_str());

    // Extract pairs from the left mappings
    const int left_start = start - (mean_insert + 3*std_dev + 2*read_length);
    const int left_end = end - (mean_insert - 3*std_dev + read_length);
    process_mates(io, tid, left_start, left_end, bloom);

    // Extract pairs from the right mappings
    const int right_start = start + (mean_insert + 3*std_dev + read_length);
    const int right_end = end + (mean_insert - 3*std_dev + read_length);
    process_mates(io, tid, right_start, right_end, bloom);

    // Output reads and count length
    find_mates(io, buffer, bloom, &reads, &seqlen, &reads_extracted);

    // Output overlapping reads
    if (insertion || !no_overlap) {
      process_region(io, tid, start, end, buffer, bloom, &reads, &seqlen, &reads_extracted);
    }
  }

  // Output unmapped reads
  if (unmapped_only || (seqlen / length) < threshold) {
    process_unmapped(io, buffer, bloom, &reads, &unmapped_extracted);
  }

  std::cout << "Extracted " << reads_extracted;
  if (unmapped_extracted > 0) {
    std::cout << "(+" << unmapped_extracted << ")";
  }
  std::cout << " out of " << num_of_reads << " reads" << std::endl;

  // Cleanup
  reads.flush();
  delete[] buffer;
  delete bloom;
  io.unload();
}
