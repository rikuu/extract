// Copyright 2016 Riku Walve

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

static const char* STR_SCAFFOLD = "-scaffold";
static const char* STR_GAP_START = "-start";
static const char* STR_GAP_END = "-end";

static const char* STR_THRESHOLD = "-unmapped";

/*****************************************************************************/

Extract::Extract() : Tool("Extract") {
  getParser()->push_front(new OptionOneParam(STR_ALIGNMENT, "Aligned BAM file", true));
  getParser()->push_front(new OptionOneParam(STR_OUTPUT, "FASTA-formatted output file", true));

  getParser()->push_front(new OptionOneParam(STR_READ_LENGTH, "Read length", true));
  getParser()->push_front(new OptionOneParam(STR_MEAN, "Mean insert size", true));
  getParser()->push_front(new OptionOneParam(STR_STD_DEV, "Insert size standard deviation", true));

  getParser()->push_front(new OptionOneParam(STR_SCAFFOLD, "Scaffold name", true));
  getParser()->push_front(new OptionOneParam(STR_GAP_START, "Gap starting position", true));
  getParser()->push_front(new OptionOneParam(STR_GAP_END, "Gap ending position", true));

  getParser()->push_front(new OptionOneParam(STR_THRESHOLD, "Threshold for using unmapped reads", false, "-1"));
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

  // std::cout << '>' << bam2string(bam) << '\n' << sequence << std::endl;

  Sequence seq(buffer);
  seq._comment = bam2string(bam);
  bank->insert(seq);
}

// Output reads that map to a region in a scaffold and not in the Bloom filter.
int Extract::process_region(const io_t io, const int tid, const int start,
    const int end, char* buffer, IBloom<std::string> *bloom, BankFasta *bank) {
  int seqlen = 0;

  sam_iterator iter(io, tid, start, end);
  while (iter.next()) {
    if (!bloom->contains(bam2string(iter.bam))) {
      print_fasta(iter.bam, buffer, bank);
      seqlen += strlen(buffer);
    }
  }

  return seqlen;
}

void Extract::process_mates(const io_t io, const int tid, const int start,
    const int end, IBloom<std::string> *bloom) {
  sam_iterator iter(io, tid, start, end);
  while (iter.next()) {
    if ((iter.bam->core.flag & BAM_FMUNMAP) != 0) {
      bloom->insert(bam2string(iter.bam));
    }
  }
}

int Extract::find_mates(const io_t io, char *buffer, IBloom<std::string> *bloom,
    BankFasta *bank) {
  int seqlen = 0;

  sam_iterator iter(io, ".");
  while (iter.next()) {
    if (bloom->contains(bam2string_mate(iter.bam))) {
      print_fasta(iter.bam, buffer, bank);
      seqlen += strlen(buffer);
    }
  }

  return seqlen;
}

// Output all unmapped reads
void Extract::process_unmapped(const io_t io, char* buffer,
    IBloom<std::string> *bloom, BankFasta *bank) {
  // Go through entire BAM file rather than the unmapped reads only.
  // BWA gives unmapped reads with a mapped mate a position, essentially making
  // it a mapped read. Outputing these reads is debateble.
  sam_iterator iter(io, ".");
  while (iter.next()) {
    if ((iter.bam->core.flag & BAM_FUNMAP) != 0 &&
          !bloom->contains(bam2string(iter.bam))) {
      print_fasta(iter.bam, buffer, bank);
    }
  }
}

// Counts the number of reads by iterating through the alignment file
int count_reads(const std::string filename) {
  io_t io(filename);
  sam_iterator iter(io, ".");

  int count = 0;
  while (iter.next()) {
    count++;
  }

  io.unload();

  return count;
}

// Execute read extraction
void Extract::execute() {
  // Get the command line arguments
  const std::string alignment = getInput()->getStr(STR_ALIGNMENT);
  const std::string output = getInput()->getStr(STR_OUTPUT);

  const int read_length = getInput()->getInt(STR_READ_LENGTH);
  const int mean_insert = getInput()->getInt(STR_MEAN);
  const int std_dev = getInput()->getInt(STR_STD_DEV);

  const std::string scaffold = getInput()->getStr(STR_SCAFFOLD);
  const int start = getInput()->getInt(STR_GAP_START);
  const int end = getInput()->getInt(STR_GAP_END);

  const bool unmapped = (getInput()->getInt(STR_THRESHOLD) != -1);
  const int threshold = getInput()->getInt(STR_THRESHOLD);

  // Load alignment file
  io_t io(alignment);
  if (!io.loaded) {
    std::cerr << "Error loading alignments" << std::endl;
    return;
  }

  // Compute scaffold id from scaffold name
  const int tid = bam_name2id(io.header, scaffold.c_str());

  // Use basic Bloom filter from GATB
  IBloom<std::string> *bloom = new BloomSynchronized<std::string>(5 * count_reads(alignment));

  // Extract pairs from the left mappings
  const int left_start = start - (mean_insert + 3*std_dev + 2*read_length);
  const int left_end = end - (mean_insert - 3*std_dev + read_length);
  process_mates(io, tid, left_start, left_end, bloom);

  // Extract pairs from the right mappings
  const int right_start = start + (mean_insert + 3*std_dev + read_length);
  const int right_end = end + (mean_insert - 3*std_dev + read_length);
  process_mates(io, tid, right_start, right_end, bloom);

  // Allocate memory for string conversions
  char *buffer = new char[read_length+1];

  // Output reads and count length
  BankFasta reads(output);
  const int seqlen =
      process_region(io, tid, start, end, buffer, bloom, &reads) +
      find_mates(io, buffer, bloom, &reads);

  // Possibly also output unmapped reads
  if (unmapped && ((seqlen / (end - start)) < threshold)) {
    process_unmapped(io, buffer, bloom, &reads);
  }

  // Cleanup
  reads.flush();
  delete[] buffer;
  delete bloom;
  io.unload();
}
