// Copyright 2016 Riku Walve

#include <vector>

#include "catch.hpp"

#include "htslib/sam.h"

#include "hash.hpp"

static char string1[] = "asd";
static char string2[] = "asd1";

inline bam1_t * construct_bam(char *string, const uint32_t flag) {
  bam1_t *bam = bam_init1();
  bam->core.flag = flag;
  bam->data = reinterpret_cast<uint8_t*>(string);

  return bam;
}

TEST_CASE("Hashes are computed", "[hash]") {
  SECTION("String hash functions") {
    REQUIRE(hash_str(string1) == hash_str(string1));
    REQUIRE(hash_str(string1) != hash_str(string2));
  }

  SECTION("Alignment hash functions") {
    bam1_t *bam1 = construct_bam(string1, BAM_FREAD1);
    bam1_t *bam2 = construct_bam(string1, BAM_FREAD2);

    REQUIRE(hash_alignment1(bam1) == hash_alignment1(bam1));
    REQUIRE(hash_alignment1(bam2) == hash_alignment1(bam2));

    REQUIRE(hash_alignment2(bam1) == hash_alignment2(bam1));
    REQUIRE(hash_alignment2(bam2) == hash_alignment2(bam2));

    REQUIRE(hash_alignment1(bam1) != hash_alignment1(bam2));
    REQUIRE(hash_alignment2(bam1) != hash_alignment2(bam2));

    REQUIRE(hash_alignment1(bam1) == hash_alignment2(bam2));
    REQUIRE(hash_alignment1(bam2) == hash_alignment2(bam1));

    bam1->data = NULL;
    bam2->data = NULL;

    bam_destroy1(bam1);
    bam_destroy1(bam2);
  }
}

TEST_CASE("Alignments are found", "[hash]") {
  bam1_t *bam1a = construct_bam(string1, BAM_FREAD1);
  bam1_t *bam1b = construct_bam(string1, BAM_FREAD2);
  bam1_t *bam2a = construct_bam(string2, BAM_FREAD1);
  bam1_t *bam2b = construct_bam(string2, BAM_FREAD2);

  bloom bloom;
  bloom.push(bam1a);

  REQUIRE(!bloom.in_alignments(bam1a));
  REQUIRE(bloom.in_alignments(bam1b));
  REQUIRE(!bloom.in_alignments(bam2a));
  REQUIRE(!bloom.in_alignments(bam2b));

  bloom.push(bam2a);

  REQUIRE(!bloom.in_alignments(bam1a));
  REQUIRE(bloom.in_alignments(bam1b));
  REQUIRE(!bloom.in_alignments(bam2a));
  REQUIRE(bloom.in_alignments(bam2b));

  bloom.clear();

  REQUIRE(!bloom.in_alignments(bam1a));
  REQUIRE(!bloom.in_alignments(bam1b));
  REQUIRE(!bloom.in_alignments(bam2a));
  REQUIRE(!bloom.in_alignments(bam2b));

  bam1a->data = NULL;
  bam1b->data = NULL;
  bam2a->data = NULL;
  bam2b->data = NULL;

  bam_destroy1(bam1a);
  bam_destroy1(bam1b);
  bam_destroy1(bam2a);
  bam_destroy1(bam2b);
}
