// Copyright 2016 Riku Walve

#include <vector>

#include "catch.hpp"

#include "htslib/sam.h"

#include "hash.hpp"

static char string1[] = "asd";
static char string2[] = "asd1";

TEST_CASE("Hashes are computed", "[hash]") {
  SECTION("String hash functions") {
    REQUIRE(hash_str(string1) == hash_str(string1));
    REQUIRE(hash_str(string1) != hash_str(string2));
  }

  SECTION("Alignment hash functions") {
    bam1_t *bam1 = bam_init1();
    bam1->core.flag = BAM_FREAD1;
    bam1->data = reinterpret_cast<uint8_t*>(string1);

    bam1_t *bam2 = bam_init1();
    bam2->core.flag = BAM_FREAD2;
    bam2->data = reinterpret_cast<uint8_t*>(string1);

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
  std::vector<size_t> vector;
  vector.push_back(hash_str(string1));

  REQUIRE(in_alignments(vector, hash_str(string1)));
  REQUIRE(!in_alignments(vector, hash_str(string2)));

  vector.push_back(hash_str(string2));
  REQUIRE(in_alignments(vector, hash_str(string2)));
}
