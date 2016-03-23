// Copyright 2016 Riku Walve

#include "catch.hpp"

#include "io.hpp"

// ATGCTGTCTN converted into htslib format
const uint8_t sequence[] = {24, 66, 132, 130, 143};
const int32_t length = 10;

TEST_CASE("Sequence querying works", "[io]") {
  SECTION("Forward strand") {
    REQUIRE(querySequence(sequence, length, 0, false) == 1);
    REQUIRE(querySequence(sequence, length, 1, false) == 8);
    REQUIRE(querySequence(sequence, length, 2, false) == 4);
    REQUIRE(querySequence(sequence, length, 3, false) == 2);
    REQUIRE(querySequence(sequence, length, 4, false) == 8);
    REQUIRE(querySequence(sequence, length, 5, false) == 4);
    REQUIRE(querySequence(sequence, length, 6, false) == 8);
    REQUIRE(querySequence(sequence, length, 7, false) == 2);
    REQUIRE(querySequence(sequence, length, 8, false) == 8);
    REQUIRE(querySequence(sequence, length, 9, false) == 15);
  }

  SECTION("Reverse strand") {
    REQUIRE(querySequence(sequence, length, 0, true) == 15);
    REQUIRE(querySequence(sequence, length, 1, true) == 1);
    REQUIRE(querySequence(sequence, length, 2, true) == 4);
    REQUIRE(querySequence(sequence, length, 3, true) == 1);
    REQUIRE(querySequence(sequence, length, 4, true) == 2);
    REQUIRE(querySequence(sequence, length, 5, true) == 1);
    REQUIRE(querySequence(sequence, length, 6, true) == 4);
    REQUIRE(querySequence(sequence, length, 7, true) == 2);
    REQUIRE(querySequence(sequence, length, 8, true) == 1);
    REQUIRE(querySequence(sequence, length, 9, true) == 8);
  }
}

TEST_CASE("Sequence converting works", "[io]") {
  REQUIRE(convertToString(sequence, length, false) == "ATGCTGTCTN");
  REQUIRE(convertToString(sequence, length, true) == "NAGACAGCAT");
}
