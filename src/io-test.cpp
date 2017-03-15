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
  char *buffer = new char[length];
  REQUIRE(convertToString(sequence, length, false, buffer) == "ATGCTGTCTN");
  REQUIRE(convertToString(sequence, length, true, buffer) == "NAGACAGCAT");
  delete[] buffer;
}
