#include "catch.hpp"

#include "hash.hpp"

TEST_CASE("Hashes are computed", "[hash]") {
  char string1[] = "asd";
  char string2[] = "asd1";

  REQUIRE(hash_str(string1) == hash_str(string1));
  REQUIRE(hash_str(string1) != hash_str(string2));
}
