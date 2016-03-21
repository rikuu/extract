#include <vector>

#include "htslib/sam.h"

#include "hash.hpp"

size_t hash_str(char *str) {
  size_t hash = 5381;

  int c;
  while ((c = *(str++)))
    hash = ((hash << 5) + hash) + c;

  return hash;
}
