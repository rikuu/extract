#ifndef EXTRACT_HPP
#define EXTRACT_HPP

#include <sys/types.h>

#include <string>
#include <vector>
#include <functional>

// GATB-core basic Bloom filter requires hash1 function for items
// TODO: Hash bam object directly
inline u_int64_t hash1(const std::string &key, u_int64_t seed=0) {
  return std::hash<std::string>{}(key);
}

#include <gatb/gatb_core.hpp>

#include "htslib/sam.h"

#include "io.hpp"

class Extract : public Tool {
 public:
    Extract();
    void execute();

    // Prints an alignment in fasta format
    void print_fasta(const bam1_t*, char*, BankFasta*);

    // Prints all alignments in a region
    int process_region(const io_t, const int, const int, const int, char*, IBloom<std::string>*, BankFasta*);
    void process_mates(const io_t, const int, const int, const int, IBloom<std::string>*);
    int find_mates(const io_t, char*, IBloom<std::string>*, BankFasta*);
    void process_unmapped(const io_t, char*, IBloom<std::string>*, BankFasta*);
};
#endif
