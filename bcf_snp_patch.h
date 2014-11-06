#ifndef _BCF_SNP_PATCH_H_
#define _BCF_SNP_PATCH_H_

#include <cstdint>
#include <utility>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <sstream>

typedef uint64_t snp_loci_t;
typedef std::unordered_map<snp_loci_t, uint32_t> snp_map_t;

#define base_to_bits(b) (((b) >> 1) & 3)

// seperator used to expand read identifier
#define READ_EXPAND_SEP "/"

// complementary bases in PREMIER's order
const char pmr_complements[4] = {'T', 'G', 'A', 'C'};
// nucleotides in samtools's order
const char sam_nucleotides[4] = {'A', 'C', 'G', 'T'};

#endif
