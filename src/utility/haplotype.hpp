#pragma once

#include <vector>
#include <utility>
#include <biovoltron/file_io/all.hpp>
#include "maf.hpp"

class Haplotype {
public:
  Haplotype() = default;
  bio::FastaRecord<false> reference;
  std::vector<bio::FastqRecord<false>> raw_reads;
  std::vector<bio::FastaRecord<false>> corrected_reads;
  std::pair<bio::VcfHeader, std::vector<bio::VcfRecord>> snp_vcf;
  std::pair<bio::VcfHeader, std::vector<bio::VcfRecord>> indel_vcf;
  std::pair<bio::MafHeader, std::vector<bio::MafRecord>> maf;
};

