#include <ranges>
#include <vector>
#include <execution>

#include <biovoltron/file_io/all.hpp>
#include "utility/all.hpp"

namespace bio = biovoltron;

auto read_preprocess(const std::vector<bio::FastaRecord<false>>& reads) {
  /**
   * 1. filter read length
   */

  const auto min_read_length = 1000ull;
  auto filtered_reads = std::vector<bio::FastaRecord<false>>{};
  /* there's no std::ranges::move_if() */
  for (auto& r : reads) {
    if (r.seq.size() >= min_read_length) {
      filtered_reads.emplace_back(std::move(r));
    }
  }
  return filtered_reads;
}

auto overlap_preprocess(const std::vector<bio::FastaRecord<false>>& raw_reads,
                        std::vector<bio::PafRecord>& overlaps) {
  /**
   * 1. filter overlap length
   * 2. filter overlap identity
   */

  assert(
      std::ranges::is_sorted(raw_reads, {}, &bio::FastaRecord<false>::name) &&
      "raw_reads must be sorted by name");
  auto valid_overlap = [&](const bio::PafRecord& overlap) {
    auto valid_read = [&](const std::string_view name) {
      return std::ranges::binary_search(raw_reads, name, {},
                                        &bio::FastaRecord<false>::name);
    };
    bool valid = valid_read(overlap.qname) && valid_read(overlap.tname);
    /* filter gogogo */
    return valid;
  };
  const auto min_overlap_length = 1000ull;
  const auto min_overlap_identity = 0.8;
  /* there's no std::ranges::move_if() */
  auto filtered_overlaps = std::vector<bio::PafRecord>{};
  for (auto& overlap : overlaps) {
    if (valid_overlap(overlap)) {
      filtered_overlaps.emplace_back(std::move(overlap));
    }
  }
  std::sort(std::execution::par, filtered_overlaps.begin(),
            filtered_overlaps.end(), [&](auto& a, auto& b) {
              return a.qname < b.qname || a.tname < b.tname;
            });
  return filtered_overlaps;
}