#pragma once

#include <optional>
#include <string_view>

class Sequence {
public:
  auto len() const noexcept -> std::size_t {
    return seq.size();
  }

  /* the given id of read where the sequence come from */
  std::size_t read_id;

  /* the left and right boundary of read */
  /* assert(right_bound - left_bound == seq.size() == qual.size()) */
  std::size_t left_bound, right_bound;

  /* the sequence, which is a subsequence of read */
  std::string_view seq;

  /**
   * the quality of `seq`, if qual.has_value() == false, then the read is come
   * from fasta file
   */
  std::optional<std::string_view> qual;

  /* forward strain(true) or reverse strain(false) */
  bool forward_strain;
};

