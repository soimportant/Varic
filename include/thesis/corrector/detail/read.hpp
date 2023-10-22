#pragma once

#include <concepts>
#include <optional>

#include <biovoltron/file_io/fasta.hpp>
#include <spdlog/spdlog.h>

#include "thesis/corrector/detail/window.hpp"

namespace bio = biovoltron;

template<class R>
  requires std::derived_from<R, bio::FastaRecord<R::encoded>>
class ReadWrapper : public R {
 public:
  ReadWrapper() = default;

  /* disable copy constructor and copy assignment */
  ReadWrapper(const ReadWrapper&) = delete;
  ReadWrapper& operator=(const ReadWrapper&) = delete;

  /* default move constructor and move assignment */
  ReadWrapper(ReadWrapper&&) = default;
  ReadWrapper& operator=(ReadWrapper&&) = default;

  auto create_rc() noexcept {
    if (rc_seq.size() != 0) {
      return true;
    }
    try {
      rc_seq = bio::Codec::rev_comp(this->seq);
      if constexpr (std::same_as<R, bio::FastqRecord<R::encoded>>) {
        rc_qual = std::string(this->qual.rbegin(), this->qual.rend());
      }
    } catch (std::bad_alloc& e) {
      spdlog::error("bad_alloc caught when creating reverse complement sequence");
      return false;
    }
    return true;
  }


  auto subview(bool forward_strain, std::size_t start, std::size_t end) const {
    if (start > this->seq.size()) {
      spdlog::error("start position {} is out of range of sequence", start);
      throw std::out_of_range("start position is out of range of sequence");
    }
    start = forward_strain ? start : this->seq.size() - start;
    end = forward_strain ? end : this->seq.size() - end;

    auto seq_view = std::string_view{};
    auto qual_view = std::optional<std::string_view>{};
    seq_view = forward_strain ? this->seq : rc_seq;
    if constexpr (std::same_as<R, bio::FastqRecord<R::encoded>>) {
      qual_view = forward_strain ? this->qual : rc_qual;
    }
    seq_view = seq_view.substr(start, end - start);
    if constexpr (std::same_as<R, bio::FastqRecord<R::encoded>>) {
      qual_view = qual_view->substr(start, end - start);
    }
    return std::make_pair(seq_view, qual_view);
  }


public:
  // TODO: coverage -> segment tree like data structure
  // ? we may add a data structure here for recording the coverage covered
  // ? by corrected_fragments, if the coverage is enough, we can assemble
  // ? the corrected read without building MSA.

  std::size_t id;

  /* reverse complement sequence */
  std::string rc_seq;

  /* reverse quality sequence */
  std::optional<std::string> rc_qual;

  /**
   * corrected sequence fragments from windows of others read, may further
   * assemble to correct version of this read
   */
  std::vector<std::string> corrected_fragments;

  /* windows of this read */
  std::vector<Window> windows;
};