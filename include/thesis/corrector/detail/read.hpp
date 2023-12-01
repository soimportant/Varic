#pragma once

#include <concepts>
#include <optional>

#include <biovoltron/file_io/fasta.hpp>
#include <spdlog/spdlog.h>

#include "thesis/corrector/detail/window.hpp"
#include "thesis/corrector/detail/sequence.hpp"
#include "thesis/algo/assemble/read_assembler.hpp"

namespace bio = biovoltron;

template <class R>
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

  ReadWrapper(R&& r) : R(std::move(r)) {};

  auto create_rc() noexcept {
    if (rc_seq.size() != 0) {
      return true;
    }
    try {
      rc_seq = bio::Codec::rev_comp(this->seq);
      if constexpr (std::same_as<R, bio::FastqRecord<R::encoded>>) {
        rev_qual = std::string(this->qual.rbegin(), this->qual.rend());
      }
    } catch (std::bad_alloc& e) {
      spdlog::error(
          "bad_alloc caught when creating reverse complement sequence");
      return false;
    }
    return true;
  }

  auto len() const noexcept {
    return this->seq.size();
  }

  /**
   * @brief return a 
   * 
   * @param forward_strain whether the position is on forward strain
   * @param start start position(0-based)
   * @param end end position(0-based)
   * @return std::pair<std::string_view, std::optional<std::string_view>>
   */
  auto subview(bool forward_strain, const std::size_t start,
               const std::size_t end) const {
    if (start > end) {
      spdlog::error("start position {} is larger than end position {}", start,
                    end);
      throw std::invalid_argument("start position is larger than end position");
    }

    if (start > this->seq.size()) {
      spdlog::error("start position {} is out of range of sequence", start);
      throw std::out_of_range("start position is out of range of sequence");
    }

    auto real_start = forward_strain ? start : this->seq.size() - end;
    auto len = end - start;

    auto seq_view = std::string_view{};
    auto qual_view = std::optional<std::string_view>{};
    seq_view = forward_strain ? this->seq : rc_seq;
    if constexpr (std::same_as<R, bio::FastqRecord<R::encoded>>) {
      qual_view = forward_strain ? this->qual : rev_qual.value();
    }
    seq_view = seq_view.substr(real_start, len);
    if constexpr (std::same_as<R, bio::FastqRecord<R::encoded>>) {
      qual_view = qual_view->substr(real_start, len);
    }
    return std::make_pair(seq_view, qual_view);
  }

  auto add_corrected_fragment(const Sequence<std::string>& fragment) {
    corrected_fragments.emplace_back(fragment);
  }

  auto operator<=>(const ReadWrapper& rhs) const {
    return this->name <=> rhs.name;
  }

private:
  // TODO: coverage -> segment tree like data structure
  // ? we may add a data structure here for recording the coverage covered
  // ? by corrected_fragments, if the coverage is enough, we can assemble
  // ? the corrected read without building MSA.

  
  /* reverse complement sequence */
  std::string rc_seq;

  /* reverse quality sequence */
  std::optional<std::string> rev_qual;

public:
  std::size_t id;
  /**
    * corrected sequence fragments from windows of others read, may further
    * assemble to correct version of this read
    */
  std::vector<Sequence<std::string>> corrected_fragments;

  /* windows of this read */
  std::vector<Window> windows;

  // /* assembler of this read */
  // ReadAssembler assembler;
 
};
