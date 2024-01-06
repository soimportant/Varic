#pragma once

#include <optional>
#include <string_view>
#include <type_traits>

template<class T>
concept StringLike = std::is_convertible_v<T, std::string_view>;

template<StringLike T = std::string_view>
class Sequence {
public:

  Sequence() = default;
  Sequence(const std::size_t read_id, const std::size_t left_bound,
           const std::size_t right_bound, T seq, std::optional<T> qual,
           const bool forward_strain)
      : read_id(read_id), left_bound(left_bound), right_bound(right_bound),
        seq(std::move(seq)), qual(std::move(qual)),
        forward_strain(forward_strain) {}

  // delete copy constructor and copy assignment
  Sequence(const Sequence&) = delete;
  Sequence& operator=(const Sequence&) = delete;

  // default move constructor and move assignment
  Sequence(Sequence&&) = default;
  Sequence& operator=(Sequence&&) = default;

  auto len() const noexcept -> std::size_t {
    return seq.size();
  }

  ~Sequence() = default;

  /* the given id of read where the sequence come from */
  std::size_t read_id;

  /* the left and right boundary of read */
  /* assert(right_bound - left_bound == seq.size() == qual.size()) */
  std::size_t left_bound, right_bound;

  /* the sequence, which is a subsequence of read */
  T seq;

  /**
   * the quality of `seq`, if qual.has_value() == false, then the read is come
   * from fasta file
   */
  std::optional<T> qual;

  /* forward strain(true) or reverse strain(false) */
  bool forward_strain;
};

