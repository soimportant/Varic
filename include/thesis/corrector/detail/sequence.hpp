#pragma once

#include <optional>
#include <string_view>
#include <type_traits>

template<class T, class CharType = T::value_type>
concept StringLike = std::is_convertible_v<T, std::basic_string_view<CharType>>;

template <StringLike T = std::string_view>
class Sequence {
 public:
  Sequence() = default;
  Sequence(const int read_id, const std::size_t left_bound,
           const std::size_t right_bound, T seq, std::optional<std::string_view> qual,
           const bool forward_strain)
      : read_id(read_id),
        left_bound(left_bound),
        right_bound(right_bound),
        seq(std::move(seq)),
        qual(std::move(qual)),
        forward_strain(forward_strain) {}

  // default copy constructor and copy assignment
  Sequence(const Sequence&) = default;
  Sequence& operator=(const Sequence&) = default;

  // default move constructor and move assignment
  Sequence(Sequence&&) = default;
  Sequence& operator=(Sequence&&) = default;

  auto len() const noexcept { return seq.size(); }

  auto empty() const noexcept { return seq.empty(); }

  ~Sequence() = default;

  /* the given id of read where the sequence come from */
  int read_id;

  /* the left and right boundary of read */
  /* assert(right_bound - left_bound == seq.size() == qual.size()) */
  std::size_t left_bound, right_bound;

  /* the sequence, which is a subsequence of read */
  T seq;

  /**
   * the quality of `seq`, the sequence doesn't have quality score if
   * qual.has_value() == false
   */
  std::optional<std::string_view> qual;

  /* forward strain(true) or reverse strain(false) */
  bool forward_strain;
};
