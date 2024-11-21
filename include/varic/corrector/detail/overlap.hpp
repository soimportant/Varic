#pragma once

#include <optional>
#include <string_view>
#include <biovoltron/file_io/cigar.hpp>

class Overlap {
public:

  using id_t = std::uint32_t;
  using pos_t = std::uint32_t;
  using len_t = pos_t;

  /* query read id */
  id_t q_id;

  /* left and right boundary of query read */
  pos_t q_idx_L, q_idx_R;

  /* query sequence length */
  len_t q_seq_len;

  /* target read id */
  id_t t_id;

  /* left and right boundary of target read */
  pos_t t_idx_L, t_idx_R;

  /* target sequence length */
  len_t t_seq_len;

  /* forward strain(true) or reverse strain(false) */
  bool forward_strain;

  auto extend(const len_t extend_len) noexcept {
    q_idx_L = q_idx_L > extend_len ? q_idx_L - extend_len : 0UL;
    q_idx_R = std::min(q_idx_R + extend_len, q_seq_len);
    t_idx_L = t_idx_L > extend_len ? t_idx_L - extend_len : 0UL;
    t_idx_R = std::min(t_idx_R + extend_len, t_seq_len);
  }

  /* three-way operator */
  auto operator<=>(const Overlap& rhs) const noexcept {
    return this->q_id <=> rhs.q_id;
  }
};


std::ostream& operator<<(std::ostream& os, const Overlap& o) {
  os << o.q_id << '\t' << o.q_idx_L << '\t' << o.q_idx_R << '\t' << o.q_seq_len
     << '\t' << o.t_id << '\t' << o.t_idx_L << '\t' << o.t_idx_R << '\t'
     << o.t_seq_len << '\t' << std::boolalpha << o.forward_strain;
  return os;
}