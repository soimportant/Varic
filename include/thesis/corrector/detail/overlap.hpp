#pragma once

#include <optional>
#include <string_view>
#include <biovoltron/file_io/cigar.hpp>

class Overlap {
public:

  /* query read id */
  std::size_t q_id;

  /* left and right boundary of query read */
  std::size_t q_idxL, q_idxR;

  /* target read id */
  std::size_t t_id;

  /* left and right boundary of target read */
  std::size_t t_idxL, t_idxR;

  /**
   * @brief string_view of query read sequence 
   * @note this not the overlap part of query read since we may extend the
   * boundary
   */
  std::string_view q_seq;

  /* string_view of query read quality */
  std::optional<std::string_view> q_qual;

  /**
   * @brief string_view of target read sequence 
   * @note this not the overlap part of target read since we may extend the
   * boundary
   */
  std::string_view t_seq;

  /* string_view of target read quality */
  std::optional<std::string_view> t_qual;

  /* forward strain(true) or reverse strain(false) */
  bool strain;

  /* cigar string of overlap part between query read and target read */
  bio::Cigar cigar;
};