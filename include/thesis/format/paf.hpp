#pragma once

#include <compare>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <algorithm>

namespace biovoltron {

/**
 * @brief Pairwise mApping Format
 * @link https://github.com/lh3/miniasm/blob/master/PAF.md
 */
class PafRecord {
 public:
  auto extend(std::size_t len) {
    len = std::min({len, q_start, q_len - q_end});
    len = std::min({len, t_start, t_len - t_end});
    q_start -= len;
    q_end += len;
    t_start -= len;
    t_end += len;
  };
  
 public:
  std::string q_name;
  std::size_t q_len;
  std::size_t q_start;
  std::size_t q_end;
  char strand;
  std::string t_name;
  std::size_t t_len;
  std::size_t t_start;
  std::size_t t_end;
  std::size_t match;
  std::size_t aln_len;
  std::size_t mapq;

  friend std::istream& operator>>(std::istream& os, PafRecord& r);
  friend std::ostream& operator<<(std::ostream& os, const PafRecord& r);
  friend auto operator<=>(const PafRecord& lhs, const PafRecord& rhs);
};

std::istream& operator>>(std::istream& os, PafRecord& r) {
  std::string line;
  std::getline(os, line);
  if (line.empty()) {
    return os;
  }
  std::istringstream iss{line};
  iss >> r.q_name >> r.q_len >> r.q_start >> r.q_end >> r.strand >> r.t_name >>
      r.t_len >> r.t_start >> r.t_end >> r.match >> r.aln_len >> r.mapq;
  return os;
}

std::ostream& operator<<(std::ostream& os, const PafRecord& r) {
  os << r.q_name << "\t" << r.q_len << "\t" << r.q_start << "\t" << r.q_end
     << "\t" << r.strand << "\t" << r.t_name << "\t" << r.t_len << "\t"
     << r.t_start << "\t" << r.t_end << "\t" << r.match << "\t" << r.aln_len
     << "\t" << r.mapq;
  return os;
}

auto operator<=>(const PafRecord& lhs, const PafRecord& rhs) {
  if (lhs.q_name != rhs.q_name) {
    return lhs.q_name <=> rhs.q_name;
  }
  if (std::pair(lhs.q_start, lhs.q_end) != std::pair(rhs.q_start, rhs.q_end)) {
    return std::pair(lhs.q_start, lhs.q_end) <=>
           std::pair(rhs.q_start, rhs.q_end);
  }
  return lhs.t_name <=> rhs.t_name;
}

}  // namespace biovoltron
