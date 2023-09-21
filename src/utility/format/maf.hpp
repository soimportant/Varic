#pragma once

#include <biovoltron/file_io/core/header.hpp>
#include <biovoltron/file_io/core/record.hpp>
#include <biovoltron/file_io/core/tuple.hpp>
#include <cassert>
#include <iomanip>
#include <spdlog/spdlog.h>
#include <vector>

namespace fs = std::filesystem;
namespace bio = biovoltron;

/**
 * @ingroup file_io
 */
namespace biovoltron {
class MafHeader : public Header {
public:
  constexpr static auto START_SYMBOLS = std::array{"#"};
};

class MafRecord {
public:
  constexpr static auto START_SYMBOLS = std::array{"a", "s"};
  enum LineType { s, q, i, e };

  /* the alignment block start with character 'a' */
  class Alignment {
  public:
    std::string name;
    std::size_t start;
    std::size_t len;
    char strand;
    std::size_t src_size;
    std::string seq;
  };

  Alignment ref;
  std::vector<Alignment> aln;
};

std::istream& operator>>(std::istream& is, MafRecord::Alignment& aln) {
  std::string line;
  std::getline(is, line);
  if (line.empty()) {
    is.clear(std::ios::eofbit);
    return is;
  }
  std::istringstream iss{line};
  iss >> aln.name >> aln.start >> aln.len >> aln.strand >> aln.src_size >>
      aln.seq;
  return is;
}

std::istream& operator>>(std::istream& is, MafRecord& r) {
  std::string line;
  std::getline(is, line);
  if (line.empty()) {
    return is;
  }
  assert(line[0] == 'a');
  /* line start with a, deal with some extra field like score=%f */
  is.get();
  is >> r.ref;
  MafRecord::Alignment aln;
  while (is.get() == 's') {
    is >> aln;
    r.aln.emplace_back(std::move(aln));
  }
  return is;
}
}; // namespace biovoltron

auto read_maf(const fs::path& path) {
  auto fin = std::ifstream(path);
  spdlog::info("Read maf records from \"{}\" ...", path.string());
  bio::MafHeader h;
  fin >> h;
  std::vector<bio::MafRecord> records;
  for (bio::MafRecord r; fin >> r;) {
    records.emplace_back(std::move(r));
  }
  return std::make_pair(std::move(h), std::move(records));
}

auto& get_maf_record(const std::vector<bio::MafRecord>& maf,
                     const std::string_view read_name) {
  auto low = (std::size_t) 0, high = maf.size();
  while (low < high) {
    auto mid = (low + (high - low) / 2);
    auto r = maf[mid];
    if (r.aln[0].name < read_name) {
      low = mid + 1;
    } else {
      high = mid;
    }
  }
  assert(low != maf.size());
  return maf[low];
}