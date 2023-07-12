#pragma once

#include <fstream>
#include <filesystem>
#include <vector>
#include <sstream>
#include <utility>
#include <string_view>
#include <map>
#include <ranges>

#include <spdlog/spdlog.h>
#include <biovoltron/file_io/vcf.hpp>

namespace fs = std::filesystem;
namespace bio = biovoltron;

auto read_vcf(const fs::path& path) {
  auto fin = std::ifstream(path);
  spdlog::info("Read vcf records from {} ...", path.string());
  bio::VcfHeader h;
  fin >> h;
  std::vector<bio::VcfRecord> records;
  for (bio::VcfRecord r; fin;) {
    std::string line;
    std::getline(fin, line);
    if (line.empty()) {
      continue;
    }
    std::stringstream ss(line);

    ss >> r.chrom >> r.pos >> r.id >> r.ref >> r.alt;
    /* supposed to be quality, but it may be '.' */
    /* use dummy r.filter instead */
    // ss >> r.qual;
    ss >> r.filter;
    ss >> r.filter >> r.info;
    records.emplace_back(std::move(r));
  }
  return std::make_pair(std::move(h), std::move(records));
}

/* write header and record into vcf file */
auto write_vcf(const fs::path& path, const bio::VcfHeader& header,
               const std::vector<bio::VcfRecord>& records) {
  auto fout = std::ofstream(path);
  assert(fout.is_open());
  fout << header << std::endl;
  for (const auto& r : records) {
    fout << r << std::endl;
  }
}

auto parse_vcf_info(const bio::VcfRecord& record) {
  std::string_view info = record.info;
  std::map<std::string_view, std::string_view> result;
  // when use std::views::split to split std::string, you should use 
  // std::string_view or character as the delimiter, if you use something like
  // ";", the function will consider the delimiter is ";\0", so you won't get
  // the proper result
  for (const auto t : std::views::split(info, ';')) {
    auto s = std::string_view(t.begin(), t.end());
    auto pos = s.find('=');
    auto key = s.substr(0, pos);
    auto value = s.substr(pos + 1);
    result[key] = value;
  }
  return result;
}

auto get_vcf_range(const std::vector<bio::VcfRecord>& vcf,
                   const std::size_t st, const std::size_t ed) {
  // std::lower_bound and std::upper_bound is too difficult to use
  // the requirement of compare function is too annoying
  // you need to handle the all combination when the type of value and iterator
  // is different
  // Also, the std::lower_bound and std::upper_bound has compare function
  // in different form, which is very annoying, either.

  const auto sz = vcf.size();
  auto vcf_info = std::vector<std::map<std::string_view, std::string_view>>(sz);
  for (std::size_t i = 0; i < sz; ++i) {
    assert(vcf[i].info.size() > 0);
    vcf_info[i] = parse_vcf_info(vcf[i]);
  }
  auto beg = 0u, end = 0u;
  {
    auto low = 0ul, high = sz;
    while (low < high) {
      auto mid = low + (high - low) / 2;
      auto v = std::stoull(std::string(vcf_info[mid]["sim_start"]));
      if (v < st) {
        low = mid + 1;
      } else {
        high = mid;
      }
    }
    beg = low;
  }
  {
    auto low = 0ul, high = sz;
    while (low < high) {
      auto mid = low + (high - low) / 2;
      auto v = std::stoull(std::string(vcf_info[mid]["sim_start"]));
      if (v <= ed) {
        low = mid + 1;
      } else {
        high = mid;
      }
    }
    end = low;
  }
  return std::ranges::subrange(vcf.begin() + beg, vcf.begin() + end);
};
