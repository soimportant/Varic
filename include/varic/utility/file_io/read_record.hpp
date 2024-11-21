#pragma once

#include <filesystem>
#include <fstream>
#include <vector>

#include <boost/type_index.hpp>
#include <spdlog/spdlog.h>

namespace fs = std::filesystem;

// template<std::derived_from<bio::Record> R>
template <class R>
auto read_records(const fs::path& path) {
  auto fin = std::ifstream(path);
  if (!fin.is_open()) {
    spdlog::error("Cannot open file: {}", path.string());
    std::exit(1);
  }
  auto records = std::vector<R>{};
  spdlog::debug("Read {} from \"{}\" ...",
                boost::typeindex::type_id_with_cvr<R>().pretty_name(),
                path.string());
  for (R r; fin >> r;) {
    records.emplace_back(std::move(r));
  }
  return records;
}