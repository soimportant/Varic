#pragma once

#include <chrono>
#include <filesystem>
#include <random>

#include <biovoltron/file_io/all.hpp>
#include <boost/describe.hpp>
#include <boost/mp11.hpp>
#include <boost/type_index.hpp>

#include "varic/format/paf.hpp"

namespace fs = std::filesystem;
namespace bpo = boost::program_options;
namespace bio = biovoltron;

using namespace std::literals;

enum FILE_FORMAT {
  FASTA = 1 << 0,
  FASTQ = 1 << 1,
  BAM = 1 << 2,
  uBAM = 1 << 3,
  GZ = 1 << 4,
  SAM = 1 << 5,
  VCF = 1 << 6,
  MAF = 1 << 7,
  PAF = 1 << 8,
  ERROR = 1 << 9
};

BOOST_DESCRIBE_ENUM(FILE_FORMAT, FASTA, FASTQ, BAM, uBAM, GZ, SAM, VCF, MAF,
                    PAF, ERROR);

int parse_file_format(fs::path p) {
  auto ext = p.extension().string();
  std::ranges::for_each(ext, [](auto& c) { c = std::tolower(c); });
  int format = 0;
  if (ext == ".gz") {
    format |= FILE_FORMAT::GZ;
    p.replace_extension("");
    ext = p.extension().string();
  }
  if (ext == ".fa" || ext == ".fasta") {
    format |= FILE_FORMAT::FASTA;
  } else if (ext == ".fq" || ext == ".fastq") {
    format |= FILE_FORMAT::FASTQ;
  } else if (ext == ".bam") {
    format |= FILE_FORMAT::BAM;
  } else if (ext == ".ubam") {
    format |= FILE_FORMAT::uBAM;
  } else if (ext == ".vcf") {
    format |= FILE_FORMAT::VCF;
  } else if (ext == ".maf") {
    format |= FILE_FORMAT::MAF;
  } else if (ext == ".paf") {
    format |= FILE_FORMAT::PAF;
  } else {
    format |= FILE_FORMAT::ERROR;
  }
  return format;
}

auto check_file_format(const std::string& opt_name, const fs::path& p,
                       const int& accept_format) {
  auto format = parse_file_format(p);
  if (format & FILE_FORMAT::ERROR || !(format & accept_format)) {
    std::string error_msg;
    boost::mp11::mp_for_each<
        boost::describe::describe_enumerators<FILE_FORMAT>>([&](auto f) {
      if (accept_format & f.value) {
        if (!error_msg.empty()) {
          error_msg += ", ";
        }
        error_msg += "."s + f.name;
      }
    });
    std::ranges::for_each(error_msg, [](auto& c) { c = std::tolower(c); });
    spdlog::error("invalid file format: \"{}\", accepted format: \"{}\"",
                  p.extension().string(), error_msg);
    throw bpo::validation_error(bpo::validation_error::invalid_option_value,
                                opt_name);
  }
  return format;
};




class ChainRecord {
 public:
  std::size_t len;
  std::size_t dq = 0;
  std::size_t dt = 0;
};

auto readChainFile(const fs::path& p) {
  auto fin = std::ifstream(p);
  if (!fin.is_open()) {
    spdlog::error("Cannot open file: {}", p.string());
    std::exit(1);
  }
  std::string line;
  while (std::getline(fin, line)) {
    if (line[0] == '#') {
      continue;
    }
    if (line.substr(0, 5) == "chain") {
      break;
    }
  }
  std::vector<ChainRecord> records;
  while (std::getline(fin, line)) {
    std::istringstream iss(line);
    if (line.empty()) {
      break;
    }
    ChainRecord record;
    iss >> record.len >> record.dq >> record.dt;
    records.push_back(record);
  }
  return records;
}
