#include <spdlog/spdlog.h>

#include <biovoltron/file_io/all.hpp>
#include <boost/program_options.hpp>
#include <execution>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <omp.h>
#include <spoa/spoa.hpp>

#include "edlib.h"
#include "utility/all.hpp"

namespace bpo = boost::program_options;
namespace fs = std::filesystem;
namespace bio = biovoltron;

auto check_arguments(bpo::variables_map& vmap) {
  bpo::notify(vmap);
  auto ploidy = vmap["ploidy"].as<std::size_t>();
  // the amount of path must match the ploidy and the format must acceptable
  auto check_amount_and_format = [&](const std::string& opt_name,
                                     const std::vector<fs::path>& paths,
                                     const int accept_format) {
    if (paths.size() != ploidy) {
      throw bpo::validation_error(
          bpo::validation_error::invalid_option_value, opt_name,
          fmt::format("Invalid number of arguments, expected {} arguments",
                      ploidy));
    }
    for (const auto& path : paths) {
      if (!fs::exists(path)) {
        throw bpo::validation_error(
            bpo::validation_error::invalid_option_value, opt_name,
            fmt::format("File {} does not exist", path.string()));
      }
      check_file_format(opt_name, path, accept_format);
    }
  };

  check_file_format("raw_read", vmap["raw_read"].as<fs::path>(),
                    FILE_FORMAT::FASTQ);
  for (auto& read : vmap["corrected_read"].as<std::vector<fs::path>>()) {
    check_file_format("corrected_read", read, FILE_FORMAT::FASTA);
  }
  check_amount_and_format("reference",
                          vmap["reference"].as<std::vector<fs::path>>(),
                          FILE_FORMAT::FASTA);
  check_amount_and_format("maf", vmap["maf"].as<std::vector<fs::path>>(),
                          FILE_FORMAT::MAF);
  check_amount_and_format("snp_vcf", 
                          vmap["snp_vcf"].as<std::vector<fs::path>>(), 
                          FILE_FORMAT::VCF);
  check_amount_and_format("indel_vcf",
                          vmap["indel_vcf"].as<std::vector<fs::path>>(),
                          FILE_FORMAT::VCF);
  auto platform = vmap["platform"].as<std::string>();
  if (platform != "ONT" && platform != "PacBio") {
    throw bpo::validation_error(bpo::validation_error::invalid_option_value,
                                "platform",
                                "Invalid platform, only accept ONT or PacBio");
  }
}

auto align(std::string_view ref, std::string_view read) {
  auto result = edlibAlign(
      ref.data(), ref.size(), read.data(), read.size(),
      edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
  auto cigar = (edlibAlignmentToCigar(result.alignment, result.alignmentLength,
                                      EDLIB_CIGAR_EXTENDED));
  edlibFreeAlignResult(result);
  auto r = bio::Cigar(cigar);
  free(cigar);
  return r;
};

/**
 * @brief Use for calculate the relative offset between reference and read, the
 * cigar string must be the alignment result from **read to reference**.
 *
 * @param cigar Cigar string between reference and read
 * @return std::vector<std::pair<std::size_t, int>>
 */
auto cal_offset(const bio::Cigar& cigar, bool rev = false) {
  std::vector<std::pair<std::size_t, std::int64_t>> offset;
  offset.emplace_back(0u, 0);
  for (auto i = 0u, haplotype_pos = 0u; i < cigar.size(); i += 1) {
    auto dx = 0;
    if (cigar[i].op == 'I') {
      dx = cigar[i].size;
    } else {
      if (cigar[i].op == 'D') {
        dx = -cigar[i].size;
      }
      haplotype_pos += cigar[i].size;
    }
    if (dx) {
      dx += offset.back().second;
      offset.emplace_back(haplotype_pos, dx);
    }
  }
  if (rev) {
    for (auto& [p, k] : offset) {
      p += k;
      k = -k;
    }
  }
  return offset;
}

/**
 * @brief Calculate the coordinate transformation from reference genome to
 * haplotype by using vcf contains only indels.
 * @param indels vcf files contains only indels
 * @param rev
 * @return std::vector<std::pair<std::size_t, int>>
 */
auto cal_offset(const std::vector<bio::VcfRecord>& indels, bool rev = false) {
  std::vector<std::pair<std::size_t, std::int64_t>> offset;
  offset.emplace_back(0u, 0);
  /**
   * NC_000913.3	8216	.	C	CGT	.	.
   * variant_type=INDEL; NC_000913.3	9190	.	A	ACA
   * .	. variant_type=INDEL;
   */
  // for (auto ref_pos = 0u; auto &v : indels) {
  //   auto dx = std::ranges::ssize(v.alt) - std::ranges::ssize(v.ref);
  //   if (rev) {
  //     auto hpos = v.pos - offset.back().second;
  //     offset.emplace_back(hpos, offset.back().second + dx);
  //   } else {
  //     offset.emplace_back(v.pos + dx, offset.back().second - dx);
  //   }
  // }
  for (auto& v : indels) {
    // ssize -> signed size()
    auto dx = std::ranges::ssize(v.alt) - std::ranges::ssize(v.ref);
    offset.emplace_back(v.pos + dx, offset.back().second - dx);
  }
  if (rev) {
    for (auto& [p, k] : offset) {
      p += k;
      k = -k;
    }
  }
  return offset;
}

auto transform_coordinate(
    const std::vector<std::pair<std::size_t, std::int64_t>>& offset,
    std::size_t old_pos) {
  auto it = std::ranges::upper_bound(
      offset, old_pos, {}, &std::pair<std::size_t, std::int64_t>::first);
  assert(it != offset.begin());
  return old_pos + std::prev(it)->second;
}

auto transform_pair_coordinate(
    const std::vector<std::pair<std::size_t, std::int64_t>>& offset,
    std::size_t st, std::size_t ed) {
  auto it_st = std::ranges::upper_bound(
      offset, st, {}, &std::pair<std::size_t, std::int64_t>::first);
  auto it_ed = std::ranges::lower_bound(
      offset, ed, {}, &std::pair<std::size_t, std::int64_t>::first);
  assert(it_st != offset.begin() && it_ed != offset.begin());
  return std::make_pair(st + std::prev(it_st)->second,
                        ed + std::prev(it_ed)->second);
}

auto eval_cigar(const bio::Cigar& cigar) {
  constexpr auto eval_len = 2000ul;
  auto clen = 0ul, ac = 0ul;
  for (const auto& [len, op] : cigar | std::views::reverse) {
    if (op == '=') {
      ac += std::min(eval_len - clen, (std::size_t) len);
    }
    clen += len;
    if (clen >= eval_len) {
      break;
    }
  }
  return ac / (double) std::min(eval_len, clen);
}

/**
 * @brief Correct result for each read
 * @todo pretty information
 */
struct CorrectedResult {

  CorrectedResult() = default;

  CorrectedResult(std::string ref, std::string raw_seq,
                  std::string corrected_seq)
      : ref(std::move(ref)), raw_seq(std::move(raw_seq)),
        corrected_seq(std::move(corrected_seq)) {}

  auto correction_eval() {
    auto alignment_engine =
        spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 5, -4, -8);
    alignment_engine->Prealloc(
        std::max({ref.size(), raw_seq.size(), corrected_seq.size()}), 128);

    // auto alignment_engine =
    //     spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -3);
    // 3 -> match score
    // -5 -> mismatch score
    // -3 -> gap score

    auto graph = spoa::Graph{};

    auto add_seq_into_msa = [&](const std::string& seq) {
      auto alignment = alignment_engine->Align(seq, graph);
      graph.AddAlignment(alignment, seq);
    };
    add_seq_into_msa(ref);
    add_seq_into_msa(corrected_seq);
    add_seq_into_msa(raw_seq);

    auto msa = graph.GenerateMultipleSequenceAlignment();
    assert(msa.size() == 3);
    auto sz = msa[0].size();

    // spdlog::debug("sz = {}, {}", ref.size(), ref.substr(0, 100));
    // spdlog::debug("sz = {}, {}", raw_seq.size(), raw_seq.substr(0, 100));
    // spdlog::debug("sz = {}, {}", corrected_seq.size(),
    //               corrected_seq.substr(0, 100));

    // for (auto &s : msa) {
    //   spdlog::debug("{}", s.substr(0, 100));
    // }

    for (auto i = 0u; i < sz; i++) {
      auto c1 = msa[0][i]; // ref
      auto c2 = msa[1][i]; // corrected
      auto c3 = msa[2][i]; // raw

      if (c1 == c3) {
        /* don't correct */
        if (c1 == c2) {
          mat.tn += 1;
        } else {
          mat.fn += 1;
        }
      } else {
        /* should be corrected */
        if (c1 == c2) {
          mat.tp += 1;
        } else {
          mat.fp += 1;
        }
      }
    }
  }

  /**
   * @brief
   *
   * @tparam R range for VcfRecord
   * @param vcf std::range::subrange for VcfRecord
   * @return auto
   * ? is it work for indel or just snp
   */
  template <std::ranges::range R> auto variant_eval(R vcf, bio::MafRecord& maf) {
    
    auto print_seq = [](std::string_view s, std::size_t pos) {
      const auto len = 25;
      auto st = std::max(0ul, pos - len);
      auto ed = std::min(s.size(), pos + len);
      spdlog::debug("{:6} : {}", st, s.substr(st, ed - st));
    };
    
    auto cigar = align(ref, corrected_seq);
    auto offset = cal_offset(cigar, true);
    total_covered_variants += std::ranges::size(vcf);
    
    // spdlog::debug("ref = {}", ref.substr(0, 100));
    // spdlog::debug("cor = {}", corrected_seq.substr(0, 100));
    // auto z = transform_coordinate(offset, 50);
    // print_seq(ref, 50);
    // spdlog::debug("{:>36}", '^');
    // print_seq(corrected_seq, z);
    // spdlog::debug("cigar = {}", std::string(cigar));

    for (const auto& v : vcf) {
      /* get variant start and end related on haplotype coordinate system */
      /* change all coordinate to 0-based */
      auto vinfo = parse_vcf_info(v);
      auto st = std::stoul(std::string(vinfo["sim_start"])) - 1 - maf.ref.start;
      auto ed = std::stoul(std::string(vinfo["sim_end"])) - 1 - maf.ref.start;
      // print_seq(ref, st);
      st = transform_coordinate(offset, st);
      ed = transform_coordinate(offset, ed);
      if (st >= corrected_seq.size()) {
        v_out_of_range += 1;
        continue;
      }
      // spdlog::debug("{:>36}", '^');
      // print_seq(corrected_seq, st);
      auto alt_in_cread = corrected_seq.substr(st, ed - st + 1);
      if (v.alt != alt_in_cread) {
        v_mismatch += 1;
      }
    }
  }

  auto print() {
    std::cout << std::fixed << std::setprecision(4) << std::endl;
    std::cout << "correction result: " << std::endl;
    std::cout << "Use elector please" << std::endl << std::endl;
    // std::cout << "accuracy: " << mat.accuracy() << std::endl;
    // std::cout << "precision: " << mat.precision() << std::endl;
    // std::cout << "senstivity: " << mat.sensitivity() << std::endl;
    // std::cout << "specificity: " << mat.specificity() << std::endl <<
    // std::endl;

    std::cout << "variant correction result: " << std::endl;
    std::cout << "v_mismatch: " << v_mismatch << std::endl;
    std::cout << "v_out_of_range: " << v_out_of_range << std::endl;
    std::cout << "total_variants: " << total_covered_variants << std::endl;
    std::cout << "variant correct successful rate: "
              << 1 - (double) (v_mismatch + v_out_of_range) /
                         total_covered_variants
              << std::endl;
  }

  CorrectedResult operator+(const CorrectedResult& rhs) const {
    CorrectedResult ret;
    ret.mat = mat + rhs.mat;
    ret.v_mismatch = v_mismatch + rhs.v_mismatch;
    ret.v_out_of_range = v_out_of_range + rhs.v_out_of_range;
    ret.total_covered_variants =
        total_covered_variants + rhs.total_covered_variants;
    for (auto i = 0u; i < MAX_STATE_LEN; ++i) {
      ret.e_insertion[i] = e_insertion[i] + rhs.e_insertion[i];
      ret.e_deletion[i] = e_deletion[i] + rhs.e_deletion[i];
      ret.e_mismatch[i] = e_mismatch[i] + rhs.e_mismatch[i];
    }
    return ret;
  }

  CorrectedResult& operator+=(const CorrectedResult& rhs) {
    mat += rhs.mat;
    v_mismatch += rhs.v_mismatch;
    v_out_of_range += rhs.v_out_of_range;
    total_covered_variants += rhs.total_covered_variants;
    for (auto i = 0u; i < MAX_STATE_LEN; ++i) {
      e_insertion[i] += rhs.e_insertion[i];
      e_deletion[i] += rhs.e_deletion[i];
      e_mismatch[i] += rhs.e_mismatch[i];
    }
    return *this;
  }

  std::string ref;
  std::string raw_seq;
  std::string corrected_seq;

  ConfusionMatrix<> mat;
  std::size_t v_mismatch = 0u;
  std::size_t v_out_of_range = 0u;
  std::size_t total_covered_variants = 0u;

  constexpr static auto MAX_STATE_LEN = 10u;
  std::array<std::size_t, MAX_STATE_LEN> e_insertion = {};
  std::array<std::size_t, MAX_STATE_LEN> e_deletion = {};
  std::array<std::size_t, MAX_STATE_LEN> e_mismatch = {};
};

auto eval(const Haplotype& haplotype) {
  auto ref = std::string_view(haplotype.reference.seq);

  /* evaluate variant */
  auto v_mismatch = 0u;
  auto v_out_of_range = 0u;
  auto v_total = 0u;

  /* evaluate correction result */
  auto total_acc = (double) 0, total_base = (double) 0;
  auto e_mismatch = std::map<std::size_t, std::atomic_int>{};
  auto e_insertion = std::map<std::size_t, std::atomic_int>{};
  auto e_deletion = std::map<std::size_t, std::atomic_int>{};

  /**
   * For each read, scan the variant that covered by the read. Check whether
   * the variant is fixed correctly or not.
   *
   * Also, evaluate the correction result in front side and rear side of read
   */

  std::vector<CorrectedResult> results(haplotype.corrected_reads.size());
  spdlog::info("Corrected read size = {}", results.size());
  auto tp = bio::make_threadpool(omp_get_max_threads());

  auto job = [&](std::size_t i) {
    auto cread = haplotype.corrected_reads[i];
    auto result = CorrectedResult{};
    auto it = std::ranges::upper_bound(haplotype.raw_reads, cread.name, {},
                                       &bio::FastqRecord<false>::name);
    it = std::ranges::prev(it);
    // spdlog::debug("corrected read name: {}", read.name);
    // spdlog::debug("raw read name: {}", it->name);

    auto maf = get_maf_record(haplotype.maf.second, it->name);
    auto ref_seq = ref.substr(maf.ref.start, maf.ref.len);
    auto corrected_seq =
        maf.aln[0].strand == '+' ? cread.seq : bio::Codec::rev_comp(cread.seq);

    auto res = CorrectedResult(std::string(ref_seq), it->seq, corrected_seq);
    auto snp_vcf_range = get_vcf_range(haplotype.snp_vcf.second, maf.ref.start,
                                       maf.ref.start + maf.ref.len - 1);
    auto indel_vcf_range =
        get_vcf_range(haplotype.indel_vcf.second, maf.ref.start,
                      maf.ref.start + maf.ref.len - 1);
    // res.correction_eval();
    res.variant_eval(snp_vcf_range, maf);
    res.variant_eval(indel_vcf_range, maf);
    // res.variant_eval(get_vcf_range(haplotype.indel_vcf.second, m.ref.start,
    //                                m.ref.start + m.ref.len - 1),
    //                  m);

    return res;
  };

  std::vector<std::future<CorrectedResult>> futures;
  for (auto i = 0u; i < results.size(); ++i) {
    auto [_, res] = tp.submit(job, i);
    futures.emplace_back(std::move(res));
  }
  CorrectedResult sum;
  for (auto& f : futures) {
    auto res = f.get();
    sum += res;
  }
  sum.print();
  auto not_corrected_reads =
      haplotype.raw_reads.size() - haplotype.corrected_reads.size();
  spdlog::info("Not corrected read = {}", not_corrected_reads);
  spdlog::info("Not corrected precent = {:.2f}%",
               not_corrected_reads / (double) haplotype.raw_reads.size() * 100);

  return sum;

  // #pragma omp parallel for reduction(+ : v_mismatch, v_out_of_range, v_total,
  //                                        total_acc, total_base)
  //   for (const auto& read : haplotype.corrected_reads) {
  //     auto m = get_maf_record(haplotype.maf.second, read.name);
  //     auto ref_seq = ref.substr(m.ref.start, m.ref.len);
  //     auto read_seq =
  //         m.aln[0].strand == '+' ? read.seq : bio::Codec::rev_comp(read.seq);

  //     auto cigar = align(read_seq, ref_seq);
  //     for (const auto [len, op] : cigar) {
  //       if (op == 'X') {
  //         e_mismatch[len] += 1;
  //       } else if (op == 'I') {
  //         e_insertion[len] += 1;
  //       } else if (op == 'D') {
  //         e_deletion[len] += 1;
  //       }
  //     }
  //     // spdlog::debug("cigar = {}", std::string(cigar));

  //     total_acc += eval_cigar(cigar);
  //     total_base += m.ref.len;

  //     /* calculate the offset based on alignment result(cigar string) */
  //     auto move_offset = cal_offset(cigar);

  //     for (auto v : get_vcf_range(haplotype.snp_vcf.second, m.ref.start,
  //                                 m.ref.start + m.ref.len - 1)) {
  //       /* get variant start and end related on haplotype coordinate system
  //       */
  //       /* change all coordinate to 0-based */
  //       auto vinfo = parse_vcf_info(v);
  //       auto st = std::stoul(std::string(vinfo["sim_start"])) - 1 -
  //       m.ref.start; auto ed = std::stoul(std::string(vinfo["sim_end"])) - 1
  //       - m.ref.start; st = transform_coordinate(move_offset, st); ed =
  //       transform_coordinate(move_offset, ed); v_total += 1; if (st >=
  //       read_seq.size()) {
  //         v_out_of_range += 1;
  //         break;
  //       }
  //       auto alt_in_cread = read_seq.substr(st, ed - st + 1);
  //       if (v.alt != alt_in_cread) {
  //         v_mismatch += 1;
  //       }
  //     }
  //   }

  //   // spdlog::info("Total {} bases", total_base);
  //   auto sum = 0ull;
  //   auto print_error = [&](const std::map<std::size_t, std::atomic_int>& e) {
  //     const int k = 10;
  //     std::vector<int> v(k + 1);
  //     for (const auto& [len, cnt] : e) {
  //       sum += 1ull * cnt * len;
  //       if (len >= k) {
  //         v[k] += cnt;
  //       } else {
  //         v[len] += cnt;
  //       }
  //     }
  //     for (int i = 1; i <= k; i++) {
  //       spdlog::info("{}\t{}", i, v[i]);
  //     }
  //   };
  //   spdlog::info("Mismatch");
  //   print_error(e_mismatch);
  //   spdlog::info("Insertion");
  //   print_error(e_insertion);
  //   spdlog::info("Deletion");
  //   print_error(e_deletion);
  //   spdlog::info("Error rate = {}", (double) sum / total_base);

  //   spdlog::info("Total {} corrected reads that belong to this haplotype",
  //                haplotype.corrected_reads.size());
  //   spdlog::info("Covered {} variants", v_total);
  //   spdlog::info(
  //       "Mismatch {} variants(corrected result not match alt in vcf file)",
  //       v_mismatch);
  //   spdlog::info("Out of range {} variants(at tail that unavilable to
  //   align)",
  //                v_out_of_range);
  //   spdlog::info("Accuracy at len 10 = {}",
  //                (double) total_acc / haplotype.corrected_reads.size());
}

int main(int argc, char* argv[]) {
  spdlog::set_level(spdlog::level::trace);
  /**
   * input:
   * 1. corrected read
   * 2. original read
   * 3. reference genome
   * 4. vcf of reference genome
   *
   *
   * TODO:
   * 1. calculate the number of different variant type
   * 2. for each variant type, calculate its precision, recall, f1-score
   *
   */

  // TODO: add threads option
  bpo::options_description opts{};
  bpo::variables_map vmap;
  try {
    opts.add_options()("help,h", "Show help message")(
        "ploidy,n", bpo::value<std::size_t>()->required())(
        "raw_read,o", bpo::value<fs::path>()->required()->notifier(
                          make_path_checker("raw_read")))(
        "corrected_read,c",
        bpo::value<std::vector<fs::path>>()->required()->multitoken())(
        "reference,r",
        bpo::value<std::vector<fs::path>>()->required()->multitoken())(
        "maf", bpo::value<std::vector<fs::path>>()->required()->multitoken())(
        "snp_vcf",
        bpo::value<std::vector<fs::path>>()->required()->multitoken())(
        "indel_vcf",
        bpo::value<std::vector<fs::path>>()->required()->multitoken())(
        "platform,p", bpo::value<std::string>()->required())(
        "threads,t", bpo::value<std::size_t>()->default_value(1));
    bpo::store(bpo::parse_command_line(argc, argv, opts), vmap);
    check_arguments(vmap);
  } catch (std::exception& ex) {
    spdlog::error("{}", ex.what());
    std::cerr << opts << std::endl;
    exit(-1);
  }

  auto ploidy = vmap["ploidy"].as<std::size_t>();
  auto haplotypes = std::vector<Haplotype>(ploidy);
  auto raw_reads = std::vector<bio::FastqRecord<false>>{};
  auto threads = vmap["threads"].as<std::size_t>();

  omp_set_num_threads(threads);

  /* --- read data --- */
  for (auto i : std::views::iota(0u, ploidy)) {
    auto ref_path = vmap["reference"].as<std::vector<fs::path>>()[i];
    auto maf_path = vmap["maf"].as<std::vector<fs::path>>()[i];
    auto snp_path = vmap["snp_vcf"].as<std::vector<fs::path>>()[i];
    auto indel_path = vmap["indel_vcf"].as<std::vector<fs::path>>()[i];

    haplotypes[i].reference =
        read_records<bio::FastaRecord<false>>(ref_path)[0];
    haplotypes[i].maf = read_maf(maf_path);
    haplotypes[i].snp_vcf = read_vcf(snp_path);
    haplotypes[i].indel_vcf = read_vcf(indel_path);
  }
  /* --- data preprocessing --- */
  /* for vechat, it will add "rr" suffix to corrected read name, holy moly!!! */
  // std::ranges::for_each(corrected_reads, [](auto &r) {
  //   r.name = r.name.substr(0, r.name.size() - 2);
  // });
  // std::ranges::sort(raw_reads, {}, &bio::FastqRecord<false>::name);
  // std::ranges::sort(corrected_reads, {}, &bio::FastaRecord<false>::name);
  std::ranges::for_each(haplotypes, [](auto& h) {
    std::ranges::sort(h.maf.second, [](const auto& lhs, const auto& rhs) {
      return lhs.aln[0].name < rhs.aln[0].name;
    });
  });

  auto raw_read_path = vmap["raw_read"].as<fs::path>();
  raw_reads = read_records<bio::FastqRecord<false>>(raw_read_path);
  std::sort(
      std::execution::par, raw_reads.begin(), raw_reads.end(),
      [](const auto& lhs, const auto& rhs) { return lhs.name < rhs.name; });
  for (auto p : std::views::iota(0u, ploidy)) {
    auto read_prefix = std::to_string(p);
    std::ranges::for_each(raw_reads, [&](auto& r) {
      if (r.name.starts_with(read_prefix)) {
        haplotypes[p].raw_reads.emplace_back(std::move(r));
      }
    });
  }
  for (const auto& path : vmap["corrected_read"].as<std::vector<fs::path>>()) {
    spdlog::info("Evalation on {}", path.string());
    auto corrected_reads = read_records<bio::FastaRecord<false>>(path);
    /* seperate raw and corrected reads according to prefix of read name */
    for (auto p : std::views::iota(0u, ploidy)) {
      auto read_prefix = fmt::format("{}", p);
      haplotypes[p].corrected_reads.clear();
      std::ranges::for_each(corrected_reads, [&](auto& r) {
        if (r.name.starts_with(read_prefix)) {
          haplotypes[p].corrected_reads.emplace_back(std::move(r));
        }
      });
    }
    std::vector<CorrectedResult> results(ploidy);
    for (const auto p : std::views::iota(0u, ploidy)) {
      spdlog::info("Evalation on haplotype {}", p);
      auto res = eval(haplotypes[p]);
      results[p] = res;
    }
    CorrectedResult sum;
    for (const auto& res : results) {
      sum += res;
    }
    sum.print();
  }
}