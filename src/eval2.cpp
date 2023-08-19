// evaluation between two haplotype

#include <bits/stdc++.h>

#include <biovoltron/file_io/all.hpp>
#include <boost/program_options.hpp>
#include <boost/process.hpp>
#include <boost/system.hpp>
#include <spdlog/spdlog.h>
#include <spoa/spoa.hpp>

#include "utility/all.hpp"

namespace bpo = boost::program_options;
namespace bio = biovoltron;
namespace fs = std::filesystem;

const auto tmp = fs::temp_directory_path();

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
  check_amount_and_format(
      "snp_vcf", vmap["snp_vcf"].as<std::vector<fs::path>>(), FILE_FORMAT::VCF);
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


auto eval(const std::vector<bio::SamRecord<false>>& sam,
          const std::vector<bio::VcfRecord>& vcf) {
  assert(std::ranges::is_sorted(vcf, {}, &bio::VcfRecord::pos));
  ConfusionMatrix mat;
  /**
   * for each alignment, check variants that located at its template region
   * if the alignment result is match at variant site, then it can be considered
   * correct wrongly, otherwise, correct successfully.
   */

  auto get_vcf_range = [&](std::size_t pos, std::size_t len) {
    auto beg = std::ranges::lower_bound(vcf, pos, {}, &bio::VcfRecord::pos);
    auto end =
        std::ranges::upper_bound(vcf, pos + len, {}, &bio::VcfRecord::pos);
    return std::ranges::subrange(beg, end);
  };

  for (const auto& s : sam) {
    // the pos in vcf and sam is both 1-based, so there's no need for adjust it.
    auto offset = cal_offset(s.cigar);
    auto cigar_seg = std::vector<std::size_t>(s.cigar.size() + 1);
    cigar_seg[0] = 0;
    for (auto i = 1u; i <= cigar_seg.size(); i++) {
      cigar_seg[i] += cigar_seg[i - 1];
    }
    auto get_ele_from_cigar = [&](std::size_t pos) {
      auto it = std::ranges::lower_bound(cigar_seg, pos);
      if (*it != pos) {
        assert(it != cigar_seg.begin());
        it = std::prev(it);
      }
      auto d = std::ranges::distance(cigar_seg.begin(), it);
      assert(d != 0);
      return d - 1;
    };

    // ? is there a direct way to check the variant is match or not?
    // ? we have alignment result from read to haplotype, and variant on
    // ? haplotype 
    // ? current method is to check the variant position on cigar
    // ? string, and check ? the state is match or not.

    // ! this method is not correct, because the variant may be located at the
    // ! middle of a cigar element, for example, a indel variant has
    // ! cigar = 10M1I10M, and the variant is located at 8-11, then the
    // ! variant is not match at 8th position, but match at 11th position.
    // ! so we need to check the variant position on cigar string, and check
    // ! the state is match or not.

    for (const auto& v : get_vcf_range(s.pos, s.tlen)) {
      // get the variant position on cigar string, check the state is match or
      // not.
      auto v_beg_idx = get_ele_from_cigar(v.pos - s.pos);
      auto v_end_idx = get_ele_from_cigar(v.pos + v.alt.size() - 1 - s.pos);
      assert(v_end_idx < s.cigar.size());
      if (v_beg_idx != v_end_idx) {
        mat.fp += 1;
      } else if (s.cigar[v_beg_idx].op == 'M') {
        mat.tp += 1;
      } else {
        // TODO: what is this case
        // print variant and sequence
        spdlog::info("v.pos = {}, v.ref = {}, v.alt = {}", v.pos, v.ref, v.alt);
        spdlog::info("s.pos = {}, s.tlen = {}, s.cigar = {}", s.pos, s.tlen,
                     std::string(s.cigar));
        spdlog::info("[{} - {}] {}", v.pos - s.pos, v.pos - s.pos + 20,
                     s.seq.substr(v.pos - s.pos, 20));
      }
    }
  }
  return mat;
}

auto eval(const Haplotype& haplotype,
          const std::vector<bio::SamRecord<false>>& sam) {
  auto snp_vcf = haplotype.snp_vcf;
  auto indel_vcf = haplotype.indel_vcf;

  auto snp_mat = eval(sam, snp_vcf.second);
  auto indel_mat = eval(sam, indel_vcf.second);


}

int main(int argc, char* argv[]) {
  // TODO: help message

  // evaluation between two haplotype
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


  auto platform = vmap["platform"].as<std::string>();
  auto threads = vmap["threads"].as<std::size_t>();
  // for each corrected reads, run evaluation
  for (const auto& path : vmap["corrected_read"].as<std::vector<fs::path>>()) {
    spdlog::info("Evalation on {}", path.string());

    auto corrected_reads = read_records<bio::FastaRecord<false>>(path);
    /* seperate corrected reads according to prefix of read name */
    {
      for (auto p : std::views::iota(0u, ploidy)) {
        auto read_prefix = fmt::format("{}", p);
        std::ranges::for_each(corrected_reads, [&](auto& r) {
          if (r.name.starts_with(read_prefix)) {
            haplotypes[p].corrected_reads.emplace_back(std::move(r));
          }
        });
      }
    }

    /* create chain file between haplotype */
    auto chain_file = std::vector<std::vector<fs::path>>(ploidy, std::vector<fs::path>(ploidy));
    {
      auto haplotypes_path = vmap["reference"].as<std::vector<fs::path>>();
      for (auto i : std::views::iota(0u, ploidy)) {
        for (auto j : std::views::iota(0u, ploidy)) {
          if (i == j) {
            continue;
          }
          auto chain_path = tmp / fmt::format("{}_{}.chain", i, j);
          chain_file[i][j] = chain_path;
          spdlog::info("start exec flo for create {}", chain_path.string());
          boost::process::system(
            boost::process::search_path("flo"),
            haplotypes_path[i].string(),
            haplotypes_path[j].string(),
            chain_path.string(),
            boost::process::std_out > boost::process::null,
            boost::process::std_err > boost::process::null
          );
          spdlog::info("finish flo");
          assert(fs::exists(chain_path));
        }
      }
    }
    

    /* align corrected reads on corresponding haplotype */
    {
      for (auto p : std::views::iota(0u, ploidy)) {
        auto cread = tmp / fmt::format("{}.fasta", p);
        auto align_res_file = tmp / fmt::format("{}.sam", p);
        {
          auto fout = std::ofstream(cread);
          for (const auto& r : haplotypes[p].corrected_reads) {
            fout << r;
          }
        }
        {
          spdlog::info("Run minimap2 on {}", cread.string());        
          boost::process::system(
            boost::process::search_path("minimap2"),
            "-ax",
            platform == "ont" ? "map-ont" : "map-pb",
            "-t",
            std::to_string(threads),
            vmap["reference"].as<std::vector<fs::path>>()[p].string(),
            cread.string(),
            boost::process::std_out > align_res_file,
            boost::process::std_err > boost::process::null
          );    
        }
        {
          for (const auto& q : std::views::iota(0u, ploidy)) {
            if (q == p) {
              continue;
            }
            auto chain_path = chain_file[p][q];
            auto liftover_file = tmp / fmt::format("{}_to_{}.sam", p, q);
             /* liftover sam file to another haplotype */
            spdlog::info("Liftover on {} to {}", align_res_file.string(), liftover_file.string());
            boost::process::system(
              boost::process::search_path("CrossMap.py"),
              "bam",
              chain_path.string(),
              align_res_file.string(),
              boost::process::std_out > liftover_file,
              boost::process::std_err > boost::process::null
            );

            auto sam_header = bio::SamHeader{};
            auto sam_records = std::vector<bio::SamRecord<false>>{};
            {
              auto fin = std::ifstream(align_res_file);
              fin >> sam_header;
              while (fin) {
                auto record = bio::SamRecord<false>{};
                fin >> record;
                sam_records.push_back(record);
              }
            }
            spdlog::info("SamRecord size = {}", sam_records.size());
            eval(haplotypes[q], sam_records);
          } 
        }
      }
    }

  }

  // auto sam_file = vmap["sam"].as<fs::path>();
  // auto vcf_snp_file = vmap["snp"].as<fs::path>();
  // auto vcf_indel_file = vmap["indel"].as<fs::path>();

  // auto sam_header = bio::SamHeader{};
  // auto sam_records = std::vector<bio::SamRecord<false>>{};
  // {
  //   auto fin = std::ifstream(sam_file);
  //   fin >> sam_header;
  //   while (fin) {
  
  //     auto record = bio::SamRecord<false>{};
  //     fin >> record;
  //     sam_records.push_back(record);
  //   }
  // }

  // auto [snp_vcf_header, snp_vcf_records] = read_vcf(vcf_snp_file);
  // auto [indel_vcf_header, indel_vcf_records] = read_vcf(vcf_indel_file);

  // auto snp_res = eval(sam_records, snp_vcf_records);
  // auto indel_res = eval(sam_records, indel_vcf_records);

  // /* try using spoa */

  // auto alignment_engine =
  // spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -3); auto
  // graph = spoa::Graph{};

  // for (auto& seq : s) {
  //   auto alignment = alignment_engine->Align(seq, graph);
  //   graph.AddAlignment(alignment, seq);
  // }

  // auto consensus = graph.GenerateConsensus();
  // spdlog::info("consensus = {}", consensus);

  // auto msa = graph.GenerateMultipleSequenceAlignment();
  // for (auto& it : msa) {
  //   spdlog::info("msa seq = {}", it);
  // }

  // auto alignment_engine = spoa::createAlignmentEngine(
  //     static_cast<spoa::AlignmentType>(spoa::AlignmentType::kSW), 2, -1, -1,
  //     0);

  // fs::path chain =
  // "/mnt/ec/ness/yolkee/thesis/data/Ecoli/K12/run/liftover.chn"; fs::path vcf
  // = "/mnt/ec/ness/yolkee/thesis/data/Ecoli/K12/h1_liftover_snp.vcf";

  // auto chain_records = readChainFile(chain);
  // auto [vcf_header, vcf_record] = read_vcf(vcf);

  // // for (auto [len, dq, dt] : chain_records) {
  // //   spdlog::info("len = {}, dq = {}, dt = {}", len, dq, dt);
  // // }

  // fs::path ref =
  // "/mnt/ec/ness/yolkee/thesis/data/Ecoli/K12/GCF_000005845.2/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fa";
  // fs::path h1 =
  // "/mnt/ec/ness/yolkee/thesis/data/Ecoli/K12/snp_ref/h1.simseq.genome.fa";

  // auto ref_seq = read_records<bio::FastaRecord<false>>(ref);
  // auto h1_seq = read_records<bio::FastaRecord<false>>(h1);

  // auto print_genome = [&](std::size_t tpos, std::size_t qpos) {
  //   const std::size_t view_len = 25;
  //   auto ref_subseq = ref_seq[0].seq.substr(tpos - view_len, view_len * 2);
  //   auto h1_subseq = h1_seq[0].seq.substr(qpos - view_len, view_len * 2);

  //   spdlog::info("st = {}, ref = {}", tpos - view_len, ref_subseq);
  //   spdlog::info("st = {}, h1  = {}", qpos - view_len, h1_subseq);
  // };

  // spdlog::info("chain_record size = {}", chain_records.size());
  // spdlog::info("vcf_record size = {}", vcf_record.size());

  // std::size_t tlen = 0, qlen = 0;
  // int dlen = 0;
  // for (auto i = 0u, cidx = 0u; i < vcf_record.size(); i++) {
  //   auto variant = vcf_record[i];
  //   auto vpos = variant.pos;
  //   auto vref = variant.ref;
  //   auto valt = variant.alt;

  //   while (vpos > tlen + chain_records[cidx].len) {
  //     dlen += chain_records[cidx].dq - chain_records[cidx].dt;
  //     tlen += chain_records[cidx].len;
  //     cidx++;
  //   }

  //   if (vref == valt) {
  //     spdlog::info("vpos = {}, dlen = {}", vpos, dlen);
  //     print_genome(vpos, vpos - dlen);
  //     exit(-1);
  //   }
  // }

  // // std::size_t tlen = 0, qlen = 0;
  // for (int i = 0; i < chain_records.size(); i++) {
  //   auto [len, dt, dq] = chain_records[i];
  //   auto variant = vcf_record[i];
  //   auto vpos = variant.pos;
  //   auto vref = variant.ref;
  //   auto valt = variant.alt;
  //   auto pref_len = std::min(vref.size(), valt.size());

  //   tlen += len;
  //   qlen += len;

  //   if (vref == valt) {
  //     spdlog::info("tlen = {}, qlen = {}", tlen, qlen);
  //     print_genome(vpos, vpos);
  //   }

  //   if (vref.size() - valt.size() != dt - dq) {
  //     spdlog::info("i = {}, vpos = {}", i, vpos);
  //     spdlog::info("tlen = {}, qlen = {}", tlen, qlen);
  //     spdlog::error("ref.size() - alt.size() != dt - dq");
  //     spdlog::error("ref = {}, alt = {}, dt = {}, dq = {}", vref, valt, dt,
  //     dq);

  //     print_genome(tlen, qlen);
  //     std::exit(1);
  //   }

  //   tlen += dt;
  //   qlen += dq;

  // }
}