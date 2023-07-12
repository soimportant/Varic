// evaluation between two haplotype

#include <bits/stdc++.h>

#include <biovoltron/file_io/all.hpp>
#include <boost/program_options.hpp>
#include <spdlog/spdlog.h>
#include <spoa/spoa.hpp>

#include "utility/all.hpp"

namespace bio = biovoltron;
namespace fs = std::filesystem;

/**
 * @brief Use for calculate the relative offset between reference and read, the
 * cigar string must be the alignment result from **read to reference**.
 *
 * @param cigar Cigar string between reference and read
 * @return std::vector<std::pair<std::size_t, int>>
 */
auto cal_offset(const bio::Cigar &cigar, bool rev = false) {
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


auto eval(std::vector<bio::SamRecord<false>>& sam, std::vector<bio::VcfRecord>& vcf) {
  assert(std::ranges::is_sorted(vcf, {}, &bio::VcfRecord::pos));
  ConfusionMatrix mat;
  /**
   * for each alignment, check variants that located at its template region
   * if the alignment result is match at variant site, then it can be considered
   * correct wrongly, otherwise, it can be considered correct successfully.
   */
  
  auto get_vcf_range = [&](std::size_t pos, std::size_t len) {
    auto beg = std::ranges::lower_bound(vcf, pos, {}, &bio::VcfRecord::pos);
    auto end = std::ranges::upper_bound(vcf, pos + len, {}, &bio::VcfRecord::pos);
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
    for (const auto& v : get_vcf_range(s.pos, s.tlen)) {
      // get the variant position on cigar string, check the state is match or not.
      auto d = get_ele_from_cigar(v.pos - s.pos);
      assert(d < s.cigar.size());
      if (s.cigar[d].op == 'M') {
        mat.tp += 1;
      } else {
        mat.fp += 1;
      }
    }
  }

}

int main(int argc, char* argv[]) {

  /**
   * argv[1] = sam file
   * argv[2] = vcf(snp)
   * argv[3] = vcf(indel)
   */

  // TODO: help message
  // evaluation between two haplotype
  bpo::options_description opts{};
  bpo::variables_map vmap;
  try {
    opts.add_options()("help,h", "Show help message")(
        "s,sam", bpo::value<fs::path>()
          ->required()->notifier(make_path_checker("sam")),
          "SAM file, which contains the alignment result of h1 reads to h1 then liftover to h2")(
        "v,snp", bpo::value<fs::path>()
          ->required()->notifier(make_path_checker("snp")),
          "VCF file contains SNP variant in h2")(
        "w,indel", bpo::value<fs::path>()->required()
          ->notifier(make_path_checker("indel")), 
          "VCF file contains INDEL variant in h2");
    bpo::store(bpo::parse_command_line(argc, argv, opts), vmap);
    bpo::notify(vmap);
    if (vmap.contains("help")) {
      exit_and_print_help(opts);
    }
  } catch (const std::exception& ex) {
    spdlog::error("{}", ex.what());
    exit_and_print_help(opts);
  }

  auto sam_file = vmap["sam"].as<fs::path>();
  auto vcf_snp_file = vmap["snp"].as<fs::path>();
  auto vcf_indel_file = vmap["indel"].as<fs::path>();

  auto sam_header = bio::SamHeader{};
  auto sam_records = std::vector<bio::SamRecord<false>>{};
  {
    auto fin = std::ifstream(sam_file);
    fin >> sam_header;
    while (fin) {
      auto record = bio::SamRecord<false>{};
      fin >> record;
      sam_records.push_back(record);
    }
  }

  auto [snp_vcf_header, snp_vcf_records] = read_vcf(vcf_snp_file);
  auto [indel_vcf_header, indel_vcf_records] = read_vcf(vcf_indel_file);


  

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