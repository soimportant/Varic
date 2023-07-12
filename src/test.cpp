#include <bits/stdc++.h>

#include <biovoltron/algo/align/inexact_match/smithwaterman_sse.hpp>
#include <biovoltron/utility/istring.hpp>
#include <biovoltron/file_io/all.hpp>

#include <spoa/spoa.hpp>
#include <spdlog/spdlog.h>

#include "utility/all.hpp"

namespace bio = biovoltron;
namespace fs = std::filesystem;

int main(int argc, char* argv[]) {

  /**
   * argv[1] = sam file
   * argv[2] = vcf(snp)
   * argv[3] = vcf(indel)
   */

  auto sam_file = "/mnt/ec/ness/yolkee/thesis/src/h2.sam";
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

  auto vcf_file = "/mnt/ec/mammoth/yolkee/thesis/data/Ecoli/K13/snp_ref/h2.refseq2simseq.SNP.vcf";
  auto vcf_header = bio::VcfHeader{};
  auto vcf_records = std::vector<bio::VcfRecord>{};
  {
    auto fin = std::ifstream(vcf_file);
    fin >> vcf_header;
    while (fin) {
      auto record = bio::VcfRecord{};
      fin >> record;
      vcf_records.push_back(record);
    }
  }

  
  









  // std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
  // int n = 10;
  // int len = 20;
  // std::vector<std::string> s(n);
  // for (int i = 0; i < n; i++) {
  //   std::string k;
  //   for (int k = 0; k < len; k++) {
  //     s[i].push_back("ACGT"[rng() % 4]);
  //   }
  // }

  // /* try using spoa */
  
  // auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -3);
  // auto graph = spoa::Graph{};

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
  //     static_cast<spoa::AlignmentType>(spoa::AlignmentType::kSW), 2, -1, -1, 0);


  // fs::path chain = "/mnt/ec/ness/yolkee/thesis/data/Ecoli/K12/run/liftover.chn";
  // fs::path vcf = "/mnt/ec/ness/yolkee/thesis/data/Ecoli/K12/h1_liftover_snp.vcf";

  // auto chain_records = readChainFile(chain);
  // auto [vcf_header, vcf_record] = read_vcf(vcf);

  // // for (auto [len, dq, dt] : chain_records) {
  // //   spdlog::info("len = {}, dq = {}, dt = {}", len, dq, dt);
  // // }
  
  // fs::path ref = "/mnt/ec/ness/yolkee/thesis/data/Ecoli/K12/GCF_000005845.2/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fa";
  // fs::path h1 = "/mnt/ec/ness/yolkee/thesis/data/Ecoli/K12/snp_ref/h1.simseq.genome.fa";

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
  //     spdlog::error("ref = {}, alt = {}, dt = {}, dq = {}", vref, valt, dt, dq);

  //     print_genome(tlen, qlen);     
  //     std::exit(1);
  //   }

  //   tlen += dt;
  //   qlen += dq;

  // }
}