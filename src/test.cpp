#include <bits/stdc++.h>
#include <execution>

#include <biovoltron/algo/align/inexact_match/smithwaterman_sse.hpp>
#include <biovoltron/utility/istring.hpp>
#include <biovoltron/file_io/all.hpp>
#include <boost/program_options.hpp>

#include <edlib.h>
#include <spoa/spoa.hpp>
#include <spdlog/spdlog.h>

#include "utility/all.hpp"
#include "graph/read_graph.hpp"

namespace bpo = boost::program_options;
namespace bio = biovoltron;
namespace fs = std::filesystem;

template<class T>
auto write_records(const fs::path& path, const std::vector<T> records) {
  auto fout = std::ofstream(path);
  if (!fout.is_open()) {
    spdlog::error("Cannot open file {}", path.string());
    exit(-1);
  }
  for (auto& it : records) {
    fout << it << '\n';
  }
}

int main(int argc, char* argv[]) {
  spdlog::set_level(spdlog::level::debug);
  
  fs::path tmp_dir = "/mnt/ec/ness/yolkee/thesis/tests/tmp";

  fs::path my_overlap = tmp_dir / "filtered_overlap.paf";
  fs::path vechat_overlap = tmp_dir / "overlap1.paf";

  auto my_overlap_records = read_records<bio::PafRecord>(my_overlap);
  auto vechat_overlap_records = read_records<bio::PafRecord>(vechat_overlap);

  spdlog::info("my_overlap_records size = {}", my_overlap_records.size());
  spdlog::info("vechat_overlap_records size = {}", vechat_overlap_records.size());
  sort(my_overlap_records.begin(), my_overlap_records.end());
  sort(vechat_overlap_records.begin(), vechat_overlap_records.end());

  fs::path a = tmp_dir / "a.paf";
  fs::path b = tmp_dir / "b.paf";

  write_records(a, my_overlap_records);
  write_records(b, vechat_overlap_records);
  

  // assert(argc == 3);
  // const int kmer = std::atoi(argv[1]);
  // const int min_occ = std::atoi(argv[2]);

  // fs::path dir = "/mnt/ec/ness/yolkee/thesis/data/Ecoli/K12/seqs";
  // for (auto d : fs::directory_iterator{dir}) {
  //   auto fin = std::ifstream(d.path());
  //   std::vector<std::string> seqs;
  //   auto len_sum = 0;
  //   for (std::string s; std::getline(fin, s);) {
  //     if (!s.empty()) {
  //       len_sum += s.size();
  //       seqs.emplace_back(std::move(s));
  //     }
  //   }    
  //   spdlog::info("Read {}, seqs size = {}, len_sum = {}", d.path().string(), seqs.size(), len_sum);
  //   ReadGraph g(kmer, min_occ);
  //   g.build(seqs);
  //   g.print();
  // }

  


  /**
   * argv[1] = sam file
   * argv[2] = vcf(snp)
   * argv[3] = vcf(indel)
   */

  // auto sam_file = "/mnt/ec/ness/yolkee/thesis/src/h2.sam";
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

  // auto vcf_file = "/mnt/ec/mammoth/yolkee/thesis/data/Ecoli/K13/snp_ref/h2.refseq2simseq.SNP.vcf";
  // auto vcf_header = bio::VcfHeader{};
  // auto vcf_records = std::vector<bio::VcfRecord>{};
  // {
  //   auto fin = std::ifstream(vcf_file);
  //   fin >> vcf_header;
  //   while (fin) {
  //     auto record = bio::VcfRecord{};
  //     fin >> record;
  //     vcf_records.push_back(record);
  //   }
  // }

 

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