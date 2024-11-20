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

// auto foo(std::string s);

auto align_overlap(std::string_view s, std::string_view t) {
  // align overlaps with edlib
  EdlibAlignResult result = edlibAlign(
    s.data(), s.size(), t.data(), t.size(),
      edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));

  bio::Cigar cigar;
  if (result.status == EDLIB_STATUS_OK) {
    char* tmp_cigar = edlibAlignmentToCigar(
        result.alignment, result.alignmentLength, EDLIB_CIGAR_EXTENDED);
    cigar = tmp_cigar;
    std::free(tmp_cigar);
  }
  edlibFreeAlignResult(result);
  return cigar;
}

int main(int argc, char* argv[]) {
  spdlog::set_level(spdlog::level::debug);
  const fs::path read_path = "/mnt/ec/ness/yolkee/thesis/data/Ecoli/K12/reads/ONT/D20/merged_reads.fastq";
  const fs::path overlap_path = "/mnt/ec/ness/yolkee/thesis/data/Ecoli/K12/reads/ONT/D20/merged_reads.overlap.paf";

  auto reads = read_records<bio::FastqRecord<false>>(read_path);
  auto overlaps = read_records<bio::PafRecord>(overlap_path);
  std::sort(std::execution::par, reads.begin(), reads.end(), [](auto& lhs, auto& rhs) {
    return lhs.name < rhs.name;
  });
  std::sort(std::execution::par, overlaps.begin(), overlaps.end(), [](auto& lhs, auto& rhs) {
    return lhs.qname < rhs.qname;
  });
  auto name2id = std::unordered_map<std::string_view, std::size_t>{};
  for (auto i = 0u; i < reads.size(); i++) {
    name2id[reads[i].name] = i;
  }
  auto get_overlap_range = [&](const std::string_view qread_name) {
    auto subrange = std::ranges::equal_range(overlaps, qread_name, std::ranges::less{}, &bio::PafRecord::qname);
    return subrange;
  };

  spdlog::debug("reads size = {}", reads.size());  
  for (const auto& read : reads | std::views::take(1)) {
    spdlog::debug("read name = {}", read.name);
    fs::path dir = "/mnt/ec/ness/yolkee/thesis/tests/tmp/window_cigar";
    for (auto& overlap : get_overlap_range(read.name)) {
      std::stringstream ss;
      ss << overlap;

      auto qlen = overlap.qend - overlap.qstart + 1;
      auto tlen = overlap.tend - overlap.tstart + 1;
      auto error = 1.0 - (double) std::min(qlen, tlen) / std::max(qlen, tlen);
      if (error > 0.3) {
        continue;
      }
      fs::path file = dir / fmt::format("{}({})x{}({}).txt", overlap.qname, qlen, overlap.tname, tlen);
      std::ofstream fout(file);

      auto qseq = reads[name2id[overlap.qname]].seq.substr(overlap.qstart, qlen);
      auto tseq = reads[name2id[overlap.tname]].seq.substr(overlap.tstart, tlen);
      if (overlap.strand == '-') {
        tseq = bio::Codec::rev_comp(tseq);
      }
      
      auto cigar = align_overlap(qseq, tseq);
      auto extended_cigar = std::string{};
      for (auto& [len, op] : cigar) {
        extended_cigar += std::string(len, op);
      }
      const std::size_t view_len = 100;
      for (int l = 0, r = 0, len = 0; l < extended_cigar.size(); l++) {
        while (r < extended_cigar.size() && len < view_len) {
          if (extended_cigar[r] != 'I') {
            len++;
          }
          r++;
        }
        int match = 0, insertion = 0, deletion = 0;
        for (int k = l; k < r; k++) {
          if (extended_cigar[k] == 'M' || extended_cigar[k] == 'X') {
            match++;
          } else if (extended_cigar[k] == 'I') {
            insertion++;
          } else if (extended_cigar[k] == 'D') {
            deletion++;
          }
        }
        // print ratio
        
        double match_ratio = match / (double)(match + insertion + deletion);
        double insertion_ratio = insertion / (double)(match + insertion + deletion);
        double deletion_ratio = deletion / (double)(match + insertion + deletion);
        fout << l << " " << match_ratio << " " << insertion_ratio << " " << deletion_ratio << "\n";
        if (extended_cigar[l] != 'I') {
          len--;
        }
      }
    }
  }

  



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