#include <bits/stdc++.h>
#include <algorithm>
#include <execution>
#include <omp.h>
#include <edlib.h>
#include <spdlog/spdlog.h>

#include <biovoltron/algo/align/inexact_match/smithwaterman_sse.hpp>
#include <biovoltron/file_io/all.hpp>
#include <biovoltron/utility/istring.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/asio.hpp>
#include <filesystem>
#include <spoa/spoa.hpp>

// #include "thesis/algo/assemble/read_assembler.hpp"
#include "thesis/corrector/detail/sequence.hpp"
#include "thesis/corrector/detail/test_read.hpp"

#include "thesis/format/maf.hpp"
#include "thesis/format/paf.hpp"
#include "thesis/utility/file_io/read_record.hpp"
#include "thesis/utility/threadpool/threadpool.hpp"
// #include "utility/all.hpp"

namespace bpo = boost::program_options;
namespace bio = biovoltron;
namespace fs = std::filesystem;

using R = bio::FastqRecord<false>;
using Read = TestReadWrapper<R>;

/**
 * @brief semi-global alignment
 *
 * @param ref subsequence of reference genome
 * @param read corrected read
 * @return bio::Cigar
 */
auto semi_global_align(std::string_view ref, std::string_view read) {
  auto result = edlibAlign(
      read.data(), read.size(), ref.data(), ref.size(),
      edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
  auto cigar = (edlibAlignmentToCigar(result.alignment, result.alignmentLength,
                                      EDLIB_CIGAR_EXTENDED));
  assert(result.numLocations == 1);
  auto start = result.startLocations[0];
  auto end = result.endLocations[0];
  edlibFreeAlignResult(result);
  auto r = bio::Cigar(cigar);
  std::free(cigar);
  return std::make_tuple(r, start, end);
};

int main(int argc, char* argv[]) {
  spdlog::set_level(spdlog::level::debug);

  bpo::options_description opts{};
  bpo::variables_map vmap;
  
  try {
    opts.add_options()("thread,t", bpo::value<int>()->default_value(1),
                       "the maximum number of threads");
    bpo::store(bpo::parse_command_line(argc, argv, opts), vmap);
    // bpo::notify(vmap);
    // if (vmap.contains("help")) {
    //   exit_and_print_help(opts);
    // }
    // check_argument(vmap);
  } catch (const std::exception& ex) {
    spdlog::error("{}", ex.what());
    // exit_and_print_help(opts);
  }


  auto threads = vmap["thread"].as<int>();
  omp_set_num_threads(threads);
  
  fs::path raw_reads_path = DATA_PATH "/Ecoli/K12_2_QS/reads/ONT/D20/merged_reads.fastq";
  fs::path perfect_reads_path = DATA_PATH "/Ecoli/K12_2_QS/reads/ONT/D20/merged_perfect_reads.fasta";
  fs::path data_dir = "/mnt/ec/ness/yolkee/thesis/tmp/fragments";
  auto paths = std::vector<fs::path>{};
  for (auto d : fs::directory_iterator{data_dir}) {
    paths.emplace_back(d.path());
  }

  std::sort(paths.begin(), paths.end(), [&](const auto p1, const auto p2) {
    auto a = p1.stem().string();
    auto b = p2.stem().string();
    auto ah = std::stoi(a.substr(0, 2));
    auto bh = std::stoi(b.substr(0, 2));
    if (ah != bh) {
      return ah < bh;
    }
    auto an = std::stoi(a.substr(3));
    auto bn = std::stoi(b.substr(3));
    if (an != bn) {
      return an < bn;
    }
    return p1 < p2;
  });

  auto raw_reads = read_records<bio::FastqRecord<false>>(raw_reads_path);
  auto perfect_reads = read_records<bio::FastaRecord<false>>(perfect_reads_path);
  sort(raw_reads.begin(), raw_reads.end(),
       [](const auto& a, const auto& b) { return a.name < b.name; });
  sort(perfect_reads.begin(), perfect_reads.end(),
       [](const auto& a, const auto& b) { return a.name < b.name; });

  auto get_seq_size = [&](const std::string_view name) -> std::size_t {
    auto it = std::ranges::lower_bound(raw_reads, name, std::less{},
                                       &bio::FastqRecord<false>::name);
    if (it == raw_reads.end() || it->name != name) {
      spdlog::error("Cannot find read with name = {}", name);
      return 0;
    }
    return it->seq.size();
  };

  auto get_ans = [&](const std::string_view name) -> std::size_t {
    auto it = std::ranges::lower_bound(perfect_reads, name, std::less{},
                                       &bio::FastaRecord<false>::name);
    if (it == perfect_reads.end() || it->name != name) {
      spdlog::error("Cannot find read with name = {}", name);
      return 0;
    }
    return it->seq.size();
  };

  auto get_ans_seq = [&](const std::string_view name) -> std::string {
    auto it = std::ranges::lower_bound(perfect_reads, name, std::less{},
                                       &bio::FastaRecord<false>::name);
    if (it == perfect_reads.end() || it->name != name) {
      spdlog::error("Cannot find read with name = {}", name);
      return "";
    }
    return it->seq;
  };

  auto fragments = std::map<std::string, std::vector<Sequence<std::string>>>{};
  auto take = 3;
  auto sample_path = std::vector<fs::path>{};
  std::ranges::sample(paths, std::back_inserter(sample_path), take,
                      std::mt19937(std::random_device{}()));
  for (auto path : sample_path) {
    auto read_name = path.stem().string();
    // if (!wrong_reads.contains(read_name)) {
    //   continue;
    // }
    auto read_len = get_seq_size(read_name);
    auto perfect_read_len = get_ans(read_name);
    auto sequences = std::vector<Sequence<std::string>>{};

    auto fin = std::ifstream(path);
    assert(fin.is_open());

    for (std::string s; std::getline(fin, s);) {
      auto ss = std::stringstream{};
      ss << s;
      auto seq = Sequence<std::string>{};
      ss >> seq.left_bound >> seq.right_bound >> seq.seq;
      sequences.emplace_back(std::move(seq));
    }
#pragma omp critical
    fragments.emplace(read_name, std::move(sequences));
  }

  spdlog::debug("Read file done");
  spdlog::debug("reads size = {}", raw_reads.size());
  spdlog::debug("perfect_reads size = {}", perfect_reads.size());
  spdlog::debug("paths size = {}", paths.size());

  auto init_reads = [&](){
    auto reads = std::vector<Read>{};
    for (auto& read : raw_reads) {
      if (fragments.contains(read.name)) {
        reads.emplace_back(std::move(read));
      }
    }
    return reads;    
  };

  auto reads = init_reads();
  auto precorrect = 0;
  auto total_fragment = 0;
  auto saved_fragment = 0;
  auto total_len = 0;
  auto corrected_len = 0;

  for (auto& read : reads) {
    const auto& fragment = fragments[read.name];
    total_fragment += fragment.size();
    
    std::vector<int> order(fragment.size());
    spdlog::debug("Read {}, len = {}, ans = {}, total fragments = {}", read.name, read.len(), get_ans(read.name), fragment.size());
    std::iota(order.begin(), order.end(), 0);
    std::ranges::shuffle(order, std::mt19937(std::random_device{}()));
    auto flag = false;
    for (auto i = 0; i < order.size(); i++) {
      int x = order[i];
      read.add_corrected_fragment(fragment[x]);

      // if (read.ready_to_assemble()) {
      //   spdlog::debug("ready when adding {} fragments", i + 1);
      //   read.assemble_corrected_seq();
      //   if (read.corrected_seq.size() > 0) {
      //     flag = true;
      //     total_len += get_ans(read.name);
      //     saved_fragment += fragment.size() - i - 1;
      //     corrected_len += read.corrected_seq.size();
      //     spdlog::debug("assemble success, size = {}", read.corrected_seq.size());
      //     break;
      //   } else {
      //     spdlog::debug("assemble failed");
      //   }
      // }
    }
    if (flag) {
      precorrect++;
    } else {
      read.assemble_corrected_seq();
      auto [cigar, st, ed] = semi_global_align(get_ans_seq(read.name), read.corrected_seq);
      if (read.corrected_seq.size() > 0) {
        spdlog::debug("assemble success, size = {}, cigar = {}", read.corrected_seq.size(), std::string(cigar));
      } else {
        spdlog::debug("assemble failed");
      }
    }
    spdlog::debug("=====================================================\n\n");
  }

  spdlog::debug("precorrect = {}", precorrect);
  spdlog::debug("correct ratio = {}", (double)corrected_len / total_len);
  spdlog::debug("reduced ratio = {}", (double)saved_fragment / total_fragment);

  // std::set<std::string> wrong_reads;
  // // fs::path wrong = "/mnt/ec/ness/yolkee/thesis/tests/wrong.txt";
  // // auto fin = std::ifstream(wrong);
  // // for (std::string s; std::getline(fin, s);) {
  // //   wrong_reads.emplace(s);
  // // }

  // auto tp = boost::asio::thread_pool(threads);  
  // std::atomic_int cnt = 0;
  // std::atomic_int ok = 0;

  // auto assemble_read = [&](const auto& name, auto& seqs) {
  //   auto read_len = get_seq_size(name);
  //   auto perfect_read_len = get_ans(name);
    
  //   // spdlog::debug("Thread {} start to assemble read {}", tid, name);
  //   // ReadAssembler assembler(read_len);
  //   // auto kmer_size = std::log2(read_len) * 2.5;

  //   ReadAssembler assembler(read_len);
  //   // spdlog::debug("Thread {} start to assemble", tid);

  //   for (const auto& seq : seqs) {
  //     assembler.add_seq(seq);
  //     // spdlog::debug("Thread {} add seq = {}", tid, seq.seq);
  //     // auto encoded_seq = seq.seq;  
  //     // spdlog::debug("encoded seq size = {}", encoded_seq.size());   
  //     // graph.add_seq(encoded_seq, seq.left_bound);
  //     // assembler.add_seq(seq.seq, seq.left_bound, seq.right_bound);
  //   }
  //   auto corrected_read = assembler.assemble();

   
  //   // auto corrected_read = assembler.assemble();

  //   spdlog::info(
  //       "Read {}, raw_read = {}, corrected_read = {}, perfect_read = {}",
  //       name, read_len, corrected_read.size(), perfect_read_len);
  //   spdlog::info("=====================================================\n\n");
  //   // cnt++;
  //   // if (read_len * 1.05 > corrected_read.size() &&
  //   //     corrected_read.size() > read_len * 0.9) {
  //   //   ok += 1;
  //   // }
  //   // if (cnt % 100 == 0) {
  //   //   spdlog::debug("ok = {}, total = {}", ok, cnt);
  //   // }

    
  //   // bio::FastaRecord<false> record;
  //   // record.name = name;
  //   // record.seq = corrected_read;
  // };

  // std::vector<std::future<void>> futures;
  // for (auto& [name, seqs] : fragments) {
  //   boost::asio::post(tp, [&]() {
  //     assemble_read(name, seqs);
  //     seqs.clear();
  //   });
  // }
  // spdlog::debug("ok = {}, total = {}", ok, wrong_reads.size());
  // tp.wait();

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

  // auto vcf_file =
  // "/mnt/ec/mammoth/yolkee/thesis/data/Ecoli/K13/snp_ref/h2.refseq2simseq.SNP.vcf";
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