#include <string>
#include <fstream>

#include <spoa/spoa.hpp>
#include <spdlog/spdlog.h>

#include <biovoltron/file_io/all.hpp>
#include <varic/utility/file_io/read_record.hpp>

namespace bio = biovoltron;

auto get_msa(const std::vector<bio::FastaRecord<false>>& records) {

  auto alignment_engine = spoa::AlignmentEngine::Create(
      spoa::AlignmentType::kNW,  // Needleman-Wunsch(global alignment)
      5,                         // match (default parameter form SPOA)
      -4,                        // mismatch
      -8,                        // gap
      -6                         // gap extension
  );

  auto sz = records.size();
  auto graph = spoa::Graph();
  for (const auto& record : records) {
    auto alignment = alignment_engine->Align(record.seq, graph);
    graph.AddAlignment(alignment, record.seq);
  }
  auto msa = graph.GenerateMultipleSequenceAlignment(true);

  std::vector<bio::FastaRecord<false>> msa_records;
  for (auto i = 0u; i < sz; i++) {
    msa_records.emplace_back(records[i].name, msa[i]);
  }
  msa_records.emplace_back("consensus", msa.back());
  
  return msa_records;
}

int main(int argc, char* argv[]) {
  auto records = read_records<bio::FastaRecord<false>>(argv[1]);
  auto msa_records = get_msa(records);

  std::ofstream out(argv[2]);
  for (const auto& record : msa_records) {
    out << record << '\n';
  }
}