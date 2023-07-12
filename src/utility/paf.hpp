#include <string>
#include <iostream>
#include <sstream>


namespace biovoltron {
  
  class PafRecord {
  public:
    std::string qname;
    std::size_t qlen;
    std::size_t qstart;
    std::size_t qend;
    char strand;
    std::string tname;
    std::size_t tlen;
    std::size_t tstart;
    std::size_t tend;
    std::size_t match;
    std::size_t aln_len;
    std::size_t mapq;  
  };

  auto& operator>> (std::istream& os, PafRecord& r) {
    std::string line;
    std::getline(os, line);
    if (line.empty()) {
      return os;
    }
    std::istringstream iss{line};
    iss >> r.qname >> r.qlen >> r.qstart >> r.qend >> r.strand >> r.tname >> r.tlen >> r.tstart >> r.tend >> r.match >> r.aln_len >> r.mapq;
    return os;
  }

  auto& operator<< (std::ostream& os, const PafRecord& r) {
    os << r.qname << "\t" << r.qlen << "\t" << r.qstart << "\t" << r.qend << "\t" << r.strand << "\t" << r.tname << "\t" << r.tlen << "\t" << r.tstart << "\t" << r.tend << "\t" << r.match << "\t" << r.aln_len << "\t" << r.mapq;
    return os;
  }

  auto operator<=> (const PafRecord& lhs, const PafRecord& rhs) {
    if (lhs.qname != rhs.qname) {
      return lhs.qname <=> rhs.qname;
    }
    return std::pair(lhs.qstart, lhs.qend) <=> std::pair(rhs.qstart, rhs.qend);
  }

} // namespace biovoltron
