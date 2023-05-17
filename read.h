#ifndef INTERVAL_H
#define INTERVAL_H

#include <string>


struct Interval {
  std::string name;
  unsigned long start, end;
  unsigned int mapq;
  Interval (std::string &n, unsigned long s, unsigned long e, unsigned int m): name(n), start(s), end(e), mapq(m) {}
  bool is_solid () {
    return ((mapq >= Globals::min_mapq_solid) || ((end - start + 1) > Globals::min_len_solid));
  }
  unsigned long int get_distance (Interval &i) {
    if (start   > i.end) return start - i.end;
    if (i.start > end)   return i.start - end;
    return 0;
  }
};

inline bool operator< (const Interval& lhs, const Interval& rhs) {
  return (lhs.start < rhs.start);
}

using Barcodes = std::unordered_map < std::string, std::vector < std::vector < Interval > > >;

#endif
