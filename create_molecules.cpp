#include <cassert>
#include <algorithm>    // std::sort
#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>     // std::cin, std::cout, std::cerr
#include <sstream>      // std::stringstream
#include <fstream>

#include "constants.h"
#include "globals.h"
#include "parse_parameters.h"
#include "sam_parser.h"
#include "create_molecules.h"



unsigned int count_n_names (std::vector < std::string > &names) {
  std::sort(names.begin(), names.end());
  auto last = std::unique(names.begin(), names.end());
  return std::distance(names.begin(), last);
}

unsigned int count_n_reads (std::vector < std::vector < Interval > > &reads) {
  unsigned int n = 0;
  std::vector < std::string > names;
  for (auto &r1: reads) {
    for (auto &r2: r1) {
      names.push_back(r2.name);
    }
  }
  return count_n_names(names);
}

void trim_barcodes (Barcodes &barcodes) {
  std::cerr << "Trimming barcodes...\n";
  auto it = barcodes.begin();
  while (it != barcodes.end()) {
    if (count_n_reads(it->second) < Globals::min_n_reads) it = barcodes.erase(it);
    else                                                ++it;
  }
}

void sort_barcodes (Barcodes &barcodes) {
  std::cerr << "Sorting barcodes...\n";
  for (auto &p: barcodes) {
    p.second.shrink_to_fit();
    for (auto &q: p.second) {
      q.shrink_to_fit();
      std::sort(q.begin(), q.end());
    }
  }
}

// Returns the index of the last read of the current split.
//   A split is made iff the distance between consecutive reads is greater than max_read_distance
unsigned int find_first_split (std::vector < Interval > &intervals, unsigned int start_id) {
  for (unsigned int prev = start_id, next = start_id + 1; next < intervals.size(); ++prev, ++next) {
    if (intervals[prev].get_distance(intervals[next]) > Globals::max_read_distance) {
      return next;
    }
  }
  return intervals.size();
}

// Returns true iff at least one reads has the given mapq
bool check_min_mapq (std::vector < Interval > &intervals, unsigned int start_id, unsigned int end_id) {
  assert(start_id <= end_id);
  assert(end_id <= intervals.size());
  unsigned int mapq = 0;
  for (unsigned int id = start_id; id < end_id; ++id) {
    mapq = std::max(mapq, intervals[id].mapq);
  }
  return (mapq >= Globals::min_mapq_solid);
}

// Returns true iff at least one solid read is found
bool find_solid_ends (std::vector < Interval > &intervals, unsigned int start_id, unsigned int end_id, unsigned int &start_solid_id, unsigned int &end_solid_id) {
  assert(start_id <= end_id);
  assert(end_id <= intervals.size());
  constexpr unsigned int no_id = -1;
  start_solid_id = end_solid_id = no_id;
  for (unsigned int id = start_id; id < end_id; ++id) {
    if (intervals[id].is_solid()) {
      if (start_solid_id == no_id) {
        start_solid_id = id;
      }
      end_solid_id = id + 1;
    }
  }
  return (start_solid_id != no_id);
}

void add_molecule (Molecules &molecules, std::vector < Interval > &intervals, unsigned int start_id, unsigned int end_id, const std::string &barcode, unsigned int chrid) {
  assert(start_id <= end_id);
  assert(end_id <= intervals.size());
  std::vector < std::string > names;
  names.reserve(end_id - start_id);
  for (unsigned int id = start_id; id < end_id; ++id) {
    names.push_back(intervals[id].name);
  }
  molecules.emplace_back(chrid, intervals[start_id].start, intervals[end_id - 1].end, barcode, count_n_names(names));
}

void _join_to_molecules (Molecules &molecules, std::vector < Interval > &intervals, const std::string &barcode, unsigned int chrid) {
  if (intervals.empty()) return;
  unsigned int start_id = 0, end_id = 0;
  unsigned int start_solid_id, end_solid_id;
  do {
    end_id = find_first_split(intervals, start_id);
    assert(end_id <= intervals.size());
    if (check_min_mapq(intervals, start_id, end_id) && find_solid_ends(intervals, start_id, end_id, start_solid_id, end_solid_id)) {
      add_molecule(molecules, intervals, start_solid_id, end_solid_id, barcode, chrid);
    }
    start_id = end_id;
  }
  while (start_id < intervals.size());
}

void join_to_molecules (Barcodes &barcodes, Molecules &molecules) {
  for (auto &p: barcodes) {
    for (unsigned int chrid = 0; chrid < Globals::chrs.size(); ++chrid) {
      _join_to_molecules(molecules, p.second[chrid], p.first, chrid);
    }
  }
}

void sort_molecules (Molecules &molecules) {
  std::cerr << "Sorting " << molecules.size() << " molecules...\n";
  molecules.shrink_to_fit();
  std::sort(molecules.begin(), molecules.end());
}

void print_molecules (Molecules &molecules) {
  std::ofstream output_file(Globals::output_molecules_file_name);
  if (! output_file) {
    std::cerr << "Cannot write output molecule file to " << Globals::output_molecules_file_name << ".\nExting." << "\n";
    exit(EXIT_FAILURE);
  }
  for (Molecule &molecule: molecules) {
    output_file << molecule;
  }
  output_file.close();
}

void make_molecules(Barcodes &barcodes, Molecules &molecules) {
  trim_barcodes(barcodes);
  sort_barcodes(barcodes);
  join_to_molecules(barcodes, molecules);
  sort_molecules(molecules);
  print_molecules(molecules);
}
