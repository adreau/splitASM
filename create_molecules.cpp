#include <cassert>
#include <algorithm>    // std::sort
#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>     // std::cin, std::cout, std::cerr
#include <sstream>      // std::stringstream


bool starts_with(const std::string &s1, const std::string &s2) {
  return s1.rfind(s2, 0) == 0;
}


struct Interval {
  unsigned long start, end;
  unsigned int mapq;
  Interval (unsigned long s, unsigned long e, unsigned int m): start(s), end(e), mapq(m) {}
};

inline bool operator< (const Interval& lhs, const Interval& rhs) {
  return (lhs.start < rhs.start);
}


struct Molecule {
  unsigned int chrid;
  unsigned long start, end;
  std::string barcode;
  unsigned int n_reads;
  Molecule (unsigned int c, unsigned long s, unsigned long e, const std::string &b, unsigned int n): chrid(c), start(s), end(e), barcode(b), n_reads(n) {}
  friend std::ostream& operator<<(std::ostream& os, const Molecule& molecule);
};

inline bool operator< (const Molecule& lhs, const Molecule& rhs) {
  if (lhs.chrid < rhs.chrid) {
    return true;
  }
  return (lhs.start < rhs.start);
}

std::vector < std::string > chrs;
std::unordered_map < std::string, std::size_t > chrids;
std::unordered_map < std::string, std::vector < std::vector < Interval > > > barcodes;
std::vector < Molecule > molecules;
const unsigned int max_read_distance = 60'000;
const unsigned int min_mapq = 30;
const unsigned int min_mapq_solid = 20;

std::ostream& operator<<(std::ostream& os, const Molecule& molecule) {
  os << chrs[molecule.chrid] << "\t" << molecule.start << "\t" << molecule.end << "\t" << molecule.barcode << "\t" << molecule.n_reads << "\n";
  return os;
}


static void show_usage(std::string name)
{
    std::cerr << "Usage: " << name << " <option(s)> < input_file.sam > output_file.tsv\n"
              << "Options:\n"
              << "  -h, --help             Show this help message\n"
              << std::endl;
}

void add_chr (std::stringstream &split_line) {
  std::string value;
  split_line >> value;
  assert(starts_with(value, "SN"));
  value = value.substr(3);
  chrids[value] = chrs.size();
  chrs.push_back(value);
}

void read_header_line (std::string &line) {
  std::stringstream split_line(line);
  std::string value;
  split_line >> value;
  if (value == "@SQ") {
    add_chr(split_line);
  }
}

void parse_main_line (std::string &line, bool &proper_pair, unsigned int &chrid, unsigned long &start, unsigned long &end, unsigned int &mapq, std::string &barcode) {
  std::stringstream split_line(line);
  std::string value;
  unsigned int flag = -1;
  for (unsigned int i = 0; split_line >> value; ++i) {
    if (i == 1) {
      flag = std::stoi(value);
      proper_pair = (flag & 4 == 0) && (flag & 8 == 0);
    }
    else if (i == 2) {
      assert(chrids.find(value) != chrids.end());
      chrid = chrids[value];
    }
    else if (i == 3) {
      start = std::stoi(value);
    }
    else if (i == 4) {
      mapq = std::stoi(value);
    }
    else if (i == 9) {
      end = value.size();
    }
    else if (i >= 11) {
      if (starts_with(value, "BX:Z:")) {
        barcode = value.substr(5);
      }
    }
  }
  assert(chrid != -1);
  assert(start != -1);
  assert(end != -1);
  assert(mapq != -1);
}

void add_barcode (unsigned int chrid, unsigned long start, unsigned long end, unsigned int mapq, std::string &barcode) {
  if (barcodes.find(barcode) == barcodes.end()) {
    barcodes[barcode] = std::vector < std::vector < Interval > > (chrs.size());
  }
  assert(chrid < chrs.size());
  if ((mapq >= min_mapq) && (! barcode.empty())) {
    barcodes[barcode][chrid].emplace_back(start, end, mapq);
  }
}

void read_main_line (std::string &line) {
  bool proper_pair = false;
  unsigned int chrid = -1;
  unsigned long start = -1, end = -1;
  unsigned int mapq = -1;
  std::string barcode;
  parse_main_line(line, proper_pair, chrid, start, end, mapq, barcode);
  if (proper_pair) {
    add_barcode(chrid, start, end, mapq, barcode);
  }
}

void read_sam() {
  unsigned long int cpt = 0;
  for (std::string line; std::getline(std::cin, line); ++cpt) {
    if (! line.empty()) {
      if (line[0] == '@') {
        read_header_line(line);
      }
      else {
        read_main_line(line);
      }
    }
    if (cpt % 1'000'000 == 0) {
      std::cerr << cpt << " lines read." << "\n";
    }
  }
  std::cerr << cpt << " lines read." << "\n";
}

void sort_barcodes() {
  std::cerr << "Sorting barcodes...\n";
  for (auto &p: barcodes) {
    for (auto &q: p.second) {
      std::sort(q.begin(), q.end());
    }
  }
}

unsigned int find_first_split (std::vector < Interval > &intervals, unsigned int start_id) {
  for (unsigned int prev = start_id, next = start_id + 1; next < intervals.size(); ++prev, ++next) {
    if (intervals[next].start - intervals[prev].end > max_read_distance) {
      return prev;
    }
  }
  return intervals.size() - 1;
}

void find_solid_ends (std::vector < Interval > &intervals, unsigned int start_id, unsigned int end_id, unsigned int &start_solid_id, unsigned int &end_solid_id) {
  start_solid_id = end_solid_id = -1;
  for (unsigned int id = start_id; id <= end_id; ++id) {
    if (intervals[id].mapq >= min_mapq_solid) {
      if (start_solid_id == -1) {
        start_solid_id = id;
      }
      end_solid_id == id;
    }
  }
}

void add_molecule (std::vector < Interval > &intervals, unsigned int start_id, unsigned int end_id, const std::string &barcode, unsigned int chrid) {
  molecules.emplace_back(chrid, intervals[start_id].start, intervals[end_id].end, barcode, end_id + start_id + 1);
}

void _join_to_molecules (std::vector < Interval > &intervals, const std::string &barcode, unsigned int chrid) {
  unsigned int start_id = 0, end_id = 0;
  unsigned int start_solid_id = 0, end_solid_id = 0;
  do {
    end_id = find_first_split(intervals, start_id);
    find_solid_ends(intervals, start_id, end_id, start_solid_id, end_solid_id);
    add_molecule(intervals, start_solid_id, end_solid_id, barcode, chrid);
    start_id = end_id + 1;
  }
  while (start_id < intervals.size());
}

void join_to_molecules() {
  for (auto &p: barcodes) {
    for (unsigned int chrid = 0; chrid < chrs.size(); ++chrid) {
      _join_to_molecules(p.second[chrid], p.first, chrid);
    }
  }
}

void sort_molecules() {
  std::cerr << "Sorting molecules...\n";
  std::sort(molecules.begin(), molecules.end());
}

void print_molecules() {
  for (Molecule &molecule: molecules) {
    std::cout << molecule;
  }
}

void make_molecules() {
  sort_barcodes();
  join_to_molecules();
  sort_molecules();
  print_molecules();
}

int main() {
  read_sam();
  make_molecules();
  return 0;
}
