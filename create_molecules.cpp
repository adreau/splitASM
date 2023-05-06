#include <cassert>
#include <algorithm>    // std::sort
#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>     // std::cin, std::cout, std::cerr
#include <sstream>      // std::stringstream

const unsigned int max_read_distance = 60'000;
const unsigned int min_mapq = 0;
const unsigned int min_mapq_solid =  30;
const unsigned int min_len_solid  = 120;
const unsigned int min_n_reads    =   2;

const std::string SEQUENCE_LINE {   "@SQ" };
const std::string SEQUENCE_NAME {   "SN:" };
const std::string BARCODE_FLAG  { "BX:Z:" };


bool starts_with(const std::string &s1, const std::string &s2) {
  return s1.rfind(s2, 0) == 0;
}


struct Interval {
  std::string name;
  unsigned long start, end;
  unsigned int mapq;
  Interval (std::string &n, unsigned long s, unsigned long e, unsigned int m): name(n), start(s), end(e), mapq(m) {}
  bool is_solid () {
    return ((mapq >= min_mapq_solid) || ((end - start + 1) > min_len_solid));
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
  if (lhs.chrid > rhs.chrid) {
    return false;
  }
  return (lhs.start < rhs.start);
}

std::vector < std::string > chrs;
std::unordered_map < std::string, std::size_t > chrids;
std::unordered_map < std::string, std::vector < std::vector < Interval > > > barcodes;
std::vector < Molecule > molecules;

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
  assert(starts_with(value, SEQUENCE_NAME));
  value = value.substr(SEQUENCE_NAME.size());
  chrids[value] = chrs.size();
  chrs.push_back(value);
}

void read_header_line (std::string &line) {
  std::stringstream split_line(line);
  std::string value;
  split_line >> value;
  if (value == SEQUENCE_LINE) {
    add_chr(split_line);
  }
}

bool is_mapped (unsigned int flag) {
  return ((flag & 4) == 0);
}

bool is_duplicate (unsigned int flag) {
  return ((flag & 1024) != 0);
}

// Return the match length
unsigned int parse_cigar (std::string &cigar) {
  unsigned int match_len = 0;
  unsigned int number = 0;
  for (char c: cigar) {
    switch (c) {
      case '0':
      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
      case '8':
      case '9':
        number = number * 10 + (c - '0');
        break;
      case 'M':
      case 'X':
      case '=':
      case 'D':
        match_len += number;
        number = 0;
        break;
      case 'I':
      case 'S':
      case 'P':
      case 'H':
        number = 0;
        break;
      default:
        assert(false);
    }
  }
  return match_len;
}

// Returns true iff:
//  - read is mapped
//  - flagged in proper pair
//  - mapq is higher than threshold
//  - barcode in BX tag is set
bool parse_main_line (std::string &line, std::string &name, unsigned int &chrid, unsigned long &start, unsigned long &end, unsigned int &mapq, std::string &barcode) {
  std::stringstream split_line(line);
  std::string value;
  unsigned int flag = -1;
  for (unsigned int i = 0; split_line >> value; ++i) {
    if (i == 0) {
      name = value;
    }
    else if (i == 1) {
      flag = std::stoi(value);
      if ((! is_mapped(flag)) || (is_duplicate(flag))) {
        return false;
      }
    }
    else if (i == 2) {
      if (value == "*") {
        return false;
      }
      assert(chrids.find(value) != chrids.end());
      chrid = chrids[value];
    }
    else if (i == 3) {
      start = std::stoi(value) - 1;
    }
    else if (i == 4) {
      mapq = std::stoi(value);
      if (mapq < min_mapq) {
        return false;
      }
    }
    else if (i == 5) {
      end = start + parse_cigar(value);
    }
    else if (i >= 11) {
      if (starts_with(value, BARCODE_FLAG)) {
        barcode = value.substr(BARCODE_FLAG.size());
      }
    }
  }
  if (barcode.empty()) {
    return false;
  }
  assert(chrid != -1);
  assert(start != -1);
  assert(end != -1);
  assert(mapq != -1);
  return true;
}

void add_barcode (std::string &name, unsigned int chrid, unsigned long start, unsigned long end, unsigned int mapq, std::string &barcode) {
  if (barcodes.find(barcode) == barcodes.end()) {
    barcodes[barcode] = std::vector < std::vector < Interval > > (chrs.size());
  }
  assert(chrid < chrs.size());
  barcodes[barcode][chrid].emplace_back(name, start, end, mapq);
}

// Returns true iff read passed the filters
bool read_main_line (std::string &line) {
  std::string name;
  unsigned int chrid = -1;
  unsigned long start = -1, end = -1;
  unsigned int mapq = -1;
  std::string barcode;
  if (parse_main_line(line, name, chrid, start, end, mapq, barcode)) {
    add_barcode(name, chrid, start, end, mapq, barcode);
    return true;
  }
  return false;
}

void read_header(std::string &line) {
  for (line; std::getline(std::cin, line);) {
    if (line[0] == '@') {
      read_header_line(line);
    }
    else {
      return;
    }
  }
}

void read_main(std::string &line) {
  unsigned long int n_reads_kept = 0;
  unsigned long int cpt = 1;
  if (read_main_line(line)) {
    ++n_reads_kept;
  }
  for (; std::getline(std::cin, line); ++cpt) {
    if (read_main_line(line)) {
      ++n_reads_kept;
    }
    if (cpt % 10'000'000 == 0) {
      std::cerr << cpt << " reads read, " << n_reads_kept << " kept in " << barcodes.size() << " barcodes." << "\n";
    }
  }
  std::cerr << cpt << " reads read, " << n_reads_kept << " kept in " << barcodes.size() << " barcodes." << "\n";
}

void read_sam() {
  std::string line;
  read_header(line);
  std::cerr << "Seen " << chrs.size() << " references.\n";
  read_main(line);
}

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

void trim_barcodes() {
  std::cerr << "Trimming barcodes...\n";
  auto it = barcodes.begin();
  while (it != barcodes.end()) {
    if (count_n_reads(it->second) < min_n_reads) it = barcodes.erase(it);
    else                                         ++it;
  }
}

void sort_barcodes() {
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
    if (intervals[prev].get_distance(intervals[next]) > max_read_distance) {
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
  return (mapq >= min_mapq_solid);
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

void add_molecule (std::vector < Interval > &intervals, unsigned int start_id, unsigned int end_id, const std::string &barcode, unsigned int chrid) {
  assert(start_id <= end_id);
  assert(end_id <= intervals.size());
  std::vector < std::string > names;
  names.reserve(end_id - start_id);
  for (unsigned int id = start_id; id < end_id; ++id) {
    names.push_back(intervals[id].name);
  }
  molecules.emplace_back(chrid, intervals[start_id].start, intervals[end_id - 1].end, barcode, count_n_names(names));
}

void _join_to_molecules (std::vector < Interval > &intervals, const std::string &barcode, unsigned int chrid) {
  if (intervals.empty()) return;
  unsigned int start_id = 0, end_id = 0;
  unsigned int start_solid_id, end_solid_id;
  do {
    end_id = find_first_split(intervals, start_id);
    assert(end_id <= intervals.size());
    if (check_min_mapq(intervals, start_id, end_id) && find_solid_ends(intervals, start_id, end_id, start_solid_id, end_solid_id)) {
      add_molecule(intervals, start_solid_id, end_solid_id, barcode, chrid);
    }
    start_id = end_id;
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
  molecules.shrink_to_fit();
  std::sort(molecules.begin(), molecules.end());
}

void print_molecules() {
  for (Molecule &molecule: molecules) {
    std::cout << molecule;
  }
}

void make_molecules() {
  trim_barcodes();
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
