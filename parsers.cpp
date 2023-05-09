#include <cassert>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

#include "globals.h"
#include "parsers.h"

void parse_contig_file () {
  std::ifstream contig_file (Globals::contigFileName);
  std::string line, contig;
  long int size;
  if (! contig_file.is_open()) {
    std::cerr << "Cannot open contig file.\n";
    exit(EXIT_FAILURE);
  }
  while(getline(contig_file, line)) {
    std::stringstream ss (line);
    ss >> contig >> size;
    Globals::chrids[contig] = Globals::chrs.size();
    Globals::chrs.push_back(contig);
    Globals::chr_sizes.push_back(size);
  }
  Globals::chrs.shrink_to_fit();
  Globals::chr_sizes.shrink_to_fit();
}

void parse_molecule_file (Molecules &molecules) {
  std::ifstream molecule_file (Globals::molecule_file);
  std::string line, contig, barcode;
  unsigned long start, end;
  unsigned int n_reads;
  if (! molecule_file.is_open()) {
    std::cerr << "Cannot open molecule file.\n";
    exit(EXIT_FAILURE);
  }
  while(getline(molecule_file, line)) {
    std::stringstream ss (line);
    ss >> contig >> start >> end >> barcode >> n_reads;
    assert(Globals::chrids.find(contig) != Globals::chrids.end());
    molecules.emplace_back(Globals::chrids[contig], start, end, barcode, n_reads);
  }
}
