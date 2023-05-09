#ifndef GLOBALS_H
#define GLOBALS_H

#include <unordered_map>
#include <vector>
#include <string>


namespace Globals {
  extern unsigned int max_read_distance;
  extern unsigned int min_mapq;
  extern unsigned int min_mapq_solid;
  extern unsigned int min_len_solid;
  extern unsigned int min_n_reads;
  extern std::string  output_molecules_file_name;
  extern float        threshold;
  extern int          window;
  extern std::string  countsFileName;
  extern std::string  scoresFileName;
  extern std::string  contigFileName;
  extern std::string  molecule_file;
  extern std::string  bedFile;
  extern long         n_sample;
  extern int          min_ctg_size;

  extern std::vector < std::string >                     chrs;
  extern std::unordered_map < std::string, std::size_t > chrids;
  extern std::vector < std::size_t >                     chr_sizes;
}


struct Molecule_stat {
  double coverage;
  double length;
  double read_density;
  double start;
  double end;
};

using Molecule_stats = std::vector < std::vector < Molecule_stat > >;



bool starts_with (const std::string &s1, const std::string &s2);


#endif
