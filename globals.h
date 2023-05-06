#ifndef GLOBALS_H
#define GLOBALS_H

#include <string>

namespace Globals {
  extern unsigned int max_read_distance;
  extern unsigned int min_mapq;
  extern unsigned int min_mapq_solid;
  extern unsigned int min_len_solid;
  extern unsigned int min_n_reads;
  extern std::string  output_molecules_file_name;
}

#endif
