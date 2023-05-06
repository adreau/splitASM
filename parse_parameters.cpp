#include <iostream>

#include "parse_parameters.h"


static void show_usage(char *name)
{
    std::cerr << "Usage: " << name << " <option(s)> < BAM_FILE \n"
              << "Options:\n"
              << "  -h, --help             Show this help message\n"
              << "  -d, --distance    INT  Max. distance to cluster reads (default: " << Globals::max_read_distance << ")\n"
              << "  -m, --mapq        INT  Min. MAPQ (default: " << Globals::min_mapq << ")\n"
              << "  -M, --mapqSolid   INT  Min. MAPQ for solid reads (default: " << Globals::min_mapq_solid << ")\n"
              << "  -l, --lengthSolid INT  Min. mapping size for solid reads (default: " << Globals::min_len_solid << ")\n"
              << "  -r, --nReads      INT  Min. #reads per barcode (default: " << Globals::min_n_reads << ")\n"
              << std::endl;
}


void parse_parameters (int argc, char* argv[])
{
  for (int i = 1; i < argc; ++i) {
      std::string arg = argv[i];
      if ((arg == "-h") || (arg == "--help")) {
          show_usage(argv[0]);
          exit(EXIT_SUCCESS);
      } else if ((arg == "-d") || (arg == "--distance")) {
          Globals::max_read_distance = std::stof(argv[++i]);
      } else if ((arg == "-m") || (arg == "--mapq")) {
          Globals::min_mapq = std::stoi(argv[++i]);
      } else if ((arg == "-M") || (arg == "--mapqSolid")) {
          Globals::min_mapq_solid = std::stoi(argv[++i]);
      } else if ((arg == "-l") || (arg == "--lengthSolid")){
          Globals::min_len_solid = std::stoi(argv[++i]);
      } else if ((arg == "-r") || (arg == "--nReads")){
          Globals::min_n_reads = std::stoi(argv[++i]);
      }
  }
}
