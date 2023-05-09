#include <iostream>

#include "parse_parameters.h"


static void show_usage(char *name)
{
    std::cerr << "Usage: " << name << " <option(s)> < BAM_FILE \n"
              << "Options:\n"
              << "  -h, --help               Show this help message\n"
              << "=== Input parameters ===\n"
              << "  -w, --window      INT    Window size for outliers detection (default 10kb) \n"
              << "  -c, --contigs     FILE   Contig size file name \n"
              << "=== Molecule creating step ===\n"
              << "  -d, --distance    INT    Max. distance to cluster reads (default: "    << Globals::max_read_distance          << ")\n"
              << "  -m, --mapq        INT    Min. MAPQ (default: "                         << Globals::min_mapq                   << ")\n"
              << "  -M, --mapqSolid   INT    Min. MAPQ for solid reads (default: "         << Globals::min_mapq_solid             << ")\n"
              << "  -l, --lengthSolid INT    Min. mapping size for solid reads (default: " << Globals::min_len_solid              << ")\n"
              << "  -r, --nReads      INT    Min. #reads per barcode (default: "           << Globals::min_n_reads                << ")\n"
              << "  -c, --ouputMol    FILE   Output file name (default: "                  << Globals::output_molecules_file_name << ")\n"
              << "=== Split ASM step ===\n"
              << "  -t, --threshold   FLOAT  Stringency threshold, higher is less stringent, (default 0.01) \n"
              << "  -o, --output      FILE   Output bed file name \n"
              << "  -s, --sampleSize  INT    Sample size for outlier detection, if zero or not specified then the sample is the input size\n"
              << "  -z, --minSize     INT    Minimum contig size (default 2 x window size)\n"
              << "  -a, --counts      FILE   write raw counts to file\n"
              << "  -A, --scores      FILE   Write scores to file\n"
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
    } else if ((arg == "-r") || (arg == "--nReads")) {
      Globals::min_n_reads = std::stoi(argv[++i]);
    } else if ((arg == "-c") || (arg == "--ouputMol")){
      Globals::output_molecules_file_name = argv[++i];
    } else if ((arg == "-t") || (arg == "--threshold")) {
      Globals::threshold = std::stof(argv[++i]);
    } else if ((arg == "-w") || (arg == "--window")) {
      Globals::window = std::stoi(argv[++i]);
    } else if ((arg == "-c") || (arg == "--contigs")){
      Globals::contigFileName = argv[++i];
    } else if ((arg == "-o") || (arg == "--output")){
      Globals::bedFile = argv[++i];
    } else if ((arg == "-s") || (arg == "--sampleSize")){
      Globals::n_sample = std::stoi(argv[++i]);
    } else if ((arg == "-z") || (arg == "--minSize")){
      Globals::min_ctg_size = std::stoi(argv[++i]);
    } else if ((arg == "-a") || (arg == "--counts")) {
      Globals::countsFileName = argv[++i];
    } else if ((arg == "-A") || (arg == "--scores")) {
      Globals::scoresFileName = argv[++i];
    }
  }
}
