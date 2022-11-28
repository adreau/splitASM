#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <experimental/filesystem>
#include <vector>
#include <thread>

#include "mol_stat.h"
#include "stat_per_interval.h"
#include "outliers_det.h"
#include "create_bed_file.h"


namespace fs = std::experimental::filesystem;
const char tab='\t';


static void show_usage(std::string name)
{
    std::cerr << "Usage: " << name << " <option(s)> MOLECULE FILE \n"
              << "Options:\n"
              << "  -h, --help             Show this help message\n"
              << "  -t, --threshold FLOAT  Stringency threshold (default 0.01, higher is less stringent) \n"
              << "  -w, --window INT       Window size for outliers detection (default 10kb) \n"
              << "  -c, --contigs FILE     Contig size file name \n"
              << "  -o, --output FILE      Output bed file name \n"
              << "  -s, --sampleSize INT   Sample size for outlier detection, if zero or not specified then the sample is the input size\n"
              << "  -m, --minSize INT      Minimum contig size (default 2 x window size)\n"
              << "  -a FILE                write raw counts to file\n"
              << "  -A FILE                write scores to file\n"
              << std::endl;
}

std::vector<std::string> split_molecule_file(std::string molecule_file, int tmp_files_size, int threads){

    std::ifstream molecules(molecule_file.c_str());
    int line_count = 0, file_count = 0;
    std::string current_tmp_file = "";
    std::string ctg="", prev_ctg="";
    std::string molecule_line;

    std::ofstream out_molecules_per_thread;

    std::vector<std::string> tmpfiles;

    while (getline(molecules,molecule_line)) {

      std::stringstream  linestream(molecule_line);
      linestream >> ctg;
      if ((line_count >= tmp_files_size && prev_ctg.compare(ctg)!=0 && file_count<threads) || (file_count == 0)) {

        out_molecules_per_thread.close();
        line_count = 0;
        file_count++;

        current_tmp_file = std::string(fs::path(molecule_file).stem()).append("_");
        current_tmp_file.append(std::to_string(file_count));
        current_tmp_file.append(".tsv");

        tmpfiles.push_back(current_tmp_file);

        out_molecules_per_thread.open(current_tmp_file.c_str());

      }

      prev_ctg = ctg;
      line_count++;
      out_molecules_per_thread << molecule_line << std::endl;

    }

    out_molecules_per_thread.close();

    return tmpfiles;
}

void split_by_outliers(std::string tmp_file, int window, long n_sample){

    std::string basename_tmp_file = std::string(fs::path(tmp_file).stem());
    std::string stat_file = basename_tmp_file + "_molStats.tsv";
    molecule_stat(tmp_file, stat_file);
    std::string stat_per_interval_file = basename_tmp_file + "_molStats_perInterval.tsv";
    statistics_per_interval(stat_file, window, stat_per_interval_file);
    std::string struct_file = basename_tmp_file + "_splits.tsv";
    std::string distance_scores_file = basename_tmp_file + "_distanceScores.tsv";
    detect_outliers(stat_per_interval_file, struct_file, distance_scores_file, n_sample);

}

int main (int argc, char* argv[])
{
  if (argc < 2) {
      show_usage(argv[0]);
      return 1;
  }

  float threshold = 0.01;
  int window = 10000;
  std::string statsFileName1;
  std::string statsFileName2;
  std::string contigFileName;
  std::string molecule_file;
  std::string bedFile;
  int lines_per_thread;
  long n_sample = 0;
  int min_ctg_size = 2 * window;

  for (int i = 1; i < argc; ++i) {
      std::string arg = argv[i];
      if ((arg == "-h") || (arg == "--help")) {
          show_usage(argv[0]);
          return 0;
      } else if ((arg == "-t") || (arg == "--threshold")) {
          threshold = std::stof(argv[++i]);
      } else if ((arg == "-w") || (arg == "--window")) {
          window = std::stoi(argv[++i]);
      } else if ((arg == "-c") || (arg == "--contigs")){
          contigFileName = argv[++i];
      } else if ((arg == "-o") || (arg == "--output")){
          bedFile = argv[++i];
      } else if ((arg == "-s") || (arg == "--sampleSize")){
          n_sample = std::stoi(argv[++i]);
      } else if ((arg == "-m") || (arg == "--minSize")){
          min_ctg_size = std::stoi(argv[++i]);
      } else if (arg == "-a") {
          statsFileName1 = argv[++i];
      } else if (arg == "-A") {
          statsFileName2 = argv[++i];
      } else {
          molecule_file = argv[i++];
      }
  }

  std::vector < std::string > chr_names;
  std::vector < long int > chr_sizes;

  std::ifstream contigFile(contigFileName);
  std::string line, contig;
  long int size;
  if (! contigFile.is_open()) {
    std::cerr << "Cannot open contig file.\n";
    exit(EXIT_FAILURE);
  }
  while(getline(contigFile, line)) {
    std::stringstream streamedLine (line);
    streamedLine >> contig >> size;
    chr_names.push_back(contig);
    chr_sizes.push_back(size);
  }
  chr_names.shrink_to_fit();
  chr_sizes.shrink_to_fit();

  std::string output_file_name ("tmp.out");
  std::vector < std::vector < double > > molecule_coverages;
  std::vector < std::vector < double > > molecule_lengths;
  std::vector < std::vector < double > > molecule_read_densities;
  std::vector < std::vector < double > > starting_molecules;
  std::vector < std::vector < double > > ending_molecules;
  molecule_stat2(
    molecule_file,
    molecule_coverages,
    molecule_lengths,
    molecule_read_densities,
    starting_molecules,
    ending_molecules,
    chr_names,
    chr_sizes,
    window,
    statsFileName1);
  detect_outliers2(
    molecule_coverages,
    molecule_lengths,
    molecule_read_densities,
    starting_molecules,
    ending_molecules,
    chr_names,
    chr_sizes,
    window,
    min_ctg_size,
    n_sample,
    bedFile,
    statsFileName2,
    threshold);


  exit(EXIT_SUCCESS);


}
