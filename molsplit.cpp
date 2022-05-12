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


using namespace std;
namespace fs = std::experimental::filesystem;
const char tab='\t';


static void show_usage(string name)
{
    cerr << "Usage: " << name << " <option(s)> MOLECULE FILE \n"
              << "Options:\n"
              << "\t-h,--help\t\tShow this help message\n"
              << "\t-t,--threads INT\tNumber of threads (default 1) \n"
              << "\t-w, --window INT\tWindow size for outliers detection (default 10kb) \n"
              << "\t-c, --contigs FILE\tContig size file name \n"
              << "\t-o, --output FILE\tOutput bed file name \n"
	            << "\t-s, --sampleSize INT\tSample size for outlier detection, if zero or not specified then the sample is the input size"
              << endl;
}

vector<string> split_molecule_file(string molecule_file, int tmp_files_size, int threads){

    ifstream molecules(molecule_file.c_str());
    int line_count = 0, file_count = 0;
    string current_tmp_file = "";
    string ctg="", prev_ctg="";
    string molecule_line;

    ofstream out_molecules_per_thread;

    vector<string> tmpfiles;

    while (getline(molecules,molecule_line)) {

      stringstream  linestream(molecule_line);
      linestream >> ctg;
      if ((line_count >= tmp_files_size && prev_ctg.compare(ctg)!=0 && file_count<threads) || (file_count == 0)) {

        out_molecules_per_thread.close();
        line_count = 0;
        file_count++;

        current_tmp_file = string(fs::path(molecule_file).stem()).append("_");
        current_tmp_file.append(to_string(file_count));
        current_tmp_file.append(".tsv");

        tmpfiles.push_back(current_tmp_file);

        out_molecules_per_thread.open(current_tmp_file.c_str());

      }

      prev_ctg = ctg;
      line_count++;
      out_molecules_per_thread << molecule_line <<endl;

    }

    out_molecules_per_thread.close();

    return tmpfiles;
}

void split_by_outliers(string tmp_file, int window, long n_sample){

    string basename_tmp_file = string(fs::path(tmp_file).stem());
    string stat_file = basename_tmp_file + "_molStats.tsv";
    molecule_stat(tmp_file, stat_file);
    string stat_per_interval_file = basename_tmp_file + "_molStats_perInterval.tsv";
    statistics_per_interval(stat_file, window, stat_per_interval_file);
    string struct_file = basename_tmp_file + "_splits.tsv";
    string distance_scores_file = basename_tmp_file + "_distanceScores.tsv";
    detect_outliers(stat_per_interval_file, struct_file, distance_scores_file, n_sample);

}

int main (int argc, char* argv[])
{
  if (argc < 2) {
      show_usage(argv[0]);
      return 1;
  }

  int threads = 1;
  int window = 10000;
  string contigFile = "";
  string molecule_file = "";
  string bedFile = "";
  int lines_per_thread;
  long n_sample = 0;

  for (int i = 1; i < argc; ++i) {
      string arg = argv[i];
      if ((arg == "-h") || (arg == "--help")) {
          show_usage(argv[0]);
          return 0;
      } else if ((arg == "-t") || (arg == "--threads")) {
            threads = stoi(argv[++i]);
          } else if ((arg == "-w") || (arg == "--window")) {
            window = stoi(argv[++i]);
          } else if ((arg == "-c") || (arg == "--contigs")){
            contigFile = argv[++i];
          } else if ((arg == "-o") || (arg == "--output")){
            bedFile = argv[++i];
          } else if ((arg == "-s") || (arg == "--sampleSize")){
	          n_sample = stoi(argv[++i]);
	        }else {
            molecule_file = argv[i++];
          }
  }

  cout << "File : " <<molecule_file <<" Window :"<<window<<" Threads:"<<threads<<endl;


  if (threads>1){

    ifstream myfile(molecule_file);

    // new lines will be skipped unless we stop it from happening:
    myfile.unsetf(ios_base::skipws);

    // count the newlines with an algorithm specialized for counting:
    unsigned line_count = count(
        istream_iterator<char>(myfile),
        istream_iterator<char>(),
        '\n');

    lines_per_thread = floor(line_count/threads);
    cout << "Lines: " << line_count << "lines_per_thread: "<<lines_per_thread<<endl;
    vector<string> tmpfiles = split_molecule_file(molecule_file,lines_per_thread,threads);

    thread t[tmpfiles.size()];

    for(int i = 0; i<tmpfiles.size();i++){

      t[i] = thread( split_by_outliers, tmpfiles[i], window , n_sample);

    }

    for(int i = 0; i<tmpfiles.size();i++){

      t[i].join();

    }




    //TODO: concatenate all struct_file
    vector<string> tmpStructfiles;
    string base_tmp_file, struct_file;

    for(int i = 0; i<tmpfiles.size();i++){

      base_tmp_file = string(fs::path(tmpfiles[i]).stem());
      struct_file = base_tmp_file + "_splits.tsv";
      tmpStructfiles.push_back(struct_file);

    }

    createBed (tmpStructfiles, contigFile, bedFile);

  }

  return 0;
}
