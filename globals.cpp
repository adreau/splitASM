#include "globals.h"

unsigned int Globals::max_read_distance          {           60'000 };
unsigned int Globals::min_mapq                   {                0 };
unsigned int Globals::min_mapq_solid             {               30 };
unsigned int Globals::min_len_solid              {              120 };
unsigned int Globals::min_n_reads                {                2 };
std::string  Globals::output_molecules_file_name { "outputmols.tsv" };
float        Globals::threshold                  {             0.01 };
int          Globals::window                     {            10000 };
std::string  Globals::countsFileName;
std::string  Globals::scoresFileName;
std::string  Globals::contigFileName;
std::string  Globals::molecule_file;
std::string  Globals::bedFile;
long         Globals::n_sample                   {                0 };
int          Globals::min_ctg_size               {            20000 };


std::vector < std::string >                     Globals::chrs;
std::unordered_map < std::string, std::size_t > Globals::chrids;
std::vector < size_t >                          Globals::chr_sizes;


bool starts_with(const std::string &s1, const std::string &s2) {
  return s1.rfind(s2, 0) == 0;
}
