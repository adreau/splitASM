#include "globals.h"

unsigned int Globals::max_read_distance                    { 60'000 };
unsigned int Globals::min_mapq                             {      0 };
unsigned int Globals::min_mapq_solid                       {     30 };
unsigned int Globals::min_len_solid                        {    120 };
unsigned int Globals::min_n_reads                          {      2 };
std::string  Globals::output_molecules_file_name { "outputmols.tsv" };
