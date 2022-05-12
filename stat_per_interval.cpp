#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include "stat_per_interval.h"

using namespace std;


const char tab='\t';



void statistics_per_interval(string stat_file, int window, string stat_per_interval_file)
{

  ifstream molecule_stat(stat_file.c_str());


  int interval = window; //resolution

  ofstream out_ctg_stat(stat_per_interval_file.c_str(), ofstream::out);

  string line;
  int count_lines = 0;

  string molecule_line, ctg;
  int pos, mol_cov, mid_mol_cov, s_mol, e_mol;
  double mean_length, read_dens;

  int mol_cov_interval = 0, mid_mol_cov_interval = 0, starting_molecules = 0, ending_molecules = 0;
  double mean_length_interval = 0, read_dens_interval = 0;

  string prev_ctg = "", current_ctg = "";

  getline(molecule_stat,molecule_line);

  while(getline(molecule_stat,molecule_line)){

    stringstream  linestream(molecule_line);
    linestream >> ctg >> pos >> mol_cov >> mid_mol_cov >> mean_length >> read_dens >> s_mol >> e_mol;

    current_ctg = ctg;

    if(current_ctg.compare(prev_ctg)!=0){

      mol_cov_interval = 0;
      mid_mol_cov_interval = 0;
      mean_length_interval = 0;
      read_dens_interval = 0;
      count_lines = 0;
      starting_molecules = 0;
      ending_molecules = 0;

      prev_ctg = ctg;

      getline(molecule_stat,molecule_line);
      stringstream  linestream(molecule_line);
      linestream >> ctg >> pos >> mol_cov >> mid_mol_cov >> mean_length >> read_dens >> s_mol >> e_mol;

    }

    mol_cov_interval+=mol_cov;
    mid_mol_cov_interval+=mid_mol_cov;
    mean_length_interval+=mean_length;
    read_dens_interval+=read_dens;
    starting_molecules += s_mol;
    ending_molecules += e_mol;

    count_lines++;

    if(count_lines == interval){

      out_ctg_stat << ctg << tab << max(pos-interval,0) << tab << pos << tab
          << (double)mol_cov_interval/interval << tab
          << (double)mid_mol_cov_interval/interval << tab
          << (double)mean_length_interval/interval << tab
          << (double)read_dens_interval/interval << tab
          << starting_molecules << tab
          << ending_molecules << endl;

      count_lines = 0;
      mol_cov_interval = 0;
      mid_mol_cov_interval = 0;
      mean_length_interval = 0;
      read_dens_interval = 0;
      starting_molecules = 0;
      ending_molecules = 0;

    }

  }

  out_ctg_stat << ctg << tab << max(pos-interval,0) << tab << pos << tab
          << (double)mol_cov_interval/count_lines << tab
          << (double)mid_mol_cov_interval/count_lines << tab
          << (double)mean_length_interval/count_lines << tab
          << (double)read_dens_interval/count_lines << tab
          << starting_molecules << tab
          << ending_molecules << endl;

  out_ctg_stat.close();


}
