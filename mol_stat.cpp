#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

#include "ring.h"
#include "mol_stat.h"


const char tab = '\t';
const int buffer_size = 300000; //maximum number of overlapping molecules

const int min_extremity_length = 1500; //the core of a molecule is at 1.5kb from extremities

void read_molecule_line_old(std::string molecule_line, std::string &ctg, int &beg, int &end,
                        std::string &bc, int &reads, int &id){

    std::stringstream  linestream(molecule_line);
    linestream >> ctg >> beg >> end >> bc >> reads >> id;

}

void read_molecule_line(std::string molecule_line, std::string &ctg, int &beg, int &end,
                        std::string &bc, double &reads){

    std::stringstream  linestream(molecule_line);
    linestream >> ctg >> beg >> end >> bc >> reads;

}

void init_ring_buffer(ring<int>& buffer){

  for (size_t i = 0; i < buffer_size; ++i) {
        buffer.push(0);
  }

}

void init_ring_buffer(ring<double>& buffer){

  for (size_t i = 0; i < buffer_size; ++i) {
        buffer.push(0);
  }

}

void molecule_stat(std::string tmp_file, std::string stat_file)
{

  std::ifstream molecule_list(tmp_file.c_str());

  std::ofstream out_mol_stat(stat_file.c_str(), std::ofstream::out);

  ring<int> molecule_coverage(buffer_size);
  init_ring_buffer(molecule_coverage);

  ring<int> middle_mol_coverage(buffer_size);
  init_ring_buffer(middle_mol_coverage);

  ring<int> molecule_length(buffer_size);
  init_ring_buffer(molecule_length);

  ring<double> molecule_read_density(buffer_size);
  init_ring_buffer(molecule_read_density);

  ring<int> starting_molecules(buffer_size);
  init_ring_buffer(starting_molecules);

  ring<int> ending_molecules(buffer_size);
  init_ring_buffer(ending_molecules);


  /*out_mol_stat << "Contig" << tab << "Position" << tab << "Molecule coverage" << tab <<
      "Middle molecules" << tab << "Mean molecule length" << tab <<
      "Mean Read Density" << endl;

  */

  std::string prev_ctg = "", current_ctg = "";

  int begin_pos_stack,pos_end_stack;
  int end_ctg;

  std::string ctg, bc;
  int beg, end, reads, id;

  std::string line;
  int mol_size = 0;
  double read_density = 0;

  while(getline(molecule_list,line)){

    //TODO: id is useless, remove it
    read_molecule_line_old(line, ctg, beg, end, bc, reads, id);

    current_ctg = ctg;


    if(current_ctg.compare(prev_ctg) != 0){ //we start a new contig

      if(prev_ctg.length()>0){ //if it's not the first contig we need to finish the previous contig

        pos_end_stack = end_ctg;

        for(int i = begin_pos_stack; i <= pos_end_stack;i++){


          out_mol_stat << prev_ctg << tab << i << tab << molecule_coverage.top() << tab <<
                middle_mol_coverage.top() << tab <<
                (double)molecule_length.top()/molecule_coverage.top() << tab <<
                (double)molecule_read_density.top()/molecule_coverage.top() << tab <<
                starting_molecules.top() << tab << ending_molecules.top() << std::endl;

            molecule_coverage.pop();
            middle_mol_coverage.pop();
            molecule_length.pop();
            molecule_read_density.pop();
            starting_molecules.pop();
            ending_molecules.pop();


        }
      }

      prev_ctg = current_ctg;
      end_ctg = end;

      for(int i = 0; i < beg; i++){ //from the beging of the new ctg to the begining of the 1st mol print 0s

        out_mol_stat << ctg << tab << i << tab << molecule_coverage.top() << tab <<
        middle_mol_coverage.top() << tab << molecule_length.top() << tab <<
        molecule_read_density.top() << tab << 
        starting_molecules.top() << tab << ending_molecules.top() << std::endl;

        molecule_coverage.push(0);
        middle_mol_coverage.push(0);
        molecule_length.push(0);
        molecule_read_density.push(0);
        starting_molecules.push(0);
        ending_molecules.push(0);

      }


      begin_pos_stack = beg;
      mol_size = end-beg+1;
      read_density = (double)reads/mol_size;

      for(int i = 0; i < mol_size; i++){

        molecule_coverage[i]++;
        molecule_length[i] += mol_size;
        molecule_read_density[i] += read_density;

      }

      starting_molecules[0]++;
      ending_molecules[mol_size-1]++;

      for(int i = min_extremity_length; i < mol_size-min_extremity_length; i++){

        middle_mol_coverage[i]++;

      }

    }else{//then we are in the middle of a contig


      end_ctg = std::max(end_ctg,end);

      //print the positions before the beging of the molecule (no next molecule will cover them)
      for(int i = begin_pos_stack; i < beg; i++){

        out_mol_stat << ctg << tab << i << tab << molecule_coverage.top() << tab <<
        middle_mol_coverage.top() << tab <<
        (double)molecule_length.top()/molecule_coverage.top() << tab <<
        (double)molecule_read_density.top()/molecule_coverage.top() << tab <<
        starting_molecules.top() << tab << ending_molecules.top() << std::endl;

        molecule_coverage.push(0);
        middle_mol_coverage.push(0);
        molecule_length.push(0);
        molecule_read_density.push(0);
        starting_molecules.push(0);
        ending_molecules.push(0);

      }

      begin_pos_stack = beg;

      //add to the positions covered by the molecule
      mol_size = end-begin_pos_stack+1;
      read_density = (double)reads/mol_size;

      for(int i = 0; i < mol_size;i++){

        molecule_coverage[i]++;
        molecule_length[i] += mol_size;
        molecule_read_density[i] += read_density;

      }

      starting_molecules[0]++;
      ending_molecules[mol_size-1]++;


      for(int i = 500; i < mol_size-500; i++){

        middle_mol_coverage[i]++;

      }



    }

  }


  pos_end_stack = end_ctg;

  for(int i = begin_pos_stack; i <= pos_end_stack;i++){


    out_mol_stat << prev_ctg << tab << i << tab << molecule_coverage.top() << tab <<
          middle_mol_coverage.top() << tab <<
          (double)molecule_length.top()/molecule_coverage.top() << tab <<
          (double)molecule_read_density.top()/molecule_coverage.top() << tab <<
          starting_molecules.top() << tab << ending_molecules.top() << std::endl;

      molecule_coverage.pop();
      middle_mol_coverage.pop();
      molecule_length.pop();
      molecule_read_density.pop();
      starting_molecules.pop();
      ending_molecules.pop();

  }

  out_mol_stat.close();

}

void molecule_stat2(
    std::string                            &input_file_name,
    std::vector < std::vector < double > > &molecule_coverage,
    std::vector < std::vector < double > > &middle_mol_coverage,
    std::vector < std::vector < double > > &molecule_length,
    std::vector < std::vector < double > > &molecule_read_density,
    std::vector < std::vector < double > > &starting_molecules,
    std::vector < std::vector < double > > &ending_molecules,
    std::vector < std::string >            &chr_names,
    std::vector < long int >               &chr_sizes,
    int                                     window,
    std::string                            &statsFileName1)
{

  size_t nchrs = chr_names.size();
  std::unordered_map < std::string, size_t > chr_map; // map chr name to ids
  molecule_coverage.resize(nchrs);
  middle_mol_coverage.resize(nchrs);
  molecule_length.resize(nchrs);
  molecule_read_density.resize(nchrs);
  starting_molecules.resize(nchrs);
  ending_molecules.resize(nchrs);

  for (size_t chrid = 0; chrid < nchrs; ++chrid) {
    chr_map[chr_names[chrid]] = chrid;
    long int size = (chr_sizes[chrid] - 1) / window + 1;
    molecule_coverage[chrid]     = std::vector < double > (size, 0);
    middle_mol_coverage[chrid]   = std::vector < double > (size, 0);
    molecule_length[chrid]       = std::vector < double > (size, 0);
    molecule_read_density[chrid] = std::vector < double > (size, 0);
    starting_molecules[chrid]    = std::vector < double > (size, 0);
    ending_molecules[chrid]      = std::vector < double > (size, 0);
  }

  std::string line;
  std::string ctg, bc;
  int beg, end;
  double reads;
  unsigned long int cpt = 0;

  std::ifstream input_file(input_file_name.c_str());
  while(getline(input_file, line)) {

    read_molecule_line(line, ctg, beg, end, bc, reads);

    // Easier to work with 0-based positions
    --beg;
    --end;

    int window_start = beg / window;
    int window_end   = end / window;
    int molecule_size = end - beg + 1;
    
    if (chr_map.find(ctg) == chr_map.end()) {
      std::cerr << "Contig " << ctg << " is not present in the contig file.\n";
      exit(EXIT_FAILURE);
    }
    size_t chrid = chr_map[ctg];
    for (int windowid = window_start; windowid <= window_end; ++windowid) {
      if (windowid >= molecule_coverage[chrid].size()) {
        std::cerr << "Molecule size " << ctg << ":" << beg << "-" << end << " out of range.\n";
        exit(EXIT_FAILURE);
      }
      int beg_window = std::max < int > (windowid * window, beg);
      int end_window = std::min < int > (beg_window + window - 1, end);
      int size = end_window - beg_window + 1;
      molecule_coverage[chrid][windowid] += size;
      molecule_length[chrid][windowid] += molecule_size * size;
      molecule_read_density[chrid][windowid] += reads / molecule_size * size;
    }
    ++starting_molecules[chrid][window_start];
    ++ending_molecules[chrid][window_end];

    // Considering inner part of molecule
    beg += min_extremity_length;
    end -= min_extremity_length;

    if (beg <= end) {
      window_start = beg / window;
      window_end   = end / window;
      
      for (int windowid = window_start; windowid <= window_end; ++windowid) {
        int beg_window = std::max < int > (windowid * window, beg);
        int end_window = std::min < int > (beg_window + window - 1, end);
        int size = end_window - beg_window + 1;
        middle_mol_coverage[chrid][windowid] += size;
      }
    }
    ++cpt;
    if (cpt % 10000000 == 0) std::cout << cpt << " lines read.\r" << std::flush;
  }
  std::cout << cpt << " lines read, done.\n";

  std::ofstream statsFile1;
  if (! statsFileName1.empty()) statsFile1.open(statsFileName1, std::ofstream::out);

  for (size_t chrid = 0; chrid < nchrs; ++chrid) {
    std::string &chr = chr_names[chrid];
    long int size = chr_sizes[chrid] / window;
    for (size_t windowid = 0; windowid <= size; ++windowid) {
      int beg_window = windowid * window;
      int end_window = std::min < int > (beg_window + window - 1, chr_sizes[chrid]);
      double size_dbl = static_cast < double > (end_window - beg_window + 1);
      // Careful in the order of normalization
      molecule_coverage[chrid][windowid]   /= size_dbl;
      middle_mol_coverage[chrid][windowid] /= size_dbl;
      if (molecule_coverage[chrid][windowid] == 0) {
        molecule_length[chrid][windowid]       = 0;
        molecule_read_density[chrid][windowid] = 0;
      }
      else {
        molecule_length[chrid][windowid]       /= molecule_coverage[chrid][windowid] * window;
        molecule_read_density[chrid][windowid] /= molecule_coverage[chrid][windowid] * window;
      }

      // Switch back to 1-based positions
      if (! statsFileName1.empty()) statsFile1 << chr << tab << windowid * window + 1 << tab << (windowid + 1) * window << tab << molecule_coverage[chrid][windowid] << tab << middle_mol_coverage[chrid][windowid] << tab << molecule_length[chrid][windowid] << tab << molecule_read_density[chrid][windowid] << tab << starting_molecules[chrid][windowid] << tab << ending_molecules[chrid][windowid] << "\n";
    }
  }
  if (! statsFileName1.empty()) statsFile1.close();
}
