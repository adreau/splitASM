#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

#include "ring.h"
#include "mol_stat.h"

using namespace std;

const char tab = '\t';
const int buffer_size = 300000; //maximum number of overlapping molecules

const int min_extremity_length = 1500; //the core of a molecule is at 1.5kb from extremities

void read_molecule_line(string molecule_line, string &ctg, int &beg, int &end,
                        string &bc, int &reads, int &id){

    stringstream  linestream(molecule_line);
    linestream >> ctg >> beg >> end >> bc >> reads >> id;

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

void molecule_stat(string tmp_file, string stat_file)
{

  ifstream molecule_list(tmp_file.c_str());

  ofstream out_mol_stat(stat_file.c_str(), ofstream::out);

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

  string prev_ctg = "", current_ctg = "";

  int begin_pos_stack,pos_end_stack;
  int end_ctg;

  string ctg, bc;
  int beg, end, reads, id;

  string line;
  int mol_size = 0;
  double read_density = 0;

  while(getline(molecule_list,line)){

    read_molecule_line(line, ctg, beg, end, bc, reads, id);

    current_ctg = ctg;


    if(current_ctg.compare(prev_ctg) != 0){ //we start a new contig

      if(prev_ctg.length()>0){ //if it's not the first contig we need to finish the previous contig

        pos_end_stack = end_ctg;

        for(int i = begin_pos_stack; i <= pos_end_stack;i++){


          out_mol_stat << prev_ctg << tab << i << tab << molecule_coverage.top() << tab <<
                middle_mol_coverage.top() << tab <<
                (double)molecule_length.top()/molecule_coverage.top() << tab <<
                (double)molecule_read_density.top()/molecule_coverage.top() << tab <<
                starting_molecules.top() << tab << ending_molecules.top() << endl;

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
        starting_molecules.top() << tab << ending_molecules.top() << endl;

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


      end_ctg = max(end_ctg,end);

      //print the positions before the beging of the molecule (no next molecule will cover them)
      for(int i = begin_pos_stack; i < beg; i++){

        out_mol_stat << ctg << tab << i << tab << molecule_coverage.top() << tab <<
        middle_mol_coverage.top() << tab <<
        (double)molecule_length.top()/molecule_coverage.top() << tab <<
        (double)molecule_read_density.top()/molecule_coverage.top() << tab <<
        starting_molecules.top() << tab << ending_molecules.top() << endl;

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
          starting_molecules.top() << tab << ending_molecules.top() << endl;

      molecule_coverage.pop();
      middle_mol_coverage.pop();
      molecule_length.pop();
      molecule_read_density.pop();
      starting_molecules.pop();
      ending_molecules.pop();

  }

  out_mol_stat.close();

}
