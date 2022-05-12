#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

const char tab='\t';



void createBed (vector<string> tmpfiles, string contigSizeFile, string bedFile)
{

  ifstream contig_size(contigSizeFile.c_str());

  ofstream split_bed_corrected(bedFile.c_str(), ofstream::out);

  string line;


  string split_line, contig_line;
  string prev_ctg = "",ctg="",ctg_off="";
  int pos_beg,pos_end,ctg_size = 0;
  int prev_pos_beg, prev_pos_end;
  int concat_beg;
  bool coordinated = false;


  for(int i = 0; i<tmpfiles.size();i++){

    ifstream split_bed(tmpfiles[i].c_str());

    while(getline(split_bed,split_line)){

      stringstream  splitstream(split_line);
      splitstream >> ctg >> pos_beg >> pos_end;

      if(prev_ctg.compare(ctg)!=0){

        if(prev_ctg.length()>0){ //last split of prev contig
          split_bed_corrected << prev_ctg << tab << concat_beg << tab << ctg_size <<endl;
        }

        concat_beg = pos_beg;
        coordinated = false;

        while(!coordinated){
          getline(contig_size,contig_line);

          stringstream contigstream(contig_line);
          contigstream >> ctg_off >> ctg_size;

          if(ctg_off.compare(ctg)!=0){//contigs not present in the split
            split_bed_corrected << ctg_off << tab << "0" << tab << ctg_size << endl;
            coordinated=false;
          }else{
            coordinated=true;
          }
        }


      }else{

          if(prev_pos_end !=0 && ( (pos_end-pos_beg>20000) || (prev_pos_end-prev_pos_beg>20000) )){
              split_bed_corrected << prev_ctg << tab << concat_beg << tab << prev_pos_end << endl;
              concat_beg = pos_beg;
          }

      }

      prev_ctg = ctg;
      prev_pos_beg = pos_beg;
      prev_pos_end = pos_end;




    }


    split_bed.close();
  }

  while(getline(contig_size,contig_line)){

    stringstream contigstream(contig_line);
    contigstream >> ctg_off >> ctg_size;
    split_bed_corrected << ctg_off << tab << "0" << tab  << ctg_size <<endl;

  }

  split_bed_corrected.close();

}

/*int main (int argc, char* argv[])
{


  vector<string> tmpfiles;
  for(int i=1; i<=20; i++){

    tmpfiles.push_back(argv[i]);
  }

  string contigSizeFile =  argv[21];
  string bedFile =  argv[22];
  createBed (tmpfiles,contigSizeFile, bedFile);

}*/
