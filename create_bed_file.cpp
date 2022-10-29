#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

const char tab='\t';



void createBed (std::vector<std::string> tmpfiles, std::string contigSizeFile, std::string bedFile)
{

    std::ifstream contig_size(contigSizeFile.c_str());

    std::ofstream split_bed_corrected(bedFile.c_str(), std::ofstream::out);

  std::string line;


  std::string split_line, contig_line;
  std::string prev_ctg = "",ctg="",ctg_off="";
  int pos_beg,pos_end,ctg_size = 0;
  int prev_pos_beg, prev_pos_end;
  int concat_beg;
  bool coordinated = false;


  for(int i = 0; i<tmpfiles.size();i++){

      std::ifstream split_bed(tmpfiles[i].c_str());

    while(getline(split_bed,split_line)){

        std::stringstream  splitstream(split_line);
      splitstream >> ctg >> pos_beg >> pos_end;

      if(prev_ctg.compare(ctg)!=0){

        if(prev_ctg.length()>0){ //last split of prev contig
          split_bed_corrected << prev_ctg << tab << concat_beg << tab << ctg_size << std::endl;
        }

        concat_beg = pos_beg;
        coordinated = false;

        while(!coordinated){
          getline(contig_size,contig_line);

          std::stringstream contigstream(contig_line);
          contigstream >> ctg_off >> ctg_size;

          if(ctg_off.compare(ctg)!=0){//contigs not present in the split
            split_bed_corrected << ctg_off << tab << "0" << tab << ctg_size << std::endl;
            coordinated=false;
          }else{
            coordinated=true;
          }
        }


      }else{

          if(prev_pos_end !=0 && ( (pos_end-pos_beg>20000) || (prev_pos_end-prev_pos_beg>20000) )){
              split_bed_corrected << prev_ctg << tab << concat_beg << tab << prev_pos_end << std::endl;
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

      std::stringstream contigstream(contig_line);
    contigstream >> ctg_off >> ctg_size;
    split_bed_corrected << ctg_off << tab << "0" << tab  << ctg_size << std::endl;

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
