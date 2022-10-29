#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>
#include <cmath>
#include <ctime>
#include <iomanip>

#include "outliers_det.h"
const char tab='\t';

long checkArray(long id_current, long n_sample, std::vector<long> &id_sample) {
  long i;

  for (i = 0; i < n_sample; i++) {
    if (id_current == id_sample[i])
      return (1);
  }
  return (0);
}
// random sampling
// if seed = 0, seed is randomely set
void sampling(long n_sample, long min, long max, long seed, std::vector<long> &id_sample) {
  long i, j = 0;
  time_t t;

  if (seed == 0)
    seed = time(&t);

  for (i = 0; i < n_sample; i++) {
    do {
      srand(seed + j++);
      id_sample[i] = min + rand() % (max - min + 1);
    } while (checkArray(id_sample[i], i, id_sample) == 1);
  }
}


// input  -> X: data matrix (array), n: # of objects, d: # of dimensions, Xs:
// sample set, sid: sample indexes, ns: # of samples
// output -> result: an array of qsps
void qsp(std::vector<double> &X, long n, long d, long n_sample, long seed, std::vector<double> &score) {
  long i, j, point;
  std::vector<long> id_sample;
  double sum, res;

  // random sampling
  for(int k=0;k<n_sample;k++){
    id_sample.push_back(0);
  }

  sampling(n_sample, 0, n - 1, seed, id_sample);

  // compute the outlierness score qsp for each data point
  for (point = 0; point < n; point++) {
    res = 0;
    for (i = 0; i < n_sample; i++) {
      if (point != id_sample[i]) {
        sum = 0;
        for (j = 0; j < d; j++)
          sum += fabs(X[point * d + j] - X[id_sample[i] * d + j]) *
                 fabs(X[point * d + j] - X[id_sample[i] * d + j]);
        if (res == 0)
          res = sum;
        else if (sum < res && sum > 0)
          res = sum;
      }
    }
    score[point] = sqrt(res);
  }
}


// normalization of data (divide by SDs for each dimension)
void normalize(std::vector<double> &X, long n, long d) {
  long i, j, flag;
  double sum;
  std::vector<double> X_means;

  for(int k=0;k<d;k++)
    X_means.push_back(0);

  for (j = 0; j < d; j++) {
    sum = 0;
    for (i = 0; i < n; i++) {
      sum += X[i * d + j];
    }
    X_means[j] = sum / n;
  }

  for (j = 0; j < d; j++) {
    flag = 0;
    sum = X[j];
    for (i = 1; i < n; i++) {
      if (sum != X[i * d + j]) {
        flag = 1;
        break;
      }
    }
    if (flag == 1) {
      sum = 0;
      for (i = 0; i < n; i++)
        sum += (X[i * d + j] - X_means[j]) * (X[i * d + j] - X_means[j]);
      sum = sqrt(sum / (n - 1)); // unbiased SD
      for (i = 0; i < n; i++)
        X[i * d + j] = X[i * d + j] / sum;
    }
  }

}

void compute_score(std::vector<double> &stat, long n, long d, long n_sample, long seed, std::vector<double> &score){

  score = std::vector < double > (stat.size(), 0);
  normalize(stat, n, d);

  if(n_sample != 0) {
    n_sample = std::min(n_sample,n);
  }else{
    n_sample = n;
  }
  qsp(stat, n, d, n_sample, seed, score);

}

bool splitCondition(double score_molCov, double score_molLength, double score_readDens, double score_molBeg, double score_molEnd){

  int countSuperOut = (score_molCov >= 0.1 ? 1 : 0) + (score_molLength >= 0.1 ? 1 : 0) + (score_readDens >= 0.1 ? 1 : 0) + (score_molBeg >= 0.1 ? 1 : 0) + (score_molEnd >= 0.1 ? 1 : 0);

  int countOut = (score_molCov >= 0.01 ? 1 : 0) + (score_molLength >= 0.01 ? 1 : 0) + (score_readDens >= 0.01 ? 1 : 0) + (score_molBeg >= 0.01 ? 1 : 0) + (score_molEnd >= 0.01 ? 1 : 0);

  bool condition1 = countSuperOut >= 2 && countOut >= 3 ;
  bool condition2 = countSuperOut >= 1 && countOut >= 4 ;

  return condition1 || condition2;
}

void detect_outliers(std::string stat_per_interval_file, std::string out_contig_struct_file, std::string out_distance_scores_file, long n_sample)
{


  std::string molecule_stat_file = stat_per_interval_file;
  std::ifstream molecule_stat(molecule_stat_file.c_str());


  std::ofstream out_contig_struct(out_contig_struct_file.c_str(), std::ofstream::out);
  std::ofstream out_distance_scores(out_distance_scores_file.c_str(), std::ofstream::out);

  std::string line;
  int count_lines = 0;

  std::string molecule_line, ctg;
  int pos_beg, pos_end;
  double mol_cov, mid_mol_cov, mean_length, read_dens, mol_beg, mol_end;

  std::string prev_ctg = "", current_ctg = "";

  std::vector<int> posBegList;
  std::vector<int> posEndList;
  std::vector<double> midMolCovList;
  std::vector<double> meanLengthList;
  std::vector<double> readDensList;
  std::vector<double> molBegList;
  std::vector<double> molEndList;

  std::vector<double> score_molCov;
  std::vector<double> score_molLength;
  std::vector<double> score_readDens;
  std::vector<double> score_molBeg;
  std::vector<double> score_molEnd;

  std::vector<double> values_multipleD;
  std::vector<double> score_multipeD;


  long n, d, seed;
// n_sample;
  int beginInterval, endInterval, count_interval;
  bool cut = false;

  while(getline(molecule_stat,molecule_line)){

      std::stringstream  linestream(molecule_line);
    linestream >> ctg >> pos_beg >> pos_end >> mol_cov >> mid_mol_cov >> mean_length >> read_dens >> mol_beg >> mol_end;
    current_ctg = ctg;

    if(current_ctg.compare(prev_ctg)!=0 && prev_ctg.length()>0){

        n = meanLengthList.size();
        d = 1;
        seed = 0;
       // n_sample = meanLengthList.size()-1;

        compute_score(midMolCovList, n, d, n_sample, seed, score_molCov);
        compute_score(meanLengthList, n, d, n_sample, seed, score_molLength);
        compute_score(readDensList, n, d, n_sample, seed, score_readDens);
        compute_score(molBegList, n, d, n_sample, seed, score_molBeg);
        compute_score(molEndList, n, d, n_sample, seed, score_molEnd);


        /*for multiple dimensions once */
        values_multipleD.clear();
        values_multipleD.insert(values_multipleD.end(), midMolCovList.begin(), midMolCovList.end());
        values_multipleD.insert(values_multipleD.end(), meanLengthList.begin(), meanLengthList.end());
        values_multipleD.insert(values_multipleD.end(), readDensList.begin(), readDensList.end());
        values_multipleD.insert(values_multipleD.end(), molBegList.begin(), molBegList.end());
        values_multipleD.insert(values_multipleD.end(), molEndList.begin(), molEndList.end());
        compute_score(values_multipleD, n, 5, n_sample, seed, score_multipeD);

        beginInterval = posBegList[0];
        endInterval = posBegList[0];
        count_interval = 0;
        cut = false;

        for (int i = 0; i < n; i++) {
          out_distance_scores << std::fixed << std::setprecision(6)<< prev_ctg << '\t'<< posBegList[i]
                              <<'\t'<< posEndList[i]<<'\t'<<score_molCov[i]<<'\t'<<score_molLength[i]
                              <<'\t'<< score_readDens[i]<<'\t'<<score_molBeg[i]<<'\t'<<score_molEnd[i]<< '\t'<<values_multipleD[i]<< std::endl;

          if( splitCondition(score_molCov[i], score_molLength[i],score_readDens[i], score_molBeg[i], score_molEnd[i]) ){

               if(!cut){

                 out_contig_struct << prev_ctg <<'\t' << beginInterval << '\t' << posBegList[i]<< std::endl;
                 beginInterval = posBegList[i];
                 cut = true;

               }
          }else{
            if(cut){

              out_contig_struct << prev_ctg << '\t' << beginInterval << '\t' << posBegList[i]<< std::endl;
              cut = false;
              beginInterval = posBegList[i];

            }
          }

        }

        out_contig_struct << prev_ctg << '\t' << beginInterval << '\t' << posEndList[n-1]<< std::endl;

        posBegList.clear();
        posEndList.clear();
        midMolCovList.clear();
        meanLengthList.clear();
        readDensList.clear();
        molBegList.clear();
        molEndList.clear();

    }

    prev_ctg = ctg;
    //cout<<mid_mol_cov<<endl;
    posBegList.push_back(pos_beg);
    posEndList.push_back(pos_end);
    midMolCovList.push_back(mid_mol_cov);
    meanLengthList.push_back(mean_length);
    readDensList.push_back(read_dens);
    molBegList.push_back(mol_beg);
    molEndList.push_back(mol_end);

  }

  //for the last contig (duplicate code => TODO: put in function)

  n = meanLengthList.size();
  d = 1;
  seed = 0;
 // n_sample = meanLengthList.size()-1;

  compute_score(midMolCovList, n, d, n_sample, seed, score_molCov);
  compute_score(meanLengthList, n, d, n_sample, seed, score_molLength);
  compute_score(readDensList, n, d, n_sample, seed, score_readDens);
  compute_score(molBegList, n, d, n_sample, seed, score_molBeg);
  compute_score(molEndList, n, d, n_sample, seed, score_molEnd);

  beginInterval = posBegList[0];
  endInterval = posBegList[0];
  count_interval = 0;
  cut = false;

  for (int i = 0; i < n; i++) {
    out_distance_scores << std::fixed << std::setprecision(6)<< ctg << '\t'<< posBegList[i]
                        <<'\t'<< posEndList[i]<<'\t'<<score_molCov[i]<<'\t'<<score_molLength[i]
                        <<'\t'<< score_readDens[i]<<'\t'<<score_molBeg[i]<<'\t'<<score_molEnd[i]<<'\t'<<values_multipleD[i]<< std::endl;

    if( splitCondition(score_molCov[i], score_molLength[i],score_readDens[i], score_molBeg[i], score_molEnd[i])){

         if(!cut){

           if (posBegList[i]>beginInterval ){
             out_contig_struct << ctg <<'\t' << beginInterval << '\t' << posBegList[i]<< std::endl;
           }
           beginInterval = posBegList[i];
           cut = true;

         }
    }else{
      if(cut){

        out_contig_struct << ctg << '\t' << beginInterval << '\t' << posBegList[i]<< std::endl;
        cut = false;
        beginInterval = posBegList[i];

      }
    }

  }

  out_contig_struct << ctg << '\t' << beginInterval << '\t' << posEndList[n-1]<< std::endl;




  out_distance_scores.close();
  out_contig_struct.close();


}


void detect_outliers2(
    std::vector < std::vector < double > >   &molecule_coverages,
    std::vector < std::vector < double > >   &mid_molecule_coverages,
    std::vector < std::vector < double > >   &mean_lengths,
    std::vector < std::vector < double > >   &read_densities,
    std::vector < std::vector < double > >   &n_starts,
    std::vector < std::vector < double > >   &n_ends,
    std::vector < std::vector < double > >   &score_molecule_coverages,
    std::vector < std::vector < double > >   &score_mid_molecule_coverages,
    std::vector < std::vector < double > >   &score_mean_lengths,
    std::vector < std::vector < double > >   &score_read_densities,
    std::vector < std::vector < double > >   &score_n_starts,
    std::vector < std::vector < double > >   &score_n_ends,
    std::vector < std::string >              &ctg_names,
    std::vector < long int >                 &ctg_sizes,
    int                                       window,
    int                                       min_ctg_size,
    long                                      n_samples,
    std::string                              &output_file_name,
    std::string                              &statsFileName2) {

  std::ofstream output_file(output_file_name, std::ofstream::out);
  std::ofstream statsFile2;
  if (! statsFileName2.empty()) statsFile2.open(statsFileName2, std::ofstream::out);

  size_t nchrs = ctg_names.size();
  int seed = 0;
  score_molecule_coverages.resize(nchrs);
  score_mid_molecule_coverages.resize(nchrs);
  score_mean_lengths.resize(nchrs);
  score_read_densities.resize(nchrs);
  score_n_starts.resize(nchrs);
  score_n_ends.resize(nchrs);
  for (size_t chrid = 0; chrid < nchrs; ++chrid) {
    int d = 1;
    size_t npos = molecule_coverages[chrid].size();
    std::string &ctg = ctg_names[chrid];
    compute_score(molecule_coverages[chrid],     npos, d, n_samples, seed, score_molecule_coverages[chrid]);
    compute_score(mid_molecule_coverages[chrid], npos, d, n_samples, seed, score_mid_molecule_coverages[chrid]);
    compute_score(mean_lengths[chrid],           npos, d, n_samples, seed, score_mean_lengths[chrid]);
    compute_score(read_densities[chrid],         npos, d, n_samples, seed, score_read_densities[chrid]);
    compute_score(n_starts[chrid],               npos, d, n_samples, seed, score_n_starts[chrid]);
    compute_score(n_ends[chrid],                 npos, d, n_samples, seed, score_n_ends[chrid]);

    bool incut          = false;
    int  bin_frag_start = 0;

    for (int pos = 0; pos < npos; ++pos) {
      if (splitCondition(
            score_molecule_coverages[chrid][pos],
            score_mean_lengths[chrid][pos],
            score_read_densities[chrid][pos],
            score_n_starts[chrid][pos],
            score_n_ends[chrid][pos])) {
        if (! incut) {
          // Check contig size
          int size = (pos - bin_frag_start) * window;
          if (size >= min_ctg_size) {
            // Do not print if we start with a cut
            // Output in BED format
            if (pos > 0) {
              output_file << ctg << '\t' << bin_frag_start * window << '\t' << pos * window << '\n';
            }
          }
        }
        incut = true;
      }
      else if (incut) {
        incut          = false;
        bin_frag_start = pos;
      }
      if (! statsFileName2.empty()) statsFile2 << std::fixed << std::setprecision(6) <<
          ctg                                      << '\t' <<
          pos * window + 1                         << '\t' <<
          (pos + 1) * window                       << '\t' <<
          score_molecule_coverages[chrid][pos]     << '\t' <<
          score_mid_molecule_coverages[chrid][pos] << '\t' <<
          score_mean_lengths[chrid][pos]           << '\t' <<
          score_read_densities[chrid][pos]         << '\t' <<
          score_n_starts[chrid][pos]               << '\t' <<
          score_n_ends[chrid][pos]                 << '\n';
    }
    std::cout << "Analyzing contig #" << chrid << "/" << nchrs << ".\r" << std::flush;
    if (! incut) {
      // Output in BED format
      output_file << ctg << '\t' << bin_frag_start * window << '\t' << ctg_sizes[chrid] << '\n';
    }
  }
  std::cout << "Analyzing contig #" << nchrs << "/" << nchrs << ".\n";

  output_file.close();

  if (! statsFileName2.empty()) statsFile2.close();
}
