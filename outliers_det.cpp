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

bool checkArray(long id_current, long n_sample, std::vector<long> &id_sample) {
  long i;

  for (i = 0; i < n_sample; i++) {
    if (id_current == id_sample[i])
      return true;
  }
  return false;
}

// random sampling
// if seed = 0, seed is randomely set
void sampling(long n_sample, long max, long seed, std::vector<long> &id_sample) {
  std::vector < bool > used (max, false);
  long i, j = 0, nused = 0;
  time_t t;

  // # samples is not less than sampling size
  if (n_sample >= max) {
    std::cout << "\tNo sampling.\n";
    for (int i = 0; i < n_sample; ++i) {
      id_sample[i] = i;
    }
    return;
  }

  std::cout << "\tSampling...\n";
  if (seed == 0)
    seed = time(&t);

  while (nused < n_sample) {
    srand(seed + j++);
    long id = rand() % max;
    //long id = rand() % (max + 1);
    if (! used[id]) {
      id_sample[nused] = id;
      used[id] = true;
      ++nused;
    }
  }
  std::cout << "\t... done.\n";
}


// input  -> X: data matrix (array), n: # of objects, d: # of dimensions, Xs:
// sample set, sid: sample indexes, ns: # of samples
// output -> result: an array of qsps
void qsp(std::vector < double > &X, long n, int d, long n_sample, long seed, std::vector < double > &score) {
  std::vector < long > id_sample (n_sample, 0);

  sampling(n_sample, n, seed, id_sample);

  // compute the outlierness score qsp for each data point
  for (long unsigned point = 0; point < n; point++) {
    double res = 0;
    for (long unsigned i = 0; i < n_sample; i++) {
      if (point != id_sample[i]) {
        double sum = 0;
        for (int j = 0; j < d; j++) {
          // std::cout << "point: " << point << ", d: " << d << ", j: " << j << " / " << X.size() << "\n";
          // std::cout << "sample: " << id_sample[i] << " / " << n << " -> " << (id_sample[i] * d + j) << "\n";
          sum += fabs(X[point * d + j] - X[id_sample[i] * d + j]) *
                 fabs(X[point * d + j] - X[id_sample[i] * d + j]);
        }
        if (res == 0)
          res = sum;
        else if (sum < res && sum > 0)
          res = sum;
      }
    }
    score[point] = sqrt(res);
    if (point % 1000 == 0) std::cout << "\t" << point << "/" << n << "\r" << std::flush;
  }
  std::cout << "\t" << n << "/" << n << "\n";
}


// normalization of data (divide by SDs for each dimension)
void normalize(std::vector < double > &X, long n, int d) {
  std::vector < double > X_means (d, 0);

  for (int j = 0; j < d; j++) {
    double sum = 0;
    for (unsigned long int i = 0; i < n; i++) {
      sum += X[i * d + j];
    }
    X_means[j] = sum / n;
  }

  for (int j = 0; j < d; j++) {
    bool flag = false;
    double sum = X[j];
    for (unsigned long int i = 1; i < n; i++) {
      if (sum != X[i * d + j]) {
        flag = true;
        break;
      }
    }
    if (flag) {
      sum = 0;
      for (unsigned long i = 0; i < n; i++)
        sum += (X[i * d + j] - X_means[j]) * (X[i * d + j] - X_means[j]);
      sum = sqrt(sum / (n - 1)); // unbiased SD
      for (unsigned long i = 0; i < n; i++)
        X[i * d + j] = X[i * d + j] / sum;
    }
  }
}

void compute_score(std::vector<double> &stat, long n, long d, long n_sample, long seed, std::vector<double> &score){

  score = std::vector < double > (stat.size(), 0);
  normalize(stat, d, n);

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
    std::vector < std::vector < double > >   &mean_lengths,
    std::vector < std::vector < double > >   &read_densities,
    std::vector < std::vector < double > >   &n_starts,
    std::vector < std::vector < double > >   &n_ends,
    std::vector < std::string >              &ctg_names,
    std::vector < long int >                 &ctg_sizes,
    int                                       window,
    int                                       min_ctg_size,
    long                                      n_samples,
    std::string                              &output_file_name,
    std::string                              &statsFileName2,
    float                                     threshold) {

  std::ofstream output_file(output_file_name, std::ofstream::out);
  std::ofstream statsFile2;
  if (! statsFileName2.empty()) statsFile2.open(statsFileName2, std::ofstream::out);

  size_t nchrs = ctg_names.size();
  size_t n_elements = 0;
  int seed = 0;
  for (size_t chrid = 0; chrid < nchrs; ++chrid) {
    n_elements += molecule_coverages[chrid].size();
  }
  int d = 5;
  std::vector < double > all_values (d * n_elements);
  std::vector < double > score_all_values (n_elements);
  /*
  for (size_t chrid = 0; chrid < nchrs; ++chrid) {
    all_values.insert(all_values.end(),     molecule_coverages[chrid].begin(),     molecule_coverages[chrid].end());
  }
  for (size_t chrid = 0; chrid < nchrs; ++chrid) {
    all_values.insert(all_values.end(),           mean_lengths[chrid].begin(),           mean_lengths[chrid].end());
  }
  for (size_t chrid = 0; chrid < nchrs; ++chrid) {
    all_values.insert(all_values.end(),         read_densities[chrid].begin(),         read_densities[chrid].end());
  }
  for (size_t chrid = 0; chrid < nchrs; ++chrid) {
    all_values.insert(all_values.end(),               n_starts[chrid].begin(),               n_starts[chrid].end());
  }
  for (size_t chrid = 0; chrid < nchrs; ++chrid) {
    all_values.insert(all_values.end(),                 n_ends[chrid].begin(),                 n_ends[chrid].end());
  }
  */
  size_t cpt = 0;
  for (size_t chrid = 0; chrid < nchrs; ++chrid) {
    for (size_t i = 0; i < molecule_coverages[chrid].size(); ++i) {
      all_values[cpt++] = molecule_coverages[chrid][i];
      all_values[cpt++] = mean_lengths      [chrid][i];
      all_values[cpt++] = read_densities    [chrid][i];
      all_values[cpt++] = n_starts          [chrid][i];
      all_values[cpt++] = n_ends            [chrid][i];
    }
  }
  if (cpt != n_elements) std::cerr << "Problem while filling the values: expected " << n_elements << ", got " << cpt << "...\n";

  compute_score(all_values, n_elements, d, n_samples, seed, score_all_values);

//std::cerr << "Detecting outliers in coverages...\n";
//compute_score(all_molecule_coverages,     n_elements, d, n_samples, seed, all_score_molecule_coverages);
//std::cerr << "Detecting outliers in middle coverages...\n";
//compute_score(all_mid_molecule_coverages, n_elements, d, n_samples, seed, all_score_mid_molecule_coverages);
//std::cerr << "Detecting outliers in molecule lengths...\n";
//compute_score(all_mean_lengths,           n_elements, d, n_samples, seed, all_score_mean_lengths);
//std::cerr << "Detecting outliers in read densities...\n";
//compute_score(all_read_densities,         n_elements, d, n_samples, seed, all_score_read_densities);
//std::cerr << "Detecting outliers in # start positions...\n";
//compute_score(all_n_starts,               n_elements, d, n_samples, seed, all_score_n_starts);
//std::cerr << "Detecting outliers in # end positions...\n";
//compute_score(all_n_ends,                 n_elements, d, n_samples, seed, all_score_n_ends);
//std::cerr << "... done.\n";

  cpt = 0;
  for (size_t chrid = 0; chrid < nchrs; ++chrid) {
    int  npos           = molecule_coverages[chrid].size();
    bool incut          = false;
    int  bin_frag_start = 0;
    std::string &ctg    = ctg_names[chrid];

    for (int pos = 0; pos < npos; ++pos, ++cpt) {
      if (score_all_values[cpt] > threshold) {
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
          ctg                   << '\t' <<
          pos * window + 1      << '\t' <<
          (pos + 1) * window    << '\t' <<
          score_all_values[cpt] << '\n';
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
