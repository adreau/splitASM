#include <string>
#include <vector>


void detect_outliers(std::string stat_per_interval_file, std::string out_contig_struct_file, std::string out_distance_scores_file, long n_sample);
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
    std::string                              &statsFileName2);
