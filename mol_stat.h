#pragma once

#include <string>

void molecule_stat(std::string tmp_file, std::string stat_file);
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
    std::string                            &statsFileName1);

