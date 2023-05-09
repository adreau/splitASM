#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

#include "mol_stat.h"

void compute_stats (Molecule_stats &molecule_stats, Molecules &molecules) {

  unsigned long cpt = 0;
  for (auto &molecule: molecules) {
    // Easier to work with 0-based positions
    --molecule.start;
    --molecule.end;

    unsigned int window_start = molecule.start / Globals::window;
    unsigned int window_end   = molecule.end   / Globals::window;
    unsigned int molecule_size = molecule.end - molecule.start + 1;
    assert(molecule.chrid < Globals::chrs.size());
    for (int windowid = window_start; windowid <= window_end; ++windowid) {
      assert(windowid < molecule_stats[molecule.chrid].size());
      unsigned int beg_window = std::max < unsigned int > (windowid * Globals::window, molecule.start);
      unsigned int end_window = std::min < unsigned int > (beg_window + Globals::window - 1, molecule.end);
      unsigned int size = end_window - beg_window + 1;
      molecule_stats[molecule.chrid][windowid].coverage     += size;
      molecule_stats[molecule.chrid][windowid].length       += molecule_size * size;
      molecule_stats[molecule.chrid][windowid].read_density += molecule.n_reads / molecule_size * size;
    }
    ++molecule_stats[molecule.chrid][window_start].start;
    ++molecule_stats[molecule.chrid][window_end].end;

    ++cpt;
    if (cpt % 10000000 == 0) std::cerr << cpt << " lines read.\r" << std::flush;
  }
  std::cerr << cpt << " lines read, done.\n";

  std::ofstream countsFile;
  if (! Globals::countsFileName.empty()) countsFile.open(Globals::countsFileName, std::ofstream::out);

  for (size_t chrid = 0; chrid < Globals::chrs.size(); ++chrid) {
    std::string &chr = Globals::chrs[chrid];
    long int size = Globals::chr_sizes[chrid] / Globals::window;
    for (size_t windowid = 0; windowid <= size; ++windowid) {
      unsigned int beg_window = windowid * Globals::window;
      unsigned int end_window = std::min < int > (beg_window + Globals::window - 1, Globals::chr_sizes[chrid]);
      double size_dbl = static_cast < double > (end_window - beg_window + 1);
      // Careful in the order of normalization
      molecule_stats[chrid][windowid].coverage   /= size_dbl;
      if (molecule_stats[chrid][windowid].coverage == 0) {
        molecule_stats[chrid][windowid].length       = 0;
        molecule_stats[chrid][windowid].read_density = 0;
      }
      else {
        molecule_stats[chrid][windowid].length       /= molecule_stats[chrid][windowid].coverage * Globals::window;
        molecule_stats[chrid][windowid].read_density /= molecule_stats[chrid][windowid].coverage * Globals::window;
      }

      // Switch back to 1-based positions
      if (! Globals::countsFileName.empty()) {
        countsFile << chr                              << tab <<
          windowid * Globals::window + 1               << tab <<
          (windowid + 1) * Globals::window             << tab <<
          molecule_stats[chrid][windowid].coverage     << tab <<
          molecule_stats[chrid][windowid].coverage     << tab <<
          molecule_stats[chrid][windowid].read_density << tab <<
          molecule_stats[chrid][windowid].start        << tab <<
          molecule_stats[chrid][windowid].end          << "\n";
      }
    }
  }
  if (! Globals::countsFileName.empty()) countsFile.close();
}
