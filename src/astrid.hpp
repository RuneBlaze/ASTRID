#ifndef ASTRID_HPP__
#define ASTRID_HPP__

#include "Args.hpp"
#include "DistanceMethods/DistanceMethods.hpp"
#include "multind.hpp"
#include "octal.hpp"
#include "phylokit/newick.hpp"


DistanceMatrix get_distance_matrix(TaxonSet &ts,
                                   const std::vector<std::string>& newicks,
                                   std::vector<double> weights,
                                   std::vector<Clade> &tree_taxa,
                                   IndSpeciesMapping *imap);

DistanceMatrix get_distance_matrix(TaxonSet &ts,
                                   const std::vector<std::string>& newicks,
                                   IndSpeciesMapping *imap);

DistanceMatrix mk_distance_matrix(TaxonSet &ts,
                                   const std::vector<std::string>& newicks);

DistanceMatrix get_distance_matrix(TaxonSet &ts,
                                   const std::vector<std::string>& newicks,
                                   std::vector<Clade> &tree_taxa,
                                   IndSpeciesMapping *imap);

void fill_in(TaxonSet &ts, DistanceMatrix &dm, const std::string& tree, bool fill_in = true);
// void finalize(TaxonSet &ts, DistanceMatrix &dm, const std::string& tree);

#endif