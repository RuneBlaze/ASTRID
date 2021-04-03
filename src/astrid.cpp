#include <fstream>
#include <glog/logging.h>
#include <iostream>
#include "astrid.hpp"

bool has_missing(TaxonSet &ts, DistanceMatrix &dm) {
  for (size_t i = 0; i < ts.size(); i++) {
    for (size_t j = i + 1; j < ts.size(); j++) {
      if (!dm.has(i, j)) {
        return true;
      }
    }
  }
  return false;
}

TaxonSet get_ts(std::vector<std::string> &newicks) {
  std::unordered_set<std::string> taxa;
  for (std::string n : newicks) {
    newick_to_ts(n, taxa);
  }
  TaxonSet ts(taxa.size());
  for (std::string t : taxa) {
    ts.add(t);
  }
  return ts;
}

DistanceMatrix get_distance_matrix(TaxonSet &ts,
                                   std::vector<std::string> newicks,
                                   std::vector<double> weights,
                                   std::vector<Clade> &tree_taxa,
                                   IndSpeciesMapping *imap) {

  DistanceMatrix result(ts);
  if (imap) {
    result = DistanceMatrix(imap->species());
  }

  for (size_t i = 0; i < newicks.size(); i++) {
    std::string newick_derooted = deroot(newicks[i]);
    double w = weights[i];

    DistanceMatrix dm(ts, newick_derooted);
    if (tree_taxa.size() > i) {
      for (Taxon t1 : ts) {
        for (Taxon t2 : ts) {
          if (!(tree_taxa[i].contains(t1) && tree_taxa[i].contains(t2))) {
            dm(t1, t2) = 0;
            dm.masked(t1, t2) = 0;
          }
        }
      }
    }

    dm *= w;

    if (imap) {
      dm = imap->average(dm);
    }

    result += dm;
  }

  for (size_t i = 0; i < ts.size(); i++) {
    for (size_t j = i; j < ts.size(); j++) {
      if (result.masked(i, j))
        result(i, j) /= result.masked(i, j);
    }
  }

  return result;
}

DistanceMatrix get_distance_matrix(TaxonSet &ts,
                                   std::vector<std::string> newicks,
                                   IndSpeciesMapping *imap) {
  std::vector<Clade> vc;
  return get_distance_matrix(ts, newicks,
                             std::vector<double>(newicks.size(), 1), vc, imap);
}

DistanceMatrix get_distance_matrix(TaxonSet &ts,
                                   std::vector<std::string> newicks,
                                   std::vector<Clade> &tree_taxa,
                                   IndSpeciesMapping *imap) {
  return get_distance_matrix(
      ts, newicks, std::vector<double>(newicks.size(), 1), tree_taxa, imap);
}

std::string run_astrid(std::vector<std::string> newicks) { return ""; }

void fill_in(TaxonSet &ts, DistanceMatrix &dm, double cval) {
  int count = 0;
  for (size_t i = 0; i < ts.size(); i++) {
    for (size_t j = i; j < ts.size(); j++) {
      if (!dm.has(i, j)) {
        dm(i, j) = cval;
        dm.masked(i, j) = 1;
        count++;
      }
    }
  }
  std::cerr << "Filled in " << count << " elements" << std::endl;
}

void fill_in(TaxonSet &ts, DistanceMatrix &dm, std::string tree, bool fill_in) {
  std::vector<std::string> trees;
  trees.push_back(tree);
  DistanceMatrix dm_tree = get_distance_matrix(ts, trees, NULL);

  for (size_t i = 0; i < ts.size(); i++) {
    for (size_t j = i; j < ts.size(); j++) {
      if (!dm.has(i, j)) {
        dm(i, j) = dm_tree(i, j);
        if (fill_in) dm.masked(i, j) = 1;
      }
    }
  }
}

void prune(TaxonSet &ts, DistanceMatrix &dm, int threshold) {
  for (size_t i = 0; i < ts.size(); i++) {
    for (size_t j = i; j < ts.size(); j++) {
      if (dm.masked(i, j) <= threshold) {
        dm(i, j) = 0;
        dm.masked(i, j) = 0;
      }
    }
  }
}

void finalize(TaxonSet &ts, DistanceMatrix &dm) {
  for (size_t i = 0; i < ts.size(); i++) {
    for (size_t j = i; j < ts.size(); j++) {
      if (!dm.has(i, j)) {
        dm.masked(i, j) = 1;
      }
    }
  }
}

void progressbar(double pct) {
  std::cerr << "[";
  for (int i = 0; i < (int)(68 * pct); i++) {
    std::cerr << "#";
  }
  for (int i = (int)(68 * pct); i < 68; i++) {
    std::cerr << "-";
  }
  std::cerr << "]\r";
}

int main(int argc, char **argv) {
  Args args(argc, argv);

  std::vector<std::string> input_trees;
  std::ifstream inf(args.infile);

  std::string buf;
  LOG(INFO) << "Reading trees..." << std::endl;
  while (!inf.eof()) {
    getline(inf, buf);
    if (buf.size() > 3)
      input_trees.push_back(buf);
  }
  LOG(INFO) << "Read " << input_trees.size() << " trees" << std::endl;

  TaxonSet ts = get_ts(input_trees);
  int iter = 1;
  std::string tree;

  // Set up multi-individual mapping
  IndSpeciesMapping *multind_mapping = nullptr;
  if (args.multindfile.size()) {
    multind_mapping = new IndSpeciesMapping(ts);
    multind_mapping->load(args.multindfile);
  }

  DistanceMatrix dm = get_distance_matrix(ts, input_trees, multind_mapping);

  std::cerr << "Estimating tree" << std::endl;
  for (std::string method : args.dms) {
    std::cerr << "Running " << method << std::endl;

    // OCTAL completion of trees with missing data
    if (tree.size() && args.octal) {

      std::vector<std::string> completed_trees;
      std::vector<Clade> tree_taxa;
      Tree T = newick_to_treeclades(tree, ts);
      for (std::string t_s : input_trees) {
        Tree t = newick_to_treeclades(t_s, ts);
        tree_taxa.push_back(t.taxa());

        octal_complete(T, t);
        std::stringstream ss;
        ss << t;
        completed_trees.push_back(ss.str());
      }
      dm = get_distance_matrix(ts, input_trees, multind_mapping);
    }

    TaxonSet *species_ts = &ts;
    if (multind_mapping) {
      species_ts = &(multind_mapping->species());
    }

    // fill in missing elements on second and later iterations
    if (iter > 1) {
      fill_in(*species_ts, dm, tree);
    } else if (args.constant != 0) {
      fill_in(*species_ts, dm, args.constant);
    }

    if (method == "auto") {
      if (has_missing(*species_ts, dm)) {
        std::cerr << "Missing entries in distance matrix, trying to run BioNJ*"
                  << std::endl;
        std::cerr << "You may have better results adding -u -s to your command "
                     "line to use UPGMA completion instead."
                  << std::endl;
        tree = BioNJStar(*species_ts, dm, args.java_opts);
      } else {
        std::cerr
            << "No missing entries in distance matrix, running FastME2+SPR"
            << std::endl;
        tree = FastME(*species_ts, dm, 1, 1);
      }
    } else if (method == "upgma") {
      tree = UPGMA(*species_ts, dm);
    } else if (method == "fastme") {
      tree = FastME(*species_ts, dm, 0, 0);
    } else if (method == "fastme_nni") {
      tree = FastME(*species_ts, dm, 1, 0);
    } else if (method == "fastme_spr") {
      tree = FastME(*species_ts, dm, 1, 1);
    } else if (method == "bionj") {
      tree = BioNJStar(*species_ts, dm, args.java_opts);
    } else if (method == "rapidnj") {
      tree = RapidNJ(*species_ts, dm);
    }

    std::ofstream outfile(args.outfile + "." + std::to_string(iter));
    outfile << tree << std::endl;
    if (args.cache) {
      std::ofstream outfile_cache(args.cachefile + "." + std::to_string(iter));
      dm.writePhylip(outfile_cache);
    }
    iter++;
  }
  std::ofstream outfile(args.outfile);
  outfile << tree << std::endl;

  std::cout << tree << std::endl;

  return 0;
}

// FIXME: I probably introduced memory leaks
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

DistanceMatrix mk_distance_matrix(TaxonSet &ts,
                                   std::vector<std::string> newicks) {
  return get_distance_matrix(ts, newicks,nullptr);
}

PYBIND11_MODULE(asterid, m) {
    m.doc() = "pybind11 astrid";
    
    py::class_<TaxonSet>(m, "TaxonSet")
        .def(py::init<int>())
        .def("size", &TaxonSet::size)
        .def("__len__", &TaxonSet::size)
        .def("__iter__", [](TaxonSet &ts){
          return py::make_iterator(ts.begin(), ts.end());
        },py::keep_alive<0, 1>())
        .def("add", &TaxonSet::add)
        .def("__getitem__", [](TaxonSet &ts, std::string name) {
          return ts[name];
        });

    py::class_<DistanceMatrix>(m, "DistanceMatrix")
        .def(py::init<const TaxonSet &>())
        .def("str", [](DistanceMatrix &m) {
          return m.str();
        })
        .def("fill_in", [](DistanceMatrix &m, TaxonSet &ts, std::string tree) {
          fill_in(ts, m, tree);
        })
        .def("has", [](DistanceMatrix &m, Taxon first, Taxon second){
          return m.has(first, second) > 0;
        })
        .def("fill_in_transient", [](DistanceMatrix &m, TaxonSet &ts, std::string tree) {
          fill_in(ts, m, tree, false);
        })
        .def("finalize",  [](DistanceMatrix &m, TaxonSet &ts){
          finalize(ts, m);
        })
        .def("prune", [](DistanceMatrix &dm, TaxonSet &ts, int threshold){
          prune(ts, dm, threshold);
        })
        .def("__getitem__", [](DistanceMatrix &m, std::pair<Taxon, Taxon> taxons){
          // this possibly UBs
          return m.get(taxons.first, taxons.second);
        })
        .def("__setitem__", [](DistanceMatrix &m, std::pair<Taxon, Taxon> taxons, int value){
          // this possibly UBs
          m(taxons.first, taxons.second) = value;
        })
        .def("getmask", [](DistanceMatrix &m, std::pair<Taxon, Taxon> taxons){
          // this possibly UBs
          return m.masked(taxons.first, taxons.second);
        })
        .def("setmask", [](DistanceMatrix &m, std::pair<Taxon, Taxon> taxons, int value){
          // this possibly UBs
          m.masked(taxons.first, taxons.second) = value;
        });
    m.def("mk_distance_matrix", &mk_distance_matrix, "making distance matrices");
    m.def("get_ts", &get_ts, "getting taxonset");
    m.def("upgma_star", &UPGMA, "UPGMA*");
    m.def("fastme_balme", &FastME, "FastME");
    m.def("rapidnj", &RapidNJ, "RapidNJ");
    m.def("deroot", &deroot, "deroot");
}