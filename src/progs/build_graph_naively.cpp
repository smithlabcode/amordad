/*
 *    Part of AMORDAD software
 *
 *    Copyright (C) 2014 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 */

#include <string>
#include <vector>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include "RegularNearestNeighborGraph.hpp"
#include "FeatureVector.hpp"

using std::string;
using std::vector;
using std::cerr;
using std::endl;
using std::cout;
using std::tr1::unordered_map;


/*****
 * This function uses naive comparisons to construct the graph.
 */
static void
construct_graph_naively(const bool VERBOSE,
                        const size_t max_degree,
                        const vector<FeatureVector> &fvs, 
                        RegularNearestNeighborGraph &nng) {
  
  for (size_t i = 0; i < fvs.size(); ++i)
    nng.add_vertex(fvs[i].get_id());
  
  const size_t total_comparisons = fvs.size()*(fvs.size() - 1)/2;
  size_t comparisons_done = 0;
  
  for (size_t i = 0; i < fvs.size(); ++i) {
    for (size_t j = 0; j < i; ++j) {
      const double w = fvs[i].compute_angle(fvs[j]);
      nng.update_vertex(fvs[i].get_id(), fvs[j].get_id(), w);
      nng.update_vertex(fvs[j].get_id(), fvs[i].get_id(), w);
    }
    
    comparisons_done += i;
    if (VERBOSE)
      cerr << '\r' << "computing edges: " 
           << percent(comparisons_done, total_comparisons) << "%\r";
  }
  if (VERBOSE)
    cerr << '\r' << "computing edges: 100%" << endl;
}


static void
load_filenames(const string &paths_file, vector<string> &filenames) { 
  std::ifstream in(paths_file.c_str());
  if (!in)
    throw SMITHLABException("could not read file: " + paths_file);
  filenames.clear();
  string tmp;
  while (in >> tmp)
    filenames.push_back(tmp);
}


/*
 * Reads the feature_vector file and fills fv_map strucute.
 */
static void
load_feature_vectors(const bool VERBOSE,
                     const vector<string> &fv_filenames, 
                     vector<FeatureVector> &fvs) {
  fvs.clear();
  vector<string> labels;
  for (size_t i = 0; i < fv_filenames.size(); ++i) {
    FeatureVector fv;
    vector<string> curr_labels;
    load_features_and_labels(fv_filenames[i], fv, curr_labels);
    if (labels.empty())
      curr_labels.swap(labels);
    else if (labels != curr_labels)
      throw SMITHLABException("incompatible labels: " + fv_filenames[i]);
    
    fvs.push_back(fv);
    if (VERBOSE)
      cerr << '\r' << "loading data: " 
           << percent(i, fv_filenames.size()) << "%\r";
  }
  if (VERBOSE)
    cerr << '\r' << "loading data: 100%" << endl;
}



int
main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    string graph_name("THE_GRAPH");
    string outfile;
    size_t max_degree;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "build m-NNG for feature vectors using "
                           "naive algorithm", "<feature-vectors-file>");
    opt_parse.add_opt("degree", 'd', "max out-degree of graph", 
                      true, max_degree);
    opt_parse.add_opt("name", 'n', "name for the graph", false, graph_name);
    opt_parse.add_opt("out", 'o', "output filename (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc < 3 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string fv_paths_file(leftover_args.front());
    /****************** END COMMAND LINE OPTIONS *****************/
    
    vector<string> fv_filenames;
    load_filenames(fv_paths_file, fv_filenames);
    
    vector<FeatureVector> fvs;
    load_feature_vectors(VERBOSE, fv_filenames, fvs);
    
    if (VERBOSE)
      cerr << "number of vectors: " << fvs.size() << endl;
    
    RegularNearestNeighborGraph naive_graph(graph_name, max_degree);
    
    construct_graph_naively(VERBOSE, max_degree, fvs, naive_graph);
    
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    out << naive_graph << endl;
  }
  catch (const SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
