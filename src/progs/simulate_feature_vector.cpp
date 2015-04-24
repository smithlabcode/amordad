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
#include <iostream>
#include <unistd.h>
#include <fstream>
#include <unordered_map>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "OptionParser.hpp"
#include "smithlab_os.hpp"

#include "FeatureVector.hpp"

using std::string;
using std::vector;
using std::cerr;
using std::endl;


int
main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    string feature_labels_file;
    string feature_label_prefix = "FEAT";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "simulate a feature vector",
			   "<name> <dimension> <outfile>");
    opt_parse.add_opt("features", 'f', "feature labels",
                      false, feature_labels_file);
    opt_parse.add_opt("prefix", 'p', "feature name prefix (default: FEAT)",
                      false, feature_label_prefix);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);

    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
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
    if (leftover_args.size() != 3) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string feature_vector_name(leftover_args.front());
    std::istringstream dim_in(leftover_args[1]);
    size_t n_dimensions = 0;
    if (!(dim_in >> n_dimensions))
      throw SMITHLABException("bad dimension number: " +
			      leftover_args[1]);
    const string outfile(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    vector<string> labels;
    if (!feature_labels_file.empty()) {
      if (VERBOSE)
        cerr << "loading labels" << endl;
      std::ifstream feat_in(feature_labels_file.c_str());
      if (!feat_in)
        throw SMITHLABException("bad file: " + feature_labels_file);
      string tmp_feat;
      while (feat_in >> tmp_feat)
        labels.push_back(tmp_feat);
    }
    else
      for (size_t i = 0; i < n_dimensions; ++i)
	labels.push_back(feature_label_prefix + toa(i));

    if (n_dimensions != labels.size())
      throw SMITHLABException("inconsistent number of labels: " +
			      toa(labels.size()));

    //setup rng
    srand(time(0) + getpid());
    gsl_rng_env_setup();
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, rand());

    if (VERBOSE)
      cerr << "generating coordinates" << endl;
    vector<double> vals(n_dimensions);
    for (size_t i = 0; i < n_dimensions; ++i)
      vals[i] = gsl_ran_gaussian(rng, 1.0);

    const FeatureVector fv(feature_vector_name, vals);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    if (!of) throw SMITHLABException("cannot write to file: " + outfile);
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    out << fv.tostring_with_labels(labels);

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
