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
#include <fstream>

#include <gsl/gsl_statistics_double.h>

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

    string outfile_suffix = ".norm";

    bool VERBOSE = false;
    string features_file;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), 
                           "normalizes a feature vector",
                           "<infile> <outfile>");
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
    if (leftover_args.size() != 2) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string input_file(leftover_args.front());
    const string outfile(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/
    
    FeatureVector fv;
    vector<string> labels;
    load_features_and_labels(input_file, fv, labels);
    
    const double mean = gsl_stats_mean(&fv[0], 1, fv.size());
    const double sd = gsl_stats_sd_m(&fv[0], 1, fv.size(), mean);
    
    transform(fv.begin(), fv.end(), fv.begin(), 
              bind2nd(std::minus<double>(), mean));
    transform(fv.begin(), fv.end(), fv.begin(), 
              bind2nd(std::divides<double>(), sd));
    
    // write the hash table to disk
    std::ofstream out(outfile.c_str());
    if (!out)
      throw SMITHLABException("bad output file: " + outfile);
    out << fv.tostring_with_labels(labels) << endl;
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
