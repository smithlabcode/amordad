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
#include <cstdlib>
#include <ctime>

#include "OptionParser.hpp"
#include "smithlab_os.hpp"

#include "LSHAngleHashFunction.hpp"

using std::string;
using std::vector;
using std::cerr;
using std::cout;
using std::endl;



int
main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    size_t n_bits = 0;
    size_t n_features = 0;
    size_t rng_key = 0;
    string id;
    string feature_set_id;
    string outfile;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), 
                           "generate lsh-angle hash function", "");
    opt_parse.add_opt("id", 'i', "id of hash function", true, id);
    opt_parse.add_opt("fs", 'f', "feature set id", true, feature_set_id);
    opt_parse.add_opt("nfeat", 'n', "number of features", true, n_features);
    opt_parse.add_opt("bits", 'b', "bits in hash value", true, n_bits);
    opt_parse.add_opt("out", 'o', "output file (default: stdout)", false, outfile);
    opt_parse.add_opt("key", 'k', "random generator key (for debugging)", 
                      false, rng_key);
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
    if (leftover_args.size() != 0) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    /****************** END COMMAND LINE OPTIONS *****************/
    
    srand((rng_key != 0) ? rng_key : time(0) + getpid());
    
    const LSHAngleHashFunction hash_function(id, feature_set_id, 
                                             n_features, n_bits);
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    if (!of) throw SMITHLABException("cannot write to file: " + outfile);
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    
    out << hash_function << endl;
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
