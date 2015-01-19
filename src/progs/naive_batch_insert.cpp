/*
 *    Part of AMORDAD software
 *
 *    Copyright (C) 2014 University of Southern California and
 *                       Andrew D. Smith
 *                       Ehsan Behnam
 *
 *    Authors: Andrew D. Smith and Ehsan Behnav
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
#include <algorithm>
#include <climits>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <iterator>
#include <queue>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include "FeatureVector.hpp"

using std::string;
using std::vector;
using std::cerr;
using std::endl;
using std::cout;
using std::pair;
using std::make_pair;

static void
get_filenames(const string &path_file, vector<string> &file_names) {
  std::ifstream in(path_file.c_str());
  if (!in)
    throw SMITHLABException("bad path file: " + path_file);

  file_names.clear();
  string line;
  while (getline(in, line))
    file_names.push_back(line);
}

static void
load_feature_vectors(const bool VERBOSE, const string &fvs_file,
                     vector<FeatureVector> &fvs) {

  vector<string> fv_files;
  get_filenames(fvs_file, fv_files);

  for(size_t i = 0; i < fv_files.size(); ++i) {
    FeatureVector fv;
    std::ifstream in(fv_files[i].c_str());
    if (!in)
      throw SMITHLABException("bad feature vector file: " + fv_files[i]);
    in >> fv;
    fvs.push_back(fv);
    if (VERBOSE)
      cerr << "\rloading feature vectors: "
           << percent(i, fv_files.size()) << "%\r";
  }
  if (VERBOSE)
    cerr << "\rloading feature vectors: 100%" << endl;
}


static bool
validate_file(const string &filename, char open_mode) {
  if (open_mode == 'r')
    return (get_filesize(filename) > 0);
  else { // if (open_mode == 'w') {
    std::ofstream check_out(filename.c_str());
    if (!check_out) return false;
    else return true;
  }
}


int
main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "batch insert a set of "
                           "feature vectors",
                           "<database-file> <insertion-file> <outfile>");
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
    const string database_file(leftover_args.front());
    const string insertions_file(leftover_args[1]);
    const string outfile(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    if (!validate_file(insertions_file, 'r'))
      throw SMITHLABException("bad insertions file: " + insertions_file);

    if (!validate_file(database_file, 'r'))
      throw SMITHLABException("bad database file: " + database_file);

    if (!validate_file(outfile, 'w'))
      throw SMITHLABException("bad output file: " + outfile);

    vector<string> fv_files, insertion_files;
    get_filenames(database_file, fv_files);
    get_filenames(insertions_file, insertion_files);
    ////////////////////////////////////////////////////////////////////////
    ///// STARTING THE INSERTION PROCESS ///////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    // reading insertions to make sure each insertion is a valid feature vector
    vector<FeatureVector> insertions;
    if (VERBOSE)
      cerr << "loading insertions" << endl;
    load_feature_vectors(VERBOSE, insertions_file, insertions);
    if (VERBOSE)
      cerr << "number of insertions: " << insertions.size() << endl;

    ////////////////////////////////////////////////////////////////////////
    ///// NOW WRITE THE OUTPUT /////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    std::ofstream out(outfile.c_str());
    if (!out)
      throw SMITHLABException("bad output file: " + outfile);

    for (size_t i = 0; i < fv_files.size(); ++i) {
      out << fv_files[i] << endl;
    }

    for (size_t i = 0; i < insertion_files.size(); ++i) {
      out << insertion_files[i] << endl;
    }
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
