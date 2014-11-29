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
#include <tr1/unordered_map>

#include "OptionParser.hpp"
#include "smithlab_os.hpp"

#include "FeatureVector.hpp"

using std::string;
using std::vector;
using std::cerr;
using std::endl;
using std::tr1::unordered_map;


static void
load_normalizers(const string &normalizers_file,
                 vector<string> &labels,
                 vector<double> &means, vector<double> &sds) {
  
  std::ifstream in(normalizers_file.c_str());
  if (!in)
    throw SMITHLABException("cannot open file: " + normalizers_file);
  
  string l;
  double m = 0.0, s = 0.0;
  while (in >> l >> m >> s) {
    labels.push_back(l);
    means.push_back(m);
    sds.push_back(s);
  }
}


int
main(int argc, const char **argv) {

  try {

    string outfile_suffix = ".norm";

    bool VERBOSE = false;
    string outfile;
    string features_file;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), 
                           "normalizes feature vectors",
                           "<normalizer-file> <vectors-path-file>");
    opt_parse.add_opt("out", 'o', "output file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("suff", 's', "output file suffic (default: norm)",
                      false, outfile_suffix);
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
    const string normalizers_file(leftover_args.front());
    const string vectors_path_file(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE)
      cerr << "loading normalizers" << endl;
    vector<double> means, sds;
    vector<string> labels;
    load_normalizers(normalizers_file, labels, means, sds);
    
    if (VERBOSE)
      cerr << "extracting feature vector paths" << endl;
    std::ifstream paths_in(vectors_path_file.c_str());
    if (!paths_in)
      throw SMITHLABException("bad feature vectors locations: " +
                              vectors_path_file);
    vector<string> feat_vec_filenames;
    string feat_vec_file;
    while (paths_in >> feat_vec_file)
      feat_vec_filenames.push_back(feat_vec_file);
    
    if (VERBOSE)
      cerr << "processing feature vectors" << endl;
    
    for (size_t i = 0; i < feat_vec_filenames.size(); ++i) {
      
      FeatureVector fv;
      vector<string> curr_labels;
      load_features_and_labels(feat_vec_filenames[i], fv, curr_labels);
      
      if (labels != curr_labels)
        throw SMITHLABException("inconsistent labels: " +
                                feat_vec_filenames[0] + "\t" + 
                                feat_vec_filenames[i]);
      
      for (size_t j = 0; j < means.size(); ++j)
        fv[j] = (fv[j] - means[j])/sds[j];
      
      if (VERBOSE)
        cerr << '\r' << percent(i, feat_vec_filenames.size()) << "%\r";
      
      // write the hash table to disk
      std::ofstream out(string(feat_vec_filenames[i] + outfile_suffix).c_str());
      out << fv.tostring_with_labels(labels) << endl;
    }
    if (VERBOSE)
      cerr << '\r' << "100%" << endl;
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
