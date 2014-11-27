/*
 *    Part of SMITHLAB software
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
#include <numeric>

#include "OptionParser.hpp"
#include "smithlab_os.hpp"

#include "FeatureVector.hpp"

using std::string;
using std::vector;
using std::cerr;
using std::endl;
using std::tr1::unordered_map;
using std::accumulate;
using std::transform;


static void
remove_tails(const double tail_fraction, vector<double> &vals) {
  const size_t offset = tail_fraction*vals.size();
  const size_t number_to_copy = vals.size() - 2*offset;
  copy(vals.begin() + offset, vals.begin() + offset + number_to_copy,
       vals.begin());
  vals.erase(vals.begin() + number_to_copy, vals.end());  
}



int
main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    string outfile;
    string features_file;

    double tail_fraction = 0.01;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), 
                           "compute mean and sd for each feature",
                           "<vectors-path-file>");
    opt_parse.add_opt("out", 'o', "output file (default: stdout)",
                      false, outfile);
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
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string vectors_path_file(leftover_args.front());
    /****************** END COMMAND LINE OPTIONS *****************/
    
    if (VERBOSE)
      cerr << "extracting feature vector paths" << endl;
    // GET FEATURE VECTOR FILE NAMES
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

    vector<vector<double> > vals;
    vector<string> labels;
    
    // ITERATE OVER EACH FEATURE VECTOR FILE AND HASH IT
    for (size_t i = 0; i < feat_vec_filenames.size(); ++i) {
      
      FeatureVector fv;
      vector<string> curr_labels;
      load_features_and_labels(feat_vec_filenames[i], fv, curr_labels);
      if (i == 0) {
        vals = vector<vector<double> >(fv.size(), vector<double>());
        labels.swap(curr_labels);
      }
      else if (labels != curr_labels) 
        throw SMITHLABException("inconsistent labels: " +
                                feat_vec_filenames[0] + "\t" +
                                feat_vec_filenames[i]);
      for (size_t j = 0; j < vals.size(); ++j)
        vals[j].push_back(fv[j]);
      
      if (VERBOSE)
        cerr << '\r' << percent(i, feat_vec_filenames.size()) << "%\r";
    }
    if (VERBOSE)
      cerr << '\r' << "100%" << endl;

    // shift the vals
    for (size_t i = 0; i < vals.size(); ++i)
      remove_tails(tail_fraction, vals[i]);
    
    vector<double> means(vals.size()), sds(vals.size());
    for (size_t j = 0; j < vals.size(); ++j)
      means[j] = accumulate(vals[j].begin(), vals[j].end(), 0.0)/vals[j].size();
    
    for (size_t j = 0; j < vals.size(); ++j) {
      transform(vals[j].begin(), vals[j].end(), vals[j].begin(),
                bind2nd(std::minus<double>(), means[j]));    
      transform(vals[j].begin(), vals[j].end(), vals[j].begin(), 
                std::ptr_fun(::fabs));
      sds[j] = accumulate(vals[j].begin(), vals[j].end(), 0.0)/vals[j].size();
    }
    
    // write the hash table to disk
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    
    for (size_t i = 0; i < means.size(); ++i)
      out << labels[i] << '\t' << means[i] << '\t' << sds[i] << endl;

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
