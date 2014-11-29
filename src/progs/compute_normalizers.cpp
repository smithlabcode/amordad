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

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_statistics_double.h>

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
remove_upper_tail(const size_t tail_size, vector<double> &vals) {
  assert(tail_size < vals.size());
  vals.erase(vals.end() - tail_size, vals.end());
}


static void
remove_lower_tail(const size_t tail_size, vector<double> &vals) {
  assert(tail_size < vals.size());
  vals.erase(copy(vals.begin() + tail_size, vals.end(), vals.begin()));
}


int
main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    string outfile;
    string features_file;

    double tail_fraction = 0.10;
    bool use_upper_tail = false;
    bool use_lower_tail = false;
    bool use_median_values = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "compute mean and sd for each feature",
                           "<vectors-path-file>");
    opt_parse.add_opt("out", 'o', "output file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("upper", 'u', "trim upper tail", false, use_upper_tail);
    opt_parse.add_opt("lower", 'l', "trim lower tail", false, use_lower_tail);
    opt_parse.add_opt("med", 'M', "use medians", false, use_median_values);
    opt_parse.add_opt("tail", 't', "size of tail to remove", false, tail_fraction);
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

    ////////////////////////////////////////////////////////////
    //// GET FEATURE VECTOR FILE NAMES
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

    ////////////////////////////////////////////////////////////
    //// READ IN THE FEATURE VECTORS
    vector<vector<double> > vals;
    vector<string> labels;
    size_t n_values = feat_vec_filenames.size();
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
        cerr << '\r' << "loading feature vectors: "
             << percent(i, n_values) << "%\r";
    }
    if (VERBOSE)
      cerr << '\r' << "loading feature vectors: 100%" << endl;

    ////////////////////////////////////////////////////////////
    //// REMOVE OUTLIER VALUES
    const size_t tail_size = tail_fraction*static_cast<double>(n_values);
    if (VERBOSE)
      cerr << "tail size: " << tail_size << endl;
    for (size_t i = 0; i < vals.size(); ++i) {
      sort(vals[i].begin(), vals[i].end());
      if (use_upper_tail) {
        remove_upper_tail(tail_size, vals[i]);
        n_values -= tail_size;
      }
      if (use_lower_tail) {
        remove_lower_tail(tail_size, vals[i]);
        n_values -= tail_size;
      }
      if (VERBOSE)
        cerr << '\r' << "trimming feature vectors: "
             << percent(i, vals.size()) << "%\r";
    }
    if (VERBOSE)
      cerr << '\r' << "trimming feature vectors: 100%" << endl;

    ////////////////////////////////////////////////////////////
    //// COMPUTE MEANS AND VARIANCES
    vector<double> means(vals.size(), 0.0);
    vector<double> medians(vals.size(), 0.0);
    for (size_t i = 0; i < vals.size(); ++i) {
      medians[i] = gsl_stats_median_from_sorted_data(&vals[i].front(), 1, n_values);
      means[i] = gsl_stats_mean(&vals[i].front(), 1, n_values);
      if (VERBOSE)
        cerr << '\r' << "computing means/medians: " << percent(i, vals.size()) << "%\r";
    }
    if (VERBOSE)
      cerr << '\r' << "computing means/medians: 100%" << endl;

    vector<double> sds(vals.size(), 0.0);
    for (size_t i = 0; i < vals.size(); ++i) {
      sds[i] = gsl_stats_sd_m(&vals[i].front(), 1, n_values, means[i]);
      if (VERBOSE)
        cerr << '\r' << "computing variances: "
             << percent(i, vals.size()) << "%\r";
    }
    if (VERBOSE)
      cerr << '\r' << "computing variances: 100%" << endl;

    ////////////////////////////////////////////////////////////
    //// WRITE THE NORMALIZERS TO DISK
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    if (!of) throw SMITHLABException("cannot write to file: " + outfile);
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    for (size_t i = 0; i < means.size(); ++i)
      out << labels[i] << '\t'
          << (use_median_values ? medians[i] : means[i]) << '\t'
          << sds[i] << endl;
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
