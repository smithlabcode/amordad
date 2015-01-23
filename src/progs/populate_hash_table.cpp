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

#include "LSHAngleHashFunction.hpp"
#include "LSHEuclideanHashFunction.hpp"
#include "FeatureVector.hpp"
#include "LSHAngleHashTable.hpp"

using std::string;
using std::vector;
using std::cerr;
using std::endl;
using std::tr1::unordered_map;


int
main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    string outfile;
    string features_file;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "build angle-lsh hash table",
                           "<vectors-path-file> <hash-fun-file>");
    opt_parse.add_opt("features", 'f', "feature set file (for validation)",
                      false, features_file);
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
    if (leftover_args.size() != 2) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string vectors_path_file(leftover_args.front());
    const string hash_function_file(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    vector<string> features;
    if (!features_file.empty()) {
      if (VERBOSE)
        cerr << "loading features" << endl;
      std::ifstream feat_in(features_file.c_str());
      if (!feat_in)
        throw SMITHLABException("bad file: " + features_file);
      string tmp_feat;
      while (feat_in >> tmp_feat)
        features.push_back(tmp_feat);
    }
    


    if (VERBOSE)
      cerr << "loading hash function" << endl;
    // LOAD THE HASH FUNCTION
    std::ifstream hash_fun_in(hash_function_file.c_str());
    if (!hash_fun_in)
      throw SMITHLABException("cannot open: " + hash_function_file);
    // LSHAngleHashFunction hash_fun;
    LSHEuclideanHashFunction hash_fun;
    hash_fun_in >> hash_fun;
    


    if (VERBOSE)
      cerr << "loading hash table" << endl;
    // INITIALIZE THE HASH TABLE
    LSHAngleHashTable hash_table(hash_fun.get_id());


    
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
      cerr << "hashing feature vectors" << endl;
    // ITERATE OVER EACH FEATURE VECTOR FILE AND HASH IT
    for (size_t i = 0; i < feat_vec_filenames.size(); ++i) {
      
      FeatureVector fv;
      if (!features_file.empty()) {
        vector<string> labels;
        load_features_and_labels(feat_vec_filenames[i], fv, labels);
        if (features != labels)
          throw SMITHLABException("inconsistent features: " + 
                                  feat_vec_filenames[i] + "/" + features_file);
      }
      else {
        std::ifstream feat_vec_in(feat_vec_filenames[i].c_str());
        if (!feat_vec_in)
          throw SMITHLABException("bad file: " + feat_vec_filenames[i]);
        feat_vec_in >> fv;
      }
      hash_table.insert(fv, hash_fun(fv));
      
      if (VERBOSE)
        cerr << '\r' << percent(i, feat_vec_filenames.size()) << "%\r";
    }

    // write the hash table to disk
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    if (!of) throw SMITHLABException("cannot write to file: " + outfile);
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    out << hash_table << endl;
    
    if (VERBOSE)
      cerr << "hash table size:" << '\t' 
           << hash_table.size() << endl
           << "max bucket load:" << '\t'
           << hash_table.max_bucket_load() << endl;
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
