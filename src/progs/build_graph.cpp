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
#include <algorithm>
#include <unordered_set>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include "RegularNearestNeighborGraph.hpp"
#include "FeatureVector.hpp"
#include "LSHAngleHashTable.hpp"

using std::string;
using std::vector;
using std::cerr;
using std::endl;
using std::cout;
using std::ifstream;
using std::unordered_map;
using std::unordered_set;


typedef unordered_map<string, unordered_set<string> > ComparedLookup;
typedef unordered_map<string, FeatureVector> FeatVecLookup;


/* read all feature vectors (fvs)
 */
static void
load_feature_vectors(const bool VERBOSE,
                     const string &feat_vecs_file, FeatVecLookup &fvs) {
  
  ifstream fv_filenames_in(feat_vecs_file.c_str());
  if (!fv_filenames_in)
    throw SMITHLABException("problem reading: " + feat_vecs_file);
  vector<string> filenames;
  string filename;
  while (fv_filenames_in >> filename)
    filenames.push_back(filename);
  
  fvs.clear();
  for (size_t i = 0; i < filenames.size(); ++i) {
    std::ifstream in(filenames[i].c_str());
    if (!in)
      throw SMITHLABException("problem reading: " + filenames[i]);
    
    FeatureVector fv;
    in >> fv;
    fvs[fv.get_id()] = fv;
    if (VERBOSE)
      cerr << '\r' << "loading data: " 
           << percent(i, filenames.size()) << "%\r";
  }
  if (VERBOSE)
    cerr << '\r' << "loading data: 100%" << endl;
}


/* reads all hash tables (hts)
 */
static void
load_hash_tables(const string &hash_tables_file, 
                 vector<LSHAngleHashTable> &hts) {
  
  ifstream ht_filenames_in(hash_tables_file.c_str());
  if (!ht_filenames_in)
    throw SMITHLABException("problem reading: " + hash_tables_file);
  
  hts.clear();
  string filename;
  while (ht_filenames_in >> filename) {
    ifstream in(filename.c_str());
    if (!in)
      throw SMITHLABException("problem reading: " + filename);
    
    LSHAngleHashTable ht;
    in >> ht;
    hts.push_back(ht);
  }
}


static void
add_relations_from_bucket(const vector<string> &bucket, 
                          const FeatVecLookup &featvecs,
                          ComparedLookup &compared, 
                          RegularNearestNeighborGraph &nng) {
  
  // iterate over bucket
  for (size_t i = 0; i < bucket.size(); ++i) {
    FeatVecLookup::const_iterator ii(featvecs.find(bucket[i]));
    assert(ii != featvecs.end());
    
    // iterate over other members of bucket
    for (size_t j = i + 1; j < bucket.size(); ++j) {
      FeatVecLookup::const_iterator jj(featvecs.find(bucket[j]));
      assert(jj != featvecs.end());
      
      // check if previously compared
      ComparedLookup::const_iterator c(compared.find(bucket[i]));
      assert(c != compared.end());
      if (c->second.find(bucket[j]) == c->second.end()) {
        
        // compare and update graph
        const double w = ii->second.compute_angle(jj->second);
        
        nng.update_vertex(bucket[j], bucket[i], w);
        nng.update_vertex(bucket[i], bucket[j], w);
        
        // add them both to "compared"
        compared[bucket[j]].insert(bucket[i]);
        compared[bucket[i]].insert(bucket[j]);
      }
    }
  }
}


int
main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;

    string graph_name("THE_GRAPH");
    string id;
    string outfile;
    size_t max_degree = 1;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "generate m-NNG graph", "<feat-vecs> <hash-tables>");
    opt_parse.add_opt("deg", 'd', "max out degree of graph", true, max_degree);
    opt_parse.add_opt("name", 'n', "name for the graph", false, graph_name);
    opt_parse.add_opt("out", 'o', "output file (default: stdout)", 
                      true, outfile);
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
    if (leftover_args.size() != 2) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string feat_vecs_filename(leftover_args.front());
    const string hash_tables_filename(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    // first load the feature vectors
    FeatVecLookup featvecs;
    load_feature_vectors(VERBOSE, feat_vecs_filename, featvecs);
    if (VERBOSE)
      cerr << "number of feature vectors: " << featvecs.size() << endl;
    
    // next load the hash tables
    vector<LSHAngleHashTable> hts;
    load_hash_tables(hash_tables_filename, hts);
    if (VERBOSE)
      cerr << "number of hash tables: " << hts.size() << endl;
    
    // now intialize the graph
    RegularNearestNeighborGraph nng(graph_name, max_degree);
    for (FeatVecLookup::const_iterator i(featvecs.begin()); 
         i != featvecs.end(); ++i)
      nng.add_vertex(i->first);
    
    // Initializing "compared": 
    // compared<FV_ID, {Chk_ID1, Chk_ID2, ...}  for each feature
    // vector "FV_ID", keeps the set of feature vectors that have been
    // already compared in searching for nearest feature vector to
    // "FV_ID".
    ComparedLookup compared;
    for (FeatVecLookup::const_iterator i(featvecs.begin()); 
         i != featvecs.end(); ++i)
      compared.insert(make_pair(i->first, unordered_set<string>()));
    
    // iterate over hash tables
    for (size_t i = 0; i < hts.size(); ++i) {
      if (VERBOSE)
        cerr << '\r' << "hashing: " << percent(i, hts.size()) << "%\r";
      // iterate over buckets
      for (BucketMap::const_iterator j(hts[i].begin()); j != hts[i].end(); ++j)
        add_relations_from_bucket(j->second, featvecs, compared, nng);
    }
    if (VERBOSE)
      cerr << '\r' << "hashing: 100%" << endl;
    
    //Finally: write the graph on the outfile.
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    if (!of) throw SMITHLABException("cannot write to file: " + outfile);
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    
    out << nng << endl;
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
