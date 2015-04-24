/*
 *    Part of AMORDAD software
 *
 *    Copyright (C) 2014 University of Southern California and
 *                       Andrew D. Smith
 *                       Wenzheng Li
 *
 *    Authors: Andrew D. Smith and Wenzheng Li
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
#include <unordered_set>
#include <cstdio>
#include <iterator>
#include <fstream>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include "RegularNearestNeighborGraph.hpp"

#include "FeatureVector.hpp"
#include "LSHAngleHashTable.hpp"
#include "LSHAngleHashFunction.hpp"

using std::string;
using std::vector;
using std::cerr;
using std::endl;
using std::cout;
using std::pair;
using std::make_pair;

using std::unordered_map;
using std::unordered_set;

typedef LSHAngleHashTable LSHTab;
typedef LSHAngleHashFunction LSHFun;

size_t comparisons = 0;

static void
execute_deletion(unordered_map<string, FeatureVector> &fvs,
                 const unordered_map<string, LSHFun> &hfs,
                 const FeatureVector &query,
                 const size_t n_neighbors,
                 unordered_map<string, LSHTab> &hts,
                 RegularNearestNeighborGraph &g) {
  
  
  // iterate over hash tables
  for (unordered_map<string, LSHTab>::iterator i(hts.begin());
       i != hts.end(); ++i) {

    unordered_map<string, LSHFun>::const_iterator hf(hfs.find(i->first));
    assert(hf != hfs.end());

    // hash the query
    const size_t bucket_number = hf->second(query);
    ++comparisons;
    unordered_map<size_t, vector<string> >::const_iterator
      bucket = i->second.find(bucket_number);

    // delete the query from each hash table
    if (bucket != i->second.end())
      i->second.remove(query, bucket_number);
  }

  g.remove_vertex(query.get_id());

  // delete the query from the feature vectors map
  fvs.erase(query.get_id());
}


/*
 * Config file has this structure:
 * feature_vectors_file: path to file containing <fv_id, fv_filename>
 * hash_functions_file: path to file containing filenames of hash functions
 * hash_tables_file: path to file containing filenames of all hash tables
 * graph_file: the filename for the m-NNG
 */
static void
read_config_file(const string &config_file,
                 string &feature_vectors_file,
                 string &hash_functions_file,
                 string &hash_tables_file,
                 string &graph_file) {

  std::ifstream in(config_file.c_str());
  if (!in)
    throw SMITHLABException("bad configuration file: " + config_file);

  getline(in, feature_vectors_file);
  getline(in, hash_functions_file);
  getline(in, hash_tables_file);
  getline(in, graph_file);
}


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
get_deletions(const string &deletions_file, vector<FeatureVector> &deletions) {

  vector<string> deletion_files;
  get_filenames(deletions_file, deletion_files);

  for(size_t i = 0; i < deletion_files.size(); ++i) {
    FeatureVector fv;
    std::ifstream in(deletion_files[i].c_str());
    if (!in)
      throw SMITHLABException("bad feature vector file: " + deletion_files[i]);
    in >> fv;
    deletions.push_back(fv);
  }
}


static void
get_database(const bool VERBOSE,
             const string &db_file, 
             unordered_map<string, FeatureVector> &db) {
  
  vector<string> fv_files;
  get_filenames(db_file, fv_files);
  for (size_t i = 0; i < fv_files.size(); ++i)
    fv_files[i] = fv_files[i].substr(fv_files[i].find(' ') + 1);

  for(size_t i = 0; i < fv_files.size(); ++i) {
    FeatureVector fv;
    std::ifstream in(fv_files[i].c_str());
    if (!in)
      throw SMITHLABException("bad feature vector file: " + fv_files[i]);
    in >> fv;
    db[fv.get_id()] = fv;

    if (VERBOSE)
      cerr << '\r' << "loading feature vectors: "
           << percent(i, fv_files.size()) << "%\r";
  }

  if (VERBOSE)
    cerr << '\r' << "loading feature vectors: 100% ("
         << fv_files.size() << ")" << endl;
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
    OptionParser opt_parse(strip_path(argv[0]), "batch deletion an amordad "
                           "database residing on disk "
                           "and add '.up' to all updated file names",
                           "<config-file> <deletion-dir>");
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
    const string database_config_file(leftover_args.front());
    const string deletions_file(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/
    
    if (!validate_file(deletions_file, 'r'))
      throw SMITHLABException("bad deletions file: " + deletions_file);
    
    ////////////////////////////////////////////////////////////////////////
    ////// READING HASH {FUNCTIONS, TABLES, FEATURE_VECTOR} FILES //////////////
    ////////////////////////////////////////////////////////////////////////

    string fv_paths_file, hf_paths_file, ht_paths_file, graph_file;
    read_config_file(database_config_file, fv_paths_file, hf_paths_file,
                     ht_paths_file, graph_file);
    
    // loading feature vectors
    unordered_map<string, FeatureVector> fv_lookup;
    get_database(VERBOSE, fv_paths_file, fv_lookup);
    
    vector<string> hash_function_files, hash_table_files;
    get_filenames(hf_paths_file, hash_function_files);
    get_filenames(ht_paths_file, hash_table_files);

    ////////////////////////////////////////////////////////////////////////
    ////// READING THE GRAPH ///////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    std::ifstream g_in(graph_file.c_str());
    if (!g_in)
      throw SMITHLABException("cannot load graph: " + graph_file);

    if (VERBOSE)
      cerr << "loading graph" << endl;

    RegularNearestNeighborGraph nng;
    g_in >> nng;

    /// THIS IS THE DEGREE OF THE GRAPH AND SHOULD BE OBTAINED
    /// FROM THE DATABASE ITSELF, AS IT IS ENCODED IN THE GRAPH FILE
    size_t n_neighbors = nng.get_maximum_degree();

    if (VERBOSE)
      cerr << "GRAPH: "
           << "[name=" << nng.get_graph_name() << "]"
           << "[vertices=" << nng.get_vertex_count() << "]"
           << "[edges=" << nng.get_edge_count() << "]"
           << "[max_degree=" << nng.get_maximum_degree() << "]" << endl;

    ////////////////////////////////////////////////////////////////////////
    ////// READING THE HASH FUNCTIONS //////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    //loading hash functions
    unordered_map<string, LSHFun> hf_lookup;
    for (size_t i = 0; i < hash_function_files.size(); ++i) {
      std::ifstream hf_in(hash_function_files[i].c_str());
      if (!hf_in)
        throw SMITHLABException("bad hash function file: " +
                                hash_function_files[i]);
      LSHFun hf;
      hf_in >> hf;
      hf_lookup[hf.get_id()] = hf;
      if (VERBOSE)
        cerr << '\r' << "load hash functions: "
             << percent(i, hash_function_files.size()) << "%\r";
    }
    if (VERBOSE)
      cerr << "load hash functions: 100% (" << hf_lookup.size() << ")" << endl;

    ////////////////////////////////////////////////////////////////////////
    ////// READING THE HASH TABLES /////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    unordered_map<string, LSHTab> ht_lookup;
    unordered_map<string, string> id_to_path_ht;
    for (size_t i = 0; i < hash_table_files.size(); ++i) {
      std::ifstream ht_in(hash_table_files[i].c_str());
      if (!ht_in)
        throw SMITHLABException("bad hash table file: " +
                                hash_table_files[i]);
      LSHTab ht;
      ht_in >> ht;
      ht_lookup[ht.get_id()] = ht;
      id_to_path_ht[ht.get_id()] = hash_table_files[i];
      if (VERBOSE)
        cerr << '\r' << "load hash tables: "
             << percent(i, hash_table_files.size()) << "%\r";
    }
    if (VERBOSE)
      cerr << "load hash tables: 100% (" << ht_lookup.size() << ")" << endl;

    if (VERBOSE)
      cerr << "database loaded" << endl;

    ////////////////////////////////////////////////////////////////////////
    ///// STARTING THE DELETIONS PROCESS ///////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    // reading deletions
    vector<FeatureVector> deletions;
    get_deletions(deletions_file, deletions);
    if (VERBOSE)
      cerr << "number of deletions: " << deletions.size() << endl;

    for (size_t i = 0; i < deletions.size(); ++i) {
      execute_deletion(fv_lookup, hf_lookup, deletions[i], n_neighbors,
                        ht_lookup, nng);
      if (VERBOSE)
        cerr << '\r' << "processing deletions: "
             << percent(i, deletions.size()) << "%\r";
    }
    if (VERBOSE)
      cerr << '\r' << "processing deletions: 100% ("
           << deletions.size() << ")" << endl;
    
    ////////////////////////////////////////////////////////////////////////
    ///// NOW WRITE THE DATABASE BACK TO DISK //////////////////////////////
    ////////////////////////////////////////////////////////////////////////


    // writing graph back to the disk
    if (VERBOSE)
      cerr << "UPDATED GRAPH: "
           << "[name=" << nng.get_graph_name() << "]"
           << "[vertices=" << nng.get_vertex_count() << "]"
           << "[edges=" << nng.get_edge_count() << "]"
           << "[max_degree=" << nng.get_maximum_degree() << "]" << endl;

    std::ofstream of_g((strip_path(graph_file) + ".up").c_str());
    of_g << nng << endl;
    if (VERBOSE)
      cerr << "writing graph back: 100%" << endl;

    // writing every hash table back
    size_t count = 0;
    for (unordered_map<string,LSHTab>::const_iterator i(ht_lookup.begin());
         i != ht_lookup.end(); ++i) {
      string path = id_to_path_ht[i->second.get_id()];
      std::ofstream of_ht((strip_path(path) + ".up").c_str());
      of_ht << i->second << endl;

      count++;
      if (VERBOSE)
        cerr << '\r' << "writing hashtables back: "
             << percent(count, ht_lookup.size()) << "%\r";
    }
    if (VERBOSE)
      cerr << '\r' << "writing hashtables back: 100% ("
        << ht_lookup.size() << ")" << endl;
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
