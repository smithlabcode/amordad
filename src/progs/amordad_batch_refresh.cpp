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
#include <tr1/unordered_set>
#include <cstdio>
#include <iterator>
#include <queue>

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

using std::tr1::unordered_map;
using std::tr1::unordered_set;

typedef LSHAngleHashTable LSHTab;
typedef LSHAngleHashFunction LSHFun;


static void
execute_refresh(const unordered_map<string, FeatureVector> &fvs,
                const unordered_map<string, LSHFun> &replaced_hfs,
                const unordered_map<string, LSHTab> &hts,
                RegularNearestNeighborGraph &g) {

  // iterate over replacing hash tables
  for(unordered_map<string, LSHFun>::const_iterator i(replaced_hfs.begin());
      i != replaced_hfs.end(); ++i) {

    unordered_map<string, LSHTab>::const_iterator ht(hts.find(i->first));
    assert(ht != hts.end());

    // check every bucket to update the graph
    for(BucketMap::const_iterator b(ht->second.begin());
        b != ht->second.end(); ++b) {

      // update every pair of fvs in the same bucket
      for(size_t j = 0; j < b->second.size(); ++j) {
        const FeatureVector fv1(fvs.find(b->second[j])->second);
        for(size_t k = j + 1; k < b->second.size(); ++k) {
          const FeatureVector fv2(fvs.find(b->second[k])->second);
          dist = fv1.compute_angle(fv2);
          g.update_vertex(b->second[j], b->second[k], dist);
          g.update_vertex(b->second[k], b->second[j], dist);
        }
      }
    }
  }
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
get_database(const bool VERBOSE, const string &db_file,
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
    size_t n_hfunc = 1;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "batch refresh an amordad "
                           "database residing on disk",
                           "<config-file>");
    opt_parse.add_opt("hfuncs", 'n', "number of hash functions to replace",
                      false, n_hfunc);
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
    const string database_config_file(leftover_args.back();
    /****************** END COMMAND LINE OPTIONS *****************/


    ////////////////////////////////////////////////////////////////////////
    ////// READING HASH {FUNCTIONS, TABLES, FEATURE_VECTOR} FILES //////////////
    ////////////////////////////////////////////////////////////////////////

    string fv_paths_file, hf_paths_file, ht_paths_file, graph_file;
    read_config_file(database_config_file, fv_paths_file, hf_paths_file,
                     ht_paths_file, graph_file);


    // reading database
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

    RegularNearestNeighborGraph nng;
    g_in >> nng;

    if (VERBOSE)
      cerr << "GRAPH: "
           << "[name=" << nng.get_graph_name() << "]"
           << "[vertices=" << nng.get_vertex_count() << "]"
           << "[edges=" << nng.get_edge_count() << "]"
           << "[max_degree=" << nng.get_maximum_degree() << "]" << endl;

    ////////////////////////////////////////////////////////////////////////
    ////// READING THE HASH FUNCTIONS //////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    // loading hash functions
    unordered_map<string, lshfun> hf_lookup;
    for (size_t i = 0; i < hash_function_files.size(); ++i) {
      std::ifstream hf_in(hash_function_files[i].c_str());
      if (!hf_in)
        throw smithlabexception("bad hash function file: " +
                                hash_function_files[i]);
      lshfun hf;
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
    unordered_map<string, string> ht_to_paths;
    for (size_t i = 0; i < hash_table_files.size(); ++i) {
      std::ifstream ht_in(hash_table_files[i].c_str());
      if (!ht_in)
        throw SMITHLABException("bad hash table file: " +
                                hash_table_files[i]);
      LSHTab ht;
      ht_in >> ht;
      ht_lookup[ht.get_id()] = ht;
      ht_to_paths[ht.get_id()] = hash_table_files[i];
      if (VERBOSE)
        cerr << '\r' << "load hash tables: "
             << percent(i, hash_table_files.size()) << "%\r";
    }
    if (VERBOSE)
      cerr << "load hash tables: 100% (" << ht_lookup.size() << ")" << endl;

    if (VERBOSE)
      cerr << "database loaded" << endl;

    ////////////////////////////////////////////////////////////////////////
    ///// STARTING THE REFRESH PROCESS ///////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    execute_refresh(fv_lookup, replaced_hf_lookup, ht_lookup, nng);


    ////////////////////////////////////////////////////////////////////////
    ///// NOW WRITE THE OUTPUT /////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    std::ofstream out(outfile.c_str());
    if (!out)
      throw SMITHLABException("bad output file: " + outfile);
    for (size_t i = 0; i < queries.size(); ++i) {
      out << queries[i].get_id() << '\t';
      copy(results[i].begin(), results[i].end(),
          std::ostream_iterator<Result>(out, "\t"));
      out << endl;
    }

    if (VERBOSE)
      cerr << comparisons << endl;
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
