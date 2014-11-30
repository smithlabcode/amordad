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

size_t comparisons = 0;

struct Result {
  Result(const string &i, const double v) : id(i), val(v) {}
  Result() : val(std::numeric_limits<double>::max()) {}
  bool operator>(const Result &other) const {return val > other.val;}
  string id;
  double val;
};


std::ostream &
operator<<(std::ostream &os, const Result &r) {
  return os << r.id << '\t' << r.val;
}


// static void
// get_feature_vector(const string &id,
//                    const unordered_map<string, string> &fv_path_lookup,
//                    FeatureVector &fv) {
//   unordered_map<string, string>::const_iterator i(fv_path_lookup.find(id));
//   if (i == fv_path_lookup.end())
//     throw SMITHLABException("cannot find file for: " + id);

//   std::ifstream in(i->second.c_str());
//   if (!in)
//     throw SMITHLABException("bad feature vector file: " + i->second);
//   in >> fv;
// }

// static void
// get_feature_vector_paths_lookup(const string &fv_paths_file,
//                                 unordered_map<string, string> &fv_paths_lookup) {
//   std::ifstream in(fv_paths_file.c_str());
//   if (!in)
//     throw SMITHLABException("bad feature vector paths file: " + fv_paths_file);

//   string fv_id, fv_path;
//   while (in >> fv_id >> fv_path)
//     fv_paths_lookup[fv_id] = fv_path;
// }

// static void
// evaluate_candidates(const unordered_map<string, string>& fv_file_lookup,
//                     const FeatureVector &query,
//                     const size_t n_neighbors,
//                     const double max_proximity_radius,
//                     const unordered_set<string> &candidates,
//                     vector<Result> &results) {

//   std::priority_queue<Result, vector<Result>, std::greater<Result> > pq;
//   double current_dist_cutoff = max_proximity_radius;
//   for (unordered_set<string>::const_iterator i(candidates.begin());
//        i != candidates.end(); ++i) {
//     FeatureVector fv;
//     get_feature_vector(*i, fv_file_lookup, fv);
//     const double dist = query.compute_angle(fv);
//     if (dist < current_dist_cutoff) {
//       if (pq.size() == n_neighbors) {
//         pq.pop();
//         current_dist_cutoff = dist;
//       }
//       pq.push(Result(*i, dist));
//     }
//   }

//   results.clear();
//   while (!pq.empty()) {
//     results.push_back(pq.top());
//     pq.pop();
//   }
//   reverse(results.begin(), results.end());
// }


static void
evaluate_candidates(const unordered_map<string, FeatureVector> &fvs,
                    const FeatureVector &query,
                    const size_t n_neighbors,
                    const double max_proximity_radius,
                    const unordered_set<string> &candidates,
                    vector<Result> &results) {

  std::priority_queue<Result, vector<Result>, std::greater<Result> > pq;
  double current_dist_cutoff = max_proximity_radius;
  for (unordered_set<string>::const_iterator i(candidates.begin());
       i != candidates.end(); ++i) {
    const FeatureVector fv(fvs.find(*i)->second);
    const double dist = query.compute_angle(fv);
    ++comparisons;
    if (dist < current_dist_cutoff) {
      if (pq.size() == n_neighbors) {
        pq.pop();
        current_dist_cutoff = dist;
      }
      pq.push(Result(*i, dist));
    }
  }

  results.clear();
  while (!pq.empty()) {
    results.push_back(pq.top());
    pq.pop();
  }
  reverse(results.begin(), results.end());
}



static void
execute_query(const unordered_map<string, FeatureVector> &fvs,
              const unordered_map<string, LSHFun> &hfs,
              const unordered_map<string, LSHTab> &hts,
              const RegularNearestNeighborGraph &g,
              const FeatureVector &query,
              const size_t n_neighbors,
              const double max_proximity_radius,
              vector<Result> &results) {

  unordered_set<string> candidates;

  // iterate over hash tables
  for (unordered_map<string, LSHTab>::const_iterator i(hts.begin());
       i != hts.end(); ++i) {

    unordered_map<string, LSHFun>::const_iterator hf(hfs.find(i->first));
    assert(hf != hfs.end());

    // hash the query
    const size_t bucket_number = hf->second(query);
    ++comparisons;
    unordered_map<size_t, vector<string> >::const_iterator
      bucket = i->second.find(bucket_number);

    if (bucket != i->second.end())
      candidates.insert(bucket->second.begin(), bucket->second.end());
  }

  // gather neighbors of candidates
  unordered_set<string> candidates_from_graph;
  for (unordered_set<string>::const_iterator i(candidates.begin());
       i != candidates.end(); ++i) {
    vector<string> neighbors;
    vector<double> neighbor_dists;
    g.get_neighbors(*i, neighbors, neighbor_dists);
    candidates_from_graph.insert(neighbors.begin(), neighbors.end());
  }

  candidates.insert(candidates_from_graph.begin(), candidates_from_graph.end());

  // evaluate_candidates(fv_files_lookup, query, n_neighbors,
  //                     max_proximity_radius, candidates, results);
  evaluate_candidates(fvs, query, n_neighbors,
                      max_proximity_radius, candidates, results);
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


/*
 * See what is inside query_dir and loads all the queries
 */
static void
get_queries(const string &queries_file, vector<FeatureVector> &queries) {

  vector<string> query_files;
  get_filenames(queries_file, query_files);

  for(size_t i = 0; i < query_files.size(); ++i) {
    FeatureVector fv;
    std::ifstream in(query_files[i].c_str());
    if (!in)
      throw SMITHLABException("bad feature vector file: " + query_files[i]);
    in >> fv;
    queries.push_back(fv);
  }
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
    size_t n_neighbors = 1;
    double max_proximity_radius = 0.75;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "batch query an amordad "
                           "database residing on disk",
                           "<config-file> <query-dir> <outfile>");
    opt_parse.add_opt("neighbors", 'n', "number of nearest neighbors to report",
                      false, n_neighbors);
    opt_parse.add_opt("mpr", 'r', "maximum proximity radius",
                      false, max_proximity_radius);
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
    const string database_config_file(leftover_args.front());
    const string queries_file(leftover_args[1]);
    const string outfile(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    if (!validate_file(queries_file, 'r'))
      throw SMITHLABException("bad queries file: " + queries_file);

    if (!validate_file(outfile, 'w'))
      throw SMITHLABException("bad output file: " + outfile);

    ////////////////////////////////////////////////////////////////////////
    ////// READING HASH {FUNCTIONS, TABLES, FEATURE_VECTOR} FILES //////////////
    ////////////////////////////////////////////////////////////////////////

    string fv_paths_file, hf_paths_file, ht_paths_file, graph_file;
    read_config_file(database_config_file, fv_paths_file, hf_paths_file,
                     ht_paths_file, graph_file);

    // unordered_map<string, string> fv_path_lookup;
    // get_feature_vector_paths_lookup(fv_paths_file, fv_path_lookup);

    // reading queries
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
    for (size_t i = 0; i < hash_table_files.size(); ++i) {
      std::ifstream ht_in(hash_table_files[i].c_str());
      if (!ht_in)
        throw SMITHLABException("bad hash table file: " +
                                hash_table_files[i]);
      LSHTab ht;
      ht_in >> ht;
      ht_lookup[ht.get_id()] = ht;
      if (VERBOSE)
        cerr << '\r' << "load hash tables: "
             << percent(i, hash_table_files.size()) << "%\r";
    }
    if (VERBOSE)
      cerr << "load hash tables: 100% (" << ht_lookup.size() << ")" << endl;

    if (VERBOSE)
      cerr << "database loaded" << endl;

    ////////////////////////////////////////////////////////////////////////
    ///// STARTING THE QUERY PROCESS ///////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    // reading queries
    vector<FeatureVector> queries;
    get_queries(queries_file, queries);
    if (VERBOSE)
      cerr << "number of queries: " << queries.size() << endl;

    // "n" query points requires a "n*t" results
    vector<vector<Result> > results(queries.size());
    for (size_t i = 0; i < queries.size(); ++i) {
      execute_query(fv_lookup, hf_lookup, ht_lookup, nng, queries[i],
                    n_neighbors, max_proximity_radius, results[i]);
      if (VERBOSE)
        cerr << '\r' << "processing queries: "
             << percent(i, queries.size()) << "%\r";
    }
    if (VERBOSE)
      cerr << '\r' << "processing queries: 100% ("
           << queries.size() << ")" << endl;

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
