/*
 *    Part of AMORDAD software
 *
 *    Copyright (C) 2014 University of Southern California and
 *                       Andrew D. Smith
 *                       Ehsan Behnam
 *                       Wenzheng Li
 *
 *    Authors: Andrew D. Smith, Ehsan Behnam and Wenzheng Li
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

#include "EngineDB.hpp"

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
typedef unordered_map<string, FeatureVector> FeatVecLookup;


static FeatureVector
get_query(const string &query_file) {

  FeatureVector fv;
  std::ifstream in(query_file.c_str());
  if (!in)
    throw SMITHLABException("bad feature vector file: " + query_file);
  in >> fv;
  return fv;
}


static void
evaluate_candidates(const unordered_map<string, FeatureVector> &fvs,
                    const FeatureVector &query,
                    const size_t n_neighbors,
                    const double max_proximity_radius,
                    const unordered_set<string> &candidates,
                    vector<Result> &results) {

  std::priority_queue<Result, vector<Result>, std::less<Result> > pq;
  double current_dist_cutoff = max_proximity_radius;
  for (unordered_set<string>::const_iterator i(candidates.begin());
       i != candidates.end(); ++i) {
    const FeatureVector fv(fvs.find(*i)->second);
    const double dist = query.compute_angle(fv);
    if (dist < current_dist_cutoff) {
      if (pq.size() == n_neighbors) 
        pq.pop();
      pq.push(Result(*i, dist));

      if(pq.size() == n_neighbors)
        current_dist_cutoff = pq.top().val;
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
              RegularNearestNeighborGraph &g,
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

  evaluate_candidates(fvs, query, n_neighbors,
                      max_proximity_radius, candidates, results);
}


static bool
execute_insertion(unordered_map<string, FeatureVector> &fvs,
                  const unordered_map<string, LSHFun> &hfs,
                  unordered_map<string, LSHTab> &hts,
                  RegularNearestNeighborGraph &g,
                  const string  &query_path,
                  const size_t n_neighbors,
                  EngineDB &eng) {

  FeatureVector query = get_query(query_path);

  /// TEST WHETHER QUERY IS ALREADY IN GRAPH
  /// IF NOT ADD QUERY AS A NEW VERTEX
  if (!g.add_vertex_if_new(query.get_id()))
    throw SMITHLABException("cannot insert existing node: " + query.get_id());

  /// INSERT QUERY INTO THE FEATURE VECTOR MAP
  fvs[query.get_id()] = query;

  unordered_set<string> candidates;
  
  // iterate over hash tables
  for (unordered_map<string, LSHTab>::iterator i(hts.begin());
       i != hts.end(); ++i) {

    unordered_map<string, LSHFun>::const_iterator hf(hfs.find(i->first));
    assert(hf != hfs.end());

    // hash the query
    const size_t bucket_number = hf->second(query);
    unordered_map<size_t, vector<string> >::const_iterator
      bucket = i->second.find(bucket_number);

    if (bucket != i->second.end())
      candidates.insert(bucket->second.begin(), bucket->second.end());
    
    // INSERT THE QUERY INTO EACH HASH TABLE
    i->second.insert(query, bucket_number);
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

  vector<Result> neighbors;
  double max = std::numeric_limits<double>::max();
  evaluate_candidates(fvs, query, n_neighbors, max, candidates, neighbors);
  
  // CONNECT THE QUERY TO EACH OF THE "RESULT" NEIGHBORS
  // IN THE GRAPH
  for (vector<Result>::const_iterator i(neighbors.begin());
       i != neighbors.end(); ++i)
    g.update_vertex(query.get_id(), i->id, i->val);

  // UPDATE THE DATABASE
  eng.process_insertion(query, query_path, hfs, neighbors);

  return true;
}


static void
execute_deletion(unordered_map<string, FeatureVector> &fvs,
                 const unordered_map<string, LSHFun> &hfs,
                 unordered_map<string, LSHTab> &hts,
                 RegularNearestNeighborGraph &g,
                 const FeatureVector &query,
                 EngineDB &eng) {
  
  // iterate over hash tables
  for (unordered_map<string, LSHTab>::iterator i(hts.begin());
       i != hts.end(); ++i) {

    unordered_map<string, LSHFun>::const_iterator hf(hfs.find(i->first));
    assert(hf != hfs.end());

    // hash the query
    const size_t bucket_number = hf->second(query);
    unordered_map<size_t, vector<string> >::const_iterator
      bucket = i->second.find(bucket_number);

    // delete the query from each hash table
    if (bucket != i->second.end())
      i->second.remove(query, bucket_number);
  }

  g.remove_vertex(query.get_id());

  // delete the query from the feature vectors map
  fvs.erase(query.get_id());

  // update the database
  eng.process_deletion(query.get_id());
}

static void
add_relations_from_bucket(const vector<string> &bucket, 
                          const FeatVecLookup &featvecs,
                          RegularNearestNeighborGraph &nng) {

  // iterate over bucket
  for (size_t i = 0; i < bucket.size(); ++i) {
    FeatVecLookup::const_iterator ii(featvecs.find(bucket[i]));
    assert(ii != featvecs.end());

    // iterate over other members of bucket
    for (size_t j = i + 1; j < bucket.size(); ++j) {
      FeatVecLookup::const_iterator jj(featvecs.find(bucket[j]));
      assert(jj != featvecs.end());


      // compare and update graph
      const double w = ii->second.compute_angle(jj->second);

      nng.update_vertex(bucket[j], bucket[i], w);
      nng.update_vertex(bucket[i], bucket[j], w);
    }
  }
}


static void
execute_refresh(const unordered_map<string, FeatureVector> &fvs,
                unordered_map<string, LSHFun> &hfs,
                unordered_map<string, LSHTab> &hts,
                RegularNearestNeighborGraph &g,
                const LSHAngleHashFunction &hash_fun) {

  // INITIALIZE THE HASH TABLE
  LSHAngleHashTable hash_table(hash_fun.get_id());
  for (unordered_map<string, FeatureVector>::const_iterator i(fvs.begin());
       i != fvs.end(); ++i)
    hash_table.insert(i->second, hash_fun(i->second));

  // iterate over buckets
  for (BucketMap::const_iterator j(hash_table.begin()); j != hash_table.end(); ++j)
    add_relations_from_bucket(j->second, fvs, g);

  // remove the oldest hash function and associated hash table
  // replaced by the new ones
  // TODO: retrive time info from database

  string oldest_hf = hfs.begin()->first;
  unordered_map<string, LSHTab>::const_iterator ht(hts.find(oldest_hf));
  hts.erase(ht->first);
  hfs.erase(oldest_hf);
  hts[hash_table.get_id()] = hash_table;
  hfs[hash_fun.get_id()] = hash_fun;
}
 
/*
 * Config file has this structure:
 * feature_vectors_file: path to file containing filenames of feature vectors
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


static void
execute_commands(const string &command_file,
                 unordered_map<string, FeatureVector> &fvs,
                 unordered_map<string, LSHFun> &hfs,
                 unordered_map<string, LSHTab> &hts,
                 RegularNearestNeighborGraph &g) {

  size_t n_neighbors = 20;
  double max_proximity_radius = 0.75;

  std::ifstream in(command_file.c_str());
  if (!in)
    throw SMITHLABException("bad command file: " + command_file);

  vector<pair<string,string> > commands;
  string line;
  while (getline(in, line)) {
    std::istringstream iss(line);
    string operation;
    string query_path;
    if(!(iss >> operation >> query_path))
      throw SMITHLABException("bad command format: " + line);

    commands.push_back(make_pair(operation, query_path));
  }

  string db = "amorgin";
  string server = "localhost";
  string user = "root";
  string pass = "580230mysql";

  EngineDB eng(db,server,user,pass);

  for(size_t i = 0; i < commands.size(); ++i) {

    // execute different functions based on the command
    if(commands[i].first == "query") {
      FeatureVector fv = get_query(commands[i].second);
      vector<Result> results;
      execute_query(fvs, hfs, hts, g, fv, n_neighbors, 
          max_proximity_radius, results);
    }
    else if(commands[i].first == "insert") {
      string query_path = commands[i].second;
      execute_insertion(fvs, hfs, hts, g, query_path, n_neighbors, eng); 
    }
    else if(commands[i].first == "delete") {
      FeatureVector fv = get_query(commands[i].second);
      execute_deletion(fvs, hfs, hts, g, fv, eng); 
    }
    else if(commands[i].first == "refresh") {
      string hash_fun_file = commands[i].second;
      std::ifstream hash_fun_in(hash_fun_file.c_str());
      if (!hash_fun_in)
        throw SMITHLABException("cannot open: " + hash_fun_file);
      LSHAngleHashFunction hf;
      hash_fun_in >> hf;
      execute_refresh(fvs, hfs, hts, g, hf); 
    }
    else
      throw SMITHLABException("unknown command");
  }
}


int
main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "amordad server supporting search, "
                           "insertion, deletion and refresh with "
                           "database residing on disk",
                           "<config-file> <command-file>");
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
    const string command_file(leftover_args[1]);
    /****************** END COMMAND LINE OPTIONS *****************/

    if (!validate_file(command_file, 'r'))
      throw SMITHLABException("bad command file: " + command_file);


    ////////////////////////////////////////////////////////////////////////
    ////// READING HASH {FUNCTIONS, TABLES, FEATURE_VECTOR} FILES //////////////
    ////////////////////////////////////////////////////////////////////////

    string fv_paths_file, hf_paths_file, ht_paths_file, graph_file;
    read_config_file(database_config_file, fv_paths_file, hf_paths_file,
                     ht_paths_file, graph_file);

    // reading samples in database
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
    ///// EXECUTE THE COMMANDS ///////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    execute_commands(command_file, fv_lookup, hf_lookup, ht_lookup, nng);
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
