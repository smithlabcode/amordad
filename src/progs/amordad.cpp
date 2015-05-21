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
#include <unordered_set>
#include <cstdio>
#include <iterator>
#include <queue>
#include <iostream>
#include <chrono>

#include <gsl/gsl_statistics_double.h>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include "RegularNearestNeighborGraph.hpp"

#include "FeatureVector.hpp"
#include "LSHAngleHashTable.hpp"
#include "LSHAngleHashFunction.hpp"

#include "EngineDB.hpp"
#include "crow.h"
#include "json.h"

using std::string;
using std::vector;
using std::cerr;
using std::endl;
using std::cout;
using std::pair;
using std::make_pair;
using std::queue;

using std::unordered_map;
using std::unordered_set;

typedef LSHAngleHashTable LSHTab;
typedef LSHAngleHashFunction LSHFun;
typedef unordered_map<string, FeatureVector> FeatVecLookup;


static FeatureVector
get_feat_vec(const string &fv_path) {

  FeatureVector fv;
  std::ifstream in(fv_path.c_str());
  if (!in)
    throw SMITHLABException("bad feature vector file: " + fv_path);
  in >> fv;

  // normalize feature vector
  const double mean = gsl_stats_mean(&fv[0], 1, fv.size());
  const double sd = gsl_stats_sd_m(&fv[0], 1, fv.size(), mean);

  transform(fv.begin(), fv.end(), fv.begin(), 
            bind2nd(std::minus<double>(), mean));
  transform(fv.begin(), fv.end(), fv.begin(), 
            bind2nd(std::divides<double>(), sd));

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


static void
execute_insertion(unordered_map<string, FeatureVector> &fvs,
                  const unordered_map<string, LSHFun> &hfs,
                  unordered_map<string, LSHTab> &hts,
                  RegularNearestNeighborGraph &g,
                  const string  &query_path,
                  EngineDB &eng) {

  FeatureVector query = get_feat_vec(query_path);

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
  size_t max_deg = g.get_maximum_degree();
  evaluate_candidates(fvs, query, max_deg, max, candidates, neighbors);
  
  // CONNECT THE QUERY TO EACH OF THE "RESULT" NEIGHBORS
  // IN THE GRAPH
  for (vector<Result>::const_iterator i(neighbors.begin());
       i != neighbors.end(); ++i)
    g.update_vertex(query.get_id(), i->id, i->val);

  // UPDATE THE DATABASE
  eng.process_insertion(query, query_path, hfs, neighbors);
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
                          RegularNearestNeighborGraph &nng,
                          vector<Edge> &added_edges) {

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

      if(nng.update_vertex(bucket[j], bucket[i], w))
        added_edges.push_back(Edge(bucket[j], bucket[i], w));
      if(nng.update_vertex(bucket[i], bucket[j], w))
        added_edges.push_back(Edge(bucket[i], bucket[j], w));
    }
  }
}


static void
execute_refresh(const unordered_map<string, FeatureVector> &fvs,
                unordered_map<string, LSHFun> &hfs,
                queue<string> &hf_queue,
                unordered_map<string, LSHTab> &hts,
                RegularNearestNeighborGraph &g,
                const string &hash_fun_file,
                EngineDB &eng) {

  // READ THE HASH FUNCTION
  std::ifstream hash_fun_in(hash_fun_file.c_str());
  if (!hash_fun_in)
    throw SMITHLABException("cannot open: " + hash_fun_file);
  LSHAngleHashFunction hash_fun;
  hash_fun_in >> hash_fun;

  if(hfs.find(hash_fun.get_id()) != hfs.end())
    throw SMITHLABException("attempt to insert an existing hash function");

  // INITIALIZE THE HASH TABLE
  LSHAngleHashTable hash_table(hash_fun.get_id());
  for (unordered_map<string, FeatureVector>::const_iterator i(fvs.begin());
       i != fvs.end(); ++i)
    hash_table.insert(i->second, hash_fun(i->second));

  vector<Edge> added_edges;
  // iterate over buckets
  for (BucketMap::const_iterator j(hash_table.begin()); 
       j != hash_table.end(); ++j)
    add_relations_from_bucket(j->second, fvs, g, added_edges);

  // remove the oldest hash function and associated hash table
  // replaced by the new ones

  string oldest_hf = hf_queue.front();
  hf_queue.pop();
  hts.erase(oldest_hf);
  hfs.erase(oldest_hf);
  hts[hash_table.get_id()] = hash_table;
  hfs[hash_fun.get_id()] = hash_fun;
  hf_queue.push(hash_fun.get_id());

  // update the database
  eng.process_refresh(hash_fun, hash_fun_file, fvs,
                      added_edges, g.get_maximum_degree());
}
 

static void
get_database(const bool VERBOSE, 
             unordered_map<string, string> &paths,
             unordered_map<string, FeatureVector> &db) {

  size_t count = 0;
  for(unordered_map<string, string>::const_iterator i(paths.begin());
      i != paths.end(); ++i) {

    FeatureVector fv = get_feat_vec(i->second.c_str());

    if(fv.get_id() != i->first)
      throw SMITHLABException("unconsistent feature vector ids");
    db[fv.get_id()] = fv;
    if (VERBOSE)
      cerr << "\rloading feature vectors: "
           << percent(count++, paths.size()) << "%\r";
  }
  if (VERBOSE)
    cerr << "\rloading feature vectors: 100%" << endl;
}


static
void add_hash_functions(size_t qsize, size_t n_bits, size_t n_features, 
                        const string &feature_set_id, const string &hf_dir, 
                        unordered_map<string, string> &hf_paths,
                        queue<string> &hash_func_queue) {

  for(size_t i = 0; i < qsize; ++i) {
    string id = toa(i);
    const LSHAngleHashFunction hash_function(id, feature_set_id,
                                             n_features, n_bits);
    std::ofstream of;
    string outfile = path_join(hf_dir,id);
    outfile = outfile + ".hf";
    if (!outfile.empty()) of.open(outfile.c_str());
    if (!of) throw SMITHLABException("cannot write to file: " + outfile);
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    out << hash_function << endl;
    hf_paths[hash_function.get_id()] = outfile;
    hash_func_queue.push(hash_function.get_id());
  }
}


static
string add_new_hash_function(size_t n_bits, size_t n_features, 
                             const string &feature_set_id, const string &hf_dir,
                             const queue<string> &hash_func_queue) {

    // find the newest hash function in queue, increase id by 1
    string newest = hash_func_queue.back();
    string id = toa(std::stoi(newest) + 1);
    const LSHAngleHashFunction hash_function(id, feature_set_id,
                                             n_features, n_bits);
    std::ofstream of;
    string outfile = path_join(hf_dir,id);
    outfile = outfile + ".hf";
    if (!outfile.empty()) of.open(outfile.c_str());
    if (!of) throw SMITHLABException("cannot write to file: " + outfile);
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    out << hash_function << endl;
    return outfile;
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

    // hash function parameters
    size_t n_bits = 0;
    size_t n_features = 0;
    string feature_set_id = "FEATURES";

    // regular NNG parameters
    string graph_name("THE_GRAPH");
    size_t max_degree = 1;

    // hash queue size and location
    size_t hf_queue_size = 0;
    string hf_dir;

    // engine database parameters
    string db;
    string pass;
    string user = "root";
    string server = "localhost";

    // server parameter
    size_t PORT = 18080;

    string init_file;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), 
                           "amordad server supporting search, "
                           "insertion, deletion and refresh with "
                           "database residing on RDBMS");

    opt_parse.add_opt("bits", 'b', "bits in hash value", true, n_bits);
    opt_parse.add_opt("nfeat", 'n', "number of features", true, n_features);
    opt_parse.add_opt("deg", 'd', "max out degree of graph", true, max_degree);
    opt_parse.add_opt("qsize", 'q', "queue size for hash functions", 
                      true, hf_queue_size);
    opt_parse.add_opt("hfdir", 'h', "folder for hash functions", true, hf_dir);
    opt_parse.add_opt("mysql", 'm', "name of the mysql database", true, db);
    opt_parse.add_opt("pass", 'p', "password for the mysql database",
                      true, pass);
    opt_parse.add_opt("user", 'u', "username for the mysql database "
                      "(Default: root)", false, user);
    opt_parse.add_opt("server", 's', "server for the mysql database "
                      "(Default: localhost)", false, server);
    opt_parse.add_opt("PORT", 'P', "Port for the server to run at "
                      "(Default: 18080)", false, PORT);
    opt_parse.add_opt("initfile", 'i', "initialize database by providing "
                      "feature paths", false, init_file);
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
    /****************** END COMMAND LINE OPTIONS *****************/

    
    ////////////////////////////////////////////////////////////////////////
    ///// READ DATA FROM ENGINE DATABASE ///////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////


    EngineDB eng(db,server,user,pass);
    unordered_map<string, string> fv_path_lookup;
    unordered_map<string, string> hf_path_lookup;
    queue<string> hash_func_queue;
    unordered_map<string, LSHTab> ht_lookup;
    RegularNearestNeighborGraph nng(graph_name, max_degree);

    if(eng.get_num_hash_functions() == 0) {

      if(VERBOSE)
        cerr << "INITIALIZING HASH FUNCTIONS" << endl;

      add_hash_functions(hf_queue_size, n_bits, n_features, 
                         feature_set_id, hf_dir, hf_path_lookup,
                         hash_func_queue);
      eng.initialize_db(fv_path_lookup, hf_path_lookup, hash_func_queue,
                        ht_lookup, nng, VERBOSE); 
    }

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    eng.read_db(fv_path_lookup, hf_path_lookup, hash_func_queue,
                ht_lookup, nng, VERBOSE);
    
    // READING SAMPLES IN DATABASE
    unordered_map<string, FeatureVector> fv_lookup;
    get_database(VERBOSE, fv_path_lookup, fv_lookup);

    ////////////////////////////////////////////////////////////////////////
    ////// READING THE HASH FUNCTIONS //////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    unordered_map<string, LSHFun> hf_lookup;
    size_t count = 0;
    for(unordered_map<string, string>::const_iterator i(hf_path_lookup.begin());
        i != hf_path_lookup.end(); ++i) {
      std::ifstream hf_in(i->second.c_str());
      if (!hf_in)
        throw SMITHLABException("bad hash function file: " +
                                i->second);
      LSHFun hf;
      hf_in >> hf;
      if(hf.get_id() != i->first) 
        throw SMITHLABException("unconsistent hash function ids");
      hf_lookup[hf.get_id()] = hf;
      if (VERBOSE)
        cerr << '\r' << "load hash functions: "
             << percent(count++, hf_path_lookup.size()) << "%\r";
    }
    if (VERBOSE)
      cerr << "load hash functions: 100% (" << hf_lookup.size() << ")" << endl;

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    if(VERBOSE)
      cerr << "load database time = " << elapsed.count() << "s\n";


    ////////////////////////////////////////////////////////////////////////
    // IF INITIALIZATION FEATURE PATHS FILE PROVIDED, INITIALIZE DATABASE ///
    // //////////////////////////////////////////////////////////////////////
    if(!init_file.empty()) {
      if(!validate_file(init_file, 'r'))
        throw SMITHLABException("bad feature paths file: " + init_file);

      std::ifstream init_in(init_file.c_str());
      if (!init_in)
        throw SMITHLABException("bad feature paths file: " + init_file);

      if(VERBOSE)
        cerr << "INITIALIZING DATABASE FROM FEATURE PATHS" << endl;

      std::chrono::time_point<std::chrono::system_clock> start, end;
      start = std::chrono::system_clock::now();

      string fv_path;
      size_t num_features_initialized = 0;
      while(getline(init_in, fv_path)) {
        execute_insertion(fv_lookup, hf_lookup, ht_lookup, 
                          nng, fv_path, eng);
        num_features_initialized++;
      }

      end = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed = end - start;
      if(VERBOSE) {
        cerr << num_features_initialized << " feature vectors inserted" << endl;
        cerr << "initialization time = " << elapsed.count() << "s\n";
      }
    }


    ////////////////////////////////////////////////////////////////////////
    ///// EXECUTE THE REQUESTS FROM URL ///////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    crow::SimpleApp app;
    CROW_ROUTE(app, "/")
    ([]() {
     return "Amordad Web Server";
     });

    CROW_ROUTE(app, "/query")
    ([&](const crow::request &req) {

      crow::json::wvalue ret;

      try {
        string fv_path = req.url_params.get("path");

        if(fv_path.empty())
          throw SMITHLABException("invalid file path");

        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();
        const size_t n_neighbors = 30;
        const double max_proximity_radius = 0.75;
        vector<Result> result;
        FeatureVector fv = get_feat_vec(fv_path);
        execute_query(fv_lookup, hf_lookup, ht_lookup, 
                      nng, fv, n_neighbors, max_proximity_radius,
                      result);
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed = end - start;

        if(VERBOSE) {
          cerr << "Wall time = " << elapsed.count() << "s\n";
          copy(result.begin(), result.end(),
               std::ostream_iterator<Result>(cerr, "\n"));
          cerr << endl;
        }

        ret["total"] = fv_lookup.size();
        ret["time"] = elapsed.count();
        ret["id"] = fv.get_id();
        for (size_t i = 0; i < result.size(); ++i)
          ret[result[i].id] = result[i].val;

        return ret;
      }
      catch (const SMITHLABException &e) {
        cerr << e.what() << endl;
        ret["error"] = e.what();
        return ret;
      }
    });

    CROW_ROUTE(app, "/insert")
    ([&](const crow::request &req) {

      crow::json::wvalue ret;

      try {
        string fv_path = req.url_params.get("path");

        if(fv_path.empty())
          throw SMITHLABException("invalid file path");

        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();
        execute_insertion(fv_lookup, hf_lookup, ht_lookup, 
                          nng, fv_path, eng);
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        if(VERBOSE)
          cerr << "Wall time = " << elapsed.count() << "s\n";

        ret["total"] = fv_lookup.size();
        ret["time"] = elapsed.count();

        return ret;
      }
      catch (const SMITHLABException &e) {
        cerr << e.what() << endl;
        ret["error"] = e.what();
        return ret;
      }
    });


    CROW_ROUTE(app, "/delete")
    ([&](const crow::request &req) {

      crow::json::wvalue ret;

      try {
        string fv_path = req.url_params.get("path");

        if(fv_path.empty())
          throw SMITHLABException("invalid file path");

        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();
        FeatureVector fv = get_feat_vec(fv_path);
        execute_deletion(fv_lookup, hf_lookup, ht_lookup, 
                         nng, fv, eng);
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        if(VERBOSE)
          cerr << "Wall time = " << elapsed.count() << "s\n";

        ret["total"] = fv_lookup.size();
        ret["time"] = elapsed.count();

        return ret;
      }
      catch (const SMITHLABException &e) {
        cerr << e.what() << endl;
        ret["error"] = e.what();
        return ret;
      }
    });


    CROW_ROUTE(app, "/refresh")
    ([&]() {

      crow::json::wvalue ret;

      try {
        string hf_path = add_new_hash_function(n_bits, n_features, feature_set_id,
                                               hf_dir, hash_func_queue);
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();
        execute_refresh(fv_lookup, hf_lookup, hash_func_queue, ht_lookup, 
                        nng, hf_path, eng);
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed = end - start;

        if(VERBOSE)
          cerr << "Wall time = " << elapsed.count() << "s\n";

        ret["total"] = fv_lookup.size();
        ret["time"] = elapsed.count();
        ret["id"] = hash_func_queue.back();

        return ret;
      }
      catch (const SMITHLABException &e) {
        cerr << e.what() << endl;
        ret["error"] = e.what();
        return ret;
      }
    });

    app.port(PORT)
       .run();
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
