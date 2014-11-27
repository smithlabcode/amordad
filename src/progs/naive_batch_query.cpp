/*
 *    Part of SMITHLAB software
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
#include <cstdio>
#include <iterator>
#include <queue>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include "FeatureVector.hpp"

using std::string;
using std::vector;
using std::cerr;
using std::endl;
using std::cout;
using std::pair;
using std::make_pair;

size_t comparisons = 0;



/******************************************************************************
 * Do the query process and returns the
 * index of the first hash table that identifies the query.
 *****************************************************************************/

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


static void
exec_query(const vector<FeatureVector> &database,
           const FeatureVector &query,
           const size_t n_neighbors,
           const double max_proximity_radius,
           vector<Result> &results) {

  std::priority_queue<Result, vector<Result>, std::greater<Result> > pq;
  double current_dist_cutoff = max_proximity_radius;
  for (vector<FeatureVector>::const_iterator i(database.begin());
       i != database.end(); ++i) {
    const double dist = query.compute_angle(*i);
    ++comparisons;
    if (dist < current_dist_cutoff) {
      if (pq.size() == n_neighbors) {
        pq.pop();
        current_dist_cutoff = dist;
      }
      pq.push(Result(i->get_id(), dist));
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
    OptionParser opt_parse(strip_path(argv[0]), "batch query a set of "
                           "feature vectors",
                           "<database-file> <query-file> <outfile>");
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
    const string database_file(leftover_args.front());
    const string queries_file(leftover_args[1]);
    const string outfile(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/
    
    if (!validate_file(queries_file, 'r'))
      throw SMITHLABException("bad queries file: " + queries_file);

    if (!validate_file(database_file, 'r'))
      throw SMITHLABException("bad database file: " + database_file);
    
    if (!validate_file(outfile, 'w'))
      throw SMITHLABException("bad output file: " + outfile);
    
    ////////////////////////////////////////////////////////////////////////
    ////// READING DATABASE ////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    // reading database
    vector<FeatureVector> database;
    get_queries(database_file, database);
    if (VERBOSE)
      cerr << "database size: " << database.size() << endl;
    
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
      exec_query(database, queries[i],
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
