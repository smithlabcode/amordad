
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
#include <sys/stat.h>

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

bool FileExists(string filename) {
    struct stat fileInfo;
    return stat(filename.c_str(), &fileInfo) == 0;
}


struct Result {
  Result(const string &i, const double v) : id(i), val(v) {}
  Result() : val(std::numeric_limits<double>::max()) {}
  bool operator<(const Result &other) const {return val < other.val;}
  string id;
  double val;
};


std::ostream &
operator<<(std::ostream &os, const Result &r) {
  return os << r.id << '\t' << r.val;
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
    ++comparisons;
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
    if (!FileExists(query_files[i]))
        throw SMITHLABException("File not found: " + query_files[i]);
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
    if (!FileExists(fv_files[i]))
        throw SMITHLABException("File not found: " + fv_files[i]);
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
