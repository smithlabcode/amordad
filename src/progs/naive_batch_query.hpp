#include <string>
#include <vector>
#include <algorithm>
#include <climits>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <iterator>
#include <queue>
#include <sys/stat.h>

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


bool FileExists(string filename) {
    struct stat fileInfo;
    return stat(filename.c_str(), &fileInfo) == 0;
}

/******************************************************************************
 * Do the query process and returns the
 * index of the first hash table that identifies the query.
 *****************************************************************************/

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
exec_query(const vector<FeatureVector> &database,
           const FeatureVector &query,
           const size_t n_neighbors,
           const double max_proximity_radius,
           vector<Result> &results) {

  std::priority_queue<Result, vector<Result>, std::less<Result> > pq;
  double current_dist_cutoff = max_proximity_radius;
  for (vector<FeatureVector>::const_iterator i(database.begin());
       i != database.end(); ++i) {
    const double dist = query.compute_angle(*i);
    ++comparisons;
    if (dist < current_dist_cutoff) {
      if (pq.size() == n_neighbors)
        pq.pop();
      pq.push(Result(i->get_id(), dist));

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
load_feature_vectors(const bool VERBOSE, const string &fvs_file,
                     vector<FeatureVector> &fvs) {

  vector<string> fv_files;
  get_filenames(fvs_file, fv_files);

  for(size_t i = 0; i < fv_files.size(); ++i) {
    FeatureVector fv;
    if (!FileExists(fv_files[i]))
        throw SMITHLABException("File not found: " + fv_files[i]);
    std::ifstream in(fv_files[i].c_str());
    if (!in)
      throw SMITHLABException("bad feature vector file: " + fv_files[i]);
    in >> fv;
    fvs.push_back(fv);
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


