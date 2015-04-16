/*
 *    Part of AMORDAD software
 *
 *    Copyright (C) 2014 University of Southern California,
 *                       Andrew D. Smith and Wenzheng Li
 *
 *    Authors: Andrew D. Smith, Wenzheng Li
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
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ENGINE_DB
#define ENGINE_DB

#include <mysql++/mysql++.h>

#include <string>
#include <vector>
#include <tr1/unordered_map>
#include <limits>


struct Result {
  Result(const std::string &i, const double v) : id(i), val(v) {}
  Result() : val(std::numeric_limits<double>::max()) {}
  bool operator<(const Result &other) const {return val < other.val;}
  std::string id;
  double val;
};

std::ostream &
operator<<(std::ostream &os, const Result &r);

struct Edge {
  Edge(const std::string &u, const std::string &v, const double d) 
    : src(u), dst(v), dist(d) {}
  Edge() : dist(std::numeric_limits<double>::max()) {}
  std::string src;
  std::string dst;
  double dist;
};

std::ostream &
operator<<(std::ostream &os, const Edge &e);


class LSHAngleHashFunction;
class LSHAngleHashTable;
class FeatureVector;
class RegularNearestNeighborGraph;
typedef std::tr1::unordered_map<std::string, LSHAngleHashFunction> HashFunLookup;
typedef std::tr1::unordered_map<std::string, LSHAngleHashTable> HashTabLookup;
typedef std::tr1::unordered_map<std::string, FeatureVector> FeatVecLookup;
typedef std::tr1::unordered_map<std::string, std::string> PathLookup;

class EngineDB {
public:
  EngineDB() {}
  EngineDB(const std::string db, const std::string server, 
           const std::string user, const std::string pass);

  bool process_deletion(const std::string &fv_id);

  bool process_insertion(const FeatureVector &fv, const std::string &path,
                         const HashFunLookup &hfs,
                         const std::vector<Result> &neighbors);

  bool process_refresh(const LSHAngleHashFunction &hf, const std::string &path,
                       const FeatVecLookup &fvs,
                       const std::vector<Edge> &added_edges);

  std::string get_oldest_hash_function();

  void initialize_db(const PathLookup &fv_paths,
                     const PathLookup &hf_paths,
                     const HashTabLookup &hts, 
                     RegularNearestNeighborGraph &g,
                     bool VERBOSE);

  void read_db(PathLookup &fv_paths, PathLookup &hf_paths, HashTabLookup &hts, 
               RegularNearestNeighborGraph &g, bool VERBOSE);

private:
  std::string db;
  std::string server;
  std::string user;
  std::string pass;
  mysqlpp::Connection conn;

  bool delete_feature_vec(const std::string &fv_id);
  bool insert_feature_vec(const std::string &fv_id, const std::string &path);
  bool insert_hash_occupant(const std::string &hf_id, const size_t hash_value, 
                            const std::string &fv_id);
  bool insert_graph_edge(const std::string &fv_id,
                         const std::string &ng_id,
                         const double dist);
  bool insert_hash_function(const std::string &hf_id, const std::string &path);
  bool delete_hash_function(const std::string &hf_id);
  bool delete_oldest_hash_function();
public:
  void get_feature_vecs(PathLookup &fv_paths);
  void get_hash_funcs(PathLookup &hf_paths);
  void get_hash_table(LSHAngleHashTable &ht);
};

#endif
