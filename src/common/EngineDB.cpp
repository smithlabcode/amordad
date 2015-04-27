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

#include "EngineDB.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <climits>
#include <numeric>

#include "LSHAngleHashFunction.hpp"
#include "FeatureVector.hpp"
#include "LSHAngleHashTable.hpp"
#include "RegularNearestNeighborGraph.hpp"

using std::string;
using std::cerr;
using std::vector;
using std::endl;

std::ostream &
operator<<(std::ostream &os, const Result &r) {
  return os << r.id << '\t' << r.val;
}


std::ostream &
operator<<(std::ostream &os, const Edge &e) {
  return os << e.src << "->" << e.dst << '\t' << e.dist;
}


EngineDB:: EngineDB(const std::string db, const std::string server, 
                    const std::string user, const std::string pass) :
  db(db), server(server), user(user), pass(pass) { 

    // conn.set_option(new mysqlpp::MultiStatementsOption(true));
    if(!conn.connect(db.c_str(), server.c_str(), user.c_str(), pass.c_str()))
      throw SMITHLABException("Cannot connect the database");
}


bool
EngineDB::process_deletion(const std::string &fv_id) {
  return delete_feature_vec(fv_id);
}
  

bool
EngineDB::process_insertion(const FeatureVector &fv, 
                            const std::string &path,
                            const HashFunLookup &hfs,
                            const std::vector<Result> &neighbors) {
  mysqlpp::Transaction trans(conn, 
      mysqlpp::Transaction::serializable,
      mysqlpp::Transaction::session);

  // insert fv_id and path to table feature_vector
  insert_feature_vec(fv.get_id(), path);

  // insert fv_id and its hash value to each hash table
  for (HashFunLookup::const_iterator i(hfs.begin());
       i != hfs.end(); ++i) {
    size_t hash_value = i->second(fv);
    insert_hash_occupant(i->first, hash_value, fv.get_id());
  }

  // insert fv_id and its neighbor to graph
  for (size_t i = 0; i < neighbors.size(); ++i)
    insert_graph_edge(fv.get_id(), neighbors[i].id, neighbors[i].val);

  trans.commit();
  return true;
}


bool
EngineDB::process_refresh(const LSHAngleHashFunction &hf, 
                          const std::string &path,
                          const FeatVecLookup &fvs,
                          const std::vector<Edge> &added_edges) {

  mysqlpp::Transaction trans(conn, 
      mysqlpp::Transaction::serializable,
      mysqlpp::Transaction::session);

  // insert new hash function and its path to hash_function
  insert_hash_function(hf.get_id(), path);

  // insert fv_id and its hash value for each fv
  for (FeatVecLookup::const_iterator i(fvs.begin());
      i != fvs.end(); ++i) {
    size_t hash_value = hf(i->second);
    insert_hash_occupant(hf.get_id(), hash_value, i->first);
  }

  // insert updated edges to graph
  for (size_t i = 0; i < added_edges.size(); ++i)
    insert_graph_edge(added_edges[i].src, added_edges[i].dst, added_edges[i].dist);

  // delete the oldest hash function
  delete_oldest_hash_function();

  trans.commit();
  return true;
}


void
EngineDB::initialize_db(const PathLookup &fv_paths,
                        const PathLookup &hf_paths,
                        const HashTabLookup &hts, 
                        RegularNearestNeighborGraph &g,
                        bool VERBOSE) {

  if(VERBOSE)
    cerr << "BEGIN INITIALIZING THE ENGINE DB" << endl;

  size_t count = 0;

  // insert feature vectors
  for (PathLookup::const_iterator i(fv_paths.begin());
       i != fv_paths.end(); ++i) {
    insert_feature_vec(i->first, i->second);
    if (VERBOSE)
      cerr << '\r' << "insert feature vectors: "
           << percent(count, fv_paths.size()) << "%\r";
    count += 1;
  }

  if (VERBOSE)
    cerr << "insert feature vectors: 100% (" << fv_paths.size() << ")" << endl;

  count = 0;
  // insert hash functions
  for (PathLookup::const_iterator i(hf_paths.begin());
       i != hf_paths.end(); ++i) {
    insert_hash_function(i->first, i->second);
    if (VERBOSE)
      cerr << '\r' << "insert hash functions: "
           << percent(count, hf_paths.size()) << "%\r";
    count += 1;
  }
  if (VERBOSE)
    cerr << "insert hash functions: 100% (" << hf_paths.size() << ")" << endl;


  count = 0;
  // insert buckets info
  for (HashTabLookup::const_iterator i(hts.begin());
       i != hts.end(); ++i) {
    for (BucketMap::const_iterator j(i->second.begin());
         j != i->second.end(); ++j)
      for (size_t k = 0; k < j->second.size(); ++k) 
        insert_hash_occupant(i->first, j->first, j->second[k]);

    if (VERBOSE)
      cerr << '\r' << "insert hash table buckets: "
           << percent(count, hts.size()) << "%\r";
    count += 1;
  }
  if (VERBOSE)
    cerr << "insert hash table buckets: 100% (" << hts.size() << ")" << endl;


  count = 0;
  // insert graph edges
  for (PathLookup::const_iterator i(fv_paths.begin());
       i != fv_paths.end(); ++i) {

    vector<string> neighbors;
    vector<double> distances;
    g.get_neighbors(i->first, neighbors, distances);
    if(neighbors.size() != distances.size())
      throw SMITHLABException("neighbors size must be equal to distances size");

    for (size_t j = 0; j < neighbors.size(); ++j) {
      insert_graph_edge(i->first, neighbors[j], distances[j]);
      if (VERBOSE)
        cerr << '\r' << "insert graph edges: "
             << percent(count, g.get_edge_count()) << "%\r";
      count += 1;
    }
  }
  if (VERBOSE)
    cerr << "insert graph edges: 100% (" << g.get_edge_count() << ")" << endl;
}


void
EngineDB::read_db(PathLookup &fv_paths, 
                  PathLookup &hf_paths,
                  HashTabLookup &hts,
                  RegularNearestNeighborGraph &g,
                  bool VERBOSE) {

  size_t count = 0;

  get_feature_vecs(fv_paths);

  if (VERBOSE)
    cerr << "read from db feature vectors: 100% (" 
         << fv_paths.size() << ")" << endl;

  get_hash_funcs(hf_paths);

  if (VERBOSE)
    cerr << "read from db hash functions: 100% (" 
         << hf_paths.size() << ")" << endl;


  for(PathLookup::const_iterator i(hf_paths.begin());
      i != hf_paths.end(); ++i) {
    LSHAngleHashTable hash_table(i->first);
    get_hash_table(hash_table);
    hts[hash_table.get_id()] = hash_table;

    if (VERBOSE)
      cerr << '\r' << "reading from db hash tables: "
           << percent(count++, hf_paths.size()) << "%\r";
  }

  if (VERBOSE)
    cerr << "read from db hash tables: 100% ("
         << hf_paths.size() << ")" << endl;

  for(PathLookup::const_iterator i(fv_paths.begin());
      i != fv_paths.end(); ++i)
    g.add_vertex(i->first);
  
  get_graph_edges(g);
  if (VERBOSE)
    cerr << "read from db graph: 100%" << endl;
}


bool
EngineDB::delete_feature_vec(const std::string &fv_id) {

  mysqlpp::Query query = conn.query();
  query << "delete from feature_vector where id=" 
        << mysqlpp::quote << fv_id;
  return query.execute();
}
 

bool
EngineDB::insert_feature_vec(const std::string &fv_id, 
                             const std::string &path) {

  mysqlpp::Query query = conn.query();
  query << "insert into feature_vector values (" 
        << mysqlpp::quote << fv_id << ","
        << mysqlpp::quote << path << ");";
  return query.execute();
}


bool 
EngineDB::insert_hash_occupant(const std::string &hf_id,
                               const size_t hash_value, 
                               const std::string &fv_id) {

  mysqlpp::Query query = conn.query();
  query << "insert into hash_table_bucket values (" 
    << mysqlpp::quote << hf_id << "," << hash_value << ","
    << mysqlpp::quote << fv_id << ");";
  return query.execute();
}


bool
EngineDB::insert_graph_edge(const std::string &fv_id,
                            const std::string &ng_id,
                            const double dist) {

  mysqlpp::Query query = conn.query();
  query << "insert into graph_edge values (" 
    << mysqlpp::quote << fv_id << ","
    << mysqlpp::quote << ng_id << ","
    << dist << ");";
  return query.execute();
}


bool 
EngineDB::insert_hash_function(const std::string &hf_id, 
                               const std::string &path) {

  mysqlpp::Query query = conn.query();
  query << "insert into hash_function (id, path) values (" 
    << mysqlpp::quote << hf_id << ","
    << mysqlpp::quote << path << ");";
  return query.execute();
}


bool 
EngineDB::delete_hash_function(const std::string &hf_id) {

  mysqlpp::Query query = conn.query();
  query << "delete from hash_function where id = " 
        << mysqlpp::quote << hf_id;
  return query.execute();
}


string
EngineDB::get_oldest_hash_function() {

  string old_hash = "";
  mysqlpp::Query query = conn.query();
  query << "select id from hash_function order by update_time asc limit 1"; 
  if(mysqlpp::StoreQueryResult res = query.store()) {
    if(res.num_rows() > 0) {
      res[0][0].to_string(old_hash);
      return old_hash;
    }
  }
  else
    throw SMITHLABException("Failed to retrive hash functions");
  return old_hash;
}


string
EngineDB::get_newest_hash_function() {

  string new_hash = "";
  mysqlpp::Query query = conn.query();
  query << "select id from hash_function order by update_time desc limit 1"; 
  if(mysqlpp::StoreQueryResult res = query.store()) {
    if(res.num_rows() > 0) {
      res[0][0].to_string(new_hash);
      return new_hash;
    }
  }
  else
    throw SMITHLABException("Failed to retrive hash functions");
  return new_hash;
}


bool 
EngineDB::delete_oldest_hash_function() {

  mysqlpp::Query query = conn.query();
  query << "delete from hash_function order by update_time asc limit 1"; 
  return query.execute();
}


void
EngineDB::get_feature_vecs(PathLookup &fv_paths) {

  mysqlpp::Query query = conn.query();
  query << "select id, path from feature_vector"; 
  if(mysqlpp::StoreQueryResult res = query.store()) {
    for(size_t i = 0; i < res.num_rows(); ++i) {
      string id = "";
      res[i][0].to_string(id);
      string path = "";
      res[i][1].to_string(path);
      fv_paths[id] = path;
    }
  }
  else
    throw SMITHLABException("Failed to retrive hash functions");
}


void
EngineDB::get_hash_funcs(PathLookup &hf_paths) {

  mysqlpp::Query query = conn.query();
  query << "select id, path from hash_function"; 
  if(mysqlpp::StoreQueryResult res = query.store()) {
    for(size_t i = 0; i < res.num_rows(); ++i) {
      string id = "";
      res[i][0].to_string(id);
      string path = "";
      res[i][1].to_string(path);
      hf_paths[id] = path;
    }
  }
  else
    throw SMITHLABException("Failed to retrive hash functions");
}


void
EngineDB::get_hash_table(LSHAngleHashTable &ht) {

  mysqlpp::Query query = conn.query();
  query << "select hash_key, occupant from hash_table_bucket where id = "
        << mysqlpp::quote << ht.get_id();
  if(mysqlpp::StoreQueryResult res = query.store()) {
    for(size_t i = 0; i < res.num_rows(); ++i) {
      size_t hash_key = res[i][0];
      string occupant = "";
      res[i][1].to_string(occupant);
      ht.insert(occupant, hash_key);
    }
  }
  else
    throw SMITHLABException("Failed to retrive hash functions");
}


void
EngineDB::get_graph_edges(RegularNearestNeighborGraph &nng) {

  mysqlpp::Query query = conn.query();
  query << "select src, dst, dist from graph_edge";
  if(mysqlpp::StoreQueryResult res = query.store()) {
    for(size_t i = 0; i < res.num_rows(); ++i) {
      string src = "";
      res[i][0].to_string(src);
      string dst = "";
      res[i][1].to_string(dst);
      double dist = res[i][2];
      nng.update_vertex(src, dst, dist);
    }
  }
  else
    throw SMITHLABException("Failed to retrive hash functions");
}


void
EngineDB::get_graph(RegularNearestNeighborGraph &nng) {
  PathLookup fv_paths;
  get_feature_vecs(fv_paths);
  for(PathLookup::const_iterator i(fv_paths.begin());
      i != fv_paths.end(); ++i)
    nng.add_vertex(i->first);

  get_graph_edges(nng);
}


size_t
EngineDB::get_num_hash_functions() {

  mysqlpp::Query query = conn.query();
  query << "select count(*) from hash_function";
  if(mysqlpp::StoreQueryResult res = query.store())
    return res[0][0];
  else
    throw SMITHLABException("Failed to retrive hash functions");
  return 0;
}
