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


// #include <iostream>
// #include <string>
// #include <iterator>
#include "EngineDB.hpp"


#include "LSHAngleHashFunction.hpp"
#include "FeatureVector.hpp"

using std::string;
using std::cout;
using std::cerr;
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
    conn.set_option(new mysqlpp::MultiStatementsOption(true));
    if(!conn.connect(db.c_str(), server.c_str(), user.c_str(), pass.c_str()))
      std::cerr << "DB connection failed: " << conn.error() << std::endl;
}


bool
EngineDB::process_deletion(const std::string &fv_id) {
  return delete_feature_vec(fv_id);
}
  

bool 
EngineDB::process_insertion( const FeatureVector &fv, 
                             const std::string &path,
                             const HashFunLookup &hfs,
                             const std::vector<Result> &neighbors) {
  mysqlpp::Transaction trans(conn, 
      mysqlpp::Transaction::serializable,
      mysqlpp::Transaction::session);

  // insert fv_id and path to table feature vector
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
  insert_hash_function(hf.get_id(), hf_path);

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
  query << "insert into hash_function values (" 
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

  mysqlpp::Query query = conn.query();
  query << "delete from hash_function where id = " 
    << mysqlpp::quote << hf_id;
  // return query.execute();
  return string("oldest_hf");
}


bool 
EngineDB::delete_oldest_hash_function() {

  mysqlpp::Query query = conn.query();
  query << "delete from hash_function order by update_time asc limit 1"; 
  return query.execute();
}
