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
#include <iostream>
#include <string>
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


EngineDB:: EngineDB(const std::string db, const std::string server, 
           const std::string user, const std::string pass) :
  db(db), server(server), user(user), pass(pass) { 
    conn.set_option(new mysqlpp::MultiStatementsOption(true));
    if(!conn.connect(db.c_str(), server.c_str(), user.c_str(), pass.c_str()))
      std::cerr << "DB connection failed: " << conn.error() << std::endl;
}


bool
EngineDB::process_deletion(const std::string &fv_id) const {
  return delete_feature_vec(fv_id);
}
  

bool 
EngineDB::process_insertion( const FeatureVector &fv, 
                             const std::string &path,
                             const HashFunLookup &hfs,
                             const std::vector<Result> &neighbors) const {
  mysqlpp::Transaction trans(conn, 
      mysqlpp::Transaction::serializable,
      mysqlpp::Transaction::session);

  // insert fv_id and path to table feature vector
  insert_feature_vec(fv.get_id(), path);

  // insert fv_id and its hash value to each hash table
  for (size_t i = 0; i < hfs.size(); ++i) {
    size_t hash_value = hfs[i].second(fv);
    insert_hash_occupant(hfs[i].first, hash_value, fv.get_id());
  }

  // insert fv_id and its neighbor to graph
  for (size_t i = 0; i < neighbors.size(); ++i)
    insert_graph_edge(fv.get_id(), neighbors[i].id, neighbors.val);

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
