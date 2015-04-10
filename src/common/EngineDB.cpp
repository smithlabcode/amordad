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

using std::string;
using std::cout;
using std::cerr;
using std::endl;

EngineDB:: EngineDB(const std::string db, const std::string server, 
           const std::string user, const std::string pass) :
  db(db), server(server), user(user), pass(pass) { 
    conn.set_option(new mysqlpp::MultiStatementsOption(true));
    if(!conn.connect(db.c_str(), server.c_str(), user.c_str(), pass.c_str()))
      std::cerr << "DB connection failed: " << conn.error() << std::endl;
}


bool
EngineDB::delete_feature_vec(const std::string &fv_id) {

  mysqlpp::Query query = conn.query();
  query << "delete from feature_vector where id=" 
        << mysqlpp::quote << fv_id;
  if(query.execute())
    return true;

  return false;
}
  

bool
EngineDB::insert_feature_vec(const std::string &fv_id, 
                             const std::string &path) {
  string query_str = "insert into feature_vector values (\"" + fv_id + "\",\"" + path + "\");";
  // string query_str_2 = "insert into feature_vector values (\"" + path + "\",\"" + fv_id + "\");";
  // mysqlpp::Query query = conn.query(query_str.c_str());
  // mysqlpp::Transaction trans(conn, 
  //                            mysqlpp::Transaction::serializable,
  //                            mysqlpp::Transaction::session);
  mysqlpp::Query query1 = conn.query();
  query1 << query_str;
  query1.execute();
  // mysqlpp::Query query2 = conn.query();
  // query2 << query_str_2;
  // query2.execute();

  // std::cout << "\nRow inserted, but not committed. Please confirm!" << std::endl;
  // int yes;
  // std::cin >> yes;
  // std::cout << "\nCommiting transaction gives us:" << std::endl;
  insert_feature_vec_t(fv_id,path);
  // trans.commit();
  return true;
}

bool
EngineDB::insert_feature_vec_t(const std::string &fv_id, 
                             const std::string &path) {
  string query_str_2 = "insert into feature_vector values (\"" + path + "\",\"" + fv_id + "\");";
  mysqlpp::Query query2 = conn.query();
  query2 << query_str_2;
  query2.execute();
  return true;
}
