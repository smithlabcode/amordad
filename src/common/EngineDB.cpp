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

bool
EngineDB::delete_feature_vec(const std::string &fv_id) {

  std::cout << fv_id << std::endl;
  mysqlpp::Connection conn(false);
  if(conn.connect(db.c_str(), server.c_str(), user.c_str(), pass.c_str())) {
    mysqlpp::Query query = conn.query("select * from sample");
    if (mysqlpp::StoreQueryResult res = query.store()) {
      std::cout << "We have:" << std::endl;
      for (size_t i = 0; i < res.num_rows(); ++i) {
        std::cout << '\t' << res[i][0] << endl;
      }
    }
    else {
      std::cerr << "Failed to get item list: " << query.error() << std::endl;
      return false;
    }
    return true;
  }
  else {
    std::cerr << "DB connection failed: " << conn.error() << std::endl;
    return false;
  }

}
  

bool
EngineDB::insert_feature_vec(const std::string &fv_id, 
                             const std::string &path) {
  std::cout << fv_id << std::endl;
  mysqlpp::Connection conn(false);
  if(conn.connect(db.c_str(), server.c_str(), user.c_str(), pass.c_str())) {
    string query_str = "insert into feature_vector values (\"" + fv_id + "\",\"" + path + "\")";
    cout << query_str << endl;
    mysqlpp::Query query = conn.query(query_str.c_str());
    if (mysqlpp::StoreQueryResult res = query.store()) {
      std::cout << "We have:" << std::endl;
      for (size_t i = 0; i < res.num_rows(); ++i) {
        std::cout << '\t' << res[i][0] << endl;
      }
    }
    else {
      std::cerr << "Failed to get item list: " << query.error() << std::endl;
      return false;
    }
    return true;
  }
  else {
    std::cerr << "DB connection failed: " << conn.error() << std::endl;
    return false;
  }
}
