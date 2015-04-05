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

#include </usr/include/mysql++/mysql++.h>

#include <string>
#include <vector>


class EngineDB {
public:
  EngineDB() {}
  EngineDB(const std::string db, const std::string server, 
           const std::string user, const std::string pass) :
    db(db), server(server), user(user), pass(pass) {}

  bool delete_feature_vec(const std::string &fv_id);

private:
  std::string db;
  std::string server;
  std::string user;
  std::string pass;
};

#endif
