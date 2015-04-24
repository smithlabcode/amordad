/*
 *    Part of AMORDAD software
 *
 *    Copyright (C) 2014 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
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

#ifndef LSHANGLEHASHTABLE_HPP
#define LSHANGLEHASHTABLE_HPP

#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>

class FeatureVector;

typedef std::unordered_map<size_t, std::vector<std::string> > BucketMap;

class LSHAngleHashTable {
public:

  // Constructors
  LSHAngleHashTable() {}
  LSHAngleHashTable(const std::string id_in) : id(id_in) {}
  LSHAngleHashTable(const std::string &i, const BucketMap &bm) : 
    id(i), buckets(bm) {}

  // Accessors
  std::string tostring() const;
  std::string get_id() const {return id;}
  size_t size() const {return buckets.size();}
  size_t max_bucket_load() const;
  bool good() const {return !buckets.empty() && !id.empty();}
  
  BucketMap::const_iterator begin() const {return buckets.begin();}
  BucketMap::const_iterator end() const {return buckets.end();}
  BucketMap::const_iterator find(const size_t &x) const {return buckets.find(x);}
  
  // Mutators
  void insert(const FeatureVector &fv, const size_t hash_value);
  void insert(const std::string &fv_id, const size_t hash_value);
  void remove(const FeatureVector &fv, const size_t hash_value);
  
private:
  std::string id;
  BucketMap buckets;
};


std::ostream&
operator<<(std::ostream &out, const LSHAngleHashTable &hash_table);

std::istream&
operator>>(std::istream &in, LSHAngleHashTable &hash_table);

void
read_non_single_buckets(std::istream &in, LSHAngleHashTable &hash_table);

#endif
