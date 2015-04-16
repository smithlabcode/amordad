/*    Part of AMORDAD software
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

#include "LSHAngleHashTable.hpp"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <limits>
#include <cassert>

#include "smithlab_utils.hpp"

#include "FeatureVector.hpp"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

size_t
LSHAngleHashTable::max_bucket_load() const {
  size_t max_load = 0;
  for (BucketMap::const_iterator i(buckets.begin()); i != buckets.end(); ++i)
    max_load = std::max(max_load, i->second.size());
  return max_load;
}


void
LSHAngleHashTable::insert(const FeatureVector &fv, const size_t hash_key) {
  const BucketMap::iterator x(buckets.find(hash_key));
  if (x == buckets.end())
    buckets[hash_key] = vector<string>(1, fv.get_id());
  else 
    x->second.push_back(fv.get_id());
}


void
LSHAngleHashTable::insert(const string &fv_id, const size_t hash_key) {
  const BucketMap::iterator x(buckets.find(hash_key));
  if (x == buckets.end())
    buckets[hash_key] = vector<string>(1, fv_id);
  else 
    x->second.push_back(fv_id);
}


void
LSHAngleHashTable::remove(const FeatureVector &fv, const size_t hash_key) {
  const BucketMap::iterator x(buckets.find(hash_key));
  if (x == buckets.end())
    throw SMITHLABException("attempt to remove from unkonwn hash key: " 
                            + toa(hash_key));
  else {
    // locate fv in the hashed bucket
    vector<string>::iterator pos = std::find(x->second.begin(),
                                             x->second.end(),
                                             fv.get_id());
    if (pos == x->second.end()) 
      throw SMITHLABException("attempt to remove unknown point: " 
                              + fv.get_id());
    else {
      x->second.erase(pos);
      if(x->second.empty())
        buckets.erase(hash_key);
    }
  }
}

string
LSHAngleHashTable::tostring() const {
  std::ostringstream oss;
  oss << id;
  for (BucketMap::const_iterator i(buckets.begin()); i != buckets.end(); ++i) {
    oss << '\n' << i->first << '\t';
    std::copy(i->second.begin(), i->second.end(),
              std::ostream_iterator<string>(oss, "\t"));
  }
  return oss.str();
}


std::ostream&
operator<<(std::ostream &out, const LSHAngleHashTable &hash_table) {
  return out << hash_table.tostring();
}


static bool
parse_bucket_line(const string &buffer, pair<size_t, vector<string> > &bucket) {
  std::istringstream iss(buffer);
  
  size_t hash_key = 0;
  if (!(iss >> hash_key))
    return false;
  
  vector<string> values;
  string val;
  while (iss >> val)
    values.push_back(val);
  
  bucket = make_pair(hash_key, values);
  return true;
}


std::istream&
operator>>(std::istream &is, LSHAngleHashTable &hash_table) {
  // read the hash table id (same as hash function id)
  string id;
  getline(is, id);
  
  // now read each bucket (line by line)
  BucketMap bm;
  pair<size_t, vector<string> > current_bucket;
  
  string buffer;
  while (getline(is, buffer)) {
    
    if (!parse_bucket_line(buffer, current_bucket))
      throw SMITHLABException("bad hash bucket line: " + buffer);
    
    bm.insert(current_bucket);
  }
  hash_table = LSHAngleHashTable(id, bm);
  return is;
}


/*
 * In many applications, a bucket containing only a single metagenome
 * does not matter. This function is a copycat of operator>> except
 * for single element buckets which are ignored here.
 */
void
read_non_single_buckets(std::istream &is, LSHAngleHashTable &hash_table) {
  // read the hash table id (same as hash function id)
  string id;
  getline(is, id);
  
  // now read each bucket (line by line)
  BucketMap bm;
  pair<size_t, vector<string> > current_bucket;
  
  string buffer;
  while (getline(is, buffer)) {
    
    if (!parse_bucket_line(buffer, current_bucket))
      throw SMITHLABException("bad hash bucket line: " + buffer);
    
    if (current_bucket.second.size() > 1)
      bm.insert(current_bucket);
  }
  hash_table = LSHAngleHashTable(id, bm);
}
