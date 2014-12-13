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
#ifndef GRAPH_TABLE_CONTAINER_HPP
#define GRAPH_TABLE_CONTAINER_HPP

#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include <tr1/unordered_map>
#include <tr1/unordered_set>
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

class FeatureVector;
class LSHAngleHashFunction;
class LSHAngleHashTable;
class RegularNearestNeighborGraph;

class GraphTableContainer {
public:
  GraphTableContainer(const std::tr1::unordered_map<std::string, FeatureVector> &fvs, 
                      const std::tr1::unordered_map<std::string, LSHAngleHashFunction> &hfs,
                      const std::tr1::unordered_map<std::string, LSHAngleHashTable> &hts,
                      const RegularNearestNeighborGraph &nng) :
    fvs(fvs), hfs(hfs), hts(hts), nng(nng) {}

  void collect_candidates_from_hts(const FeatureVector &query, 
                                   std::tr1::unordered_set<std::string> &candidates) const;
  void collect_candidates_from_hts_and_nng(const FeatureVector &query,
                                   std::tr1::unordered_set<std::string> &candidates) const;

  void find_nearest_neighbors(const FeatureVector &query, 
                              const size_t n_neighbors,
                              const double max_proximity_radius,
                              std::vector<Result> &results) const;

private:
  const std::tr1::unordered_map<std::string, FeatureVector> &fvs;
  const std::tr1::unordered_map<std::string, LSHAngleHashFunction> &hfs;
  const std::tr1::unordered_map<std::string, LSHAngleHashTable> &hts;
  const RegularNearestNeighborGraph &nng;
};

#endif
