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
 */

#include "GraphTableContainer.hpp"

#include "FeatureVector.hpp"
#include "LSHAngleHashFunction.hpp"
#include "LSHAngleHashTable.hpp"
#include "RegularNearestNeighborGraph.hpp"

using std::string;
using std::vector;
using std::cerr;
using std::endl;
using std::cout;
using std::pair;
using std::make_pair;

using std::tr1::unordered_map;
using std::tr1::unordered_set;

typedef LSHAngleHashTable LSHTab;
typedef LSHAngleHashFunction LSHFun;


std::ostream &
operator<<(std::ostream &os, const Result &r) {
  return os << r.id << '\t' << r.val;
}

void
GraphTableContainer::collect_candidates_from_hts(const FeatureVector &query,
                                                 unordered_set<string> 
                                                 &candidates) const {

  // iterate over hash tables
  for (unordered_map<string, LSHTab>::const_iterator i(hts.begin());
      i != hts.end(); ++i) {

    unordered_map<string, LSHFun>::const_iterator hf(hfs.find(i->first));
    assert(hf != hfs.end());

    // hash the query
    const size_t bucket_number = hf->second(query);
    unordered_map<size_t, vector<string> >::const_iterator
      bucket = i->second.find(bucket_number);

    if (bucket != i->second.end())
      candidates.insert(bucket->second.begin(), bucket->second.end());
  }
}

void
GraphTableContainer::collect_candidates_from_hts_and_nng(const 
                                                         FeatureVector &query,
                                                         unordered_set<string> 
                                                         &candidates) const {

  collect_candidates_from_hts(query, candidates);

  // gather neighbors of candidates from graph
  unordered_set<string> candidates_from_graph;
  for (unordered_set<string>::const_iterator i(candidates.begin());
      i != candidates.end(); ++i) {
    vector<string> neighbors;
    vector<double> neighbor_dists;
    nng.get_neighbors(*i, neighbors, neighbor_dists);
    candidates_from_graph.insert(neighbors.begin(), neighbors.end());
  }

  candidates.insert(candidates_from_graph.begin(), candidates_from_graph.end());
}

void
GraphTableContainer::find_nearest_neighbors(const FeatureVector &query, 
                                            const size_t n_neighbors,
                                            const double max_proximity_radius,
                                            vector<Result> &results) const {

  unordered_set<string> candidates;
  collect_candidates_from_hts_and_nng(query, candidates);

  std::priority_queue<Result, vector<Result>, std::less<Result> > pq;
  double current_dist_cutoff = max_proximity_radius;
  for (unordered_set<string>::const_iterator i(candidates.begin());
      i != candidates.end(); ++i) {
    const FeatureVector fv(fvs.find(*i)->second);
    const double dist = query.compute_angle(fv);
    if (dist < current_dist_cutoff) {
      if (pq.size() == n_neighbors) 
        pq.pop();
      pq.push(Result(*i, dist));

      if(pq.size() == n_neighbors)
        current_dist_cutoff = pq.top().val;
    }
  }

  results.clear();
  while (!pq.empty()) {
    results.push_back(pq.top());
    pq.pop();
  }
  reverse(results.begin(), results.end());
}
