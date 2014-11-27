/*
 *    Part of SMITHLAB software
 *
 *    Copyright (C) 2014 University of Southern California,
 *                       Andrew D. Smith and
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
#ifndef REGULAR_NEAREST_NEIGHBOR_GRAPH_HPP
#define REGULAR_NEAREST_NEIGHBOR_GRAPH_HPP

#include <fstream>
#include <string>
#include <vector>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <tr1/unordered_map>

typedef 
boost::adjacency_list<boost::setS, boost::vecS, 
                      boost::directedS,
                      boost::no_property,
                      boost::property<boost::edge_weight_t, double>
                      > internal_graph;

typedef
boost::graph_traits<internal_graph>::vertex_descriptor nng_vertex;

class RegularNearestNeighborGraph{
public:
  RegularNearestNeighborGraph() {}
  RegularNearestNeighborGraph(const std::string &gn, const size_t deg) : 
    graph_name(gn), maximum_degree(deg) {}
  RegularNearestNeighborGraph(const std::string &gn, 
                              const std::vector<std::string> &ids,
                              const size_t deg) : 
    graph_name(gn), maximum_degree(deg) {add_vertices(ids);}
  
  // graph accessors
  std::string get_graph_name() const {return graph_name;}
  size_t get_maximum_degree() const {return maximum_degree;}
  size_t get_vertex_count() const;
  size_t get_edge_count() const;

  /* vertex accessors */
  double get_distance(const nng_vertex &u, const nng_vertex &v) const;
  double get_distance(const std::string &u, const std::string &v) const;
  double get_distance(const boost::graph_traits<internal_graph>::edge_descriptor
                      &the_edge) const;
  void get_neighbors(const std::string &query,
                     std::vector<std::string> &neighbors,
                     std::vector<double> &distances) const;
  
  // mutators
  void add_edge(const nng_vertex &u, const nng_vertex &v, const double &w);
  void add_edge(const std::string &u, const std::string &v, const double &w);

  void add_vertex(const std::string &id);
  void add_vertices(const std::vector<std::string> &id_list);
  bool add_vertex_if_new(const std::string &id);
  
  bool update_vertex(const nng_vertex &u, const nng_vertex &v, const double &w);
  bool update_vertex(const std::string &u, const std::string &v, 
                     const double &w);
  
  std::string tostring() const;
  
private:

  std::string graph_name;
  internal_graph the_graph;
  std::tr1::unordered_map<std::string, size_t> name_to_index;
  std::tr1::unordered_map<size_t, std::string> index_to_name;
  size_t maximum_degree;

  // lookups
  size_t convert_name_to_index(const std::string &name) const;
  std::string convert_index_to_name(const size_t &index) const;

  void
  get_most_distant_neighbor(const nng_vertex &query,
                            nng_vertex &result, double &max_distance) const;
  
};

std::ostream&
operator<<(std::ostream &is, const RegularNearestNeighborGraph &g);

std::istream&
operator>>(std::istream &in, RegularNearestNeighborGraph &g);

#endif
