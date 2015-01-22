/*
 *    Part of AMORDAD software
 *
 *    Copyright (C) 2014 University of Southern California,
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

#include <string>
#include <vector>
#include <sstream>
#include <iterator>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <numeric>

#include "smithlab_utils.hpp"

#include "LSHEuclideanHashFunction.hpp"
#include "FeatureVector.hpp"

using std::string;
using std::vector;
using std::pair;


static double
get_uniform(const double start, const double end) {
  assert(start <= end);
  const double r = static_cast<double>((rand() % RAND_MAX)) / RAND_MAX;
  const double u = start + ((end - start) * r);
  assert(u >= start && u <= end);
  return u;
}


static double
get_gaussian(const double mu, const double sd) {
  const double u1 = static_cast<double>((rand() % RAND_MAX)) / RAND_MAX;
  const double u2 = static_cast<double>((rand() % RAND_MAX)) / RAND_MAX;
  return (((rand() % 2) == 0) ?
          sd*sqrt(-2 * log(u1)) * cos(2 * M_PI * u2) :
          sd*sqrt(-2 * log(u1)) * sin(2 * M_PI * u2)) + mu;
}


static void 
generate_random_parameter(const size_t dim, const double w, Parameter &ep) {
  ep.rand_vec.resize(dim, 0.0);
  for (size_t i = 0; i < dim; ++i)
    ep.rand_vec[i] = get_gaussian(0.0, 1.0);
  ep.rand_uniform = get_uniform(0.0, w);
}


// CONSTRUCTORS
LSHEuclideanHashFunction::LSHEuclideanHashFunction(const string &id_in,
                                                   const string &fsi,
                                                   const size_t n_features,
                                                   const size_t n_bits,
                                                   const double w) :
  id(id_in), feature_set_id(fsi), uniform_seed(w) {
  parameters.resize(n_bits);
  for (size_t i = 0; i < n_bits; ++i) {
    generate_random_parameter(n_features, w, parameters[i]);
    assert(parameters[i].rand_vec.size() == n_features);
  }
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


// INPUT
std::istream&
operator>>(std::istream &in, LSHEuclideanHashFunction &hf) {
  
  // first read the id for this hash function
  string hf_id;
  getline(in, hf_id);
  
  // second line is for the feature set
  string fs_id;
  getline(in, fs_id);
  
  // third line is for the uniform seed
  string line;
  getline(in, line);
  std::istringstream w_iss(line);
  double w = 0.0;
  w_iss >> w;

  size_t n_dimensions = 0;
  
  vector<Parameter> parameters;
  while (getline(in, line)) {
    
    vector<double> current(n_dimensions);
    
    std::istringstream iss(line);
    
    double single_coord = 0.0;
    
    if (n_dimensions == 0) {
      // dimensions unknown for first unit vector
      while (iss >> single_coord)
        current.push_back(single_coord);
      n_dimensions = current.size();
    }
    else {
      size_t i = 0;
      while (iss >> single_coord)
        current[i++] = single_coord;
    }
    
    if (n_dimensions != current.size())
      throw SMITHLABException("inconsistent hash function lines");
    
    Parameter ep;
    ep.rand_vec.assign(current.begin(), current.end()-1);
    ep.rand_uniform = current.back();
    parameters.push_back(ep);
  }
  
  hf = LSHEuclideanHashFunction(hf_id, fs_id, parameters, w);
  return in;  
}


// OUTPUT
std::ostream&
operator<<(std::ostream &os, const LSHEuclideanHashFunction &hf) {
  return os << hf.tostring();
}


string
LSHEuclideanHashFunction::tostring() const {
  std::ostringstream oss;
  oss << id << '\n' << feature_set_id;
  oss << '\n' << uniform_seed;
  for (size_t i = 0; i < parameters.size(); ++i) {
    oss << '\n';
    copy(parameters[i].rand_vec.begin(), 
         parameters[i].rand_vec.end(), 
         std::ostream_iterator<double>(oss, "\t"));
    oss << parameters[i].rand_uniform << "\t";
  }
  return oss.str();
}


size_t 
LSHEuclideanHashFunction::operator()(const FeatureVector &fv) const {
  size_t value = 0;
  for (size_t i = 0; i < parameters.size(); ++i) {
    double inner = inner_product(parameters[i].rand_vec.begin(),
                                 parameters[i].rand_vec.end(),
                                 fv.begin(), 0.0);
    const double hash_value = floor((inner + parameters[i].rand_uniform)
                                    / uniform_seed);
  }
  return value;
}
