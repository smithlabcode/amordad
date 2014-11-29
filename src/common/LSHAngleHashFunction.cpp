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

#include "LSHAngleHashFunction.hpp"
#include "FeatureVector.hpp"

using std::string;
using std::vector;
using std::pair;


static double
get_gaussian(const double mu, const double sd) {
  const double u1 = static_cast<double>((rand() % RAND_MAX)) / RAND_MAX;
  const double u2 = static_cast<double>((rand() % RAND_MAX)) / RAND_MAX;
  return (((rand() % 2) == 0) ?
          sd*sqrt(-2 * log(u1)) * cos(2 * M_PI * u2) :
          sd*sqrt(-2 * log(u1)) * sin(2 * M_PI * u2)) + mu;
}


static void 
generate_random_unit_vec(const size_t dim, vector<double> &v) {
  //if x1, x2, ..., xn are normal(0, 1)
  //then normalized (x1, x2, .., xn) is a random point on unit hypersphere.
  double r = 0.0;
  for (size_t i = 0; i < dim; ++i) {
    const double n = get_gaussian(0.0, 1.0);
    v[i] = n;
    r += n*n;
  }
  for (vector<double>::iterator i = v.begin(); i != v.end(); ++i)
    (*i) /= std::sqrt(r);
}


// CONSTRUCTORS
LSHAngleHashFunction::LSHAngleHashFunction(const string &id_in,
                                           const string &fsi,
                                           const size_t n_features,
                                           const size_t n_bits) :
  id(id_in), feature_set_id(fsi) {
  // generate the hyperplanes
  unit_vecs.resize(n_bits, vector<double>(n_features, 0.0));
  for (size_t i = 0; i < n_bits; ++i) {
    generate_random_unit_vec(n_features, unit_vecs[i]);
    assert(unit_vecs[i].size() == n_features);
  }
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


// INPUT
std::istream&
operator>>(std::istream &in, LSHAngleHashFunction &hf) {
  
  // first read the id for this hash function
  string hf_id;
  getline(in, hf_id);
  
  // second line is for the feature set
  string fs_id;
  getline(in, fs_id);
  
  size_t n_dimensions = 0;
  
  vector<vector<double> > uvs;
  string line;
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
    
    uvs.push_back(vector<double>());
    uvs.back().swap(current);
  }
  
  hf = LSHAngleHashFunction(hf_id, fs_id, uvs);
  return in;  
}


// OUTPUT
std::ostream&
operator<<(std::ostream &os, const LSHAngleHashFunction &hf) {
  return os << hf.tostring();
}


string
LSHAngleHashFunction::tostring() const {
  std::ostringstream oss;
  oss << id << '\n' << feature_set_id;
  for (size_t i = 0; i < unit_vecs.size(); ++i) {
    oss << '\n';
    copy(unit_vecs[i].begin(), unit_vecs[i].end(), 
         std::ostream_iterator<double>(oss, "\t"));
  }
  return oss.str();
}


size_t 
LSHAngleHashFunction::operator()(const FeatureVector &fv) const {
  size_t value = 0;
  for (size_t i = 0; i < unit_vecs.size(); ++i) {
    value <<= 1ul;
    value += (inner_product(unit_vecs[i].begin(),
                            unit_vecs[i].end(), fv.begin(), 0.0) >= 0);
  }
  return value;
}
