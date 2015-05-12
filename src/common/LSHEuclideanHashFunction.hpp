/*
 *    Part of AMORDAD software
 *
 *    Copyright (C) 2014 University of Southern California and
 *                       Andrew D. Smith, Wenzheng Li
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

#ifndef LSHEUCLIDEANHASHFUNCTION_HPP
#define LSHEUCLIDEANHASHFUNCTION_HPP

#include <string>
#include <vector>

class FeatureVector;

struct Parameter {
  // random vector a entries chosen from Guassian dist
  std::vector<double> rand_vec;
  // random value b chosen uniformly from [0,w]
  double rand_uniform;
};

class LSHEuclideanHashFunction {
public:
  LSHEuclideanHashFunction() {}
  LSHEuclideanHashFunction(const std::string &id_in, const std::string &fsi,
                       const size_t n_features, const size_t n_bits,
                       const double w);
  LSHEuclideanHashFunction(const std::string &id_in, const std::string &fsi,
                       const std::vector<Parameter> &ps, const double w,
                       const std::vector<double> avs) :
    id(id_in), feature_set_id(fsi), parameters(ps),
    uniform_seed(w), assist_vals(avs) {}
 
  size_t operator()(const FeatureVector &fv) const;
 
  std::string tostring() const;
  size_t size() const {return parameters.size();};
  std::string get_id() const {return id;}
  std::string get_feature_set_id() const {return feature_set_id;}
 
private:
  std::string id;
  std::string feature_set_id;
  std::vector<Parameter> parameters;
  // parameter w
  double uniform_seed;
  // n_features-d vector assisting to generate final hash value
  std::vector<double> assist_vals;
};

std::ostream&
operator<<(std::ostream &os, const LSHEuclideanHashFunction &hf);

std::istream&
operator>>(std::istream &in, LSHEuclideanHashFunction &hf);

#endif
