/*
 *    Part of SMITHLAB software
 *
 *    Copyright (C) 2013 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Ehsan Behnam, Andrew D. Smith
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

#ifndef FEATUREVECTOR_HPP
#define FEATUREVECTOR_HPP

#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <numeric>

class FeatureVector {
public:

  FeatureVector() {}
  FeatureVector(const std::string &id_in,
                const std::vector<double>& fv) : 
    id(id_in), values(fv), 
    norm(std::sqrt(std::inner_product(values.begin(), values.end(), 
                                      values.begin(), 0.0))) {}
  
  const double& operator[](const size_t i) const {return values[i];}
  double& operator[](const size_t i) {return values[i];}
  
  std::vector<double>::const_iterator begin() const {return values.begin();}
  std::vector<double>::const_iterator end() const {return values.end();}
  std::vector<double>::iterator begin() {return values.begin();}
  std::vector<double>::iterator end() {return values.end();}
  
  std::string get_id() const {return id;}
  size_t size() const {return values.size();}
  
  std::string tostring() const;
  
  std::string
  tostring_with_labels(const std::vector<std::string> &labels) const;
  
  // Functions for comparing two FeatureVectors
  double compute_angle(const FeatureVector &other) const;
  
private:
  std::string id;
  std::vector<double> values;
  double norm;
};


std::istream&
operator>>(std::istream &, FeatureVector &);

std::ostream&
operator<<(std::ostream &, const FeatureVector &);

void 
load_features_and_labels(const std::string &filename,
                         FeatureVector &fv,
                         std::vector<std::string> &labels);

#endif
