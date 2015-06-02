/*
 *    Part of AMORDAD software
 *
 *    Copyright (C) 2014 University of Southern California and
 *                       Andrew D. Smith
 *                       Ehsan Behnam
 *
 *    Authors: Andrew D. Smith and Ehsan Behnam
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
#include "FeatureVector.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <climits>
#include <numeric>

using std::string;
using std::vector;
using std::endl;
using std::inner_product;
using std::sqrt;
using std::acos;
using std::copy;

////////////////////////////////////////////////////////////////////////
///////   FUNCTIONS FOR I/O

std::ostream &
operator<<(std::ostream &out, const FeatureVector &fv) {
  return out << fv.tostring();
}


string
FeatureVector::tostring() const {
  std::ostringstream oss;
  oss << id << endl;
  copy(values.begin(), values.end(), 
       std::ostream_iterator<double>(oss, "\n"));
  return oss.str();
}


string
FeatureVector::tostring_with_labels(const vector<string> &labels) const {
  if (values.size() != labels.size())
    throw SMITHLABException("feature vector values size is not "
        "equal to labels size");

  std::ostringstream oss;
  oss << id;
  for (size_t i = 0; i < labels.size(); ++i)
    oss << endl << labels[i] << '\t' << values[i];
  return oss.str();
}


std::istream&
operator>>(std::istream &in, FeatureVector &fv) {
  string id;
  if(!getline(in, id))
    throw SMITHLABException("bad feature vector: empty file");
  
  string label, value, line;
  vector<double> values;
  while (getline(in, line)) {
    std::istringstream iss(line);
    if (!(iss >> label >> value))
      throw SMITHLABException("bad feature vector line: " + line);
    values.push_back(strtod(value.c_str(), 0));
  }

  if(values.empty())
    throw SMITHLABException("bad feature vector: only with ID");

  fv = FeatureVector(id, values);
  return in;
}


void
load_features_and_labels(const string &filename,
                         FeatureVector &fv, vector<string> &labels) {
  
  std::ifstream in(filename.c_str());
  if (!in)
    throw SMITHLABException("bad file: " + filename);
  
  string id;
  getline(in, id);
  
  string label, value, line;
  vector<double> values;
  while (getline(in, line)) {
    std::istringstream iss(line);
    if (!(iss >> label >> value))
      throw SMITHLABException("bad feature vector line: " + line);
    values.push_back(strtod(value.c_str(), 0));
    labels.push_back(label);
  }
  fv = FeatureVector(id, values);
}


double
FeatureVector::compute_angle(const FeatureVector &other) const {
  if (values.size() != other.values.size())
    throw SMITHLABException("cannot compute angle: different feature"
        "vector size: (" + id + ',' + other.get_id() + ")");
  double angle = inner_product(values.begin(), values.end(), 
                               other.values.begin(), 0.0)/(norm*other.norm);
  angle = std::max(-1.0, std::min(1.0, angle));
  assert(abs(angle) <= 1);
  return acos(angle); // scale factor 180/M_PI for "degree" conversion
}
