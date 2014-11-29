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

#ifndef LSHANGLEHASHFUNCTION_HPP
#define LSHANGLEHASHFUNCTION_HPP

#include <string>
#include <vector>

class FeatureVector;

class LSHAngleHashFunction {
public:
  LSHAngleHashFunction() {}
  LSHAngleHashFunction(const std::string &id_in, const std::string &fsi,
                       const size_t n_features, const size_t n_bits);
  LSHAngleHashFunction(const std::string &id_in, const std::string &fsi,
                       const std::vector<std::vector<double> > &uvs) :
    id(id_in), feature_set_id(fsi), unit_vecs(uvs) {}
  
  size_t operator()(const FeatureVector &fv) const;
  
  std::string tostring() const;
  size_t size() const {return unit_vecs.size();};
  std::string get_id() const {return id;}
  std::string get_feature_set_id() const {return feature_set_id;}
  
private:
  std::string id;
  std::string feature_set_id;
  std::vector<std::vector<double> > unit_vecs;
};

std::ostream&
operator<<(std::ostream &os, const LSHAngleHashFunction &hf);

std::istream&
operator>>(std::istream &in, LSHAngleHashFunction &hf);

#endif
