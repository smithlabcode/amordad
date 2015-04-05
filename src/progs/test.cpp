#include <iostream>
#include <iomanip>
#include <string>
#include "EngineDB.hpp"

using namespace std;

int main() {
  string db = "amorgin";
  string server = "localhost";
  string user = "root";
  string pass = "580230mysql";

  EngineDB eng(db,server,user,pass);
  eng.insert_feature_vec("test_fv_1", "path1");
}
