#include <iostream>
#include <iomanip>
#include <string>
#include "EngineDB.hpp"

using namespace std;

int main() {
  string db = "amordad";
  string server = "localhost";
  string user = "root";
  string pass = "580230mysql";

  EngineDB eng(db,server,user,pass);
  eng.delete_feature_vec("haah");
}
