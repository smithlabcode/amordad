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
  string id, path;
  std::cout << "input:" << std::endl;
  while(cin >> id >> path) {
    eng.insert_feature_vec(id, path);
    std::cout << "input:" << std::endl;
  }
}
