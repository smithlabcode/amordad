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

  // test get_oldest_hash_function
  string oldest_hf = eng.get_oldest_hash_function();
  if(oldest_hf != "")
    cout << oldest_hf << endl;
  else
    cout << "No hash function in the database" << endl;
}
