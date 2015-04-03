#include </usr/include/mysql++/mysql++.h>
#include "/usr/include/mysql++/cmdline.h"
#include <iostream>
#include <iomanip>

using namespace std;

int main() {
  const char* db = "amordad";
  const char* server = "localhost";
  const char* user = "root";
  const char* pass = "580230mysql";
  mysqlpp::Connection conn(false);
  if (conn.connect(db, server, user, pass)) {
    mysqlpp::Query query = conn.query("select * from sample");
    if (mysqlpp::StoreQueryResult res = query.store()) {
      cout << "We have:" << endl;
      for (size_t i = 0; i < res.num_rows(); ++i) {
        cout << '\t' << res[i][0] << endl;
      }
    }
    else {
      cerr << "Failed to get item list: " << query.error() << endl;
      return 1;
    }
    return 0;
  }
  else {
    cerr << "DB connection failed: " << conn.error() << endl;
    return 1;
  }
}
