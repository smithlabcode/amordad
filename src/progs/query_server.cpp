/*
 *    Part of SMITHLAB software
 *
 *    Copyright (C) 2014 University of Southern California and
 *                       Andrew D. Smith
 *                       Saket Choudhary
 *
 *    Authors: Andrew D. Smith and Saket Choudhary
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


#include "crow_all.h"
#include "naive_batch_query.hpp"


#include <sstream>

/*
 * Load all database files in memory for once and all
 *
 */
vector<FeatureVector> init_amorad(){

    bool VERBOSE = true;
    string database_filepath = "/media/data1/Development_Version_Controlled/Research/Amordad/src/progs/mgrast_k5_paths.txt";
    const string database_file(database_filepath);
    if (!validate_file(database_file, 'r'))
      throw SMITHLABException("bad database file: " + database_file);
    vector<FeatureVector> database;
    if (VERBOSE)
      cerr << "loading database" << endl;
    load_feature_vectors(VERBOSE, database_file, database);
    if (VERBOSE)
      cerr << "database size: " << database.size() << endl;
    return database;
}

int run_naive_batch_query(vector<FeatureVector> database, string query_folder){

    bool VERBOSE = true;
    //Create a output file with the same names as the query_folder appending a '_output' to it.
    string out_filepath = query_folder+"_output";
    //out_filepath.replace(query_folder.end()-1, query_folder.end(),"_output");
    size_t n_neighbors = 1;
    double max_proximity_radius = 0.75;
    cout << "QUERY_FOLDEr" << query_folder << endl;
    const string queries_file(query_folder);
    const string outfile(out_filepath);


    ////////////////////////////////////////////////////////////////////////
    ////// READING DATABASE ////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////
    ///// STARTING THE QUERY PROCESS ///////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    // reading queries
    vector<FeatureVector> queries;
    if (VERBOSE)
      cerr << "loading queries" << endl;
    load_feature_vectors(VERBOSE, queries_file, queries);
    if (VERBOSE)
      cerr << "number of queries: " << queries.size() << endl;

    // "n" query points requires a "n*t" results
    vector<vector<Result> > results(queries.size());
    for (size_t i = 0; i < queries.size(); ++i) {
      exec_query(database, queries[i],
                 n_neighbors, max_proximity_radius, results[i]);
      if (VERBOSE)
        cerr << '\r' << "processing queries: "
             << percent(i, queries.size()) << "%\r";
    }
    if (VERBOSE)
      cerr << '\r' << "processing queries: 100% ("
           << queries.size() << ")" << endl;

    ////////////////////////////////////////////////////////////////////////
    ///// NOW WRITE THE OUTPUT /////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    std::ofstream out(outfile.c_str());
    if (!out)
      throw SMITHLABException("bad output file: " + outfile);
    for (size_t i = 0; i < queries.size(); ++i) {
      out << queries[i].get_id() << '\t';
      copy(results[i].begin(), results[i].end(),
           std::ostream_iterator<Result>(out, "\t"));
      out << endl;
    }

    if (VERBOSE)
      cerr << comparisons << endl;
    return 0;
}


int main()
{
    crow::SimpleApp app;
    vector<FeatureVector> database = init_amorad();
    CROW_ROUTE(app, "/about")
    ([](){
        return "About Crow example.";
    });

    // simple json response
    CROW_ROUTE(app, "/json")
    ([]{
        crow::json::wvalue x;
        x["message"] = "Hello, World!";
        return x;
    });

    // argument
    CROW_ROUTE(app,"/naivequery/<string>")
    ([&](string folderpath){
        folderpath = "/media/data1/Development_Version_Controlled/Research/Amordad/src/progs/"+folderpath;
        return run_naive_batch_query(database, folderpath);
    });

    // Compile error with message "Handler type is mismatched with URL paramters"
    //CROW_ROUTE(app,"/another/<int>")
    //([](int a, int b){
        //return crow::response(500);
    //});

    // more json example
    CROW_ROUTE(app, "/add_json")
    ([](const crow::request& req){
        auto x = crow::json::load(req.body);
        if (!x)
            return crow::response(400);
        int sum = x["a"].i()+x["b"].i();
        std::ostringstream os;
        os << sum;
        return crow::response{os.str()};
    });

    app.port(18080)
        .multithreaded()
        .run();
}
