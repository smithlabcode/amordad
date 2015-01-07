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

#define PORT 18080 // Port to run the server at
// Set this to absolute location of the database file
// TODO Document the database file format
#define DATABASE_FILEPATH "/media/data1/Development_Version_Controlled/Research/Amordad/src/progs/mgrast_k5_paths.txt"
// Set this to location where the PHP code creates the upload directories
#define UPLOAD_DIR_LOCATION "/media/data1/Development_Version_Controlled/Research/Amordad/src/progs/"
#include <sstream>

/*
 * Load all database files in memory for once and all
 * @return database A vector of Featurevector
 */
vector<FeatureVector> init_amorad(){

    bool VERBOSE = true;
    string database_filepath = DATABASE_FILEPATH;
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

/**
 * Walk the directory and create an example query file. 'query_paths.query'
 * This query_paths.query contains paths to all the query files.
 * These query files in turn are identified by using the read_dir method to return all
 * files with a '.cv' extension
 * @param query_dir Absolute path to the directory where *.cv files are located
 * @output outfile An explicitly 'query_paths.query' named file that contains nothing but
 *                 the absolute paths to all *.cv files inside query_dir
 * TODO Maybe this is an overkill?
 **/
string walk_dir_and_create_query_file(string query_dir){
    static const string filename_suffix = "cv";
    string outfile = query_dir + "/" + "query_paths.query";
    vector<string> query_files;
    read_dir(query_dir, filename_suffix, query_files);
    std::ofstream out(outfile.c_str());
    if (!out)
        throw SMITHLABException("bad output file: " + outfile);
    for (size_t i = 0; i < query_files.size(); ++i) {
        out << query_files[i].c_str();
        out << endl;
    }
    return outfile;
}

/**
 * Method to run naive_batch_query
 * @param database A vector containing the featurevectors coming from the amordad database
 * @param query_dir The absolute location of directory where all query files are located
 * @output exit
 * NOTE: Ideally this step should not be modified, and will be consistent throughout
 * for all other query methods
 **/
int run_naive_batch_query(vector<FeatureVector> database, string query_dir){

    bool VERBOSE = true; //Set it to true so everything shows up in the server logs
    /**
     * TODO: these should ideally come from the query itself?
     * E.g http://localhost/naivequery/test?n_neighbors=1*max_proximity_radius=0.9
     **/
    size_t n_neighbors = 1;
    double max_proximity_radius = 0.75;
    // Inside the query directory create a ;query_filepath.query' file
    // @see walk_dir_and_create_query_file()
    string query_filepath = walk_dir_and_create_query_file(query_dir);
    //Create a output file with the same names as the query_filepath appending a '.output' to it.
    string out_filepath = query_filepath+".output";
    const string queries_file(query_filepath);
    const string outfile(out_filepath);


    // Read Queries
    vector<FeatureVector> queries;
    if (VERBOSE)
        cerr << "Loading queries" << endl;
    load_feature_vectors(VERBOSE, queries_file, queries);
    if (VERBOSE)
      cerr << "Number of queries: " << queries.size() << endl;

    // "n" query points requires a "n*t" results
    vector<vector<Result> > results(queries.size());
    for (size_t i = 0; i < queries.size(); ++i) {
      exec_query(database, queries[i],
                 n_neighbors,
                 max_proximity_radius,
                 results[i]);
      if (VERBOSE)
        cerr << '\r' << "processing queries: "
             << percent(i, queries.size()) << "%\r";
    }
    if (VERBOSE)
      cerr << '\r' << "processing queries: 100% ("
           << queries.size() << ")" << endl;

    /**
     * Output Results
     *
     **/
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
    /**
     * Load database in memory
     * This runs at the start of server for once and all
     **/
    vector<FeatureVector> database = init_amorad();


    /**
     * The section below takes care of the url mappings.
     * Adding a new query method is as easy as creating just a new mapping
     *
     **/

    // Homepage
    CROW_ROUTE(app, "/")
    ([](){
        return "Amordad Web Server";
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
    ([&](string dirpath){
        dirpath = "/media/data1/Development_Version_Controlled/Research/Amordad/src/progs/"+dirpath;
        return run_naive_batch_query(database, dirpath);
    });


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


    // Main app run
    app.port(PORT)
        .multithreaded()
        .run();
}
