/*
 *    Part of AMORDAD software
 *
 *    Copyright (C) 2014 University of Southern California and
 *                       Andrew D. Smith
 *                       Ehsan Behnam
 *
 *    Authors: Andrew D. Smith and Ehsan Behnav
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
#include "naive_batch_query.hpp"
int
main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    size_t n_neighbors = 1;
    double max_proximity_radius = 0.75;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "batch query a set of "
                           "feature vectors",
                           "<database-file> <query-folder-path> <outfile>");
    opt_parse.add_opt("neighbors", 'n', "number of nearest neighbors to report",
                      false, n_neighbors);
    opt_parse.add_opt("mpr", 'r', "maximum proximity radius",
                      false, max_proximity_radius);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() != 3) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string database_file(leftover_args.front());
    const string queries_file(leftover_args[1]);
    const string outfile(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    if (!validate_file(queries_file, 'r'))
      throw SMITHLABException("bad queries file: " + queries_file);

    if (!validate_file(database_file, 'r'))
      throw SMITHLABException("bad database file: " + database_file);

    if (!validate_file(outfile, 'w'))
      throw SMITHLABException("bad output file: " + outfile);

    ////////////////////////////////////////////////////////////////////////
    ////// READING DATABASE ////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    // reading database
    vector<FeatureVector> database;
    if (VERBOSE)
      cerr << "loading database" << endl;
    load_feature_vectors(VERBOSE, database_file, database);
    if (VERBOSE)
      cerr << "database size: " << database.size() << endl;

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
  }
  catch (const SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
