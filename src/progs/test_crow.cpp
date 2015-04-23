// ex1.cpp
#include "crow_all.hpp"

int main()
{
  crow::SimpleApp app;

  CROW_ROUTE(app, "/")
  ([]{
    return "Hello, world!";
  });

  app
    .port(18080)
    .run();
}
