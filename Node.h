#ifndef NODE_H_
#define NODE_H_

#include <iostream>
#include <vector>

#include "Edge.h"

using namespace std;

struct Edge;

struct Node {
  string name;  // Useful name accessible in repr(), used as key in maps.
  int index;  // Topological ordering index. Technically not used since nodes
              // are consistently stored in vectors...
  vector<Edge *> parent_edges;
  vector<Edge *> child_edges;
  Node(string name, int index) {
    this->name = name;
    this->index = index;
  }
  string repr() {
    return this->name;
  }
};

#endif
