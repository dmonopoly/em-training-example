#ifndef NODE_H_
#define NODE_H_

#include <iostream>
#include <vector>

#include "Edge.h"

using namespace std;

struct Edge;

struct Node {
  string name;
  int index;  // Topological ordering index.
  vector<Edge *> parent_edges;
  vector<Edge *> child_edges;
  Node(string name, int index) {
    this->name = name;
    this->index = index;
  }
  string repr() {
    return this->name;
  }
  void set_name(const string &other) {
    this->name = other;
  }
};

#endif