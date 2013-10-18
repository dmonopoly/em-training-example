// not used yet, in progress
#ifndef TRELLIS_H_
#define TRELLIS_H_

#include <cassert>
#include <cstdlib>
#include <map>
#include <sstream>
#include <vector>

#include "BasicHelper.h"
#include "Notation.h"
#include "NotationConstants.h"
#include "Node.h"
#include "Edge.h"

using namespace std;

struct Node;
struct Edge;

namespace TrellisAid {

  // WARNING: Creates data on heap. Call DestroyTrellis when done.  Post:
  // 'nodes' points to a vector where front() is the start node, back() is the
  // end, and the vector lists the nodes in topological order. Each node has a
  // unique name. 'edges' points to a vector of corresponding edges with
  // representations like P(A|X). **Nodes and edges are in topological order!**.
  void BuildTrellis(vector<Node *> *nodes, vector<Edge *> *select_edges,
                    vector<Edge *> *all_edges, const vector<string>
                    &observed_data, const vector<string> &tag_list);
  void DestroyTrellis(vector<Node *> *nodes, vector<Edge *> *all_edges);
  // TODO
//   void Viterbi(const map<string, double> &data, const vector<Node *> &nodes);
//   void ForwardBackwardAndViterbi(Notation n, const vector<Node *> &nodes,
//                                  const vector<Edge *> &select_edges,
//                                  const vector<Edge *> &all_edges,
//                                  map<string, double> *data);
}

#endif  // TRELLIS_H_
