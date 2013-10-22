#ifndef TRELLIS_H_
#define TRELLIS_H_

#include <cassert>
#include <cstdlib>
#include <map>
#include <unordered_map>
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
  // WARNING: Creates data on heap. Call DestroyTrellis when done. Post:
  // 'nodes' points to a vector where front() is the start node, back() is the
  // end, and the vector lists the nodes in topological order. Each node has a
  // unique name. 'edges' points to a vector of corresponding edges with
  // representations like P(A|X). **Nodes and edges are in topological order!**.
  void BuildTrellis(vector<Node *> *nodes, vector<Edge *> *select_edges,
                    vector<Edge *> *all_edges, const vector<string>
                    &observed_data, const vector<string> &tag_list);
  void DestroyTrellis(vector<Node *> *nodes, vector<Edge *> *all_edges);

  // Traverses the trellis with the Viterbi algorithm to find the best matching
  // tag sequence and prints the results. If very_small_data_set is true, we
  // also update saved_obs_seq_probs with newer P(observation sequences).
  void Viterbi(const map<string, double> &data, const vector<Node *> &nodes,
               const vector<string> observed_data, bool very_small_data_set,
               const vector<double> saved_obs_seq_probs);
  // Runs ForwardBackward 'num_iterations' times to determine best probabilities
  // and then calls Viterbi(). Updates data's count keys (e.g., C(X,A)) in the
  // process. If very_small_data_set is true, we also print organized data rows
  // used all parameters that follow this bool. Otherwise we don't care for any
  // of the params that follow this bool. 'select_edges' are the edges like
  // P(A|X) that we want to update.
  void ForwardBackwardAndViterbi(const int num_iterations,
                                 const vector<Node *> &nodes,
                                 const vector<Edge *> &select_edges,
                                 const vector<Edge *> &all_edges,
                                 map<string, double> *data,
                                 bool very_small_data_set=false,
                                 Notation n=NotationConstants::p1,
                                 const vector<string> observed_data={},
                                 const vector<string> tag_list={},
                                 vector<double> *saved_obs_seq_probs={});
}

#endif  // TRELLIS_H_
