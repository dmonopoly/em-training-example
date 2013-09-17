#include <cassert>
#include <ctime>
#include <iostream>
#include <map>
#include <vector>

#include "Helper.h"
#include "Notation.h"

// The number of iterations to do the EM training.
#define NUMBER_ITERATIONS 5

using namespace std;

const string X = "X";
const string Y = "Y";
const string A = "A";
const string B = "B";
const vector<string> TAG_LIST{X, Y};
const vector<string> OBSERVED_DATA{A, B, A};

// Known probabilities:
Notation pX("P", {X}, {});  // "probability of x"
Notation pY("P", {Y}, {});
Notation pXGivenX("P", {X}, {X});
Notation pYGivenX("P", {Y}, {X});
Notation pXGivenY("P", {X}, {Y});
Notation pYGivenY("P", {Y}, {Y});
// Objectives:
Notation pABA("P", {A,B,A}, {});
Notation pAGivenX("P", {A}, {X});
Notation pAGivenY("P", {A}, {Y});
Notation pBGivenX("P", {B}, {X});
Notation pBGivenY("P", {B}, {Y});
Notation cXA("C", {X, A}, {});  // "count of x intersected with a"
Notation cXB("C", {X, B}, {});

void PrepareInitialData(map<string, double> *data) {
  // Given data.
  data->emplace(pX.repr(), .6);
  data->emplace(pY.repr(), .4);
  data->emplace(pXGivenX.repr(), .6);
  data->emplace(pYGivenX.repr(), .4);
  data->emplace(pXGivenY.repr(), .9);
  data->emplace(pYGivenY.repr(), .1);

  // Initial value for unknowns. We improve upon these.
  double initVal = .5;
  data->emplace(pAGivenX.repr(), initVal);
  data->emplace(pAGivenY.repr(), initVal);
  data->emplace(pBGivenX.repr(), initVal);
  data->emplace(pBGivenY.repr(), initVal);
}

void BruteForceCompute(map<string, double> *data) {
  // Enumerate all possible tag sequences.
  vector<string> tagSequences{
    X+X+X, X+X+Y, X+Y+X, X+Y+Y,
    Y+X+X, Y+X+Y, Y+Y+X, Y+Y+Y};
  for (int i = 0; i < NUMBER_ITERATIONS; ++i) {
    cout << "#" << i+1 << ":\n";
    // Get norm P(t,w) and counts.
    for (string seq : tagSequences) {
      vector<string> tags = NotationHelper::Individualize(seq);
      Notation pTW("P", OBSERVED_DATA, tags);
      double normalizedProb = Calculator::ComputeNormalizedProbability(pTW,
          *data, TAG_LIST.size(), OBSERVED_DATA.size());
      data->emplace(pTW.repr(), normalizedProb);
      // Get counts.
      (*data)[cXA.repr()] += Calculator::NormProbFactor(normalizedProb, pTW,
          cXA); (*data)[cXB.repr()] +=
        Calculator::NormProbFactor(normalizedProb, pTW, cXB);
    }
    // Update the unknown probabilities that we want to find. Use them in the
    // next iteration.
    (*data)[pAGivenX.repr()] = (*data)[cXA.repr()]/( (*data)[cXA.repr()] +
        (*data)[cXB.repr()] );
    (*data)[pBGivenX.repr()] = (*data)[cXB.repr()]/( (*data)[cXB.repr()] +
        (*data)[cXA.repr()] );

    // The ultimate value we want to maximize. This should increase with each
    // iteration.
    Calculator::UpdateProbOfObsDataSeq(pABA, data, tagSequences);
    cout << cXA << ": " << (*data)[cXA.repr()] << endl;
    cout << cXB << ": " << (*data)[cXB.repr()] << endl;
    cout << pABA << ": " << (*data)[pABA.repr()] << endl;
  }
}

struct Node {
  string name;
  vector<Edge> parent_edges;
  vector<Edge> child_edges;
  Node(string name) {
    this->name = name;
  }
  string repr() {
    return this->name;
  }
  void set_name(const string &other) {
    this->name = other;
  }
};

struct Edge {
  Notation notation;
  Node src, dest;
  Edge(const Notation &n, const Node &src, const Node &dest) {
    this->notation = n;
    this->src = src;
    this->dest = dest;
  }
  string repr() {
    return notation.repr();
  }
}

void LinkNodeAndEdge(Node *node, Edge *edge) {

}

// Warning: Creates data on heap. Call DestroyTrellis after done.
// Post: 'nodes' points to a vector where [0] is the start node, back() is the
// end. 'edges' points to a vector of corresponding edges.
void BuildTrellis(vector<Node> *nodes, Vector<Edge> *edges) {
  Node start_node = new Node("start");
  nodes->push_back(start_node);
  Node xa1first = new Node("xa1first");
  nodes->push_back(xa1first);
  Node xa2first = new Node("xa2first");
  nodes->push_back(xa2first);
  Node xb1second = new Node("xb1second");
  nodes->push_back(xb1second);
  Node xb2second = new Node("xb2second");
  nodes->push_back(xb2second);
  Node xa1third = new Node("xa1third");
  nodes->push_back(xa1third);
  Node xa2third = new Node("xa2third");
  nodes->push_back(xa2third);
  Node ya1first = new Node("ya1first");
  nodes->push_back(ya1first);
  Node ya2first = new Node("ya2first");
  nodes->push_back(ya2first);
  Node yb1second = new Node("yb1second");
  nodes->push_back(yb1second);
  Node yb2second = new Node("yb2second");
  nodes->push_back(yb2second);
  Node ya1third = new Node("ya1third");
  nodes->push_back(ya1third);
  Node ya2third = new Node("ya2third");
  nodes->push_back(ya2third);
  Node end_node = new Node("end");
  nodes->push_back(end_node);

  // TODONOW
}

void ForwardBackwardCompute(const vector<Node> &nodes, const Node &start_node,
                            const Node &end_node, const vector<Edge> &edges,
                            map<string, double> *data) {
  map<string, double> alpha;  // Sum of all paths from start state to this node.
  map<string, double> beta;  // Sum of all paths from this node to final state.

  alpha[start_node.repr()] = 1;

  // Forward pass.
  for (int i = 0; i < nodes.size(); ++i) {
    if (nodes[i] != start_node) {
      double sum = 0;
      for (Edge e : nodes[i].parent_edges) {
        sum += alpha[e.src.repr()] * data->at(e.repr());
      }
      alpha[nodes[i].repr()] = sum;
    }
  }

  // Backward pass. TODO.

  // Counting pass.
  for (int i = 0; i < edges.size(); ++i) {

  }

}

void RunBruteForceEM() {
  map<string, double> data;  // Storage for probabilities and counts.
  PrepareInitialData(&data);

  clock_t t;
  t = clock();
  BruteForceCompute(&data);
  t = clock() - t;
  
  cout << "--Results--\n";
  cout << cXA << ": " << data[cXA.repr()] << endl;
  cout << cXB << ": " << data[cXB.repr()] << endl;
  cout << pABA << ": " << data[pABA.repr()] << endl;
  printf("It took me %lu clicks (%f seconds).\n", t, ((float)t)/CLOCKS_PER_SEC);
}

void RunEfficientEM() {
  map<string, double> data;  // Storage for probabilities and counts.
  vector<Node> nodes;
  vector<Edge> edges;

  PrepareInitialData(&data);
  BuildTrellis(&nodes, &edges);

  clock_t t;
  t = clock();
  ForwardBackwardCompute(nodes, nodes[0], nodes.back(), edges, &data);
  t = clock() - t;
  
  cout << "--Results--\n";
  cout << cXA << ": " << data[cXA.repr()] << endl;
  cout << cXB << ": " << data[cXB.repr()] << endl;
  cout << pABA << ": " << data[pABA.repr()] << endl;
  printf("It took me %lu clicks (%f seconds).\n", t, ((float)t)/CLOCKS_PER_SEC);

  DestroyTrellis(&nodes, &edges);
}

int main() {
  RunBruteForceEM();
  RunEfficientEM();
  return 0;
}
