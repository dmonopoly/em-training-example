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
  vector<pair<Node, Edge> > parents;  // Incoming edges.
  vector<pair<Node, Edge> > children;  // Outgoing edges.
  Node(string name) {
    this->name = name;
  }
  string repr() {
    return name;
  }
};

struct Edge {
  Notation notation;
  Edge(Notation n) {
    this->notation = n;
  }
  string repr() {
    return notation.repr();
  }
}

void PrepareTrellis(vector<Node> *graph) {

}

void ForwardBackwardCompute(const vector<Node> &nodes, const Node &startNode,
    const Node &endNode, map<string, double> *data) {
  map<string, double> alpha;  // Sum of all paths from start state to this node.
  map<string, double> beta;  // Sum of all paths from this node to final state.

  alpha[startNode.repr()] = 1;

  // Forward pass.
  for (int i = 0; i < nodes.size(); ++i) {
    if (nodes[i] != startNode) {
      double sum = 0;
      for (Node parent : nodes[i].parents) {
        sum += alpha[parent.first.repr()] * data->at(parent.second.repr());
      }
      alpha[nodes[i].repr()] = sum;
    }
  }

  // Backward pass. TODO.

  // Counting pass.


}

int main() {
  map<string, double> data;  // Storage for probabilities and counts.
  PrepareInitialData(&data);
  BruteForceCompute(&data);
  
  cout << "--Results--\n";
  cout << cXA << ": " << data[cXA.repr()] << endl;
  cout << cXB << ": " << data[cXB.repr()] << endl;
  cout << pABA << ": " << data[pABA.repr()] << endl;
  return 0;
}
