// Helper functions for brute force and efficient EM + Viterbi.
#ifndef HELPER_H_
#define HELPER_H_

#include <cassert>
#include <cstdlib>
#include <iomanip>
#include <map>
#include <sstream>
#include <vector>

#include "Notation.h"
#include "Node.h"
#include "Edge.h"

using namespace std;

struct Node;
struct Edge;

namespace OutputHelper {
  void PrintHeader(const vector<Notation> &nots);
  void PrintDataRow(int iteration, const vector<Notation> &nots,
                    const map<string, double> &data);
}

namespace NotationHelper {
  // Returns a string that is the concatentation of all strings in v.
  string Combine(const vector<string> &v);
  string SurroundWithParentheses(const string &predicate, const string &target);
}

namespace GraphHelper {
  void LinkNodeAndEdge(Node *src_node, Edge &edge, Node *dest_node);
}

// Notation Calculator methods that use the map of calculations.
namespace Calculator {
  // Part of brute force method:
  // Pre: Notation object's 'first' and 'second' values have same length.
  // 'first' represents the observed data; 'second', the proposed data
  // completion. Essentially, P(t, w), with AND.
  // Post: Normalized probability for the data completion represented by
  // Notation n - i.e., P(t, w).
  double ComputeUnnormalizedProbability(const Notation &n, const map<string,
      double> &data);

  // Part of brute force method:
  // Pre: pn's 'first' and 'second' values have same length. pn is like
  // P(ABA|XYX), so pn.first is ABA and pn.second is XYX (basically, probability
  // of observed data | tags). cn is like C(X, A), so cn.first has size 2, and
  // cn.second has size 0.  Uses pn (prob notation) and cn (count notation) to
  // check if factor needed. 
  // Post: Returns the normalized probability*number of matches. e.g., .216*2. 
  double NormProbFactor(const double &normalizedProb, const Notation &pn,
                        const Notation &cn);

  // Usable by both brute force and efficient EM:
  // TODO: But this goes over all tag sequences! BAD.
  // Pre: Notation is like P(ABA) (empty second list), and data has
  // appropriate key-value pairs set.
  // Post: Updates data to have computed probability, which is \sum_{t1,t2,t3}
  // P(t1,t2,t3)P(ABA|t1,t2,t3).
  void UpdateProbOfObsDataSeq(const Notation &observedNotation,
                              map<string, double> *data,
                              const vector<vector<string> > &tagSequences);
}

namespace TagHandler {
  // Generates all possible tag sequences for the brute force method.
  vector<vector<string> > GenerateTagSequences(const vector<string> &tags,
                                               int size);
  vector<string> RecGenerateTagSequences(const vector<string> &tags, int size);
}

#endif
