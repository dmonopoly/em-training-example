#ifndef HELPER_H_
#define HELPER_H_

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <map>
#include <sstream>
#include <vector>

#include "Notation.h"

using namespace std;

namespace OutputHelper {
  void PrintHeader(const vector<Notation> &nots);
  void PrintDataRow(int iteration, const vector<Notation> &nots,
                    const map<string, double> &data);
}

namespace NotationHelper {
  // Returns a vector of strings of length 1 representing each character in s.
  vector<string> Individualize(const string &s);
  // Returns a string that is the concatentation of all strings in v.
  string Combine(const vector<string> &v);
  string SurroundWithParentheses(const string &predicate, const string &target);
}

// Notation Calculator methods that use the map of calculations.
namespace Calculator {
  // Pre: Notation object's 'first' and 'second' values have same length.
  // 'first' represents the observed data; 'second', the proposed data
  // completion. Essentially, P(t, w), with AND.
  // Post: Normalized probability for the data completion represented by
  // Notation n - i.e., P(t, w).
  double ComputeNormalizedProbability(const Notation &n, const map<string,
      double> &data, const int &tag_list_size, const int &observed_data_size);

  // Pre: pn's 'first' and 'second' values have same length. pn is like
  // P(ABA|XYX), so pn.first is ABA and pn.second is XYX (basically, probability
  // of observed data | tags). cn is like C(X, A), so cn.first has size 2, and
  // cn.second has size 0.  Uses pn (prob notation) and cn (count notation) to
  // check if factor needed. 
  // Post: Returns the normalized probability*number of matches. e.g., .216*2. 
  double NormProbFactor(const double &normalizedProb, const Notation &pn, const
      Notation &cn);

  // Pre: Notation is like P(ABA) (empty second list), and data has
  // appropriate key-value pairs set.
  // Post: Updates data to have computed probability, which is \sum_{t1,t2,t3}
  // P(t1,t2,t3)P(ABA|t1,t2,t3).
  void UpdateProbOfObsDataSeq(const Notation &observedNotation, map<string,
      double> *data, const vector<string> &tagSequences);
}

namespace TagHandler {
  // Generates all possible tag sequences for the brute force method.
  vector<string> GenerateTagSequences(const vector<string> &tags, int size);
}

#endif
