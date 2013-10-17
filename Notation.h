// Models first|second, as in P(first|second) or C(first|second).
#ifndef NOTATION_H_
#define NOTATION_H_

#include <initializer_list>
#include <iostream>
#include <vector>
#include <sstream>

using namespace std;

class Notation {
 public:
  static const string GIVEN_DELIM;
  static const string AND_DELIM;
  static const string SEQ_DELIM;
  static const string NULL_DELIM;

  static const Notation p1;  // Auto-probability 1.

  string predicate;  // P or C. TODO: include SIGMA for sums...
  // Denotes the delimiter separating 'first' and 'second'.
  string delimiter;

  // Delimiters within each list. This specificity allows P(ABA|t1,t2,t3).
  string first_delimiter, second_delimiter;

  vector<string> first;  // Each string can be used as a key to the map.
  vector<string> second;
  Notation();
  Notation(string predicate, initializer_list<string> first_args);
  Notation(string predicate, vector<string> first);
  Notation(string predicate, initializer_list<string> first_args, string first_delim);
  Notation(string predicate, vector<string> first, string first_delim);
  Notation(string predicate, initializer_list<string> first_args,
           string delimiter, initializer_list<string> second_args);
  Notation(string predicate, vector<string> first, string delimiter,
           vector<string> second);

  void set_first(const vector<string>& other) {
    this->first = other;
  }
  void set_second(const vector<string>& other) {
    this->second = other;
  }
  void set_predicate(const string new_pred) {
    this->predicate = new_pred;
  }
  string repr() const;
};
ostream& operator<<(ostream& out, const Notation& n);

#endif
