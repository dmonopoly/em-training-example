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

  string predicate;  // P or C
  // Denotes the delimiter separating 'first' and 'second'.
  string delimiter;

  // Delimiters within each list. This specificity allows P(ABA|t1,t2,t3).
  string first_delimiter, second_delimiter;

  vector<string> first;  // Each string can be used as a key to the map.
  vector<string> second;
  Notation(string predicate, initializer_list<string> first_args,
           string delimiter, initializer_list<string> second_args);
  Notation(string predicate, vector<string> first, string delimiter,
           vector<string> second);
  string repr() const;
};
ostream& operator<<(ostream& out, const Notation& n);

#endif
