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

  string predicate;  // P or C
  vector<string> first;  // Each string can be used as a key to the map.
  vector<string> second;
  Notation();
  Notation(string predicate);
  Notation(string predicate, initializer_list<string> first_args,
      initializer_list<string> second_args);
  Notation(string predicate, vector<string> first, vector<string> second);

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
