#include "Notation.h"

const string Notation::GIVEN_DELIM = "|";
const string Notation::AND_DELIM = ",";

Notation::Notation(string predicate) {
  this->predicate = predicate;
}

Notation::Notation(string predicate, initializer_list<string> first_args,
    initializer_list<string> second_args) {
  this->predicate = predicate;
  for (string s : first_args)
    this->first.push_back(s);
  for (string s : second_args)
    this->second.push_back(s);
}

Notation::Notation(string predicate, vector<string> first, vector<string> second) {
  this->predicate = predicate;
  this->first = first;
  this->second = second;
}

string Notation::repr() const {
  // TODO: Clarify notation to be more flexible. Comma strictly means AND, and
  // sometimes one list may want that, e.g. P(aba|t1,t2,t3). Also, this assumes
  // GIVEN_DELIM is the intended meaning; new param in constructor to allow
  // AND_DELIM instead?
  stringstream ss;
  ss << this->predicate << "(";
  for (int i = 0; i < this->first.size()-1; ++i)
    ss << this->first[i] << ",";
  if (!this->first.empty())
    ss << this->first.back();
  if (!this->second.empty()) { 
    ss << GIVEN_DELIM;
    for (int i = 0; i < this->second.size()-1; ++i)
      ss << this->second[i] << ",";
    ss << this->second.back();
  }
  ss << ")";
  return ss.str();
}

ostream& operator<<(ostream& out, const Notation& n) {
    return out << n.repr();
}

