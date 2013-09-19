#include "Notation.h"

const string Notation::GIVEN_DELIM = "|";

Notation::Notation() {

}

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

