#include "Notation.h"

const string Notation::GIVEN_DELIM = "|";
const string Notation::AND_DELIM = ",";
const string Notation::SEQ_DELIM = " ";
const string Notation::NULL_DELIM = "NULL";  // Signifies delim never used.
// const Notation p1("P", {"1"});  // Auto-probability 1.

Notation::Notation() {
}
Notation::Notation(string predicate, initializer_list<string> first_args) {
  this->predicate = predicate;
  for (string s : first_args)
    this->first.push_back(s);
  this->delimiter = Notation::NULL_DELIM;
  this->first_delimiter = Notation::SEQ_DELIM;
  this->second_delimiter = Notation::NULL_DELIM;

  this->repr_ = this->repr();
}
Notation::Notation(string predicate, vector<string> first) {
  this->predicate = predicate;
  this->first = first;
  this->delimiter = Notation::NULL_DELIM;
  this->first_delimiter = Notation::SEQ_DELIM;
  this->second_delimiter = Notation::NULL_DELIM;

  this->repr_ = this->repr();
}
Notation::Notation(string predicate, initializer_list<string> first_args, string first_delim) {
  this->predicate = predicate;
  for (string s : first_args)
    this->first.push_back(s);
  this->delimiter = Notation::NULL_DELIM;
  this->first_delimiter = first_delim;
  this->second_delimiter = Notation::NULL_DELIM;

  this->repr_ = this->repr();
}
Notation::Notation(string predicate, vector<string> first, string first_delim) {
  this->predicate = predicate;
  this->first = first;
  this->delimiter = Notation::NULL_DELIM;
  this->first_delimiter = first_delim;
  this->second_delimiter = Notation::NULL_DELIM;

  this->repr_ = this->repr();
}
Notation::Notation(string predicate, initializer_list<string> first_args,
           string delimiter, initializer_list<string> second_args) {
  this->predicate = predicate;
  for (string s : first_args)
    this->first.push_back(s);
  this->delimiter = delimiter;
  for (string s : second_args)
    this->second.push_back(s);
  this->first_delimiter = Notation::SEQ_DELIM;
  this->second_delimiter = Notation::SEQ_DELIM;

  this->repr_ = this->repr();
}
Notation::Notation(string predicate, vector<string> first, string delimiter,
                   vector<string> second) {
  this->predicate = predicate;
  this->first = first;
  this->delimiter = delimiter;
  this->second = second;
  this->first_delimiter = Notation::SEQ_DELIM;
  this->second_delimiter = Notation::SEQ_DELIM;

  this->repr_ = this->repr();
}

string Notation::repr() const {
  // TODO: could optimize with repr_ instance field, but then must update it
  // appropriately, overload = assignment, etc.... have an update_repr method.
  stringstream ss;
  ss << this->predicate << "(";
  for (int i = 0; i < this->first.size()-1; ++i)
    ss << this->first[i] << this->first_delimiter;
  if (!this->first.empty())
    ss << this->first.back();
  if (!this->second.empty()) { 
    ss << this->delimiter;
    for (int i = 0; i < this->second.size()-1; ++i)
      ss << this->second[i] << this->second_delimiter;
    ss << this->second.back();
  }
  ss << ")";
  return ss.str();
}

// static void Notation ExtractFromString(const string &s) {
//   string predicate = s.substr(0, 1);
//   vector<string> first, second;
//   string first_delim, second_delim;
//   size_t delim_pos = s.find(Notation::GIVEN_DELIM);
// // const string Notation::GIVEN_DELIM = "|";
// // const string Notation::AND_DELIM = ",";
// // const string Notation::SEQ_DELIM = " ";
// // const string Notation::NULL_DELIM = "NULL";  // Signifies delim never used.
//   if (delim_pos != string::npos) {
// //     vector<string> &split(const string &s, char delim, vector<string> &elems);
//     string first_part = s.substr(2, delim_pos - 2);
//     if (first_part.size() > 1)
//     char delim;
//     stringstream ss;
//     ss << Notation::
//     BasicHelper::Split( )
//   }
//   Notation n(predicate, first, first_delim, second, second_delim);
//   return n;
// }

ostream& operator<<(ostream& out, const Notation& n) {
    return out << n.repr();
}

