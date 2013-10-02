// Putting variables above main is not well-defined behavior...
#include <cassert>
#include <iostream>
#include <map>
#include <vector>

//#include "Constants.h"
#include "Notation.h"

using namespace std;

const string X = "X";
const string Y = "Y";
const string A = "A";
const string B = "B";

Notation pXGivenX("P", {X}, Notation::GIVEN_DELIM, {X});

string s = Notation::GIVEN_DELIM;

int main() {
  cout << pXGivenX << endl; // GIVEN_DELIM is "" for some reason!
  Notation pXGivenX2("P", {X}, Notation::GIVEN_DELIM, {X});
  cout << pXGivenX2 << endl; // GIVEN_DELIM is properly "|"

  cout << "--" << endl;
  cout << s << endl;
  cout << Notation::GIVEN_DELIM << endl;

  return 0;
}
