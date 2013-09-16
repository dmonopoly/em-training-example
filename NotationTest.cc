#include <assert.h>

// Includes various libraries as well and sets namespace to std.
#include "Notation.h"

void TestXGivenY() {
  Notation n("P");
  vector<string> test;
  test.push_back("a");
  test.push_back("b");
  n.set_first(test);
  n.set_second(test);

  stringstream ss;  // We simulate 'cout <<' with this.

  ss << n;
  cout << ss.str() << endl;
  assert(ss.str() == "P(a,b|a,b)");
  ss.clear();
}

void TestJustX() {
  Notation n("P");
  vector<string> test;
  test.push_back("a");
  test.push_back("b");
  n.set_first(test);
  //n.set_second(test);  // This time, second is never set.

  stringstream ss;  // We simulate 'cout <<' with this.

  ss << n;
  cout << ss.str() << endl;
  assert(ss.str() == "P(a,b)" && "Nonexistent second vector didn't work.");
}

void TestJustXAgain() {
  Notation n("C", {"A","B","C"}, {});
  stringstream ss;
  ss << n;
  cout << ss.str() << endl;
  assert(ss.str() == "C(A,B,C)" && "No good for JustXAgain.");
}

void TestNotation() {
  cout << "Testing Notation..." << endl;
  TestXGivenY();
  TestJustX();
  TestJustXAgain();
  cout << "Test(s) passed." << endl;
}

int main() {
  TestNotation();
  return 0;
}

