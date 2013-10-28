#include "BasicHelper.h"

namespace Basic {

string Tab(int n) {
    string tab = "";
    for (int i=0; i<n; i++) {
        tab += "\t";
    }
    return tab;
}

void Enter(int n) {
    for (int i=0; i < n; i++) {
        cout << endl;
    }
}

bool Contains(string s1, string s2) {
    if (std::string::npos != s1.find(s2))
        return true;
    else
        return false;
}

bool Contains(long num1, long num2) {
	stringstream ss;
	ss << num1;
	string s1 = ss.str();
	ss.str("");
	ss << num2;
	string s2 = ss.str();
	return Contains(s1, s2);
}

bool AskAndDecide(string question) {
  char response;
  cout << question;
  cin >> response;
  while (response != 'n' && response != 'N' && response != 'y' && response != 'Y') {
      cout << "Improper response.\n";
      cout << question;
      cin >> response;
  }
  switch(response) {
      case 'n':
      case 'N': return false;
      case 'y':
      case 'Y': return true;
  }
  cerr << "Reached spot it shouldn't have... (AskAndDecide)" << endl;
  return 1;
}

void ToLowerCase(string &s) {
  for (int i=0; i < s.size(); i++) {
    s[i] = tolower(s[i]);
  }
}

void ReplaceAll(string& str, const string& from, const string& to) {
  if(from.empty())
    return;
  size_t start_pos = 0;
  while((start_pos = str.find(from, start_pos)) != std::string::npos) {
    str.replace(start_pos, from.length(), to);
    // In case 'to' contains 'from', like replacing 'x' with 'yx'.
    start_pos += to.length();
  }
}

vector<string>& Split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

vector<string> Split(const string &s, char delim) {
    vector<string> elems;
    Split(s, delim, elems);
    return elems;
}
}  // namespace Basic
