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

}  // namespace Basic
