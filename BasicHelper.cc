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

// Key idea for these two methods:
// add(log x, log y)
// if log x = -infinity, then output log y
// if log y = -infinity, then output log x
// if log x - log y > 16 (for floats, 32 double) then output log x
//  (y is minor anyway)
// if log x > log y then output log x + log(1 + exp(log y - log x))
// if log y - log x > 16 (for floats, 32 double) then output log y
//  (x is minor anyway)
// if log y > log x then output log y + log(1 + exp(log x - log y))
double TakeLogsAndSum(double x, double y) {
  if (x < 0 || y < 0) {
    cerr << "TakeLogsAndSum Domain Error: x or y was negative." << endl;
    cerr << "x: " << x << ", y: " << y << endl;
    return -1;
  } else if (x == 0 && y == 0) {
    cerr << "TakeLogsAndSum Error: Adding log(0) and log(0)." << endl;
    return log(x) + log(y);
  } 
  if (x == 0 && y > 0) {
    return log(y);
  } else if (x > 0 && y == 0) {
    return log(x);
  }
  double logx = log(x);
  double logy = log(y);
  return AddLogs(logx, logy);
}
double AddLogs(double logx, double logy) {
//   // add(log x, log y)
//   // if log x - log y > 16 (for floats, 32 double) then output log x
//   //  (y is minor anyway)
//   // if log x > log y then output log x + log(1 + exp(log y - log x))
//   // if log y - log x > 16 (for floats, 32 double) then output log y
//   //  (x is minor anyway)
//   // if log y > log x then output log y + log(1 + exp(log x - log y))
  if (logx - logy > 32) {
    return logx;
  } else if (logy - logx > 32) {
    return logy;
  } else if (logx > logy) {
    return logx + log(1 + exp(logy - logx));
  } else {
//     cout << "Case 4" << endl;
//     cout << "logy: " << logy << "; " << "logx: " << logx << endl;
//     cout << "logx - logy: " << logx - logy << endl;
//     cout << "exp(logx - logy): " << exp(logx - logy) << endl;
//     cout << "log(1 + exp(logx - logy)): " << log(1 + exp(logx - logy)) << endl;
//     cout << "logy + log(1 + exp(logx - logy)): " << logy + log(1 + exp(logx - logy)) << endl;
    return logy + log(1 + exp(logx - logy));
  }

  // Online method:
//   // 1. make x the max
//   if (logy > logx) {
//     double temp = logx;
//     logx = logy;
//     logy = temp;
//   }
//   // 2. now x is bigger
//   if (logx == -DBL_MAX) {
//     return logx;
//   }
//   // 3. how far "down" (think decibels) is logy from logx?
//   //    if it's really small (20 orders of magnitude smaller), then ignore
//   double negDiff = logy - logx;
//   if (negDiff < -20) {
//     return logx;
//   }
//   // 4. otherwise use some nice algebra to stay in the log domain
//   //    (except for negDiff)
//   return logx + log(1.0 + exp(negDiff));
}
double SubtractLogs(double logx, double logy) {
  // Idea is same as AddLogs, except use log(x) + log(1 - exp(logy - logx)).
  // http://d.pr/i/rGSl
  if (logx - logy > 32) {
    return logx;
  } else if (logy - logx > 32) {
    cerr << "SubtractLogs Error (>32): tried computing log(x-y) for y>x.\n";
    return logy;
  } else if (logx > logy) {
    return logx + log(1 - exp(logy - logx));
  } else if (logx == logy) {
    return -DBL_MAX;
  } else {
    cerr << "SubtractLogs Error: tried computing log(x-y) for y>x.\n";
    return logy + log(1 - exp(logx - logy));
  }
}

}  // namespace Basic
