#include "Helper.h"

namespace NotationHelper {
  vector<string> Individualize(const string &s) {
    stringstream ss;
    vector<string> result;
    for (int i = 0; i < s.size(); ++i) {
      ss << s[i];
      result.push_back(ss.str());
      ss.str("");
      ss.clear();
    }
    return result;
  }

  string SurroundWithParentheses(const string &predicate, const string &target) {
    return predicate + "(" + target + ")";
  }
}

namespace Calculator {
  double ComputeNormalizedProbability(const Notation &n, const map<string,
      double> &data, const int &tag_list_size, const int &observed_data_size) {
    double d = 1;
    for (int i = 0; i < n.second.size(); ++i) {
      string tagKey;
      if (i == 0)
        tagKey = NotationHelper::SurroundWithParentheses("P", n.second[i]);
      else {
        tagKey = NotationHelper::SurroundWithParentheses("P", n.second[i] + Notation::GIVEN_DELIM +
            n.second[i-1]); 
      }
      string wordTagKey = NotationHelper::SurroundWithParentheses("P", n.first[i] +
          Notation::GIVEN_DELIM + n.second[i]);

      auto it = data.find(tagKey);
      assert(it != data.end() && "tagKey was not found");
      double tagProb = it->second;

      it = data.find(wordTagKey);
      assert(it != data.end() && "wordTagKey was not found");
      double wordTagProb = it->second;

      //cout << tagKey << ": " << tagProb << ", " << wordTagKey << ": " << wordTagProb << endl;
      d *= tagProb;
      d *= wordTagProb;
    }
    return d*pow(tag_list_size, observed_data_size);
  }

  double NormProbFactor(const double &normalizedProb, const Notation &pn, const
      Notation &cn) {
    string one = cn.first[0];  // "X"
    string two = cn.first[1];  // "A"
    int factor = 0;
    for (int i = 0; i < pn.second.size(); ++i) {
      // could allow for "X" and "A" swapped in cn...
      if (pn.second[i] == one && pn.first[i] == two) {
        ++factor;
      }
    }
    return normalizedProb*factor;
  }

  void UpdateProbOfObsDataSeq(const Notation &observedNotation, map<string, double>
      *data, const vector<string> &tagSequences) {
    double sum = 0;
    for (string tagSeq : tagSequences) {
      // Compute P(t1,t2,t3). Compute P(aba|t1,t2,t3) = P(a|t1)P(b|t2)...
      double probOfTagSeq = 1;
      double probOfObservedGivenTagSeq = 1;
      assert(tagSeq.size() == observedNotation.first.size() && "Tag sequence "
          "and observed data sequence are not the same size.");
      for (int i = 0; i < tagSeq.size(); ++i) {
        string currTag = string(1, tagSeq[i]);
        string tagKey = NotationHelper::SurroundWithParentheses("P", currTag);
        probOfTagSeq *= data->at(tagKey);

        string obsGivenTagKey = NotationHelper::SurroundWithParentheses("P",
            observedNotation.first[i] + Notation::GIVEN_DELIM + currTag);
        probOfObservedGivenTagSeq *= data->at(obsGivenTagKey);
      }
      sum += probOfTagSeq*probOfObservedGivenTagSeq;
    }
    cout << "Setting " << observedNotation.repr() << " to " << sum << endl;
    (*data)[observedNotation.repr()] = sum;
  }
}
