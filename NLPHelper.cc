// Functions used primarily for the brute force method in TwoTagExample.cc.
// Can decouple TagHandler from TrellisAid to make this not needed at all by
// TrellisAid.
#include "NLPHelper.h"
#include "BasicHelper.h"

#define EXTRA_PRINTING false

namespace OutputHelper {
  void PrintHeader(const vector<Notation> &nots) {
    vector<string> header;
    header.push_back("Iteration");
    for (Notation n : nots) {
      header.push_back(n.repr());
    }
    Basic::PrintRow(header);
  }

  void PrintDataRow(int iteration, const vector<Notation> &nots,
                    const map<Notation, double> &data) {
    vector<double> values;
    values.push_back((double) iteration);
    for (int i = 0; i < nots.size(); ++i) {
      try {
        values.push_back(data.at(nots[i]));
      } catch (exception &e) {
        cerr << "No key: " << nots[i] << endl;
      }
    }
    Basic::PrintRow(values);
  }
}

namespace NotationHelper {
  string Combine(const vector<string> &v) {
    stringstream ss;
    for (auto it = v.begin(); it != v.end(); ++it) {
      ss << *it;
    }
    return ss.str();
  }

  string SurroundWithParentheses(const string &predicate, const string &target) {
    return predicate + "(" + target + ")";
  }

  void ReplaceSymbol(const string &old_s, const string &new_s, Notation *n) {
    vector<string> changed;
    for (int i = 0; i < n->first.size(); ++i) {
      if (n->first[i] == old_s) {
        if (EXTRA_PRINTING)
          cout << "REPLACING..." << n->first[i] << " w/ " << new_s << endl;
        changed.push_back(new_s);
      } else
        changed.push_back(n->first[i]);
    }
    n->set_first(changed);
    changed.clear();
    for (int i = 0; i != n->second.size(); ++i) {
      if (n->second[i] == old_s) {
        if (EXTRA_PRINTING)
          cout << "REPLACING..." << n->second[i] << " w/ " << new_s << endl;
        changed.push_back(new_s);
      } else
        changed.push_back(n->second[i]);
    }
    n->set_second(changed);
  }
}

namespace Calculator {
  double ComputeUnnormalizedProbability(const Notation &n,
                                        const map<Notation, double> &data) {
    double d = 1;
    for (int i = 0; i < n.second.size(); ++i) {
      Notation tagKey;
//       string tagKey;
      if (i == 0){
        Notation tmp("P", {n.second[i]});
//         tagKey = NotationHelper::SurroundWithParentheses("P", n.second[i]);
        tagKey = tmp;
      } else {
        Notation tmp("P", {n.second[i]}, Notation::GIVEN_DELIM,
                     {n.second[i-1]});
//         tagKey = NotationHelper::SurroundWithParentheses("P", n.second[i] +
//             Notation::GIVEN_DELIM + n.second[i-1]); 
        tagKey = tmp;
      }
      Notation wordTagKey("P", {n.first[i]}, Notation::GIVEN_DELIM,
                          {n.second[i]});
//       string wordTagKey = NotationHelper::SurroundWithParentheses("P", n.first[i] +
//           Notation::GIVEN_DELIM + n.second[i]);

      auto it = data.find(tagKey);
      assert(it != data.end() && "tagKey was not found");
      double tagProb = it->second;

      it = data.find(wordTagKey);
      assert(it != data.end() && "wordTagKey was not found");
      double wordTagProb = it->second;

      d *= tagProb;
      d *= wordTagProb;
    }
    assert(d < 1 && "Unnormalized probability is too high.");
    return d;
  }

  double NormProbFactor(const double &normalizedProb, const Notation &pn,
                        const Notation &cn) {
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

  void UpdateProbOfObsDataSeq(const Notation &observedNotation,
                              map<Notation, double> *data,
                              const vector<vector<string> > &tagSequences) {
    double sum = 0;
    for (vector<string> tagSeq : tagSequences) {
      // Compute P(t1,t2,t3) = P(t1)P(t2|t1)P(t3|t2).
      // Compute P(aba|t1,t2,t3) = P(a|t1)P(b|t2)...
      double probOfTagSeq = 1;
      double probOfObservedGivenTagSeq = 1;
      assert(tagSeq.size() == observedNotation.first.size() && "Tag sequence "
          "and observed data sequence are not the same size.");
      string prevTag;
      for (int i = 0; i < tagSeq.size(); ++i) {
        string currTag = tagSeq[i];
//         string tagKey;
        Notation theTagKey;
        if (i == 0) {
          Notation tagKey("P", {currTag});
//           tagKey = NotationHelper::SurroundWithParentheses("P", currTag);
          probOfTagSeq *= data->at(tagKey);
          theTagKey = tagKey;
        } else {
          Notation tagKey("P", {currTag}, Notation::GIVEN_DELIM, {prevTag});
//           tagKey = NotationHelper::SurroundWithParentheses("P", currTag +
//               Notation::GIVEN_DELIM + prevTag);
          probOfTagSeq *= data->at(tagKey);
          theTagKey = tagKey;
        }
        prevTag = currTag;
        if (EXTRA_PRINTING)
          cout << theTagKey << ": " << data->at(theTagKey) << endl;

        Notation obsGivenTagKey("P", {observedNotation.first[i]},
            Notation::GIVEN_DELIM, {currTag});
//         string obsGivenTagKey = NotationHelper::SurroundWithParentheses("P",
//             observedNotation.first[i] + Notation::GIVEN_DELIM + currTag);
        if (EXTRA_PRINTING)
          cout << obsGivenTagKey << ": " << data->at(obsGivenTagKey) << endl;
        probOfObservedGivenTagSeq *= data->at(obsGivenTagKey);
      }
      sum += probOfTagSeq*probOfObservedGivenTagSeq;
    }
    if (EXTRA_PRINTING)
      cout << "Setting " << observedNotation.repr() << " to " << sum << endl;
    (*data)[observedNotation] = sum;
  }
}

namespace TagHandler {
  // Used only in GenerateTagSequences.
  vector<string> RecGenerateTagSequences(const vector<string> &tags, int size) {
    vector<string> list;
    if (size == 1) {
      for (string t : tags) {
        list.push_back(t);
      }
    } else {
      vector<string> tmp = TagHandler::RecGenerateTagSequences(tags, size - 1);
      for (string t1 : tags) {
        for (string t2 : tmp) {
          list.push_back(t1 + "~!~" + t2);
        }
      }
    }
    return list;
  }
  vector<vector<string> > GenerateTagSequences(const vector<string> &tags, int size) {
    vector<string> unsplit_tag_seqs = TagHandler::RecGenerateTagSequences(tags, size);
    vector<vector<string> > result;
    for (string unsplit : unsplit_tag_seqs) {
      int spot = unsplit.find("~!~");  // Arbitrary length-3 delimiter.
      vector<string> new_vec;
      while (spot != string::npos) {
        string s = unsplit.substr(0, spot);
        new_vec.push_back(s);
        unsplit = unsplit.replace(0, spot + 3, "");
        spot = unsplit.find("~!~");
      }
      new_vec.push_back(unsplit);
      result.push_back(new_vec);
    }
    return result;
  }
}
