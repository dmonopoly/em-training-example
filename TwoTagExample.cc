#include <cassert>
#include <iostream>
#include <map>
#include <vector>

#include "Helper.h"
#include "Notation.h"

// The number of iterations to do the EM training.
#define NUMBER_ITERATIONS 100

// Two ways to run this program: with a short or a long observed sequence.
#define DO_SHORT_SEQ false

// TODO: Reorganize and use Notation::GIVEN_DELIM. http://bit.ly/15rbAom
#define GIVEN_DELIM "|"
#define AND_DELIM ","
#define SEQ_DELIM ""

using namespace std;

const string X = "X";
const string Y = "Y";
const string A = "A";
const string B = "B";
const vector<string> TAG_LIST{X, Y};
#if DO_SHORT_SEQ
const vector<string> OBSERVED_DATA{A, B, A};
#else
const vector<string> OBSERVED_DATA{A,A,A,B,A,A,B,A,A};
#endif

// For output.
vector<double> saved_results;

// Known probabilities:
Notation pX("P", {X});  // "probability of x"
Notation pY("P", {Y});
Notation pXGivenX("P", {X}, GIVEN_DELIM, {X});
Notation pYGivenX("P", {Y}, GIVEN_DELIM, {X});
Notation pXGivenY("P", {X}, GIVEN_DELIM, {Y});
Notation pYGivenY("P", {Y}, GIVEN_DELIM, {Y});
// Objectives:
// Short seq type.
Notation pABA("P", {A,B,A}, SEQ_DELIM);
Notation pAGivenX("P", {A}, GIVEN_DELIM, {X});
Notation pAGivenY("P", {A}, GIVEN_DELIM, {Y});
Notation pBGivenX("P", {B}, GIVEN_DELIM, {X});
Notation pBGivenY("P", {B}, GIVEN_DELIM, {Y});
Notation cXA("C", {X, A}, AND_DELIM);  // "count of x and a"
Notation cXB("C", {X, B}, AND_DELIM);
Notation cYA("C", {Y, A}, AND_DELIM);
Notation cYB("C", {Y, B}, AND_DELIM);
// Long seq type.
Notation pLong("P", {A,A,A,B,A,A,B,A,A}, SEQ_DELIM);

void PrepareInitialData(map<string, double> *data) {
  // Given data.
  data->emplace(pX.repr(), .6);
  data->emplace(pY.repr(), .4);
  data->emplace(pXGivenX.repr(), .6);
  data->emplace(pYGivenX.repr(), .4);
  data->emplace(pXGivenY.repr(), .9);
  data->emplace(pYGivenY.repr(), .1);

  // Initial value for unknowns. We improve upon these.
  double initVal = .5;
  data->emplace(pAGivenX.repr(), initVal);
  data->emplace(pAGivenY.repr(), initVal);
  data->emplace(pBGivenX.repr(), initVal);
  data->emplace(pBGivenY.repr(), initVal);
}

void ComputeDataWithBruteForce(map<string, double> *data, const Notation &n,
                               const vector<string> &tag_sequences) {
  saved_results.push_back((*data)[n.repr()]); // push back initial 0
  cout << "Initially: \n";
  cout << cXA << ": " << (*data)[cXA.repr()] << endl;
  cout << cXB << ": " << (*data)[cXB.repr()] << endl;
  cout << cYA << ": " << (*data)[cYA.repr()] << endl;
  cout << cYB << ": " << (*data)[cYB.repr()] << endl;
  cout << n << ": " << (*data)[n.repr()] << endl << endl;

  for (int i = 0; i < NUMBER_ITERATIONS; ++i) {
    cout << "#" << i+1 << ":\n";
    // Reset counts to zero.
    (*data)[cXA.repr()] = 0;
    (*data)[cXB.repr()] = 0;
    (*data)[cYA.repr()] = 0;
    (*data)[cYB.repr()] = 0;

    // Get norm P(t,w) and counts.
    for (string seq : tag_sequences) {
      vector<string> tags = NotationHelper::Individualize(seq);
      Notation pTW("P", OBSERVED_DATA, AND_DELIM, tags);
      double normalizedProb = Calculator::ComputeNormalizedProbability(pTW,
          *data, TAG_LIST.size(), OBSERVED_DATA.size());
      (*data)[pTW.repr()] = normalizedProb;

      // Get counts.
      (*data)[cXA.repr()] += Calculator::NormProbFactor(normalizedProb, pTW,
          cXA);
      (*data)[cXB.repr()] +=
        Calculator::NormProbFactor(normalizedProb, pTW, cXB);
      (*data)[cYA.repr()] += Calculator::NormProbFactor(normalizedProb, pTW,
          cYA);
      (*data)[cYB.repr()] +=
        Calculator::NormProbFactor(normalizedProb, pTW, cYB);
    }
    // Update the unknown probabilities that we want to find. Use them in the
    // next iteration.
    (*data)[pAGivenX.repr()] = (*data)[cXA.repr()]/( (*data)[cXA.repr()] +
        (*data)[cXB.repr()] );
    (*data)[pBGivenX.repr()] = (*data)[cXB.repr()]/( (*data)[cXB.repr()] +
        (*data)[cXA.repr()] );
    (*data)[pAGivenY.repr()] = (*data)[cYA.repr()]/( (*data)[cYA.repr()] +
        (*data)[cYB.repr()] );
    (*data)[pBGivenY.repr()] = (*data)[cYB.repr()]/( (*data)[cYB.repr()] +
        (*data)[cYA.repr()] );

    // The ultimate value we want to maximize. This should increase with each
    // iteration.
    Calculator::UpdateProbOfObsDataSeq(n, data, tag_sequences);
    cout << "--Summary of iteration " << i+1 << "--\n";
    cout << cXA << ": " << (*data)[cXA.repr()] << endl;
    cout << cXB << ": " << (*data)[cXB.repr()] << endl;
    cout << cYA << ": " << (*data)[cYA.repr()] << endl;
    cout << cYB << ": " << (*data)[cYB.repr()] << endl;
    cout << pAGivenX << ": " << (*data)[pAGivenX.repr()] << endl;
    cout << pBGivenX << ": " << (*data)[pBGivenX.repr()] << endl;
    cout << pAGivenY << ": " << (*data)[pAGivenY.repr()] << endl;
    cout << pBGivenY << ": " << (*data)[pBGivenY.repr()] << endl;
    cout << n << ": " << (*data)[n.repr()] << endl;
    cout << endl;
    saved_results.push_back((*data)[n.repr()]);
  }
}

void OutputResults(map<string, double> &data, Notation n, const vector<string>
                   &tag_sequences) {
  cout << "--Results based on " << NUMBER_ITERATIONS << " iterations--\n";
  cout << n << ": ";
  for (int i = 0; i < saved_results.size(); ++i) {
    cout << saved_results[i] << " ";
  }
  cout << endl << endl;
  cout << "Final " << n << ": " << data[n.repr()] << endl;
  cout << "Final " << pAGivenX << ": " << data[pAGivenX.repr()] << endl;
  cout << "Final " << pBGivenX << ": " << data[pBGivenX.repr()] << endl;
  cout << "Final " << pAGivenY << ": " << data[pAGivenY.repr()] << endl;
  cout << "Final " << pBGivenY << ": " << data[pBGivenY.repr()] << endl << endl;

  cout << "Determining the best matching tag sequence:\n";
  vector<string> tags = NotationHelper::Individualize(tag_sequences.at(0));
  Notation pTW_first("P", OBSERVED_DATA, AND_DELIM, tags);
  Notation *best_pTGivenW = NULL;
  string best_match_string_repr = pTW_first.repr();
  for (string seq : tag_sequences) {
    vector<string> tags = NotationHelper::Individualize(seq);
    Notation pTW("P", OBSERVED_DATA, AND_DELIM, tags);

    // Compute P(t|w). Technically not used because divided values seem to
    // incorrectly yield >1 (decimals too small, possibly).
    Notation pTGivenW("P", tags, GIVEN_DELIM, OBSERVED_DATA);
    data[pTGivenW.repr()] = data[pTW.repr()] / data[n.repr()];
    cout << pTW << ": " << data[pTW.repr()] << endl; // ", " << pTGivenW << ": " <<
    //  data[pTGivenW.repr()] << endl;
    if (data[pTW.repr()] > data[best_match_string_repr]) {
      best_match_string_repr = pTW.repr();
      delete best_pTGivenW;
      // Same as pTGivenW.
      best_pTGivenW = new Notation("P", tags, GIVEN_DELIM, OBSERVED_DATA);
    }
  }
  string pTAndWRepr = best_match_string_repr;
  cout << "The highest probability found belongs to " << pTAndWRepr << ": " <<
    data[pTAndWRepr] << endl;
  cout << "The best matching tag sequence is " <<
    NotationHelper::Combine(best_pTGivenW->first) << endl;
  delete best_pTGivenW;
}

int main() {
  map<string, double> data;
  PrepareInitialData(&data);
  vector<string> tag_sequences = TagHandler::GenerateTagSequences(TAG_LIST,
      OBSERVED_DATA.size());
  if (DO_SHORT_SEQ) {
    ComputeDataWithBruteForce(&data, pABA, tag_sequences);
    OutputResults(data, pABA, tag_sequences);
  } else {
    ComputeDataWithBruteForce(&data, pLong, tag_sequences);
    OutputResults(data, pLong, tag_sequences);
  }
  return 0;
}
