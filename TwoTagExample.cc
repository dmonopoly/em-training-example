#include <cassert>
#include <iostream>
#include <map>
#include <vector>

#include "Helper.h"
#include "Notation.h"

// The number of iterations to do the EM training.
#define NUMBER_ITERATIONS 100

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
const vector<string> OBSERVED_DATA{A, B, A};
// Enumerate all possible tag sequences for the brute force method.
const vector<string> TAG_SEQUENCES{X+X+X, X+X+Y, X+Y+X, X+Y+Y,
                                   Y+X+X, Y+X+Y, Y+Y+X, Y+Y+Y};

// For output.
vector<double> saved_pABA_results;

// Known probabilities:
Notation pX("P", {X});  // "probability of x"
Notation pY("P", {Y});
Notation pXGivenX("P", {X}, GIVEN_DELIM, {X});
Notation pYGivenX("P", {Y}, GIVEN_DELIM, {X});
Notation pXGivenY("P", {X}, GIVEN_DELIM, {Y});
Notation pYGivenY("P", {Y}, GIVEN_DELIM, {Y});
// Objectives:
Notation pABA("P", {A,B,A}, SEQ_DELIM);
Notation pAGivenX("P", {A}, GIVEN_DELIM, {X});
Notation pAGivenY("P", {A}, GIVEN_DELIM, {Y});
Notation pBGivenX("P", {B}, GIVEN_DELIM, {X});
Notation pBGivenY("P", {B}, GIVEN_DELIM, {Y});
Notation cXA("C", {X, A}, AND_DELIM);  // "count of x and a"
Notation cXB("C", {X, B}, AND_DELIM);
Notation cYA("C", {Y, A}, AND_DELIM);
Notation cYB("C", {Y, B}, AND_DELIM);

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

void ComputeDataWithBruteForce(map<string, double> *data) {
  saved_pABA_results.push_back((*data)[pABA.repr()]);
  cout << "Initially: \n";
  cout << cXA << ": " << (*data)[cXA.repr()] << endl;
  cout << cXB << ": " << (*data)[cXB.repr()] << endl;
  cout << cYA << ": " << (*data)[cYA.repr()] << endl;
  cout << cYB << ": " << (*data)[cYB.repr()] << endl;
  cout << pABA << ": " << (*data)[pABA.repr()] << endl << endl;
  
  for (int i = 0; i < NUMBER_ITERATIONS; ++i) {
    cout << "#" << i+1 << ":\n";
    // Reset counts to zero.
    (*data)[cXA.repr()] = 0;
    (*data)[cXB.repr()] = 0;
    (*data)[cYA.repr()] = 0;
    (*data)[cYB.repr()] = 0;

    // Get norm P(t,w) and counts.
    for (string seq : TAG_SEQUENCES) {
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
    Calculator::UpdateProbOfObsDataSeq(pABA, data, TAG_SEQUENCES);
    cout << "--Summary of iteration " << i+1 << "--\n";
    cout << cXA << ": " << (*data)[cXA.repr()] << endl;
    cout << cXB << ": " << (*data)[cXB.repr()] << endl;
    cout << cYA << ": " << (*data)[cYA.repr()] << endl;
    cout << cYB << ": " << (*data)[cYB.repr()] << endl;
    cout << pAGivenX << ": " << (*data)[pAGivenX.repr()] << endl;
    cout << pBGivenX << ": " << (*data)[pBGivenX.repr()] << endl;
    cout << pAGivenY << ": " << (*data)[pAGivenY.repr()] << endl;
    cout << pBGivenY << ": " << (*data)[pBGivenY.repr()] << endl;
    cout << pABA << ": " << (*data)[pABA.repr()] << endl;
    cout << endl;
    saved_pABA_results.push_back((*data)[pABA.repr()]);
  }
}

int main() {
  map<string, double> data;
  PrepareInitialData(&data);
  ComputeDataWithBruteForce(&data);
  
  // Goal:
  cout << "--Results based on " << NUMBER_ITERATIONS << " iterations--\n";
  cout << pABA << ": ";
  for (int i = 0; i < saved_pABA_results.size(); ++i) {
    cout << saved_pABA_results[i] << " ";
  }
  cout << endl;
  cout << "Final " << pABA << ": " << data[pABA.repr()] << endl << endl;

  cout << "Determining the best matching tag sequence:\n";
  vector<string> tags = NotationHelper::Individualize(TAG_SEQUENCES[0]);
  Notation pTW_first("P", OBSERVED_DATA, AND_DELIM, tags);
  Notation *best_pTW = NULL;
  string best_match_string_repr = pTW_first.repr();
  for (string seq : TAG_SEQUENCES) {
    vector<string> tags = NotationHelper::Individualize(seq);
    Notation pTW("P", OBSERVED_DATA, AND_DELIM, tags);
    cout << pTW << ": " << data[pTW.repr()] << endl;
    if (data[pTW.repr()] > data[best_match_string_repr]) {
      best_match_string_repr = pTW.repr();
      delete best_pTW;
      best_pTW = new Notation("P", OBSERVED_DATA, AND_DELIM, tags);
    }
  }
  cout << "The highest probability found belongs to " << best_match_string_repr
    << ": " << data[best_match_string_repr] << endl;
  cout << "The best matching tag sequence is " <<
    NotationHelper::Combine(best_pTW->second) << endl;
  delete best_pTW;

  return 0;
}
