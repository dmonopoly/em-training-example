#include <cassert>
#include <iostream>
#include <map>
#include <vector>

#include "Helper.h"
#include "Notation.h"

// The number of iterations to do the EM training.
#define NUMBER_ITERATIONS 100

using namespace std;

const string X = "X";
const string Y = "Y";
const string A = "A";
const string B = "B";
const vector<string> TAG_LIST{X, Y};
const vector<string> OBSERVED_DATA{A, B, A};
// Enumerate all possible tag sequences for the brute force method.
const vector<string> TAG_SEQUENCES{
  X+X+X, X+X+Y, X+Y+X, X+Y+Y,
  Y+X+X, Y+X+Y, Y+Y+X, Y+Y+Y};

// For output.
vector<double> saved_pABA_results;

// Known probabilities:
Notation pX("P", {X}, {});  // "probability of x"
Notation pY("P", {Y}, {});
Notation pXGivenX("P", {X}, {X});
Notation pYGivenX("P", {Y}, {X});
Notation pXGivenY("P", {X}, {Y});
Notation pYGivenY("P", {Y}, {Y});
// Objectives:
Notation pABA("P", {A,B,A}, {});
Notation pAGivenX("P", {A}, {X});
Notation pAGivenY("P", {A}, {Y});
Notation pBGivenX("P", {B}, {X});
Notation pBGivenY("P", {B}, {Y});
Notation cXA("C", {X, A}, {});  // "count of x intersected with a"
Notation cXB("C", {X, B}, {});
Notation cYA("C", {Y, A}, {});
Notation cYB("C", {Y, B}, {});


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
    // Get norm P(t,w) and counts.
    for (string seq : TAG_SEQUENCES) {
      vector<string> tags = NotationHelper::Individualize(seq);
      Notation pTW("P", OBSERVED_DATA, tags);
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
  cout << "--Results based on " << NUMBER_ITERATIONS << "iterations--\n";
  cout << pABA << ": ";
  for (int i = 0; i < saved_pABA_results.size(); ++i) {
    cout << saved_pABA_results[i] << " ";
  }
  cout << endl;
  cout << "Final " << pABA << ": " << data[pABA.repr()] << endl << endl;

  cout << "Determining the best matching tag sequence:\n";
  vector<string> tags = NotationHelper::Individualize(TAG_SEQUENCES[0]);
  Notation pTW_first("P", OBSERVED_DATA, tags);
  string best_match_string_repr = pTW_first.repr();
  for (string seq : TAG_SEQUENCES) {
    vector<string> tags = NotationHelper::Individualize(seq);
    Notation pTW("P", OBSERVED_DATA, tags);
    cout << pTW << ": " << data[pTW.repr()] << endl;
    if (data[pTW.repr()] > data[best_match_string_repr])
      best_match_string_repr = pTW.repr();
  }
  cout << "The highest probability found belongs to " << best_match_string_repr
    << ": " << data[best_match_string_repr] << endl;
  cout << "The best matching tag sequence is " <<
    best_match_string_repr.substr(8, 5) << endl;

  return 0;
}
