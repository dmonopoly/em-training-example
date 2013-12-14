#include <algorithm>
#include <cassert>
#include <ctime>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>

#include "BasicHelper.h"
#include "NLPHelper.h"
#include "Notation.h"
#include "Node.h"
#include "Edge.h"
#include "TrellisAid.h"

/*  SETTINGS  */
#define DO_SHORT_SEQ true
#define NUMBER_ITERATIONS 50

// Initial values.
// #define INIT_VAL_pAGivenX .7  // Best case for long seq: .7
// #define INIT_VAL_pAGivenY .3  // Best case for long seq: .3
#define INIT_VAL_pAGivenX .5
#define INIT_VAL_pAGivenY .5

#define EXTRA_PRINTING false
/*  END SETTINGS  */

// TODO: Reorganize and use Notation::GIVEN_DELIM. Requires moving all the
// Notation objects below. http://bit.ly/15rbAom
#define GIVEN_DELIM "|"
#define AND_DELIM ","
#define SEQ_DELIM ""

using namespace std;

const string X = "X";
const string Y = "Y";
const string A = "A";
const string B = "B";
const vector<string> TAG_LIST{X, Y};

// Only the brute force functions here use these. We have decoupled these for
// efficient EM / Viterbi.
#if DO_SHORT_SEQ
const vector<string> OBSERVED_DATA{A, B, A};
#else
const vector<string> OBSERVED_DATA{A,A,A,B,A,A,B,A,A};
#endif

// For output. Saves the changing values of pABA (or the longer seq).
vector<double> saved_obs_seq_probs;

Notation p1("P", {"1"});  // Auto-probability 1.
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
Notation cXX("C", {X, X}, AND_DELIM);
Notation cXY("C", {X, Y}, AND_DELIM);
Notation cYX("C", {Y, X}, AND_DELIM);
Notation cYY("C", {Y, Y}, AND_DELIM);
// Long seq type.
Notation pLong("P", {A,A,A,B,A,A,B,A,A}, SEQ_DELIM);

void PrepareInitialData(map<Notation, double> *data) {
  // Given data.
  data->emplace(p1, 1);
  data->emplace(pX, .6);
  data->emplace(pY, .4);
  data->emplace(pXGivenX, .6);
  data->emplace(pYGivenX, .4);
  data->emplace(pXGivenY, .9);
  data->emplace(pYGivenY, .1);

  // Initial value for unknowns. We improve upon these.
  data->emplace(pAGivenX, INIT_VAL_pAGivenX);
  data->emplace(pBGivenX, 1 - INIT_VAL_pAGivenX);
  data->emplace(pAGivenY, INIT_VAL_pAGivenY);
  data->emplace(pBGivenY, 1 - INIT_VAL_pAGivenY);

  // Initial counts can be set to 0.
  data->emplace(cXA, 0);
  data->emplace(cYA, 0);
  data->emplace(cXB, 0);
  data->emplace(cYB, 0);

  // Also make sure LM counts are 0. Actually doesn't matter if select_edges
  // does not include LM edges.
//   data->emplace(cXX, 0);
//   data->emplace(cXY, 0);
//   data->emplace(cYX, 0);
//   data->emplace(cYY, 0);
}

void ComputeDataWithBruteForce(map<Notation, double> *data, const Notation &n,
                               const vector<vector<string> > &tag_sequences) {
  saved_obs_seq_probs.push_back((*data)[n]); // push back initial 0
  vector<Notation> rowOfNots{cXA, cXB, pAGivenX, pBGivenX, cYA, cYB, pAGivenY,
    pBGivenY, n};
  OutputHelper::PrintHeader(rowOfNots);
  OutputHelper::PrintDataRow(0, rowOfNots, *data);

  for (int iter_count = 0; iter_count < NUMBER_ITERATIONS; ++iter_count) {
    // Reset counts to zero.
    (*data)[cXA] = 0;
    (*data)[cXB] = 0;
    (*data)[cYA] = 0;
    (*data)[cYB] = 0;

    // Get norm P(t,w) and counts.
    double sum_of_all_pTW = 0;  // Use this as divisor in normalization.
    for (vector<string> tags : tag_sequences) {
      Notation pTW("P", OBSERVED_DATA, AND_DELIM, tags);
      double unnormalized_prob = Calculator::ComputeUnnormalizedProbability(pTW,
          *data);
      sum_of_all_pTW += unnormalized_prob;
      (*data)[pTW] = unnormalized_prob;
    }
    for (vector<string> tags : tag_sequences) {
      Notation pTW("P", OBSERVED_DATA, AND_DELIM, tags);

      // Update counts with *normalized* values. We can also, while we have
      // access to these values, store P(tag sequence|observation seq) (like
      // P(YXY|ABA)), which is P(obs seq \cap tag seq) / P(obs seq). Note that
      // P(obs seq) by total probability is sum_of_all_pTW.
      double normalized_prob = (*data)[pTW]/sum_of_all_pTW;
      Notation pTGivenW("P", tags, GIVEN_DELIM, OBSERVED_DATA);
      (*data)[pTGivenW] = normalized_prob;

      // Get counts.
      (*data)[cXA] += Calculator::NormProbFactor(normalized_prob, pTW,
          cXA);
      (*data)[cXB] +=
        Calculator::NormProbFactor(normalized_prob, pTW, cXB);
      (*data)[cYA] += Calculator::NormProbFactor(normalized_prob, pTW,
          cYA);
      (*data)[cYB] += Calculator::NormProbFactor(normalized_prob, pTW,
          cYB);
    }
    // Update the unknown probabilities that we want to find. Use them in the
    // next iteration.
    (*data)[pAGivenX] = (*data)[cXA]/( (*data)[cXA] +
        (*data)[cXB] );
    (*data)[pBGivenX] = (*data)[cXB]/( (*data)[cXB] +
        (*data)[cXA] );
    (*data)[pAGivenY] = (*data)[cYA]/( (*data)[cYA] +
        (*data)[cYB] );
    (*data)[pBGivenY] = (*data)[cYB]/( (*data)[cYB] +
        (*data)[cYA] );

    // The ultimate value we want to maximize. This should increase with each
    // iteration.
    Calculator::UpdateProbOfObsDataSeq(n, data, tag_sequences);
    saved_obs_seq_probs.push_back((*data)[n]);
    OutputHelper::PrintDataRow(iter_count + 1, rowOfNots, *data);
  }
}

void OutputResultsForBruteForce(map<Notation, double> &data, Notation n,
                                const vector<vector<string> > &tag_sequences) {
  cout << "\n--Results based on " << NUMBER_ITERATIONS << " iterations--\n";
  ofstream fout("observed_data_probabilities.txt");
  for (int i = 0; i < saved_obs_seq_probs.size(); ++i) {
    fout << saved_obs_seq_probs[i] << endl;
  }
  cout << "Values of " << n << " have been written to "
    "observed_data_probabilities.txt." << endl << endl;

  cout << "Final " << n << ": " << data[n] << endl;
  cout << "Final " << pAGivenX << ": " << data[pAGivenX] << endl;
  cout << "Final " << pBGivenX << ": " << data[pBGivenX] << endl;
  cout << "Final " << pAGivenY << ": " << data[pAGivenY] << endl;
  cout << "Final " << pBGivenY << ": " << data[pBGivenY] << endl << endl;

  cout << "Determining the best matching tag sequence:\n";
  vector<string> tags = tag_sequences.at(0);
  Notation pTW_first("P", OBSERVED_DATA, AND_DELIM, tags);
  Notation *best_pTGivenW = NULL;
  Notation best_match_pTAndW_key = pTW_first;
//   string best_match_pTAndW_key = pTW_first.repr();
  Notation best_match_pTGivenW_key;
//   string best_match_pTGivenW_key;
  for (vector<string> tags : tag_sequences) {
    Notation pTW("P", OBSERVED_DATA, AND_DELIM, tags);
    Notation pTGivenW("P", tags, GIVEN_DELIM, OBSERVED_DATA);

    if (DO_SHORT_SEQ) { // Only print for short seq; long seq has too many.
      cout << pTW << ": " << data[pTW] << ", " << pTGivenW << ": " <<
        data[pTGivenW] << endl;
    }
    if (data[pTW] > data[best_match_pTAndW_key]) {
      best_match_pTAndW_key = pTW;
      best_match_pTGivenW_key = pTGivenW;
      delete best_pTGivenW;
      // Same as pTGivenW. Saved for future reference.
      best_pTGivenW = new Notation("P", tags, GIVEN_DELIM, OBSERVED_DATA);
    }
  }
  
  cout << "The highest probability found belongs to " << best_match_pTAndW_key
    << ": " << data[best_match_pTAndW_key] << ", " << best_match_pTGivenW_key <<
    ": " << data[best_match_pTGivenW_key] << endl;
  cout << "The best matching tag sequence is " <<
    NotationHelper::Combine(best_pTGivenW->first) << endl;
  delete best_pTGivenW;
}

void RunBruteForceEM() {
  map<Notation, double> data;  // Storage for probabilities and counts.
  PrepareInitialData(&data);

  vector<vector<string> > tag_sequences =
    TagHandler::GenerateTagSequences(TAG_LIST, OBSERVED_DATA.size());

  clock_t t;
  t = clock();
  if (DO_SHORT_SEQ) {
    ComputeDataWithBruteForce(&data, pABA, tag_sequences);
    OutputResultsForBruteForce(data, pABA, tag_sequences);
  } else {
    ComputeDataWithBruteForce(&data, pLong, tag_sequences);
    OutputResultsForBruteForce(data, pLong, tag_sequences);
  }
  t = clock() - t;
  cout << "\n--Timing Results--\n";
  printf("It took me %lu clicks (%f seconds).\n", t, ((float)t)/CLOCKS_PER_SEC);
}

void RunForwardBackwardAndViterbi(vector<string> observed_data, vector<string> tag_list) {
  map<Notation, double> data;  // Storage for probabilities and counts.
  vector<Node *> nodes;
  vector<Edge *> edges_to_update;
  vector<Edge *> all_edges; // for deletion later

  PrepareInitialData(&data);
  TrellisAid::BuildTrellis(&nodes, &edges_to_update, &all_edges, observed_data,
      tag_list);

  clock_t t;
  t = clock();
  bool very_small_data_set = true; // tmp
  if (DO_SHORT_SEQ) {
    if (EXTRA_PRINTING)
      cout << "Short sequence: " << endl;
    TrellisAid::ForwardBackwardAndViterbi(NUMBER_ITERATIONS, nodes,
        edges_to_update, all_edges, &data, observed_data, pABA,
        &saved_obs_seq_probs, very_small_data_set);
  } else {
    if (EXTRA_PRINTING)
      cout << "Long sequence: " << endl;
    TrellisAid::ForwardBackwardAndViterbi(NUMBER_ITERATIONS, nodes,
        edges_to_update, all_edges, &data, observed_data, pLong,
        &saved_obs_seq_probs, very_small_data_set);
  }
  t = clock() - t;
  cout << "\n--Timing Results--\n";
  printf("It took me %lu clicks (%f seconds).\n", t, ((float)t)/CLOCKS_PER_SEC);

  TrellisAid::DestroyTrellis(&nodes, &all_edges);
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    cout << "Running Forward-backward and Viterbi." << endl << endl;
    RunForwardBackwardAndViterbi(OBSERVED_DATA, TAG_LIST);
  } else if (argc >= 2) {
    cout << "Running brute force." << endl << endl;
    RunBruteForceEM();
  }
  return 0;
}
