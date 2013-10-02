#include <cassert>
#include <ctime>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>

#include "NLPHelper.h"
#include "Notation.h"
#include "Node.h"
#include "Edge.h"

/*  SETTINGS  */
// The number of iterations to do the EM training.
#define NUMBER_ITERATIONS 100

// Two ways to run this program: with a short or a long observed sequence.
// Applies only to brute force method.
#define DO_SHORT_SEQ false

// Initial values.
#define INIT_VAL_pAGivenX .3
#define INIT_VAL_pAGivenY .7
#define INIT_VAL_pBGivenX .3
#define INIT_VAL_pBGivenY .7
/*  END SETTINGS  */

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
// Long seq type.
Notation pLong("P", {A,A,A,B,A,A,B,A,A}, SEQ_DELIM);

void PrepareInitialData(map<string, double> *data) {
  // Given data.
  data->emplace(p1.repr(), 1);
  data->emplace(pX.repr(), .6);
  data->emplace(pY.repr(), .4);
  data->emplace(pXGivenX.repr(), .6);
  data->emplace(pYGivenX.repr(), .4);
  data->emplace(pXGivenY.repr(), .9);
  data->emplace(pYGivenY.repr(), .1);

  // Initial value for unknowns. We improve upon these.
  data->emplace(pAGivenX.repr(), INIT_VAL_pAGivenX);
  data->emplace(pAGivenY.repr(), INIT_VAL_pAGivenY);
  data->emplace(pBGivenX.repr(), INIT_VAL_pBGivenX);
  data->emplace(pBGivenY.repr(), INIT_VAL_pBGivenY);

  // Initial counts can be set to 0.
  data->emplace(cXA.repr(), 0);
  data->emplace(cYA.repr(), 0);
  data->emplace(cXB.repr(), 0);
  data->emplace(cYB.repr(), 0);
}

void ComputeDataWithBruteForce(map<string, double> *data, const Notation &n,
                               const vector<string> &tag_sequences) {
  saved_obs_seq_probs.push_back((*data)[n.repr()]); // push back initial 0

  vector<Notation> rowOfNots{cXA, cXB, pAGivenX, pBGivenX, cYA, cYB, pAGivenY,
    pBGivenY, n};
  OutputHelper::PrintHeader(rowOfNots);
  OutputHelper::PrintDataRow(0, rowOfNots, *data);

  for (int i = 0; i < NUMBER_ITERATIONS; ++i) {
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
    saved_obs_seq_probs.push_back((*data)[n.repr()]);
    OutputHelper::PrintDataRow(i + 1, rowOfNots, *data);
  }
}

void OutputResults(map<string, double> &data, Notation n, const vector<string>
                   &tag_sequences) {
  cout << "\n--Results based on " << NUMBER_ITERATIONS << " iterations--\n";
  ofstream fout("observed_data_probabilities.txt");
  for (int i = 0; i < saved_obs_seq_probs.size(); ++i) {
    fout << saved_obs_seq_probs[i] << endl;
  }
  cout << "Values of " << n << " have been written to "
    "observed_data_probabilities.txt." << endl << endl;

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
    Notation pTGivenW("P", tags, GIVEN_DELIM, OBSERVED_DATA);
    // Compute P(t|w). Technically not used because divided values seem to
    // incorrectly yield >1 (decimals too small, possibly).
    data[pTGivenW.repr()] = data[pTW.repr()] / data[n.repr()];
    if (DO_SHORT_SEQ) { // Only print for short seq; long seq has too many.
      cout << pTW << ": " << data[pTW.repr()] << endl;
    }
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

// Warning: Creates data on heap. Call DestroyTrellis after done.
// Post: 'nodes' points to a vector where [0] is the start node, back() is the
// end, and the vector lists the nodes in topological order. 'edges' points to a
// vector of corresponding edges.
void BuildTrellis(vector<Node *> *nodes, vector<Edge *> *select_edges, vector<Edge *> *all_edges) {
  Node *start_node = new Node("start", 0);
  nodes->push_back(start_node);
  Node *xa1first = new Node("xa1first", 1);
  nodes->push_back(xa1first);
  Node *ya1first = new Node("ya1first", 1);
  nodes->push_back(ya1first);
  Node *xa2first = new Node("xa2first", 2);
  nodes->push_back(xa2first);
  Node *ya2first = new Node("ya2first", 2);
  nodes->push_back(ya2first);
  Node *xb1second = new Node("xb1second", 3);
  nodes->push_back(xb1second);
  Node *yb1second = new Node("yb1second", 3);
  nodes->push_back(yb1second);
  Node *xb2second = new Node("xb2second", 4);
  nodes->push_back(xb2second);
  Node *yb2second = new Node("yb2second", 4);
  nodes->push_back(yb2second);
  Node *xa1third = new Node("xa1third", 5);
  nodes->push_back(xa1third);
  Node *ya1third = new Node("ya1third", 5);
  nodes->push_back(ya1third);
  Node *xa2third = new Node("xa2third", 6);
  nodes->push_back(xa2third);
  Node *ya2third = new Node("ya2third", 6);
  nodes->push_back(ya2third);
  Node *end_node = new Node("end", 7);
  nodes->push_back(end_node);

  // From the start point.
  Edge *pX_edge = new Edge(pX, start_node, xa1first);
  all_edges->push_back(pX_edge);
  Edge *pY_edge = new Edge(pY, start_node, ya1first);
  all_edges->push_back(pY_edge);

  // The diagonals.
  Edge *pYGivenX_edge = new Edge(pYGivenX, xa2first, yb1second);
  all_edges->push_back(pYGivenX_edge);
  Edge *pXGivenY_edge = new Edge(pXGivenY, ya2first, xb1second);
  all_edges->push_back(pXGivenY_edge);
  Edge *pYGivenX_edge2 = new Edge(pYGivenX, xb2second, ya1third);
  all_edges->push_back(pYGivenX_edge2);
  Edge *pXGivenY_edge2 = new Edge(pXGivenY, yb2second, xa1third);
  all_edges->push_back(pXGivenY_edge2);

  // Across the top.
  Edge *pAGivenX_edge = new Edge(pAGivenX, xa1first, xa2first);
  all_edges->push_back(pAGivenX_edge);
  Edge *pXGivenX_edge = new Edge(pXGivenX, xa2first, xb1second);
  all_edges->push_back(pXGivenX_edge);
  Edge *pBGivenX_edge = new Edge(pBGivenX, xb1second, xb2second);
  all_edges->push_back(pBGivenX_edge);
  Edge *pXGivenX_edge2 = new Edge(pXGivenX, xb2second, xa1third);
  all_edges->push_back(pXGivenX_edge2);
  Edge *pAGivenX_edge2 = new Edge(pAGivenX, xa1third, xa2third);
  all_edges->push_back(pAGivenX_edge2);

  // Across the bottom.
  Edge *pAGivenY_edge = new Edge(pAGivenY, ya1first, ya2first);
  all_edges->push_back(pAGivenY_edge);
  Edge *pYGivenY_edge = new Edge(pYGivenY, ya2first, yb1second);
  all_edges->push_back(pYGivenY_edge);
  Edge *pBGivenY_edge = new Edge(pBGivenY, yb1second, yb2second);
  all_edges->push_back(pBGivenY_edge);
  Edge *pYGivenY_edge2 = new Edge(pYGivenY, yb2second, ya1third);
  all_edges->push_back(pYGivenY_edge2);
  Edge *pAGivenY_edge2 = new Edge(pAGivenY, ya1third, ya2third);
  all_edges->push_back(pAGivenY_edge2);

  // To the end point.
  Edge *x_last_edge = new Edge(p1, xa2third, end_node);
  all_edges->push_back(x_last_edge);
  Edge *y_last_edge = new Edge(p1, ya2third, end_node);
  all_edges->push_back(y_last_edge);

  select_edges->push_back(pAGivenX_edge);
  select_edges->push_back(pAGivenY_edge);
  select_edges->push_back(pBGivenX_edge);
  select_edges->push_back(pBGivenY_edge);
}

void DestroyTrellis(vector<Node *> *nodes, vector<Edge *> *all_edges) {
  // Deletes nodes and edges.
  for (Node *n : *nodes) {
    delete n;
  }
  for (Edge *e : *all_edges) {
    delete e;
  }
}

void ForwardBackwardCompute(const vector<Node *> &nodes, 
                            const vector<Edge *> &select_edges,
                            map<string, double> *data) {
  map<string, double> alpha;  // Sum of all paths from start state to this node.
  map<string, double> beta;  // Sum of all paths from this node to final state.

  assert(OBSERVED_DATA.size() == 3 && "Generating wrong sized tag sequence for "
      "forward-backward. Should only be the small ABA sequence.");
  vector<string> tag_sequences = TagHandler::GenerateTagSequences(TAG_LIST,
      OBSERVED_DATA.size());

  // Set start node alpha value to 1.
  alpha[nodes.at(0)->repr()] = 1;

  // Forward pass. Assumes start node is at i = 0.
  for (int i = 1; i < nodes.size(); ++i) {
    double sum = 0;
    for (Edge *e : nodes[i]->parent_edges) {
      sum += alpha[e->src->repr()] * data->at(e->repr());
    }
    alpha[nodes[i]->repr()] = sum;
  }

  // Backward pass. TODO.

  // Counting pass. First reset and then update the counts.
  (*data)[cXA.repr()] = 0;
  (*data)[cXB.repr()] = 0;
  // need cAY too?
  for (int i = 0; i < select_edges.size(); ++i) {
    Edge *e = select_edges[i];
    string count_key = NotationHelper::GetCountKey(e->repr());
    cout << "Getting count key " << count_key << " from " << e->repr() << endl;
    (*data)[count_key] += (alpha[e->src->repr()] * data->at(e->repr())
                           * beta[e->dest->repr()]) / alpha[nodes.back()->repr()];
  }

  // Update the unknown probabilities that we want to find. Use them in the
  // next iteration.
  (*data)[pAGivenX.repr()] = (*data)[cXA.repr()]/( (*data)[cXA.repr()] +
      (*data)[cXB.repr()] );
  (*data)[pBGivenX.repr()] = (*data)[cXB.repr()]/( (*data)[cXB.repr()] +
      (*data)[cXA.repr()] );

  // The ultimate value we want to maximize. This should increase with each
  // iteration.
  Calculator::UpdateProbOfObsDataSeq(pABA, data, tag_sequences);
  cout << cXA << ": " << (*data)[cXA.repr()] << endl;
  cout << cXB << ": " << (*data)[cXB.repr()] << endl;
  cout << pABA << ": " << (*data)[pABA.repr()] << endl;
}

void RunBruteForceEM() {
  map<string, double> data;  // Storage for probabilities and counts.
  PrepareInitialData(&data);

  vector<string> tag_sequences = TagHandler::GenerateTagSequences(TAG_LIST,
      OBSERVED_DATA.size());

  clock_t t;
  t = clock();
  if (DO_SHORT_SEQ) {
    ComputeDataWithBruteForce(&data, pABA, tag_sequences);
    OutputResults(data, pABA, tag_sequences);
  } else {
    ComputeDataWithBruteForce(&data, pLong, tag_sequences);
    OutputResults(data, pLong, tag_sequences);
  }
  t = clock() - t;
  cout << "--Timing Results--\n";
  printf("It took me %lu clicks (%f seconds).\n", t, ((float)t)/CLOCKS_PER_SEC);
}

void RunEfficientEM() {
  map<string, double> data;  // Storage for probabilities and counts.
  vector<Node *> nodes;
  vector<Edge *> edges_to_update;
  vector<Edge *> all_edges; // for deletion later

  PrepareInitialData(&data);
  BuildTrellis(&nodes, &edges_to_update, &all_edges);

  clock_t t;
  t = clock();
  ForwardBackwardCompute(nodes, edges_to_update, &data);
  t = clock() - t;
  
  cout << "--Timing Results--\n";
  printf("It took me %lu clicks (%f seconds).\n", t, ((float)t)/CLOCKS_PER_SEC);
  DestroyTrellis(&nodes, &all_edges);
}

int main() {
  RunBruteForceEM();
  //RunEfficientEM();
  return 0;
}
