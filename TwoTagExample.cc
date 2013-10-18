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
#define DO_SHORT_SEQ false
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
  data->emplace(pBGivenX.repr(), 1 - INIT_VAL_pAGivenX);
  data->emplace(pAGivenY.repr(), INIT_VAL_pAGivenY);
  data->emplace(pBGivenY.repr(), 1 - INIT_VAL_pAGivenY);

  // Initial counts can be set to 0.
  data->emplace(cXA.repr(), 0);
  data->emplace(cYA.repr(), 0);
  data->emplace(cXB.repr(), 0);
  data->emplace(cYB.repr(), 0);
}

void ComputeDataWithBruteForce(map<string, double> *data, const Notation &n,
                               const vector<vector<string> > &tag_sequences) {
  saved_obs_seq_probs.push_back((*data)[n.repr()]); // push back initial 0
  vector<Notation> rowOfNots{cXA, cXB, pAGivenX, pBGivenX, cYA, cYB, pAGivenY,
    pBGivenY, n};
  OutputHelper::PrintHeader(rowOfNots);
  OutputHelper::PrintDataRow(0, rowOfNots, *data);

  for (int iter_count = 0; iter_count < NUMBER_ITERATIONS; ++iter_count) {
    // Reset counts to zero.
    (*data)[cXA.repr()] = 0;
    (*data)[cXB.repr()] = 0;
    (*data)[cYA.repr()] = 0;
    (*data)[cYB.repr()] = 0;

    // Get norm P(t,w) and counts.
    double sum_of_all_pTW = 0;  // Use this as divisor in normalization.
    for (vector<string> tags : tag_sequences) {
      Notation pTW("P", OBSERVED_DATA, AND_DELIM, tags);
      double unnormalized_prob = Calculator::ComputeUnnormalizedProbability(pTW,
          *data);
      sum_of_all_pTW += unnormalized_prob;
      (*data)[pTW.repr()] = unnormalized_prob;
    }
    for (vector<string> tags : tag_sequences) {
      Notation pTW("P", OBSERVED_DATA, AND_DELIM, tags);

      // Update counts with *normalized* values. We can also, while we have
      // access to these values, store P(tag sequence|observation seq) (like
      // P(YXY|ABA)), which is P(obs seq \cap tag seq) / P(obs seq). Note that
      // P(obs seq) by total probability is sum_of_all_pTW.
      double normalized_prob = (*data)[pTW.repr()]/sum_of_all_pTW;
      Notation pTGivenW("P", tags, GIVEN_DELIM, OBSERVED_DATA);
      (*data)[pTGivenW.repr()] = normalized_prob;

      // Get counts.
      (*data)[cXA.repr()] += Calculator::NormProbFactor(normalized_prob, pTW,
          cXA);
      (*data)[cXB.repr()] +=
        Calculator::NormProbFactor(normalized_prob, pTW, cXB);
      (*data)[cYA.repr()] += Calculator::NormProbFactor(normalized_prob, pTW,
          cYA);
      (*data)[cYB.repr()] += Calculator::NormProbFactor(normalized_prob, pTW,
          cYB);
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
    OutputHelper::PrintDataRow(iter_count + 1, rowOfNots, *data);
  }
}

void OutputResultsForBruteForce(map<string, double> &data, Notation n,
                                const vector<vector<string> > &tag_sequences) {
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
  vector<string> tags = tag_sequences.at(0);
  Notation pTW_first("P", OBSERVED_DATA, AND_DELIM, tags);
  Notation *best_pTGivenW = NULL;
  string best_match_pTAndW_key = pTW_first.repr();
  string best_match_pTGivenW_key;
  for (vector<string> tags : tag_sequences) {
    Notation pTW("P", OBSERVED_DATA, AND_DELIM, tags);
    Notation pTGivenW("P", tags, GIVEN_DELIM, OBSERVED_DATA);

    if (DO_SHORT_SEQ) { // Only print for short seq; long seq has too many.
      cout << pTW << ": " << data[pTW.repr()] << ", " << pTGivenW << ": " <<
        data[pTGivenW.repr()] << endl;
    }
    if (data[pTW.repr()] > data[best_match_pTAndW_key]) {
      best_match_pTAndW_key = pTW.repr();
      best_match_pTGivenW_key = pTGivenW.repr();
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

void Viterbi(const map<string, double> &data, const vector<Node *> &nodes,
             const vector<string> observed_data) {
  // Key: string representation of node; Value: best value of P(t, w) so far to
  // that node. Best P(t, w) is stored in opt.at(last node).
  map<string, double> opt;
  // Key: string representation of node; Value: previous node string repr.
  map<string, string> best_path;

  // Set all nodes except start to default low value.
  for (Node *n : nodes) {
    opt.emplace(n->repr(), 0);
  }
  opt[nodes.front()->repr()] = 1;

  // Run through trellis. Topological order assumed.
  for (int i = 0; i < nodes.size(); ++i) {
    Node *current_node = nodes.at(i);
    for (Edge *e : current_node->child_edges) {
      Node *next = e->dest;
      double new_val = opt.at(current_node->repr()) * data.at(e->repr());
      if (new_val > opt.at(next->repr())) {
        opt[next->repr()] = new_val;
        best_path[next->repr()] = current_node->repr();
      }
    }
  }

  // Output best result following backpointer map.
  vector<string> best_tag_seq;
  vector<string> assoc_word_seq;
  string next_node_repr = nodes.back()->repr();
  for (int i = 0; i < observed_data.size(); ++i) {
    string name = best_path.at(next_node_repr);
    string tag = name.substr(0, 1);
    string word = name.substr(1, 1);
    best_tag_seq.push_back(tag);
    assoc_word_seq.push_back(word);
    next_node_repr = best_path.at(name); // Skip sister node.
  }
  double best_prob_pTAndW = opt.at(nodes.back()->repr());
  Notation n_best_match_pTAndW("P", assoc_word_seq, AND_DELIM,
      best_tag_seq);
  Notation n_best_match_pTGivenW("P", best_tag_seq, GIVEN_DELIM,
      assoc_word_seq);
  cout << "\n--Viterbi results--\n";
  stringstream ss;
  for (int i = best_tag_seq.size() - 1; i >= 0; --i) {
    ss << best_tag_seq.at(i);
  }
  string best_match_pTAndW_str = ss.str();
  cout << "The highest probability found belongs to " << n_best_match_pTAndW <<
    ": " << best_prob_pTAndW << ", " << n_best_match_pTGivenW << ": " <<
    best_prob_pTAndW/saved_obs_seq_probs.back() << endl;
  cout << "Best matching tag sequence: " << best_match_pTAndW_str << endl;
}

void ForwardBackwardAndViterbi(Notation n, const vector<Node *> &nodes,
                               const vector<Edge *> &select_edges,
                               const vector<Edge *> &all_edges,
                               map<string, double> *data,
                               bool very_small_data_set,
                               const vector<string> observed_data,
                               const vector<string> tag_list) {
  if (EXTRA_PRINTING)
    cout << "Beginning Forward-Backward." << endl;
  saved_obs_seq_probs.push_back((*data)[n.repr()]); // push back initial 0
  vector<Notation> rowOfNots{cXA, cXB, pAGivenX, pBGivenX, cYA, cYB, pAGivenY,
    pBGivenY, n};
  OutputHelper::PrintHeader(rowOfNots);
  OutputHelper::PrintDataRow(0, rowOfNots, *data);

  // PRECONDITION: The order of nodes/edges is already in topological
  // order.
  map<string, double> alpha;  // Sum of all paths from start state to this node.
  map<string, double> beta;  // Sum of all paths from this node to final state.
  alpha[nodes.at(0)->repr()] = 1;
  beta[nodes.at(nodes.size() - 1)->repr()] = 1;
  for (int iter_count = 0; iter_count < NUMBER_ITERATIONS; ++iter_count) {
    // Forward pass. Assumes start node is at i = 0.
    for (int i = 1; i < nodes.size(); ++i) {
      double sum = 0;
      for (Edge *e : nodes[i]->parent_edges) {
        sum += alpha[e->src->repr()] * data->at(e->repr());
      }
      if (EXTRA_PRINTING){
        cout << Basic::Tab(1) << "Alpha value for " << nodes[i]->repr() << ": "
          << sum << endl;
      }
      alpha[nodes[i]->repr()] = sum;
    }
    if (EXTRA_PRINTING)
      cout << endl;

    // Backward pass. Assumes end node is at i = size - 1.
    for (int i = nodes.size() - 2; i >= 0; --i) {
      double sum = 0;
      for (Edge *e : nodes[i]->child_edges) {
        sum += beta[e->dest->repr()] * data->at(e->repr());
      }
      if (EXTRA_PRINTING) {
        cout << Basic::Tab(1) << "Beta value for " << nodes[i]->repr() << ": "
          << sum << endl;
      }
      beta[nodes[i]->repr()] = sum;
    }

    if (EXTRA_PRINTING)
      cout << endl;

    // Counting pass. First reset and then update the counts.
    (*data)[cXA.repr()] = 0;
    (*data)[cXB.repr()] = 0;
    (*data)[cYA.repr()] = 0;
    (*data)[cYB.repr()] = 0;
    for (int i = 0; i < select_edges.size(); ++i) {
      Edge *e = select_edges[i];
      // The count key can be determined by looking at any node this select edge
      // is incident on. We take that node's 'tag' and 'word' fields.
      Notation n_count_key("C", {e->dest->tag, e->dest->word}, AND_DELIM);
      string count_key = n_count_key.repr();
      if (EXTRA_PRINTING) {
        cout << Basic::Tab(1) << "Getting count key: " << count_key << endl;
        cout << Basic::Tab(1) << alpha[e->src->repr()] << "<  >" << data->at(e->repr()) << endl;
        cout << Basic::Tab(1) << beta[e->dest->repr()] << "[ ]" << alpha[nodes.back()->repr()] << endl;
      }

      (*data)[count_key] += (alpha[e->src->repr()] * data->at(e->repr())
                             * beta[e->dest->repr()]) / alpha[nodes.back()->repr()];
    }
    if (EXTRA_PRINTING) {
      cout << endl << Basic::Tab(1) << "'Given' X probabilities: " << endl <<
        Basic::Tab(1) << (*data)[cXA.repr()] << "/" << ( (*data)[cXA.repr()] +
        (*data)[cXB.repr()] ) << 
        Basic::Tab(1) << (*data)[cXB.repr()] << "/" << ( (*data)[cXB.repr()] +
        (*data)[cXB.repr()] );
      cout << endl << Basic::Tab(1) << "'Given' Y probabilities: " << endl <<
        Basic::Tab(1) << (*data)[cYA.repr()] << "/" << ( (*data)[cYA.repr()] +
        (*data)[cYB.repr()] ) << 
        Basic::Tab(1) << (*data)[cYB.repr()] << "/" << ( (*data)[cYB.repr()] +
        (*data)[cYB.repr()] );
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

    // Do not execute the following unless we are dealing with a very small
    // tag set and vocabulary. For real large texts, this would take too long -
    // but it's nice to have this to test on small self-defined texts.
    if (very_small_data_set) {
      vector<vector<string> > tag_sequences =
        TagHandler::GenerateTagSequences(tag_list, observed_data.size());

      // The ultimate value we want to maximize. This should increase with each
      // iteration.
      if (EXTRA_PRINTING)
        cout << endl << Basic::Tab(1) << "Updating probabilities." << endl;
      Calculator::UpdateProbOfObsDataSeq(n, data, tag_sequences);
      saved_obs_seq_probs.push_back((*data)[n.repr()]);
      OutputHelper::PrintDataRow(iter_count + 1, rowOfNots, *data);
    }
  }
  if (EXTRA_PRINTING)
    cout << "Done with Forward-Backward. Proceeding to Viterbi." << endl;

  // By this point, P(A|X), P(A|Y), etc. have been maximized thanks to alpha and
  // beta values in the forward-backward passes. The counting pass collected
  // fractional counts of e.g. C(X|A), which were then used to update the
  // "Given" probabilities (P(A|X), P(A|Y), etc.). Now we use Viterbi to find
  // the highest-probability path based on the collected probabilities.
  Viterbi(*data, nodes, observed_data);
}

void RunBruteForceEM() {
  map<string, double> data;  // Storage for probabilities and counts.
  PrepareInitialData(&data);

  vector<vector<string> > tag_sequences = TagHandler::GenerateTagSequences(TAG_LIST,
      OBSERVED_DATA.size());

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
  map<string, double> data;  // Storage for probabilities and counts.
  vector<Node *> nodes;
  vector<Edge *> edges_to_update;
  vector<Edge *> all_edges; // for deletion later

  vector<vector<string> > tag_sequences =
      TagHandler::GenerateTagSequences(tag_list, observed_data.size());

  PrepareInitialData(&data);
  TrellisAid::BuildTrellis(&nodes, &edges_to_update, &all_edges, observed_data,
      tag_list);

  clock_t t;
  t = clock();
  bool very_small_data_set = true; // TODO issue when this is false
  if (DO_SHORT_SEQ) {
    if (EXTRA_PRINTING)
      cout << "Short sequence: " << endl;
    ForwardBackwardAndViterbi(pABA, nodes, edges_to_update, all_edges, &data,
        very_small_data_set, observed_data, tag_list);
  } else {
    if (EXTRA_PRINTING)
      cout << "Long sequence: " << endl;
    ForwardBackwardAndViterbi(pLong, nodes, edges_to_update, all_edges, &data,
        very_small_data_set, observed_data, tag_list);
  }
  t = clock() - t;
  cout << "\n--Timing Results--\n";
  printf("It took me %lu clicks (%f seconds).\n", t, ((float)t)/CLOCKS_PER_SEC);

  TrellisAid::DestroyTrellis(&nodes, &all_edges);
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    cout << "Running forward-backward and viterbi.\n" << endl;
    RunForwardBackwardAndViterbi(OBSERVED_DATA, TAG_LIST);
  } else if (argc >= 2) {
    cout << "Running brute force.\n" << endl;
    RunBruteForceEM();
  }
  return 0;
}
