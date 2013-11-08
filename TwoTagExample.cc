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

// WARNING: Creates data on heap. Call DestroyTrellis when done.
// Post: 'nodes' points to a vector where front() is the start node, back() is
// the end, and the vector lists the nodes in topological order. Each node has a
// unique name. 'edges' points to a vector of corresponding edges with
// representations like P(A|X). **Nodes and edges are in topological order!**.
void BuildTrellis(vector<Node *> *nodes, vector<Edge *> *select_edges,
                  vector<Edge *> *all_edges) {
  if (EXTRA_PRINTING)
    cout << "Building trellis." << endl;
  // Store the last column of nodes to add links to in prev_nodes. Accumulate
  // the next set of prev_nodes in fugure_prev_nodes.
  vector<Node *> prev_nodes, future_prev_nodes;

  Node *start_node = new Node("START_NODE", 0);
  nodes->push_back(start_node);

  int topol_index = 1;
  for (int i = 0; i < OBSERVED_DATA.size(); ++i)  {
    future_prev_nodes.clear();
    if (EXTRA_PRINTING)
      cout << "obs " << OBSERVED_DATA[i] << "--\n";
    for (int j = 0; j < TAG_LIST.size(); ++j) {
      if (EXTRA_PRINTING)
        cout << Basic::Tab(1) << "tag " << TAG_LIST[j] << "--\n";
      // Encode unique names for each node. Note that 'first' and 'second' are
      // useful for retrieving the 'sister' node if you only have one.
      string the_tag = TAG_LIST[j];
      string the_word = OBSERVED_DATA[i];
      stringstream ss;
      ss << the_tag << the_word << i << j << "first";
      Node *n1 = new Node(ss.str(), topol_index, the_tag, the_word);
      ss.clear(); ss.str("");
      ss << the_tag << the_word << i << j << "second";
      Node *n2 = new Node(ss.str(), topol_index + 1, the_tag, the_word);
      if (EXTRA_PRINTING) {
        cout << Basic::Tab(2) << "node name: " << n1->repr() << endl;
        cout << Basic::Tab(2) << "node name: " << n2->repr() << endl;
      }
      nodes->push_back(n1);
      nodes->push_back(n2);
      if (topol_index == 1) {
        Notation notation_obj("P", {the_tag});
        Edge *e = new Edge(notation_obj, start_node, n1);
        all_edges->push_back(e);
      } else {
        for (Node *p : prev_nodes) {
          // P(t2|t1)
          Notation notation_obj("P", {the_tag}, GIVEN_DELIM, {p->tag});
          if (EXTRA_PRINTING)
            cout << Basic::Tab(2) << "new edge: " << notation_obj << endl;
          Edge *e = new Edge(notation_obj, p, n1);
          all_edges->push_back(e);
        }
      }
      future_prev_nodes.push_back(n2);
      Notation notation_obj("P", {the_word}, GIVEN_DELIM,
          {the_tag});
      Edge *e = new Edge(notation_obj, n1, n2);
      select_edges->push_back(e);
      all_edges->push_back(e);
    }
    topol_index += 2;
    prev_nodes = future_prev_nodes;
  }

  // To the end point. It is important that these have probability 1 for the
  Node *end_node = new Node("END_NODE", topol_index);
  nodes->push_back(end_node);
  for (Node *p : prev_nodes) {
    Edge *e = new Edge(p1, p, end_node);
    all_edges->push_back(e);
  }

  if (DO_SHORT_SEQ) {
    assert(all_edges->size() == 18 && "Incorrect number of edges");
    assert(select_edges->size() == 6 && "Incorrect number of edges");
    assert(nodes->size() == 14 && "Incorrect number of nodes");
  }
  if (EXTRA_PRINTING)
    cout << "Done building trellis." << endl;
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

void Viterbi(const map<string, double> &data, const vector<Node *> &nodes) {
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
  for (int i = 0; i < OBSERVED_DATA.size(); ++i) {
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
                            map<string, double> *data) {
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

    vector<vector<string> > tag_sequences = TagHandler::GenerateTagSequences(TAG_LIST,
        OBSERVED_DATA.size());

    // The ultimate value we want to maximize. This should increase with each
    // iteration.
    if (EXTRA_PRINTING)
      cout << endl << Basic::Tab(1) << "Updating probabilities." << endl;
    Calculator::UpdateProbOfObsDataSeq(n, data, tag_sequences);
    saved_obs_seq_probs.push_back((*data)[n.repr()]);
    OutputHelper::PrintDataRow(iter_count + 1, rowOfNots, *data);
  }
  if (EXTRA_PRINTING)
    cout << "Done with Forward-Backward. Proceeding to Viterbi." << endl;

  // By this point, P(A|X), P(A|Y), etc. have been maximized thanks to alpha and
  // beta values in the forward-backward passes. The counting pass collected
  // fractional counts of e.g. C(X|A), which were then used to update the
  // "Given" probabilities (P(A|X), P(A|Y), etc.). Now we use Viterbi to find
  // the highest-probability path based on the collected probabilities.
  Viterbi(*data, nodes);
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

void RunForwardBackwardAndViterbi() {
  map<string, double> data;  // Storage for probabilities and counts.
  vector<Node *> nodes;
  vector<Edge *> edges_to_update;
  vector<Edge *> all_edges; // for deletion later

  vector<vector<string> > tag_sequences = TagHandler::GenerateTagSequences(TAG_LIST,
      OBSERVED_DATA.size());

  PrepareInitialData(&data);
  BuildTrellis(&nodes, &edges_to_update, &all_edges);

  clock_t t;
  t = clock();
  if (DO_SHORT_SEQ) {
    if (EXTRA_PRINTING)
      cout << "Short sequence: " << endl;
    ForwardBackwardAndViterbi(pABA, nodes, edges_to_update, all_edges, &data);
  } else {
    if (EXTRA_PRINTING)
      cout << "Long sequence: " << endl;
    ForwardBackwardAndViterbi(pLong, nodes, edges_to_update, all_edges, &data);
  }
  t = clock() - t;
  cout << "\n--Timing Results--\n";
  printf("It took me %lu clicks (%f seconds).\n", t, ((float)t)/CLOCKS_PER_SEC);

  DestroyTrellis(&nodes, &all_edges);
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    cout << "Running forward-backward and viterbi.\n" << endl;
    RunForwardBackwardAndViterbi();
  } else if (argc >= 2) {
    cout << "Running brute force.\n" << endl;
    RunBruteForceEM();
  }
  return 0;
}
