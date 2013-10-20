#include "TrellisAid.h"

#define EXTRA_PRINTING true

namespace TrellisAid {
  void BuildTrellis(vector<Node *> *nodes, vector<Edge *> *select_edges,
                    vector<Edge *> *all_edges,
                    const vector<string> &observed_data,
                    const vector<string> &tag_list) {
    if (EXTRA_PRINTING)
      cout << "Building trellis." << endl;
    // Store the last column of nodes to add links to in prev_nodes. Accumulate
    // the next set of prev_nodes in future_prev_nodes.
    vector<Node *> prev_nodes, future_prev_nodes;

    Node *start_node = new Node("START_NODE", 0);
    nodes->push_back(start_node);

    int topol_index = 1;
    for (int i = 0; i < observed_data.size(); ++i)  {
      future_prev_nodes.clear();
      if (EXTRA_PRINTING)
        cout << "obs " << observed_data.at(i) << "--\n";
      for (int j = 0; j < tag_list.size(); ++j) {
        if (EXTRA_PRINTING)
          cout << Basic::Tab(1) << "tag " << tag_list.at(j) << "--\n";
        // Encode unique names for each node. Note that 'first' and 'second' are
        // useful for retrieving the 'sister' node if you only have one.
        string the_tag = tag_list.at(j);
        string the_word = observed_data.at(i);
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
            Notation notation_obj("P", {the_tag}, Notation::GIVEN_DELIM, {p->tag});
            if (EXTRA_PRINTING)
              cout << Basic::Tab(2) << "new edge: " << notation_obj << endl;
            Edge *e = new Edge(notation_obj, p, n1);
            all_edges->push_back(e);
          }
        }
        future_prev_nodes.push_back(n2);
        Notation notation_obj("P", {the_word}, Notation::GIVEN_DELIM,
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
      Edge *e = new Edge(NotationConstants::p1, p, end_node);
      all_edges->push_back(e);
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

  void Viterbi(const map<string, double> &data, const vector<Node *> &nodes,
               const vector<string> observed_data, bool very_small_data_set,
               const vector<double> saved_obs_seq_probs) {
    // Key: string representation of node; Value: best value of P(t, w) so far
    // to that node. Best P(t, w) is stored in opt.at(last node).
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
    // Reverse these because they were retrieved backwards.
    reverse(best_tag_seq.begin(), best_tag_seq.end());
    reverse(assoc_word_seq.begin(), assoc_word_seq.end());

    double best_prob_pTAndW = opt.at(nodes.back()->repr());
    Notation n_best_match_pTAndW("P", assoc_word_seq, Notation::AND_DELIM,
        best_tag_seq);
    Notation n_best_match_pTGivenW("P", best_tag_seq, Notation::GIVEN_DELIM,
        assoc_word_seq);
    cout << "\n--Viterbi results--\n";
    stringstream ss;
    for (int i = 0; i < best_tag_seq.size(); ++i) {
      ss << best_tag_seq.at(i);
    }
    string best_match_pTAndW_str = ss.str();
    cout << "The highest probability found belongs to " << n_best_match_pTAndW
        << ": " << best_prob_pTAndW;
    if (very_small_data_set) {
      cout << ", " << n_best_match_pTGivenW << ": " <<
          best_prob_pTAndW/saved_obs_seq_probs.back() << endl;
    } else {
      cout << endl;
    }
    cout << "Best matching tag sequence: " << best_match_pTAndW_str << endl;
  } // End Viterbi

  void ForwardBackwardAndViterbi(const int num_iterations,
                                 const vector<Node *> &nodes,
                                 const vector<Edge *> &select_edges,
                                 const vector<Edge *> &all_edges,
                                 map<string, double> *data,
                                 bool very_small_data_set, Notation n,
                                 const vector<string> observed_data,
                                 const vector<string> tag_list,
                                 vector<double> *saved_obs_seq_probs) {
    if (EXTRA_PRINTING)
      cout << "Beginning Forward-Backward." << endl;
    if (very_small_data_set) {
      saved_obs_seq_probs->push_back((*data)[n.repr()]); // push back initial 0
    }

    vector<Notation> rowOfNots;
    for (int i = 0; i < select_edges.size(); ++i) {
      Edge *e = select_edges[i];
      Notation n_count_key("C", {e->dest->tag, e->dest->word},
                           Notation::AND_DELIM);
      rowOfNots.push_back(n_count_key);
      rowOfNots.push_back(e->notation);
    }
//     vector<Notation> rowOfNots{NotationConstants::cXA, NotationConstants::cXB,
//                                NotationConstants::pAGivenX,
//                                NotationConstants::pBGivenX,
//                                NotationConstants::cYA, NotationConstants::cYB,
//                                NotationConstants::pAGivenY,
//                                NotationConstants::pBGivenY, n};
    if (very_small_data_set) {
      OutputHelper::PrintHeader(rowOfNots);
      OutputHelper::PrintDataRow(0, rowOfNots, *data);
    }

    // PRECONDITION: The order of nodes/edges is already in topological
    // order.
    map<string, double> alpha;  // Sum of all paths from start to this node.
    map<string, double> beta;  // Sum of all paths from this node to final.
    alpha[nodes.at(0)->repr()] = 1;
    beta[nodes.at(nodes.size() - 1)->repr()] = 1;
    for (int iter_count = 0; iter_count < num_iterations; ++iter_count) {
      // Forward pass. Assumes start node is at i = 0.
      for (int i = 1; i < nodes.size(); ++i) {
        double sum = 0;
        for (Edge *e : nodes[i]->parent_edges) {
          sum += alpha[e->src->repr()] * data->at(e->repr());
        }
        if (EXTRA_PRINTING){
          cout << Basic::Tab(1) << "Alpha value for " << nodes[i]->repr() <<
              ": " << sum << endl;
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
      // The count key can be determined by looking at any node this select
      // edge is incident on. We take that node's 'tag' and 'word' fields. For
      // count_keys, we follow the convention of C(tag, word) (e.g., C(X,A)).
      map<string, double> total_fract_counts ; // Key: tag
      for (int i = 0; i < select_edges.size(); ++i) {
        Edge *e = select_edges[i];
        Notation n_count_key("C", {e->dest->tag, e->dest->word}, Notation::AND_DELIM);
        (*data)[n_count_key.repr()] = 0;
        total_fract_counts[e->dest->tag] = 0;
      }
//       (*data)[cXA.repr()] = 0;
//       (*data)[cXB.repr()] = 0;
//       (*data)[cYA.repr()] = 0;
//       (*data)[cYB.repr()] = 0;
      // Iterate over select edges to update count_keys, used for updating
      // probabilities later.
      for (int i = 0; i < select_edges.size(); ++i) {
        Edge *e = select_edges[i];
        Notation n_count_key("C", {e->dest->tag, e->dest->word}, Notation::AND_DELIM);
        string count_key = n_count_key.repr();
        if (EXTRA_PRINTING) {
          cout << Basic::Tab(1) << "Getting count key: " << count_key << endl;
          cout << Basic::Tab(1) << alpha[e->src->repr()] <<
            "<-alpha edge prob->" << data->at(e->repr()) << endl;
          cout << Basic::Tab(1) << beta[e->dest->repr()] << "<-beta end->" <<
              alpha[nodes.back()->repr()] << endl;
        }
        (*data)[count_key] += (alpha[e->src->repr()] * data->at(e->repr())
                               * beta[e->dest->repr()]) /
                               alpha[nodes.back()->repr()];
        total_fract_counts[e->dest->tag] += (*data)[count_key];
        // TODO: why prob not right... .6 should not be there
        if (count_key == "C(X,A)" || count_key == "C(X,B)") {
          cout << "--adding this to " << e->dest->tag << "," << e->dest->word << ": " <<
            (*data)[count_key] << ". new value: " <<
            total_fract_counts[e->dest->tag] << endl;
        }
      }
//       if (EXTRA_PRINTING) {
//         cout << endl << Basic::Tab(1) << "'Given' X probabilities: " << endl <<
//           Basic::Tab(1) << (*data)[cXA.repr()] << "/" << ( (*data)[cXA.repr()] +
//           (*data)[cXB.repr()] ) << 
//           Basic::Tab(1) << (*data)[cXB.repr()] << "/" << ( (*data)[cXB.repr()] +
//           (*data)[cXB.repr()] );
//         cout << endl << Basic::Tab(1) << "'Given' Y probabilities: " << endl <<
//           Basic::Tab(1) << (*data)[cYA.repr()] << "/" << ( (*data)[cYA.repr()] +
//           (*data)[cYB.repr()] ) << 
//           Basic::Tab(1) << (*data)[cYB.repr()] << "/" << ( (*data)[cYB.repr()] +
//           (*data)[cYB.repr()] );
//       }
      // Update the unknown probabilities that we want to find. Use them in the
      // next iteration.
      for (int i = 0; i < select_edges.size(); ++i) {
        Edge *e = select_edges[i];
        Notation n_count_key("C", {e->dest->tag, e->dest->word}, Notation::AND_DELIM);
        (*data)[e->repr()] = (*data)[n_count_key.repr()] /
                             total_fract_counts.at(e->dest->tag);
        if (n_count_key.repr() == "C(X,A)" || n_count_key.repr() == "C(X,B)")
          cout << "tfc: " << e->dest->tag << ": " << total_fract_counts.at(e->dest->tag) << endl;
      }
//       (*data)[pAGivenX.repr()] = (*data)[cXA.repr()]/( (*data)[cXA.repr()] +
//           (*data)[cXB.repr()] );
//       (*data)[pBGivenX.repr()] = (*data)[cXB.repr()]/( (*data)[cXB.repr()] +
//           (*data)[cXA.repr()] );
//       (*data)[pAGivenY.repr()] = (*data)[cYA.repr()]/( (*data)[cYA.repr()] +
//           (*data)[cYB.repr()] );
//       (*data)[pBGivenY.repr()] = (*data)[cYB.repr()]/( (*data)[cYB.repr()] +
//           (*data)[cYA.repr()] );

      // Do not execute the following unless we are dealing with a very small
      // tag set and vocabulary. For real large texts, this would take too long
      // - but it's nice to have this to test on small self-defined texts.
      if (very_small_data_set) {
        vector<vector<string> > tag_sequences =
          TagHandler::GenerateTagSequences(tag_list, observed_data.size());

        // The ultimate value we want to maximize. This should increase with
        // each iteration.
        if (EXTRA_PRINTING)
          cout << endl << Basic::Tab(1) << "Updating probabilities." << endl;
        Calculator::UpdateProbOfObsDataSeq(n, data, tag_sequences);
        saved_obs_seq_probs->push_back((*data)[n.repr()]);
        OutputHelper::PrintDataRow(iter_count + 1, rowOfNots, *data);
      }
    }
    if (EXTRA_PRINTING)
      cout << "Done with Forward-Backward. Proceeding to Viterbi." << endl;

    // By this point, P(A|X), P(A|Y), etc. have been maximized thanks to alpha
    // and beta values in the forward-backward passes. The counting pass
    // collected fractional counts of e.g. C(X|A), which were then used to
    // update the "Given" probabilities (P(A|X), P(A|Y), etc.). Now we use
    // Viterbi to find the highest-probability path based on the collected
    // probabilities.
    TrellisAid::Viterbi(*data, nodes, observed_data, very_small_data_set,
                        *saved_obs_seq_probs);
  }
} // end namespace TrellisAid
