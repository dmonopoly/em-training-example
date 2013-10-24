#include "TrellisAid.h"
#include "NLPHelper.h"

#define EXTRA_PRINTING false

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
    if (EXTRA_PRINTING) {
      cout << "Initializing optimal values." << endl;
    }
    for (Node *n : nodes) {
      opt.emplace(n->repr(), 0);
    }
    opt[nodes.front()->repr()] = 1;

    // Run through trellis. Topological order assumed.
    if (EXTRA_PRINTING) {
      cout << "About to run through trellis." << endl;
    }
    for (int i = 0; i < nodes.size(); ++i) {
      Node *current_node = nodes.at(i);
      for (Edge *e : current_node->child_edges) {
        Node *next = e->dest;
        try {
          double new_val = opt.at(current_node->repr()) * data.at(e->repr());
          if (new_val > opt.at(next->repr())) {
            opt[next->repr()] = new_val;
            best_path[next->repr()] = current_node->repr();
          }
        } catch (out_of_range &e) {
          cerr << "Out of range error in Viterbi while going through trellis: "
            << e.what() << endl;
          exit(0);
        }
      }
    }
    if (EXTRA_PRINTING) {
      cout << "Done running through trellis." << endl;
    }

    // Output best result following backpointer map.
    if (EXTRA_PRINTING) {
      cout << "Printing result from backpointer map." << endl;
    }
    vector<string> best_tag_seq;
    vector<string> assoc_word_seq;
    string next_node_repr = nodes.back()->repr();
    for (int i = 0; i < observed_data.size(); ++i) {
      string name;
      try {
        name = best_path.at(next_node_repr);
      } catch (out_of_range &e) {
        cerr << "Out of range error in Viterbi while getting name: " <<
          e.what() << endl;
      }
      string tag = name.substr(0, 1);
      string word = name.substr(1, 1);
      best_tag_seq.push_back(tag);
      assoc_word_seq.push_back(word);
      try {
        next_node_repr = best_path.at(name); // Skip sister node.
      } catch (out_of_range &e) {
        cerr << "Out of range error in Viterbi while getting next " <<
          "node from best path: " << e.what() << endl;
      }
    }
    // Reverse these because they were retrieved backwards.
    reverse(best_tag_seq.begin(), best_tag_seq.end());
    reverse(assoc_word_seq.begin(), assoc_word_seq.end());

    double best_prob_pTAndW;
    best_prob_pTAndW = opt.at(nodes.back()->repr());
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
                                 const vector<string> observed_data,
                                 Notation nObsSeq,
                                 vector<double> *saved_obs_seq_probs,
                                 bool very_small_data_set) {
    if (EXTRA_PRINTING)
      cout << "Beginning Forward-Backward." << endl;
    // Key: count_key repr, like "C(X,A)". Value: true if already used. Main
    // purpose is for checking while accumulating fractional counts in counting
    // pass, but also used in output setup.
    unordered_map<string, bool> already_used;

    // Push back initial 0.
    saved_obs_seq_probs->push_back((*data)[nObsSeq.repr()]);
    vector<Notation> rowOfNots;
    if (very_small_data_set) {
      // Prepare row of Notation strings for nice column-organized output.
      for (int i = 0; i < select_edges.size(); ++i) {
        Edge *e = select_edges[i];
        Notation n_count_key("C", {e->dest->tag, e->dest->word},
            Notation::AND_DELIM);
        if (!already_used[n_count_key.repr()]) {
          rowOfNots.push_back(n_count_key);
          rowOfNots.push_back(e->notation);
          already_used[n_count_key.repr()] = true;
        }
      }
      rowOfNots.push_back(nObsSeq);
      OutputHelper::PrintHeader(rowOfNots);
      OutputHelper::PrintDataRow(0, rowOfNots, *data);
      already_used.clear();
    }

    // PRECONDITION: The order of nodes/edges is already in topological order.
    map<string, double> alpha;  // Sum of all paths from start to this node.
    map<string, double> beta;  // Sum of all paths from this node to final.
    alpha[nodes.at(0)->repr()] = 1;
    beta[nodes.at(nodes.size() - 1)->repr()] = 1;
    for (int iter_count = 0; iter_count < num_iterations; ++iter_count) {
      if (EXTRA_PRINTING)
        cout << "Forward pass... ";
      // Forward pass. Assumes start node is at i = 0.
      try {
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
      } catch (out_of_range &e) {
        cerr << "Out of range error in forward pass: " << e.what() << endl;
      } catch (exception &e) {
        cerr << "Issue in forward pass: " << e.what() << endl;
      }

      if (EXTRA_PRINTING) {
        cout << "Backward pass... ";
      }
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

      if (EXTRA_PRINTING) {
        cout << "Counting pass... " << endl;
      }

      // Counting pass.  The count key can be determined by looking at any node
      // this select edge is incident on. We take that node's 'tag' and 'word'
      // fields. For count_keys, we follow the convention of C(tag, word) (e.g.,
      // C(X,A)).

      // First reset the counts.
      for (int i = 0; i < select_edges.size(); ++i) {
        Edge *e = select_edges[i];
        Notation n_count_key("C", {e->dest->tag, e->dest->word},
                             Notation::AND_DELIM);
        (*data)[n_count_key.repr()] = 0;
      }

      // Key: tag. Value: total fractional count associated with that tag.
      unordered_map<string, double> total_fract_counts;

      // Iterate over select edges to update count_keys, used for updating
      // probabilities later.
      for (int i = 0; i < select_edges.size(); ++i) {
        Edge *e = select_edges[i];
        Notation n_count_key("C", {e->dest->tag, e->dest->word},
                             Notation::AND_DELIM);
        string count_key = n_count_key.repr();
        if (EXTRA_PRINTING) {
          cout << Basic::Tab(1) << "Getting count key from edge " << e->repr()
              << ": " << count_key << endl;
          cout << Basic::Tab(1) << alpha[e->src->repr()] <<
              "<-alpha edge prob->" << data->at(e->repr()) << endl;
          cout << Basic::Tab(1) << beta[e->dest->repr()] << "<-beta end->" <<
              alpha[nodes.back()->repr()] << endl;
        }
        // If we already saw this count key, subtract all previously accumulated
        // values for that count_key from total_fract_counts before updating
        // (*data)[count_key]. This compensates for adding too-early cXA's when
        // computing the total fractional counts for C(X,w_i).
        if (already_used[count_key])
          total_fract_counts[e->dest->tag] -= (*data)[count_key];
        (*data)[count_key] += (alpha[e->src->repr()] * data->at(e->repr())
                               * beta[e->dest->repr()]) /
                               alpha[nodes.back()->repr()];
        already_used[count_key] = true;
        total_fract_counts[e->dest->tag] += (*data)[count_key];
      }

      // Update the unknown probabilities that we want to find. Use them in the
      // next iteration.
      for (int i = 0; i < select_edges.size(); ++i) {
        Edge *e = select_edges[i];
        Notation n_count_key("C", {e->dest->tag, e->dest->word}, Notation::AND_DELIM);
        (*data)[e->repr()] = (*data)[n_count_key.repr()] /
                             total_fract_counts.at(e->dest->tag);
      }

      // Update probability of observed data sequence. This should increase
      // with each iteration.
      (*data)[nObsSeq.repr()] = alpha[nodes.back()->repr()];
      saved_obs_seq_probs->push_back((*data)[nObsSeq.repr()]);

      // If we have a small data set, print neatly organized columns of output.
      // Do not execute the following unless we are dealing with a very small
      // tag set and vocabulary. For real (large) texts, the resulting output
      // table would have too many columns.
      if (very_small_data_set) {
        OutputHelper::PrintDataRow(iter_count + 1, rowOfNots, *data);
      } else {
        // Print viterbi.
        if (saved_obs_seq_probs != NULL)
          TrellisAid::Viterbi(*data, nodes, observed_data, very_small_data_set,
              *saved_obs_seq_probs);
        else {
          TrellisAid::Viterbi(*data, nodes, observed_data, very_small_data_set,
              {});
        }
        // Print P(obs).
        cout << alpha[nodes.back()->repr()] << endl;
      }
    }
    if (EXTRA_PRINTING) {
      cout << "Done with Forward-Backward. Proceeding to Viterbi." << endl;
    }

    // By this point, P(A|X), P(A|Y), etc. have been maximized thanks to alpha
    // and beta values in the forward-backward passes. The counting pass
    // collected fractional counts of e.g. C(X|A), which were then used to
    // update the "Given" probabilities (P(A|X), P(A|Y), etc.). Now we use
    // Viterbi to find the highest-probability path based on the collected
    // probabilities and print the result.
    if (saved_obs_seq_probs != NULL)
      TrellisAid::Viterbi(*data, nodes, observed_data, very_small_data_set,
          *saved_obs_seq_probs);
    else {
      TrellisAid::Viterbi(*data, nodes, observed_data, very_small_data_set,
          {});
    }
  }
} // end namespace TrellisAid

