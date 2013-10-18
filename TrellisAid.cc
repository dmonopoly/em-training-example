#include "TrellisAid.h"

#define EXTRA_PRINTING false

namespace TrellisAid {
  void BuildTrellis(vector<Node *> *nodes, vector<Edge *> *select_edges,
                    vector<Edge *> *all_edges,
                    const vector<string> &observed_data,
                    const vector<string> &tag_list) {
    if (EXTRA_PRINTING)
      cout << "Building trellis." << endl;
    // Store the last column of nodes to add links to in prev_nodes. Accumulate
    // the next set of prev_nodes in fugure_prev_nodes.
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
}
