#include "Edge.h"

Edge::Edge(Notation n, Node *src, Node *dest) {
  this->notation = n;
  this->src = src;
  this->dest = dest;
  GraphHelper::LinkNodeAndEdge(src, *this, dest);
}
