#ifndef NOTATION_H_
#define NOTATION_H_

#include "Notation.h"

using namespace std;

class Notation;


Notation pXGivenX("P", {X}, Notation::GIVEN_DELIM, {X});
//const Notation pYGivenX("P", {Y}, Notation::GIVEN_DELIM, {X});
//const Notation pXGivenY("P", {X}, Notation::GIVEN_DELIM, {Y});
//const Notation pYGivenY("P", {Y}, Notation::GIVEN_DELIM, {Y});

#endif
