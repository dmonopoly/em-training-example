#ifndef NOTATION_CONSTANTS_TWO_TAG_H_
#define NOTATION_CONSTANTS_TWO_TAG_H_

#include <initializer_list>
#include <iostream>
#include <vector>

namespace NotationConstants {
  const Notation pAGivenX("P", {A}, GIVEN_DELIM, {X});
  const Notation pAGivenY("P", {A}, GIVEN_DELIM, {Y});
  const Notation pBGivenX("P", {B}, GIVEN_DELIM, {X});
  const Notation pBGivenY("P", {B}, GIVEN_DELIM, {Y});
  const Notation cXA("C", {X, A}, AND_DELIM);  // "count of x and a"
  const Notation cXB("C", {X, B}, AND_DELIM);
  const Notation cYA("C", {Y, A}, AND_DELIM);
  const Notation cYB("C", {Y, B}, AND_DELIM);
}

#endif
