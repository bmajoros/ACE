/****************************************************************
 MatchType.H
 bmajoros@duke.edu
 ****************************************************************/

#ifndef INCL_MatchType_H
#define INCL_MatchType_H

enum MatchType {
  MATCH,            // not necessarily identical!
  FIRST_UNMATCHED,  // insertion
  SECOND_UNMATCHED  // deletion
};

#endif
