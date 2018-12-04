/* -*- mode: C++; indent-tabs-mode: nil; -*-
 *
 * This file is a part of LEMON, a generic C++ optimization library.
 *
 * Copyright (C) 2003-2013
 * Egervary Jeno Kombinatorikus Optimalizalasi Kutatocsoport
 * (Egervary Research Group on Combinatorial Optimization, EGRES).
 *
 * Permission to use, modify and distribute this software is granted
 * provided that this copyright notice appears in all copies. For
 * precise terms see the accompanying LICENSE file.
 *
 * This software is provided "AS IS" with no warranty of any kind,
 * express or implied, and with no claim as to its suitability for any
 * purpose.
 *
 */

#include <lemon/unionfind.h>

#include "test_tools.h"
#include "lemon/maps.h"

using namespace lemon;
using namespace std;

int main() {

  typedef RangeMap<string> RangeMap;
  typedef HeapUnionFind<int, RangeMap> HeapUnionFind;

  RangeMap map = RangeMap(4);
  map[0] = "zero";
  map[1] = "one";
  map[2] = "two";
  map[3] = "three";

  
  std::cerr << "done" << std::endl;
}
