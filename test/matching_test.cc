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

#include <iostream>
#include <sstream>
#include <vector>
#include <queue>
#include <cstdlib>

#include <lemon/matching.h>
#include <lemon/smart_graph.h>
#include <lemon/concepts/graph.h>
#include <lemon/concepts/maps.h>
#include <lemon/lgf_reader.h>
#include <lemon/math.h>

#include "test_tools.h"

using namespace std;
using namespace lemon;

GRAPH_TYPEDEFS(SmartGraph);

const std::string lgf0 =
"@nodes\n"
"label\n"
"0\n"
"1\n"
"2\n"
"3\n"
"4\n"
"5\n"
"6\n"
"7\n"
"@edges\n"
"     label  weight\n"
"7 4  0      984\n"
"0 7  1      73\n"
"7 1  2      204\n"
"2 3  3      583\n"
"2 7  4      565\n"
"2 1  5      582\n"
"0 4  6      551\n"
"2 5  7      385\n"
"1 5  8      561\n"
"5 3  9      484\n"
"7 5  10     904\n"
"3 6  11     47\n"
"7 6  12     888\n"
"3 0  13     747\n"
"6 1  14     310\n";


const std::string lgf1 =
"@nodes\n"
"label\n"
"0\n"
"1\n"
"2\n"
"3\n"
"4\n"
"5\n"
"6\n"
"7\n"
"@edges\n"
"     label  weight\n"
"2 5  0      710\n"
"0 5  1      241\n"
"2 4  2      856\n"
"2 6  3      762\n"
"4 1  4      747\n"
"6 1  5      962\n"
"4 7  6      723\n"
"1 7  7      661\n"
"2 3  8      376\n"
"1 0  9      416\n"
"6 7  10     391\n";

const std::string lgf2 =
"@nodes\n"
"label\n"
"0\n"
"1\n"
"2\n"
"3\n"
"4\n"
"5\n"
"6\n"
"7\n"
"@edges\n"
"     label  weight\n"
"6 2  0      553\n"
"0 7  1      653\n"
"6 3  2      22\n"
"4 7  3      846\n"
"7 2  4      981\n"
"7 6  5      250\n"
"5 2  6      539\n";

const std::string graph4 =
"@nodes\n"
"label\n"
"0\n"
"1\n"
"2\n"
"3\n"
"@edges\n"
"     label  weight\n"
"0 1  0      10\n"
"1 2  1      20\n"
"2 3  2      30\n"
"3 0  3      40\n"
"0 2  4      50\n"
;

const int lgfn = 3;
const std::string lgf[lgfn] = {
  lgf0,
  lgf1,
  lgf2
};

void checkWeightedPerfectMatching(const SmartGraph& graph,
                          const SmartGraph::EdgeMap<int>& weight,
                          const MaxWeightedPerfectMatching<SmartGraph>& mwpm) {
  int nBlossoms = mwpm.blossomNum();
  //std::cerr << nBlossoms << " blossoms";
  for (SmartGraph::EdgeIt e(graph); e != INVALID; ++e) {
    if (graph.u(e) == graph.v(e)) continue;
    int reducedWeight = mwpm.nodeValue(graph.u(e)) + mwpm.nodeValue(graph.v(e));
    if (mwpm.matching(e)) {
      std::cerr << "IN MATCHING: ";
    }
    std::cerr << "weight: " << weight[e] << " reduced: " << reducedWeight << std::endl;
    for (int i = 0; i < mwpm.blossomNum(); ++i) {
      bool s = false, t = false;
      for (MaxWeightedPerfectMatching<SmartGraph>::BlossomIt n(mwpm, i);
           n != INVALID; ++n) {
        if (graph.u(e) == n) s = true;
        if (graph.v(e) == n) t = true;
      }
      if (s == true && t == true) {
        int blossomValue = mwpm.blossomValue(i);
        //std::cerr << "   adding blossom value " + blossomValue << std::endl;
        reducedWeight += blossomValue;
      }
    }
    reducedWeight -= weight[e] * mwpm.dualScale;
    //std::cerr << "  final: " << reducedWeight << std::endl;
    check(reducedWeight >= 0, "Negative reduced weight");
    check(reducedWeight == 0 || !mwpm.matching(e),
          "Non-zero reduced weight on matching edge");
  }

  int pv = 0;
  for (SmartGraph::NodeIt n(graph); n != INVALID; ++n) {
    check(mwpm.matching(n) != INVALID, "Non perfect");
    pv += weight[mwpm.matching(n)];
    SmartGraph::Node o = graph.target(mwpm.matching(n));
    check(mwpm.mate(n) == o, "Invalid matching");
    check(mwpm.matching(n) == graph.oppositeArc(mwpm.matching(o)),
          "Invalid matching");
  }

  int dv = 0;
  for (SmartGraph::NodeIt n(graph); n != INVALID; ++n) {
    dv += mwpm.nodeValue(n);
  }

  for (int i = 0; i < mwpm.blossomNum(); ++i) {
    check(mwpm.blossomValue(i) >= 0, "Invalid blossom value");
    check(mwpm.blossomSize(i) % 2 == 1, "Even blossom size");
    dv += mwpm.blossomValue(i) * ((mwpm.blossomSize(i) - 1) / 2);
  }

  check(pv * mwpm.dualScale == dv * 2, "Wrong duality");

  return;
}


int main() {

  for (int i = 0; i < 3; ++i) {
    SmartGraph graph;
    SmartGraph::EdgeMap<int> weight(graph);

    istringstream lgfs(lgf[i]);
    graphReader(graph, lgfs).
      edgeMap("weight", weight).run();

    bool perfect = (i!=2);

    {
      MaxWeightedPerfectMatching<SmartGraph> mwpm(graph, weight);
      bool result = mwpm.run();

      check(result == perfect, "Perfect matching found");
      if (perfect) {
        checkWeightedPerfectMatching(graph, weight, mwpm);
      }
    }

    {
      MaxWeightedPerfectMatching<SmartGraph> mwpm(graph, weight);
      mwpm.init();
      bool result = mwpm.start();

      check(result == perfect, "Perfect matching found");
      if (perfect) {
        checkWeightedPerfectMatching(graph, weight, mwpm);
      }
    }
  }

  return 0;
}
