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

#ifndef LEMON_FRACTIONAL_MATCHING_H
#define LEMON_FRACTIONAL_MATCHING_H

#include <vector>
#include <queue>
#include <set>
#include <limits>

#include <lemon/core.h>
#include <lemon/unionfind.h>
#include <lemon/bin_heap.h>
#include <lemon/maps.h>
#include <lemon/assert.h>
#include <lemon/elevator.h>

 ///\ingroup matching
 ///\file
 ///\brief Fractional matching algorithms in general graphs.

namespace lemon {

  /// \ingroup matching
  ///
  /// \brief Weighted fractional perfect matching in general graphs
  ///
  /// This class provides an efficient implementation of fractional
  /// matching algorithm. The implementation uses priority queues and
  /// provides \f$O(nm\log n)\f$ time complexity.
  ///
  /// The maximum weighted fractional perfect matching is a relaxation
  /// of the maximum weighted perfect matching problem where the odd
  /// set constraints are omitted.
  /// It can be formulated with the following linear program.
  /// \f[ \sum_{e \in \delta(u)}x_e = 1 \quad \forall u\in V\f]
  /// \f[x_e \ge 0\quad \forall e\in E\f]
  /// \f[\max \sum_{e\in E}x_ew_e\f]
  /// where \f$\delta(X)\f$ is the set of edges incident to a node in
  /// \f$X\f$. The result must be the union of a matching with one
  /// value edges and a set of odd length cycles with half value edges.
  ///
  /// The algorithm calculates an optimal fractional matching and a
  /// proof of the optimality. The solution of the dual problem can be
  /// used to check the result of the algorithm. The dual linear
  /// problem is the following.
  /// \f[ y_u + y_v \ge w_{uv} \quad \forall uv\in E\f]
  /// \f[\min \sum_{u \in V}y_u \f]
  ///
  /// The algorithm can be executed with the run() function.
  /// After it the matching (the primal solution) and the dual solution
  /// can be obtained using the query functions.
  ///
  /// The primal solution is multiplied by
  /// \ref MaxWeightedPerfectFractionalMatching::primalScale "2".
  /// If the value type is integer, then the dual
  /// solution is scaled by
  /// \ref MaxWeightedPerfectFractionalMatching::dualScale "4".
  ///
  /// \tparam GR The undirected graph type the algorithm runs on.
  /// \tparam WM The type edge weight map. The default type is
  /// \ref concepts::Graph::EdgeMap "GR::EdgeMap<int>".
#ifdef DOXYGEN
  template <typename GR, typename WM>
#else
  template <typename GR,
    typename WM = typename GR::template EdgeMap<int> >
#endif
    class MaxWeightedPerfectFractionalMatching {
    public:

      /// The graph type of the algorithm
      typedef GR Graph;
      /// The type of the edge weight map
      typedef WM WeightMap;
      /// The value type of the edge weights
      typedef typename WeightMap::Value Value;

      /// The type of the matching map
      typedef typename Graph::template NodeMap<typename Graph::Arc>
        MatchingMap;

      /// \brief Scaling factor for primal solution
      ///
      /// Scaling factor for primal solution.
      static const int primalScale = 2;

      /// \brief Scaling factor for dual solution
      ///
      /// Scaling factor for dual solution. It is equal to 4 or 1
      /// according to the value type.
      static const int dualScale =
        std::numeric_limits<Value>::is_integer ? 4 : 1;

    private:

      TEMPLATE_GRAPH_TYPEDEFS(Graph);

      typedef typename Graph::template NodeMap<Value> NodePotential;

      const Graph& _graph;
      const WeightMap& _weight;

      MatchingMap* _matching;
      NodePotential* _node_potential;

      int _node_num;
      bool _allow_loops;

      enum Status {
        EVEN = -1, MATCHED = 0, ODD = 1
      };

      typedef typename Graph::template NodeMap<Status> StatusMap;
      StatusMap* _status;

      typedef typename Graph::template NodeMap<Arc> PredMap;
      PredMap* _pred;

      typedef ExtendFindEnum<IntNodeMap> TreeSet;

      IntNodeMap *_tree_set_index;
      TreeSet *_tree_set;

      IntNodeMap *_delta2_index;
      BinHeap<Value, IntNodeMap> *_delta2;

      IntEdgeMap *_delta3_index;
      BinHeap<Value, IntEdgeMap> *_delta3;

      Value _delta_sum;

      void createStructures() {
        _node_num = countNodes(_graph);

        if (!_matching) {
          _matching = new MatchingMap(_graph);
        }
        if (!_node_potential) {
          _node_potential = new NodePotential(_graph);
        }
        if (!_status) {
          _status = new StatusMap(_graph);
        }
        if (!_pred) {
          _pred = new PredMap(_graph);
        }
        if (!_tree_set) {
          _tree_set_index = new IntNodeMap(_graph);
          _tree_set = new TreeSet(*_tree_set_index);
        }
        if (!_delta2) {
          _delta2_index = new IntNodeMap(_graph);
          _delta2 = new BinHeap<Value, IntNodeMap>(*_delta2_index);
        }
        if (!_delta3) {
          _delta3_index = new IntEdgeMap(_graph);
          _delta3 = new BinHeap<Value, IntEdgeMap>(*_delta3_index);
        }
      }

      void destroyStructures() {
        if (_matching) {
          delete _matching;
        }
        if (_node_potential) {
          delete _node_potential;
        }
        if (_status) {
          delete _status;
        }
        if (_pred) {
          delete _pred;
        }
        if (_tree_set) {
          delete _tree_set_index;
          delete _tree_set;
        }
        if (_delta2) {
          delete _delta2_index;
          delete _delta2;
        }
        if (_delta3) {
          delete _delta3_index;
          delete _delta3;
        }
      }

      void matchedToEven(Node node, int tree) {
        _tree_set->insert(node, tree);
        _node_potential->set(node, (*_node_potential)[node] + _delta_sum);

        if (_delta2->state(node) == _delta2->IN_HEAP) {
          std::cerr << "delta2 erasing " << node.Id() << std::endl;
          _delta2->erase(node);
        }

        for (InArcIt a(_graph, node); a != INVALID; ++a) {
          Node v = _graph.source(a);
          int weight = _weight[a];
          Value rw = (*_node_potential)[node] + (*_node_potential)[v] -
            dualScale * weight;
          if (node == v) {
            if (_allow_loops && _graph.direction(a)) {
              _delta3->push(a, rw / 2);
            }
          }
          else if ((*_status)[v] == EVEN) {
            _delta3->push(a, rw / 2);
          }
          else if ((*_status)[v] == MATCHED) {
            if (_delta2->state(v) != _delta2->IN_HEAP) {
 //             std::cerr << "setting pred " << v.Id() << " " << _graph.ArcTarget(a.Id()) << std::endl;

              _pred->set(v, a);
              std::cerr << "delta2 pushing " << v.Id() << rw;
              _delta2->push(v, rw);
            }
            else if ((*_delta2)[v] > rw) {
  //            std::cerr << "setting pred " << v.Id() << " " << _graph.ArcTarget(a.Id()) << std::endl;

              _pred->set(v, a);
              std::cerr << "delta2 decreasing " << v.Id() << " " << rw << std::endl;
              _delta2->decrease(v, rw);
            }
          }
        }
      }

      void matchedToOdd(Node node, int tree) {
        _tree_set->insert(node, tree);
        _node_potential->set(node, (*_node_potential)[node] - _delta_sum);

        if (_delta2->state(node) == _delta2->IN_HEAP) {
          std::cerr << "delta2 erasing " << node.Id() << std::endl;
          _delta2->erase(node);
        }
      }

      void evenToMatched(Node node, int tree) {
   //     std::cerr << "etm " << node.Id() << " " << tree << std::endl;
        _node_potential->set(node, (*_node_potential)[node] - _delta_sum);
        Arc min = INVALID;
        Value minrw = std::numeric_limits<Value>::max();
        for (InArcIt a(_graph, node); a != INVALID; ++a) {
          Node v = _graph.source(a);
          Value rw = (*_node_potential)[node] + (*_node_potential)[v] -
            dualScale * _weight[a];

          if (node == v) {
            if (_allow_loops && _graph.direction(a)) {
              _delta3->erase(a);
            }
          }
          else if ((*_status)[v] == EVEN) {
            int weight = _weight[a];
    //        std::cerr << "erasing arc "<< weight << std::endl;
            _delta3->erase(a);
            if (minrw > rw) {
              min = _graph.oppositeArc(a);
              minrw = rw;
            }
          }
          else if ((*_status)[v] == MATCHED) {
            if ((*_pred)[v] == a) {
              Arc mina = INVALID;
              Value minrwa = std::numeric_limits<Value>::max();
              for (OutArcIt aa(_graph, v); aa != INVALID; ++aa) {
                Node va = _graph.target(aa);
                if ((*_status)[va] != EVEN ||
                  _tree_set->find(va) == tree) continue;
                Value rwa = (*_node_potential)[v] + (*_node_potential)[va] -
                  dualScale * _weight[aa];
                if (minrwa > rwa) {
                  minrwa = rwa;
                  mina = aa;
                }
              }
              if (mina != INVALID) {
       //         std::cerr << "setting pred " << v.Id() << " " << _graph.ArcTarget(mina.Id()) << std::endl;
                _pred->set(v, mina);
                std::cerr << "delta2 increasing " << v.Id() << " " << minrwa << std::endl;
                _delta2->increase(v, minrwa);
              }
              else {
      //          std::cerr << "setting pred " << v.Id() << " " << "null" << std::endl;
                _pred->set(v, INVALID);
                std::cerr << "delta2 erasing " << v.Id() << std::endl;
                _delta2->erase(v);
              }
            }
          }
        }
        if (min != INVALID) {
          int minId = min.Id();
    //      std::cerr << "setting pred " << node.Id() << " " << _graph.ArcTarget(minId) << std::endl;
          _pred->set(node, min);
          std::cerr << "delta2 pushing " << node.Id() << " " << minrw << std::endl;
          _delta2->push(node, minrw);
        }
        else {
   //      std::cerr << "setting pred " << node.Id() << " null" << std::endl;
          _pred->set(node, INVALID);
        }
      }

      void oddToMatched(Node node) {
        _node_potential->set(node, (*_node_potential)[node] + _delta_sum);
        Arc min = INVALID;
        Value minrw = std::numeric_limits<Value>::max();
        for (InArcIt a(_graph, node); a != INVALID; ++a) {
          Node v = _graph.source(a);
          if ((*_status)[v] != EVEN) continue;
          Value rw = (*_node_potential)[node] + (*_node_potential)[v] -
            dualScale * _weight[a];

          if (minrw > rw) {
            min = _graph.oppositeArc(a);
            minrw = rw;
          }
        }
        if (min != INVALID) {
     //     std::cerr << "setting pred " << node.Id() << _graph.ArcTarget(min.Id()) << std::endl;
          _pred->set(node, min);
          std::cerr << "delta2 pushing " << node.Id() << " " << minrw << std::endl;
          _delta2->push(node, minrw);
        }
        else {
  //        std::cerr << "setting pred " << node.Id() << " null" << std::endl;

          _pred->set(node, INVALID);
        }
      }

      void alternatePath(Node even, int tree) {
        Node odd;

        _status->set(even, MATCHED);
        evenToMatched(even, tree);

        Arc prev = (*_matching)[even];
        while (prev != INVALID) {
          odd = _graph.target(prev);
          even = _graph.target((*_pred)[odd]);
          std::cerr << "setting matching " << odd.Id() << " " << even.Id() << std::endl;
          _matching->set(odd, (*_pred)[odd]);
          _status->set(odd, MATCHED);
          oddToMatched(odd);

          prev = (*_matching)[even];
          _status->set(even, MATCHED);
          std::cerr << "setting matching " << even.Id() << " " << odd.Id() << std::endl;
          _matching->set(even, _graph.oppositeArc((*_matching)[odd]));
          evenToMatched(even, tree);
        }
      }

      void destroyTree(int tree) {
        for (typename TreeSet::ItemIt n(*_tree_set, tree); n != INVALID; ++n) {
          if ((*_status)[n] == EVEN) {
            _status->set(n, MATCHED);
            evenToMatched(n, tree);
          }
          else if ((*_status)[n] == ODD) {
            _status->set(n, MATCHED);
            oddToMatched(n);
          }
        }
        _tree_set->eraseClass(tree);
      }

      void augmentOnEdge(const Edge& edge) {
        Node left = _graph.u(edge);
        int left_tree = _tree_set->find(left);

        alternatePath(left, left_tree);
        destroyTree(left_tree);
        Arc arc = _graph.direct(edge, true);
        Node target = _graph.target(arc);
        std::cerr << "setting matching " << left.Id() << " " << target.Id() << std::endl;
        _matching->set(left, arc);

        Node right = _graph.v(edge);
        int right_tree = _tree_set->find(right);

        alternatePath(right, right_tree);
        destroyTree(right_tree);
        Arc rightTarget = _graph.direct(edge, false);
        std::cerr << "setting matching " << right.Id() << " " << _graph.target(rightTarget).Id() << std::endl;

        _matching->set(right, rightTarget);
      }

      void augmentOnArc(const Arc& arc) {
        Node left = _graph.source(arc);
        _status->set(left, MATCHED);
        std::cerr << "setting matching " << left.Id() << " " << _graph.target(arc).Id() << std::endl;

        _matching->set(left, arc);
  //      std::cerr << "setting pred " << left.Id() << " " << _graph.ArcTarget(arc.Id()) << std::endl;
        _pred->set(left, arc);

        Node right = _graph.target(arc);
        int right_tree = _tree_set->find(right);

        alternatePath(right, right_tree);
        destroyTree(right_tree);
        std::cerr << "setting matching " << right.Id() << " " << _graph.source(arc).Id() << std::endl;

        _matching->set(right, _graph.oppositeArc(arc));
      }

      void extendOnArc(const Arc& arc) {
        Node base = _graph.target(arc);
        int tree = _tree_set->find(base);

        Node odd = _graph.source(arc);
        _tree_set->insert(odd, tree);
        _status->set(odd, ODD);
        matchedToOdd(odd, tree);

        int arcId = arc.Id();
    //    std::cerr << "setting pred " << odd.Id() << " " << _graph.ArcTarget(arcId) << std::endl;
        _pred->set(odd, arc);

        Node even = _graph.target((*_matching)[odd]);
        _tree_set->insert(even, tree);
        _status->set(even, EVEN);
        matchedToEven(even, tree);
      }

      void cycleOnEdge(const Edge& edge, int tree) {
   //     std::cerr << "cycleOnEdge" << " " << tree << std::endl;
        Node nca = INVALID;
        std::vector<Node> left_path, right_path;

        {
          std::set<Node> left_set, right_set;
          Node left = _graph.u(edge);
      //    std::cerr << "left path adding " << left.Id() << std::endl;
          left_path.push_back(left);
          left_set.insert(left);

          Node right = _graph.v(edge);
          right_path.push_back(right);
          right_set.insert(right);

          while (true) {

            if (left_set.find(right) != left_set.end()) {
              nca = right;
              break;
            }

            if ((*_matching)[left] == INVALID) break;

            left = _graph.target((*_matching)[left]);
      //      std::cerr << "left path adding " << left.Id() << std::endl;
            left_path.push_back(left);
            left = _graph.target((*_pred)[left]);
       //     std::cerr << "left path adding " << left.Id() << std::endl;
            left_path.push_back(left);

            left_set.insert(left);

            if (right_set.find(left) != right_set.end()) {
              nca = left;
              break;
            }

            if ((*_matching)[right] == INVALID) break;

            right = _graph.target((*_matching)[right]);
            right_path.push_back(right);
            right = _graph.target((*_pred)[right]);
            right_path.push_back(right);

            right_set.insert(right);

          }

          if (nca == INVALID) {
            if ((*_matching)[left] == INVALID) {
              nca = right;
              while (left_set.find(nca) == left_set.end()) {
                nca = _graph.target((*_matching)[nca]);
                right_path.push_back(nca);
                nca = _graph.target((*_pred)[nca]);
                right_path.push_back(nca);
              }
            }
            else {
              nca = left;
              while (right_set.find(nca) == right_set.end()) {
                nca = _graph.target((*_matching)[nca]);
         //       std::cerr << "left path adding " << nca.Id() << std::endl;
                left_path.push_back(nca);
         //       std::cerr << "left path adding " << nca.Id() << std::endl;
                nca = _graph.target((*_pred)[nca]);
                left_path.push_back(nca);
              }
            }
          }
        }

        alternatePath(nca, tree);
        Arc prev;

        prev = _graph.direct(edge, true);
        for (int i = 0; left_path[i] != nca; i += 2) {
          std::cerr << "setting matching " << left_path[i].Id() << " " << _graph.target(prev).Id() << std::endl;

          _matching->set(left_path[i], prev);
          _status->set(left_path[i], MATCHED);
          evenToMatched(left_path[i], tree);
          Node vert = left_path[i + 1];
          prev = _graph.oppositeArc((*_pred)[vert]);
          _status->set(left_path[i + 1], MATCHED);
          oddToMatched(left_path[i + 1]);
        }
        std::cerr << "setting matching " << nca.Id() << " " << _graph.target(prev).Id() << std::endl;

        _matching->set(nca, prev);

        for (int i = 0; right_path[i] != nca; i += 2) {
          _status->set(right_path[i], MATCHED);
          evenToMatched(right_path[i], tree);
          std::cerr << "setting matching " << right_path[i+1].Id() << " " << (*_pred)[right_path[i + 1]].Id() << std::endl;

          _matching->set(right_path[i + 1], (*_pred)[right_path[i + 1]]);
          _status->set(right_path[i + 1], MATCHED);
          oddToMatched(right_path[i + 1]);
        }

        destroyTree(tree);
      }

      void extractCycle(const Arc &arc) {
        Node left = _graph.source(arc);
        Node odd = _graph.target((*_matching)[left]);
        Arc prev;
        while (odd != left) {
          Node even = _graph.target((*_matching)[odd]);
          prev = (*_matching)[odd];
          odd = _graph.target((*_matching)[even]);
          std::cerr << "setting matching " << even.Id() << " " << _graph.source(prev).Id() << std::endl;

          _matching->set(even, _graph.oppositeArc(prev));
        }
        std::cerr << "setting matching " << left.Id() << " " << _graph.target(arc).Id() << std::endl;

        _matching->set(left, arc);

        Node right = _graph.target(arc);
        int right_tree = _tree_set->find(right);
        alternatePath(right, right_tree);
        destroyTree(right_tree);
        std::cerr << "setting matching " << right.Id() << " " << _graph.source(arc).Id() << std::endl;

        _matching->set(right, _graph.oppositeArc(arc));
      }

    public:

      /// \brief Constructor
      ///
      /// Constructor.
      MaxWeightedPerfectFractionalMatching(const Graph& graph,
        const WeightMap& weight,
        bool allow_loops = true)
        : _graph(graph), _weight(weight), _matching(0),
        _node_potential(0), _node_num(0), _allow_loops(allow_loops),
        _status(0), _pred(0),
        _tree_set_index(0), _tree_set(0),

        _delta2_index(0), _delta2(0),
        _delta3_index(0), _delta3(0),

        _delta_sum() {}

      ~MaxWeightedPerfectFractionalMatching() {
        destroyStructures();
      }

      /// \name Execution Control
      /// The simplest way to execute the algorithm is to use the
      /// \ref run() member function.

      ///@{

      /// \brief Initialize the algorithm
      ///
      /// This function initializes the algorithm.
      void init() {
     //   std::cerr << "MaxWeightedPerfectFractionalMatching init" << std::endl;
        createStructures();

        for (NodeIt n(_graph); n != INVALID; ++n) {
          (*_delta2_index)[n] = _delta2->PRE_HEAP;
        }
        for (EdgeIt e(_graph); e != INVALID; ++e) {
          (*_delta3_index)[e] = _delta3->PRE_HEAP;
        }

        _delta2->clear();
        _delta3->clear();
        _tree_set->clear();

        for (NodeIt n(_graph); n != INVALID; ++n) {
          Value max = -std::numeric_limits<Value>::max();
    
          // This calculates the initial node potentials
          for (OutArcIt e(_graph, n); e != INVALID; ++e) {
            if (_graph.target(e) == n && !_allow_loops) continue;
            if ((dualScale * _weight[e]) / 2 > max) {
              max = (dualScale * _weight[e]) / 2;
            }
          }
          _node_potential->set(n, max);

          _tree_set->insert(n);
          
          _matching->set(n, INVALID);
          Arc foo = (*_matching)[n];
          _status->set(n, EVEN);
        }

        for (EdgeIt e(_graph); e != INVALID; ++e) {
          Node left = _graph.u(e);
          Node right = _graph.v(e);
          if (left == right && !_allow_loops) continue;
          int edgeDelta = ((*_node_potential)[left] +
            (*_node_potential)[right] -
            dualScale * _weight[e]) / 2;
          int weight = _weight[e];
   //       std::cerr << "edgeDelta for edge with weight " << weight << " is " << edgeDelta << std::endl;
          _delta3->push(e, edgeDelta);
        }
      }

      /// \brief Start the algorithm
      ///
      /// This function starts the algorithm.
      ///
      /// \pre \ref init() must be called before using this function.
      bool start() {
        enum OpType {
          D2, D3
        };

        int unmatched = _node_num;
        while (unmatched > 0) {
          Value d2 = !_delta2->empty() ?
            _delta2->prio() : std::numeric_limits<Value>::max();

          Value d3 = !_delta3->empty() ?
            _delta3->prio() : std::numeric_limits<Value>::max();

          _delta_sum = d3; OpType ot = D3;
    //      std::cerr << "d2 " << d2 << " d3 " << d3 << "\n";
          if (d2 < _delta_sum) { _delta_sum = d2; ot = D2; }

          if (_delta_sum == std::numeric_limits<Value>::max()) {
            return false;
          }

          switch (ot) {
          case D2:
          {
            std::cerr << "D2" << std::endl;
            Node n = _delta2->top();
            Arc a = (*_pred)[n];
            if ((*_matching)[n] == INVALID) {
              augmentOnArc(a);
              --unmatched;
            }
            else {
              Arc nArc = (*_matching)[n];
              Node v = _graph.target(nArc);
              Arc matchV = (*_matching)[v];
              if (nArc !=
                _graph.oppositeArc(matchV)) {
                extractCycle(a);
                --unmatched;
              }
              else {
                extendOnArc(a);
              }
            }
          } break;
          case D3:
          {
     //       std::cerr << "D3" << std::endl;
            Edge e = _delta3->top();

            Node left = _graph.u(e);
            Node right = _graph.v(e);

            int left_tree = _tree_set->find(left);
            int right_tree = _tree_set->find(right);

      //      std::cerr << left.Id() << " " << left_tree << " " << right.Id() << " " << right_tree << "\n";

            if (left_tree == right_tree) {
              cycleOnEdge(e, left_tree);
              --unmatched;
            }
            else {
              augmentOnEdge(e);
              unmatched -= 2;
            }
          } break;
          }
        }
        return true;
      }

      /// \brief Run the algorithm.
      ///
      /// This method runs the \c %MaxWeightedPerfectFractionalMatching
      /// algorithm.
      ///
      /// \note mwfm.run() is just a shortcut of the following code.
      /// \code
      ///   mwpfm.init();
      ///   mwpfm.start();
      /// \endcode
      bool run() {
        init();
        return start();
      }

      /// @}

      /// \name Primal Solution
      /// Functions to get the primal solution, i.e. the maximum weighted
      /// matching.\n
      /// Either \ref run() or \ref start() function should be called before
      /// using them.

      /// @{

      /// \brief Return the weight of the matching.
      ///
      /// This function returns the weight of the found matching. This
      /// value is scaled by \ref primalScale "primal scale".
      ///
      /// \pre Either run() or start() must be called before using this function.
      Value matchingWeight() const {
        Value sum = 0;
        for (NodeIt n(_graph); n != INVALID; ++n) {
          if ((*_matching)[n] != INVALID) {
            sum += _weight[(*_matching)[n]];
          }
        }
        return sum * primalScale / 2;
      }

      /// \brief Return the number of covered nodes in the matching.
      ///
      /// This function returns the number of covered nodes in the matching.
      ///
      /// \pre Either run() or start() must be called before using this function.
      int matchingSize() const {
        int num = 0;
        for (NodeIt n(_graph); n != INVALID; ++n) {
          if ((*_matching)[n] != INVALID) {
            ++num;
          }
        }
        return num;
      }

      /// \brief Return \c true if the given edge is in the matching.
      ///
      /// This function returns \c true if the given edge is in the
      /// found matching. The result is scaled by \ref primalScale
      /// "primal scale".
      ///
      /// \pre Either run() or start() must be called before using this function.
      int matching(const Edge& edge) const {
        return (edge == (*_matching)[_graph.u(edge)] ? 1 : 0)
          + (edge == (*_matching)[_graph.v(edge)] ? 1 : 0);
      }

      /// \brief Return the fractional matching arc (or edge) incident
      /// to the given node.
      ///
      /// This function returns one of the fractional matching arc (or
      /// edge) incident to the given node in the found matching or \c
      /// INVALID if the node is not covered by the matching or if the
      /// node is on an odd length cycle then it is the successor edge
      /// on the cycle.
      ///
      /// \pre Either run() or start() must be called before using this function.
      Arc matching(const Node& node) const {
        return (*_matching)[node];
      }

      /// \brief Return a const reference to the matching map.
      ///
      /// This function returns a const reference to a node map that stores
      /// the matching arc (or edge) incident to each node.
      const MatchingMap& matchingMap() const {
        return *_matching;
      }

      /// @}

      /// \name Dual Solution
      /// Functions to get the dual solution.\n
      /// Either \ref run() or \ref start() function should be called before
      /// using them.

      /// @{

      /// \brief Return the value of the dual solution.
      ///
      /// This function returns the value of the dual solution.
      /// It should be equal to the primal value scaled by \ref dualScale
      /// "dual scale".
      ///
      /// \pre Either run() or start() must be called before using this function.
      Value dualValue() const {
        Value sum = 0;
        for (NodeIt n(_graph); n != INVALID; ++n) {
          sum += nodeValue(n);
        }
        return sum;
      }

      /// \brief Return the dual value (potential) of the given node.
      ///
      /// This function returns the dual value (potential) of the given node.
      ///
      /// \pre Either run() or start() must be called before using this function.
      Value nodeValue(const Node& n) const {
        return (*_node_potential)[n];
      }

      /// @}

  };

} //END OF NAMESPACE LEMON

#endif //LEMON_FRACTIONAL_MATCHING_H
