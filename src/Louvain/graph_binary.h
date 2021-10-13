// File: graph_binary.h
// -- graph handling header file
//-----------------------------------------------------------------------------
// Community detection
// Based on the article "Fast unfolding of community hierarchies in large networks"
// Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
//
// And based on the article
// Copyright (C) 2013 R. Campigotto, P. Conde Céspedes, J.-L. Guillaume
//
// This file is part of Louvain algorithm.
//
// Louvain algorithm is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Louvain algorithm is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Louvain algorithm.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------
// Author   : E. Lefebvre, adapted by J.-L. Guillaume and R. Campigotto
// Email    : jean-loup.guillaume@lip6.fr
// Location : Paris, France
// Time	    : July 2013
//-----------------------------------------------------------------------------
// see README.txt for more details


#ifndef GRAPH_H
#define GRAPH_H

#include <assert.h>
#include <vector>
#include <algorithm>
// #include "stdafx.h"
#include "Snap.h"

#define WEIGHTED   0
#define UNWEIGHTED 1

typedef unsigned int uint;
typedef unsigned long ulint;


typedef TVec< THash<TInt, TInt> > WeightVHM;


class Graph {
 public:
  uint nb_nodes;
  ulint nb_links;

  ulint total_weight;
  uint sum_nodes_w;

  std::vector<ulint> degrees;
  std::vector<uint> links;
  std::vector<bool> nodes_auth;
  std::vector<uint> weights;

  std::vector<uint> nodes_w;

  Graph();

  // binary file format is
  // 4 bytes for the number of nodes in the graph
  // 8*(nb_nodes) bytes for the cumulative degree for each node:
  //    deg(0)=degrees[0]
  //    deg(k)=degrees[k]-degrees[k-1]
  // 4*(sum_degrees) bytes for the links
  // IF WEIGHTED, 10*(sum_degrees) bytes for the weights in a separate file
  Graph(char *filename, char *filename_w, int type);
  Graph(const WeightVHM SnapGraph);
  Graph(char *EdgeListIn);

  // return the biggest weight of links in the graph

  // assign a weight to a node (needed after the first level)
  void assign_weight(uint node, uint weight);

  void show();


  // return the number of neighbors (degree) of the node
  inline uint nb_neighbors(uint node);

  // return the number of self loops of the node
  inline uint nb_selfloops(uint node);

  // return the weighted degree of the node
  inline ulint weighted_degree(uint node);

  // return pointers to the first neighbor and first weight of the node
  inline std::pair<std::vector<uint>::iterator, std::vector<uint>::iterator > neighbors(uint node);
};


inline uint Graph::nb_neighbors(uint node) {
  assert(node>=0 && node<nb_nodes);

  if (node==0)
    return degrees[0];
  else
    return (uint)(degrees[node]-degrees[node-1]);
}

inline uint Graph::nb_selfloops(uint node) {
  assert(node>=0 && node<nb_nodes);

  std::pair<std::vector<uint>::iterator, std::vector<uint>::iterator > p = neighbors(node);
  for (uint i=0 ; i<nb_neighbors(node) ; i++) {
    if (*(p.first+i)==node) {
      if (weights.size()!=0) {
	      return (uint)*(p.second+i);
      }
      else
	return 1.0L;
    }
  }
  return 0.0L;
}

inline ulint Graph::weighted_degree(uint node) {
  assert(node>=0 && node<nb_nodes);

  if (weights.size()==0)
    return (ulint)nb_neighbors(node);
  else {
    std::pair<std::vector<uint>::iterator, std::vector<uint>::iterator > p = neighbors(node);
    ulint res = 0L;
    for (uint i=0 ; i<nb_neighbors(node) ; i++) {
      res += (ulint)*(p.second+i);
    }
    return res;
  }
}

inline std::pair<std::vector<uint>::iterator, std::vector<uint>::iterator > Graph::neighbors(uint node) {
  assert(node>=0 && node<nb_nodes);

  if (node==0)
    return make_pair(links.begin(), weights.begin());
  else if (weights.size()!=0)
    return make_pair(links.begin()+degrees[node-1], weights.begin()+degrees[node-1]);
  else
    return make_pair(links.begin()+degrees[node-1], weights.begin());
}


#endif // GRAPH_H
