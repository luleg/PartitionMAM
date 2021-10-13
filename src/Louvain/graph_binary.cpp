// File: graph_binary.cpp
// -- graph handling source
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


#include <stdio.h>
#include "graph_binary.h"


Graph::Graph() {
  nb_nodes = 0;
  nb_links = 0UL;

  total_weight = 0UL;
  sum_nodes_w = 0;
}

Graph::Graph(char *filename, char *filename_w, int type) {

  FILE* finputF = fopen(filename,"rb");
  if (finputF==NULL) {
    fprintf(stderr,"The file %s does not exist.\n", filename);
    exit(EXIT_FAILURE);
  }
  FILE* finputF_w = fopen(filename_w,"rb");
  if (finputF==NULL) {
    fprintf(stderr,"The file %s does not exist.\n", filename_w);
    exit(EXIT_FAILURE);
  }

  // Read number of nodes on 4 bytes
  size_t toComp = fread((char *)&nb_nodes,sizeof(nb_nodes), 1, finputF);

  if ( toComp != 1) {
    fprintf(stderr,"The file %s does not contain a valid graph.\n", filename);
    exit(EXIT_FAILURE);
  }

  // Read cumulative degree sequence: 8 bytes for each node
  // cum_degree[0]=degree(0); cum_degree[1]=degree(0)+degree(1), etc.
  degrees.resize(nb_nodes);
  ulint degOld=0;
  total_weight = 0UL;
  for (uint i=0;i<nb_nodes;i++){

    uint deg;
    fread((char *)&deg,sizeof(uint),1, finputF);
    degOld+=deg;
    degrees[i]=degOld;
    for (uint k=0;k<deg;k++){
      uint neighb,weight;
      fread((char *)&neighb,sizeof(neighb),1, finputF);
      fread((char *)&weight,sizeof(weight),1, finputF_w);
      links.push_back(neighb);
      weights.push_back(weight);
      total_weight+=weight;
    }
  }
  fclose(finputF);
  fclose(finputF_w);
  /***************************************************************************/

  /***************************************************************************/

  nb_links=degOld;
  nodes_w.assign(nb_nodes, 1);
  sum_nodes_w = nb_nodes;
}


Graph::Graph(const WeightVHM SnapGraph){

  nb_nodes = SnapGraph.Len();
  degrees.resize(nb_nodes);
  total_weight = 0UL;
  ulint degOld=0;
  for (uint node=0;node<nb_nodes;node++){
    const THash<TInt, TInt>& edge_list = SnapGraph[node];
    degOld+=edge_list.Len();
    degrees[node] = degOld;
    for (THash<TInt, TInt>::TIter it = edge_list.BegI(); it < edge_list.EndI();it++) {
      uint nei = (uint)it->Key;
      uint weight = (uint)it->Dat;
      links.push_back(nei);
      weights.push_back(weight);
      total_weight+=weight;
    }
  }
  nb_links=degOld;
  nodes_w.assign(nb_nodes, 1);
  sum_nodes_w = nb_nodes;
}


Graph::Graph(char *EdgeListIn){
  WeightVHM TmpGraph;
  nb_nodes = 0;

  FILE * finput = fopen(EdgeListIn,"r");
  if (finput == NULL){
    fprintf(stderr,"The file %s does not exist\n", EdgeListIn);
    exit(EXIT_FAILURE);
  }
  while (!feof(finput)) {
    uint src,tgt,weight;
    int nbRead = fscanf(finput,"%d %d %d\n",&src,&tgt,&weight);
    if (nbRead == 3){
      uint mx_crt = std::max(src,tgt)+1;
      if (mx_crt > nb_nodes){
        nb_nodes = mx_crt;
        TmpGraph.Reserve(nb_nodes,nb_nodes);
      }
      THash<TInt, TInt>& NbrNode = TmpGraph[src];
      NbrNode(tgt) += weight;
    }
  }
  total_weight = 0UL;
  degrees.resize(nb_nodes);
  ulint degOld=0;
  for (uint node=0;node<nb_nodes;node++){
    const THash<TInt, TInt>& edge_list = TmpGraph[node];
    degOld+=edge_list.Len();
    degrees[node] = degOld;
    for (THash<TInt, TInt>::TIter it = edge_list.BegI(); it < edge_list.EndI();it++) {
      uint nei = (uint)it->Key;
      uint weight = (uint)it->Dat;
      links.push_back(nei);
      weights.push_back(weight);
      total_weight+=weight;
    }
  }
  nb_links=degOld;
  nodes_w.assign(nb_nodes, 1);
  sum_nodes_w = nb_nodes;
  TmpGraph.Clr();
}

void Graph::assign_weight(uint node, uint weight) {
  sum_nodes_w -= nodes_w[node];

  nodes_w[node] = weight;

  sum_nodes_w += weight;
}

void Graph::show(){
  int deg =0;
  for (int node = 0; node < nb_nodes ; node++){
    fprintf(stderr,"Node %d of degree %d : \n",node,degrees[node]-deg);
    for (int nei = deg; nei < degrees[node];nei ++ ){
      fprintf(stderr,"\t (%d,%d)\n",links[nei],weights[nei]);
    }
    deg = degrees[node];
  }
}
