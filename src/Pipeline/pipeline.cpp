// #include "stdafx.h"
#include "BensonGraphQuad.h"
#include "BensonGraphTriangle.h"
#include "BensonGraphEdge.h"

#include "louvain.h"
#include "modularity.h"

#include <cstring>
#include <stdio.h>

TStr InGraph = ""; // Input graph
TStr InPart  = ""; // Input Partition
TStr InMAM   = ""; // Input MAM
TStr OutMAM  = ""; // Output MAM
TStr OutPart = ""; // Output final Partition
TStr OutPTmp = ""; // Output temporary partition (without postproc)

TStr Motif = "E2"; // Motif (default none)
TInt nbProc = -1; // number of threads

long double alpha = 1.0L;
long double beta  = 0.001L;
unsigned short id_qual = 0;
int kmin = 1;
int display_level = -2;

TBool help      = false;
TBool verbose   = false;
TBool buildMAM  = false;
TBool partition = false;
TBool postproc  = false;
TBool indir     = true;

/* For the "Quality measure" to be used in Louvain */
Quality *qual = NULL;
long double precision = 0.0000001L;
bool initial_part = false;

/* To get insight in the code efficiency */
TSecTm TStart;
TExeTm ExeTm;

int readFromFile(char* file){
  FILE * Finput = fopen(file,"r");
  if (Finput == NULL){
    printf("ERROR :: FILE %S NOT FOUND\n",file);
    help = true;
    return 1;
  }
  long lSize;
  char * buffer;
  size_t result;
  fseek (Finput , 0 , SEEK_END);
  lSize = ftell (Finput);
  rewind (Finput);
  buffer = (char*) malloc (sizeof(char)*lSize);
  result = fread (buffer,1,lSize,Finput);
  if (result != lSize) {
    printf ("ERROR :: PROBLEM WHILE READING FILE %S\n",file);
    help = true;
    return 1;
  }
  fclose(Finput);

  char * oneArg;
  oneArg = strtok(buffer,"- ?");
  while (oneArg != NULL)
  {
    if (strcmp(oneArg, "ma") == 0){ // Computing the MAM
      buildMAM = true;
      oneArg = strtok (NULL, "- ?");
    }
    else if (strcmp(oneArg, "pa") == 0){ // Partitioning a graph
      partition = true;
      oneArg = strtok (NULL, "- ?");
    }
    else if (strcmp(oneArg, "po") == 0){ // Postprocessing a partitioning
      postproc = true;
      oneArg = strtok (NULL, "- ?");
    }
    else if (strcmp(oneArg, "mapa") == 0){ // Computing and partitioning the MAM
      buildMAM = true;
      partition = true;
      oneArg = strtok (NULL, "- ?");
    }
    else if (strcmp(oneArg, "papo") == 0){ // Partitioning a MAM and postprocessing the partition
      partition = true;
      postproc = true;
      oneArg = strtok (NULL, "- ?");
    }

    // INPUT FILES
    else if (strcmp(oneArg, "igraph") == 0){ // File with the edgelist of the directed graph
      oneArg = strtok (NULL, "- ?");
      InGraph = TStr(oneArg);
      oneArg = strtok (NULL, "- ?");
    }
    else if (strcmp(oneArg, "isym") == 0){ // File with the edgelist of the symmetrised graph
      oneArg = strtok (NULL, "- ?");
      InGraph = TStr(oneArg);
      indir = false;
      oneArg = strtok (NULL, "- ?");
    }
    else if (strcmp(oneArg, "imam") == 0){ // File with the edgelist of the MAM
      oneArg = strtok (NULL, "- ?");
      InMAM = TStr(oneArg);
      oneArg = strtok (NULL, "- ?");
    }
    else if (strcmp(oneArg, "ipart") == 0){ // File with the partitioning obtained on MAM
      oneArg = strtok (NULL, "- ?");
      InPart = TStr(oneArg);
      oneArg = strtok (NULL, "- ?");
    }
    // OUTPUT FILES
    else if (strcmp(oneArg, "omam") == 0){ // File in which the MAM is saved
      oneArg = strtok (NULL, "- ?");
      OutMAM = TStr(oneArg);
      oneArg = strtok (NULL, "- ?");
    }
    else if (strcmp(oneArg, "opart") == 0){ // File in which the (partial) partitioning is saved
      oneArg = strtok (NULL, "- ?");
      OutPart = TStr(oneArg);
      oneArg = strtok (NULL, "- ?");
    }
    else if (strcmp(oneArg, "oppart") == 0){ // File in which the partial partitioning is saved
      oneArg = strtok (NULL, "- ?");
      OutPTmp = TStr(oneArg);
      oneArg = strtok (NULL, "- ?");
    }
    // Other arguments
    /*For building MAM*/
    else if (strcmp(oneArg, "m") == 0){ // Motif used to build the MAM
      oneArg = strtok (NULL, "- ?");
      Motif = TStr(oneArg);
      oneArg = strtok (NULL, "- ?");
    }
    else if (strcmp(oneArg, "nth") == 0){ // Number of threads
      oneArg = strtok (NULL, "- ?");
      nbProc = atoi(oneArg);
      oneArg = strtok (NULL, "- ?");
    }
    /*For Partitioning*/
    else if (strcmp(oneArg, "c") == 0){ // Resolution parameter
      oneArg = strtok (NULL, "- ?");
      alpha = atof(oneArg);
      oneArg = strtok (NULL, "- ?");
    }
    else if (strcmp(oneArg, "l") == 0){ // Level (fine or coarse)
      oneArg = strtok (NULL, "- ?");
      display_level = atoi(oneArg);
      oneArg = strtok (NULL, "- ?");
    }
    /*For Postproc*/
    else if (strcmp(oneArg, "cc") == 0){ // Resolution Param for the postprocessing
      oneArg = strtok (NULL, "- ?");
      beta = atof(oneArg);
      oneArg = strtok (NULL, "- ?");
    }
    else if (strcmp(oneArg, "k") == 0){ // minimum size to forbid merging
      oneArg = strtok (NULL, "- ?");
      kmin = atoi(oneArg);
      oneArg = strtok (NULL, "- ?");
    }
    else if (strcmp(oneArg, "h") == 0){ // Postprocessing a partitioning
      help = true;
      oneArg = strtok (NULL, "- ?");
    }
    else if (strcmp(oneArg, "v") == 0){ // Postprocessing a partitioning
      verbose = true;
      oneArg = strtok (NULL, "- ?");
    }
    oneArg = strtok (NULL, "- ?");
  }
  if (not buildMAM && not partition && not postproc){
    buildMAM  = true;
    partition = true;
    postproc  = true;
  }

  free (buffer);
  return 0;
}

int ReadArgs(int argc, char **argv){
  int num_arg = 1;
  while( num_arg < argc){
    if (strcmp(argv[num_arg], "-f") == 0){ // output files to store the MAM (binary format)
      num_arg++;
      readFromFile(argv[num_arg]);
      return 0;
    }
    else if (strcmp(argv[num_arg], "-ma") == 0){ // Computing the MAM
      buildMAM = true;
    }
    else if (strcmp(argv[num_arg], "-pa") == 0){ // Partitioning a graph
      partition = true;
    }
    else if (strcmp(argv[num_arg], "-po") == 0){ // Postprocessing a partitioning
      postproc = true;
    }
    else if (strcmp(argv[num_arg], "-mapa") == 0){ // Computing and partitioning the MAM
      buildMAM = true;
      partition = true;
    }
    else if (strcmp(argv[num_arg], "-papo") == 0){ // Partitioning a MAM and postprocessing the partition
      partition = true;
      postproc = true;
    }
    // INPUT FILES
    else if (strcmp(argv[num_arg], "-igraph") == 0){ // File with the edgelist of the directed graph
      num_arg++;
      InGraph = argv[num_arg];
    }
    else if (strcmp(argv[num_arg], "-isym") == 0){ // File with the edgelist of the symmetrised graph
      num_arg++;
      InGraph = argv[num_arg];
      indir = false;
    }
    else if (strcmp(argv[num_arg], "-imam") == 0){ // File with the edgelist of the MAM
      num_arg++;
      InMAM = argv[num_arg];
    }
    else if (strcmp(argv[num_arg], "-ipart") == 0){ // File with the partitioning obtained on MAM
      num_arg++;
      InPart = argv[num_arg];
    }
    // OUTPUT FILES
    else if (strcmp(argv[num_arg], "-omam") == 0){ // File in which the MAM is saved
      num_arg++;
      OutMAM = argv[num_arg];
    }
    else if (strcmp(argv[num_arg], "-opart") == 0){ // File in which the (partial) partitioning is saved
      num_arg++;
      OutPart = argv[num_arg];
    }
    else if (strcmp(argv[num_arg], "-oppart") == 0){ // File in which the partial partitioning is saved
      num_arg++;
      OutPTmp = argv[num_arg];
    }
    // Other arguments
    /*For building MAM*/
    else if (strcmp(argv[num_arg], "-m") == 0){ // Motif used to build the MAM
      num_arg++;
      Motif = argv[num_arg];
    }
    else if (strcmp(argv[num_arg], "-nth") == 0){ // Number of threads
      num_arg++;
      nbProc = atoi(argv[num_arg]);
    }
    /*For Partitioning*/
    else if (strcmp(argv[num_arg], "-c") == 0){ // Resolution parameter
      num_arg++;
      alpha = atof(argv[num_arg]);
    }
    else if (strcmp(argv[num_arg], "-l") == 0){ // Level (fine or coarse)
      num_arg++;
      display_level = atoi(argv[num_arg]);
    }
    /*For Postproc*/
    else if (strcmp(argv[num_arg], "-cc") == 0){ // Resolution Param for the postprocessing
      num_arg++;
      beta = atof(argv[num_arg]);
    }
    else if (strcmp(argv[num_arg], "-k") == 0){ // minimum size to forbid merging
      num_arg++;
      kmin = atoi(argv[num_arg]);
    }
    else if (strcmp(argv[num_arg], "-h") == 0){ // Postprocessing a partitioning
      help = true;
    }
    else if (strcmp(argv[num_arg], "-v") == 0){ // Postprocessing a partitioning
      verbose = true;
    }
    num_arg++;
  }
  if (not buildMAM && not partition && not postproc){
    buildMAM  = true;
    partition = true;
    postproc  = true;
  }
  return 0;
}

void init_quality(Graph *One_graph, unsigned short nbc,long double alpha) {

  if (nbc > 0){
    delete qual;
  }

  switch (id_qual) {
  case 0:
    if (alpha <= 0.){
      alpha = 1.0L;
    }
    qual = new Modularity(*One_graph,alpha);
    break;
  default:
    if (alpha <= 0.){
      alpha = 1.0L;
    }
    qual = new Modularity(*One_graph,alpha);
    break;
  }
}

void writeFinalPart(TStr FilePart){
  std::vector<std::vector<int> >levels;

  FILE * finput = fopen(FilePart.CStr(),"r");
  if (finput == NULL){
    fprintf(stderr,"The file %s does not exist\n", FilePart.CStr());
    exit(EXIT_FAILURE);
  }

  int l = -1;
  while (!feof(finput)) {
    int node, nodecomm;
    int nbRead = fscanf(finput,"%d %d\n",&node,&nodecomm);

    if (nbRead==2) {
      if (node==0) {
	       l++;
	       levels.resize(l+1);
       }
      levels[l].push_back(nodecomm);
    }
  }
  fclose(finput);
  display_level = l+1;
  std::vector<int> n2c(levels[0].size());

  for (unsigned int i=0 ; i<levels[0].size() ; i++){
    n2c[i] = i;
  }
  for (l=0 ; l<display_level ; l++){
    for (unsigned int node=0 ; node<levels[0].size() ; node++){
      n2c[node] = levels[l][n2c[node]];
    }
  }

  // FILE *FOut  = fopen("Dmy","w");
  FILE *FOut  = fopen(FilePart.CStr(),"w");

  for (unsigned int node=0 ; node<levels[0].size() ; node++){
    fprintf(FOut,"%d %d\n",node,n2c[node]);
  }
  fclose(FOut);

}

void BuildMAM(TStr InGraph,TStr Motif, Graph& MAMGraph, TStr OutMAM){
  char type_mot = Motif.GetLc()[0]; //  What kind of motif is asked to build the MAM ?

  ExeTm.Tick();
  PNGraph SGraph = TSnap::LoadEdgeList<PNGraph>(InGraph);

  TStart = TSecTm::GetCurTm()-TStart;
  if (verbose) {
    fprintf(stderr,":::CPU Time to load the graph\t:%s\n",ExeTm.GetTmStr());
  }

  BensonGraph<TNGraph> *Benson(0);
  // Decide which operation should be performed, provided the specified graphlet/motif
  ////////// QUADRANGLE
  if (type_mot == 'q'){// The Motif is a quadrangle::
    MotifTypeQuad QuadMot(Motif);

    /// MULTIPROC
    if (nbProc>0){ // If Several processors
      TStart = TSecTm::GetCurTm();
      ExeTm.Tick();
      Benson = new BensonGraphQuadMP(SGraph,QuadMot,nbProc, verbose);
      TStart = TSecTm::GetCurTm()-TStart;
      if (verbose){
        fprintf(stderr,":::CPU Time to build the MAM\t:%s\n(Ellapsed Time\t\t:%s)\n",
                ExeTm.GetTmStr(),TStart.GetTmStr().CStr());
      }
    }
    /// SEQUENTIAL
    else{ // Otherwise
      ExeTm.Tick();
      Benson = new BensonGraphQuad(SGraph,QuadMot,true, verbose);
      TStart = TSecTm::GetCurTm()-TStart;
      if (verbose){
        fprintf(stderr,":::CPU Time to build the MAM\t:%s\n",ExeTm.GetTmStr());
      }
    }
  }

  ////////// TRIANGLE
  else if (type_mot == 't'){// The Motif is a triangle
    MotifTypeTriangle TriMot(Motif);

    /// MULTIPROC
    if (nbProc>0){
      TStart = TSecTm::GetCurTm();
      ExeTm.Tick();
      Benson = new BensonGraphTriangleMP(SGraph,TriMot,nbProc, verbose);
      TStart = TSecTm::GetCurTm()-TStart;
      if (verbose){
        fprintf(stderr,":::CPU Time to build the MAM\t:%s\n(Ellapsed Time\t\t:%s)\n",
                ExeTm.GetTmStr(),TStart.GetTmStr().CStr());
      }
    }
    /// SEQUENTIAL
    else{
      ExeTm.Tick();
      Benson = new BensonGraphTriangle(SGraph,TriMot,true, verbose);
      TStart = TSecTm::GetCurTm()-TStart;
      if (verbose){
        fprintf(stderr,":::CPU Time to build the MAM\t:%s\n",ExeTm.GetTmStr());
      }
    }
  }

  ////////// EDGE (Directed or Symmetrised)
  else if (type_mot == 'e'){
    MotifTypeEdge MEdge(Motif);

    /// MULTIPROC
    if (nbProc>0){ // If Several processors
      TStart = TSecTm::GetCurTm();
      ExeTm.Tick();
      Benson = new BensonGraphEdgeMP(SGraph,MEdge,nbProc, verbose);
      TStart = TSecTm::GetCurTm()-TStart;
      if (verbose){
        fprintf(stderr,":::CPU Time to build the MAM\t:%s\n(Ellapsed Time\t\t:%s)\n",
                ExeTm.GetTmStr(),TStart.GetTmStr().CStr());
      }
    }
    /// SEQUENTIAL
    else{ // Otherwise
      ExeTm.Tick();
      Benson = new BensonGraphEdge(SGraph,MEdge,true, verbose);
      TStart = TSecTm::GetCurTm()-TStart;
      if (verbose){
        fprintf(stderr,":::CPU Time to build the MAM\t:%s\n",ExeTm.GetTmStr());
      }
    }
  }
  ////////// UNKNOWN GRAPHLET
  else{
    TExcept::Throw("Unhandled kind of Motif.");
  }
  if (not OutMAM.Empty()){ // If the MAM is requiored to be saved in a text file
    ExeTm.Tick();
    Benson->SaveAsTxt(OutMAM);
    TStart = TSecTm::GetCurTm()-TStart;
    if (verbose){
      fprintf(stderr,":::CPU Time to write the MAM\t:%s\n",ExeTm.GetTmStr());
    }
    fprintf(stderr,"MAM saved in %s\n",OutMAM.CStr());
  }
  ExeTm.Tick();
  int singles;
  MAMGraph = Graph(Benson->GetWeight(singles));

  if (verbose){
    fprintf(stderr,":::CPU Time to convert the MAM into Louvain format\t:%s\n",ExeTm.GetTmStr());
    fprintf(stderr,":::Number of single nodes in the MAM:\t\t%d over %d\n",singles,SGraph->GetNodes());
  }

  delete Benson;
}

void Partition(Graph& AGraph, TStr OutPart){
  // Clear the file if not empty
  FILE *FOut  = fopen(OutPart.CStr(),"w");
  fclose(FOut);

  unsigned short nb_calls = 0;
  init_quality(&AGraph, nb_calls,alpha);
  nb_calls++;
  Louvain Proc(-1, precision, qual);
  bool improvement = true;
  //
  long double quality = (Proc.qual)->quality();
  long double new_qual;

  int level = 0;
  ExeTm.Tick();
  do {
    if (verbose) {
      fprintf(stderr,"--------\nlevel %d : \n\tstart computation :: %s\n", level,TSecTm::GetCurTm().GetTmStr().CStr());
      fprintf(stderr,"\tnetwork size:\n\t\t%d nodes, %d links, %d weight\n",
              (Proc.qual)->g.nb_nodes, (Proc.qual)->g.nb_links, (Proc.qual)->g.total_weight );
    }

    improvement = Proc.one_level(false,verbose);
    new_qual = (Proc.qual)->quality();
    //

    if ((display_level < 0) or (level < display_level)){
      Proc.display_partition(OutPart);
    }
    else if( level == display_level){
      improvement = false;
    }
    level++;
    //
    AGraph = Proc.partition2graph_binary(false);
    /*************************************************************************
    g.show();
    /*************************************************************************/
    init_quality(&AGraph, nb_calls,alpha);
    nb_calls++;
    Proc = Louvain(-1, precision, qual);
    if (verbose){
      fprintf(stderr,"==> Quality increased from %.4LF to %.4LF\n",quality, new_qual);
    }
    quality = new_qual;

  } while(improvement);
  if (verbose){
    fprintf(stderr,"--------\nEnd. Number of communities: %d\n",(Proc.qual)->g.nb_nodes);
    fprintf(stderr,":::CPU Time to partition the MAM\t:%s\n",ExeTm.GetTmStr());
  }

  delete qual;

  writeFinalPart(OutPart);
  if (verbose){
    fprintf(stderr,"Partial partition saved in file %s\n",OutPart.CStr());
  }
}

void PostProc(Graph& AGraph, TStr InPart){
  unsigned short nb_calls = 0;
  init_quality(&AGraph, nb_calls,beta);
  nb_calls++;

  Louvain Proc(-1, precision, qual);
  Proc.init_partition(InPart);
  if (verbose){
    fprintf(stderr,"Initial network \t%d nodes and %d edges\n",
            (Proc.qual)->g.nb_nodes, (Proc.qual)->g.nb_links);
  }

  // Clear the file if not empty
  FILE *FOut  = fopen(OutPart.CStr(),"w");
  fclose(FOut);

  AGraph = Proc.partition2graph_binary(true,true,kmin);
  Proc.display_partition(OutPart.CStr());
  init_quality(&AGraph, nb_calls,beta);
  Proc = Louvain(-1, precision, qual);
  if (verbose){
    fprintf(stderr,"Hyper Graph \t%d nodes and %d edges\n",
            (Proc.qual)->g.nb_nodes, (Proc.qual)->g.nb_links);
  }
  bool improvement = true;
  //
  long double quality = (Proc.qual)->quality();
  long double new_qual;

  int level = 0;
  do {
    if (verbose) {
      fprintf(stderr,"----------------------------------------------------------\nlevel %d : \n\tstart computation :: %s\n", level,TSecTm::GetCurTm().GetTmStr().CStr());
      fprintf(stderr,"\tnetwork size:\n\t\t%d nodes, %d links, %d weight\n", (Proc.qual)->g.nb_nodes, (Proc.qual)->g.nb_links, (Proc.qual)->g.total_weight );
    }

    improvement = Proc.one_level(true,verbose);
    new_qual = (Proc.qual)->quality();
    //
    Proc.display_partition(OutPart);
    level++;
    //
    AGraph = Proc.partition2graph_binary(true);

    init_quality(&AGraph, nb_calls,beta);
    nb_calls++;
    Proc = Louvain(-1, precision, qual);
    if (verbose){
      fprintf(stderr,"==> Quality increased from %.4LF to %.4LF\n",quality, new_qual);
    }
    quality = new_qual;

    if (level==1){
      improvement=true;
    }
  } while(improvement);
  if (verbose){
    fprintf(stderr,"--------\nEnd. Number of communities: %d\n",(Proc.qual)->g.nb_nodes);
    fprintf(stderr,":::CPU Time to partition the MAM\t:%s\n",ExeTm.GetTmStr());
  }

  delete qual;

  writeFinalPart(OutPart);

  if (verbose){
    fprintf(stderr,"Final partition saved in file %s\n",OutPart.CStr());
  }
}

void usage(){
  printf("PURPOSE :: \n\tPartitioning a directed network by applying Louvain on a Motif Adjacency Matrix (MAM) of the network.\n");
  printf("POSSIBLE USAGES :: SEE\n");
  printf("\t./PartMAM -h\t\t\t:Build the MAM, partition it, postprocess single nodes.\n");
  printf("\t./PartMAM -ma -h\t\t:Build the MAM.\n");
  printf("\t./PartMAM -pa -h\t\t:Partition a MAM.\n");
  printf("\t./PartMAM -po -h\t\t:Postprocess single nodes.\n");
  printf("\t./PartMAM -mapa -h\t\t:Build a MAM and partition it.\n");
  printf("\t./PartMAM -papo -h\t\t:Partition a MAM and postprocess single nodes.\n");
  printf("==========================================================================================\n");
  if (buildMAM && partition && postproc){
    printf("USAGE FOR THE WHOLE PIPELINE ::\n");
    printf("\t./PartMAM -igraph InEdgeList -opart OutPart [-omam OutMAM -oppart OutPPart -m Motif -nth NbTh -c alpha -l level -cc beta -k kmin]\n");
    printf("-igraph InEdgeList\t:File containing the network as an edgelist with integer nodes.\n");
    printf("\tNodes should be labelled with successive integers, starting with 0.\n");
    printf("-opart OutPart\t\t:File where the final partitioning is stored.\n");
    printf("\tEach line is of the form \"node community\".\n");
    printf("-omam OutMAM\t\t:File where the MAM is written as a weighted edgelist (Default: MAM not written).\n");
    printf("\tEach line corresponds to an edge : \"source target weight\".\n");
    printf("\tEach edge appears twice (\"source target weight\" and \"target source weight\").\n");
    printf("-oppart OutPPart\t:File where the partial partitioning is stored (Default: not stored).\n");
    printf("\tEach line is of the form \"node community\".\n");
    printf("\tNodes disconnected in the MAM have their own community.\n");
    printf("-m Motif\t\t:Motif used to build the MAM (Default: E2).\n");
    printf("\tIdentifiers of all admissible motifs listed in GraphletIdentifiers.pdf.\n");
    printf("-nth NbTh\t\t:Number of threads to use to build the MAM (Default: sequential).\n");
    printf("\tWith only this program running on a laptop, best choice is the CPU number of threads (often 8).\n");
    printf("-c alpha\t\t:Resolution parameter for the modularity to use in the Louvain algorithm (Default: 1).\n");
    printf("\tGenerally, the higher the value of alpha, the smaller the communities.\n");
    printf("-l level\t\t:Level returned by the Louvain algorithm (fine or coarse) (Default: coarse).\n");
    printf("\tlevel=1 for the fine (larger number of communities), level=-2 for the coarse.\n");
    printf("-cc beta\t\t:Resolution parameter for the modularity to use in the postprocessing (Default: 1e-3).\n");
    printf("-k kmin\t\t\t:Maximum size of (meta-)nodes authorised to be merged in the postprocessing (Default: 1).\n");
  }
  else if (buildMAM && partition){
    printf("USAGE FOR BUILDING THE MAM AND PARTITIONING IT::\n");
    printf("\t./PartMAM -mapa -igraph InEdgeList -opart OutPart [-omam OutMAM -m Motif -nth NbTh -c alpha -l level]\n");
    printf("-igraph InEdgeList\t:File containing the network as an edgelist with integer nodes.\n");
    printf("\tNodes should be labelled with successive integers, starting with 0.\n");
    printf("-opart OutPart\t\t:File where the partial partitioning is stored.\n");
    printf("\tEach line is of the form \"node community\".\n");
    printf("-omam OutMAM\t\t:File where the MAM is written as a weighted edgelist (Default: MAM not written).\n");
    printf("\tEach line corresponds to an edge : \"source target weight\".\n");
    printf("\tEach edge appears twice (\"source target weight\" and \"target source weight\").\n");
    printf("-m Motif\t\t:Motif used to build the MAM (Default: E2).\n");
    printf("\tIdentifiers of all admissible motifs listed in GraphletIdentifiers.pdf.\n");
    printf("-nth NbTh\t\t:Number of threads to use to build the MAM (Default: sequential).\n");
    printf("\tWith only this program running on a laptop, best choice is the CPU number of threads (often 8).\n");
    printf("-c alpha\t\t:Resolution parameter for the modularity to use in the Louvain algorithm (Default: 1).\n");
    printf("\tGenerally, the higher the value of alpha, the smaller the communities.\n");
    printf("-l level\t\t:Level returned by the Louvain algorithm (fine or coarse) (Default: coarse).\n");
    printf("\tlevel=1 for the fine (larger number of communities), level=-2 for the coarse.\n");
  }
  else if (partition && postproc){
    printf("USAGE FOR PARTITIONING A MAM AND POSTPROCESSING THE PARTITION::\n");
    printf("\t./PartMAM -papo -imam InMAM -opart OutPart [-igraph InDir -isym InSym -oppart OutPPart -c alpha -l level -cc beta -k kmin]\n");
    printf("-imam InMAM\t\t:File containing the MAM as returned by the software.\n");
    printf("-opart OutPart\t\t:File where the final partitioning is stored.\n");
    printf("\tEach line is of the form \"node community\".\n");
    printf("-igraph InDir\t\t:File containing the directed network as an edgelist with integer nodes.\n");
    printf("\tNodes should be labelled with successive integers, starting with 0.\n");
    printf("-isym InSym\t\t:File containing the symmetrised version of the network (computed with motif E2).\n");
    printf("/!\\/!\\/!\\\t\t -igraph or -isym must be filled. \t\t/!\\/!\\/!\\\n");
    printf("-oppart OutPPart\t:File where the partial partitioning is stored.\n");
    printf("\tEach line is of the form \"node community\".\n");
    printf("-c alpha\t\t:Resolution parameter for the modularity to use in the Louvain algorithm (Default: 1).\n");
    printf("\tGenerally, the higher the value of alpha, the smaller the communities.\n");
    printf("-l level\t\t:Level returned by the Louvain algorithm (fine or coarse) (Default: coarse).\n");
    printf("\tlevel=1 for the fine (larger number of communities), level=-2 for the coarse.\n");
    printf("-cc beta\t\t:Resolution parameter for the modularity to use in the postprocessing (Default: 1e-3).\n");
    printf("-k kmin\t\t\t:Maximum size of (meta-)nodes authorised to be merged in the postprocessing (Default: 1).\n");
  }
  else if (buildMAM){
    printf("USAGE FOR BUILDING THE MAM ::\n");
    printf("\t./PartMAM -ma -igraph InEdgeList -omam OutMAM [-m Motif -nth NbTh]\n");
    printf("-igraph InEdgeList\t:File containing the network as an edgelist with integer nodes.\n");
    printf("\tNodes should be labelled with successive integers, starting with 0.\n");
    printf("-omam OutMAM\t\t:File where the MAM is written as a weighted edgelist.\n");
    printf("\tEach line corresponds to an edge : \"source target weight\".\n");
    printf("\tEach edge appears twice (\"source target weight\" and \"target source weight\").\n");
    printf("-m Motif\t\t:Motif used to build the MAM (Default: E2).\n");
    printf("\tIdentifiers of all admissible motifs listed in GraphletIdentifiers.pdf.\n");
    printf("-nth NbTh\t\t:Number of threads to use to build the MAM (Default: sequential).\n");
    printf("\tWith only this program running on a laptop, best choice is the CPU number of threads (often 8).\n");
  }
  else if (partition){
    printf("USAGE FOR PARTITIONING A GRAPH/A MAM ::\n");
    printf("\t./PartMAM -pa -opart OutPart [-igraph InDir -isym InMAM -c alpha -l level]\n");
    printf("-opart OutPart\t:File where the partial partitioning is stored.\n");
    printf("\tEach line is of the form \"node community\".\n");
    printf("-igraph InDir\t:File containing the directed network as an edgelist with integer nodes.\n");
    printf("\tNodes should be labelled with successive integers, starting with 0.\n");
    printf("-isym InMAM\t:File containing the MAM as returned by the software.\n");
    printf("/!\\/!\\/!\\\t\t -igraph or -isym must be filled. \t\t/!\\/!\\/!\\\n");
    printf("-c alpha\t:Resolution parameter for the modularity to use in the Louvain algorithm (Default: 1).\n");
    printf("\tGenerally, the higher the value of alpha, the smaller the communities.\n");
    printf("-l level\t:Level returned by the Louvain algorithm (fine or coarse) (Default: coarse).\n");
    printf("\tlevel=1 for the fine (larger number of communities), level=-2 for the coarse.\n");
  }
  else if (postproc){
    printf("USAGE FOR POSTPROCESSING A PARTITION ::\n");
    printf("\t./PartMAM -po -ipart InPart -opart OutPart [-igraph InDir -isym InSym -cc beta -k kmin]\n");
    printf("-inpart InPart\t:File containing the partial partitioning to use, as returned by the software.\n");
    printf("-opart OutPart\t:File where the final partitioning is stored.\n");
    printf("\tEach line is of the form \"node community\".\n");
    printf("-igraph InDir\t:File containing the directed network as an edgelist with integer nodes.\n");
    printf("\tNodes should be labelled with successive integers, starting with 0.\n");
    printf("-isym InSym\t:File containing the symmetrised version of the network (computed with motif E2).\n");
    printf("/!\\/!\\/!\\\t\t -igraph or -isym must be filled. \t\t/!\\/!\\/!\\\n");
    printf("-cc beta\t:Resolution parameter for the modularity to use in the postprocessing (Default: 1e-3).\n");
    printf("-k kmin\t\t:Maximum size of (meta)-nodes authorised to be merged in the postprocessing (Default: 1).\n");
  }
}


int main(int argc, char **argv){

  srand(time(NULL)+getpid());
  ReadArgs(argc,argv);
  if (help){
    usage();
    return 1;
  }

  /****************************************************************************/
  /*                              Whole pipeline                              */
  if (buildMAM && partition && postproc){
    fprintf(stderr,"************************************************************\n");
    fprintf(stderr,"WHOLE PIPELINE\n");
    fprintf(stderr,"************************************************************\n");
    if (InGraph.Empty() || OutPart.Empty()){
      fprintf(stderr,"ERROR::Whole pipeline:\n");
      fprintf(stderr,"\t\t-igraph and -opart must be filled\n");
      usage();
      return -1;
    }
    Graph MAMGraph,SGraph;
    fprintf(stderr,"=============== Build the MAM...\n");
    BuildMAM(InGraph,Motif,MAMGraph,OutMAM);
    fprintf(stderr,"...MAM Built ===============\n");
    fprintf(stderr,"=============== Partition the MAM...\n");
    if (OutPTmp.Empty()){
      OutPTmp = OutPart;
    }
    Partition(MAMGraph,OutPTmp);
    fprintf(stderr,"...MAM Partitioned ===============\n");
    fprintf(stderr,"=============== Postprocess the partition...\n");
    BuildMAM(InGraph,"E2",SGraph,"");
    PostProc(SGraph,OutPTmp);
    fprintf(stderr,"...Partition Postprocessed ===============\n");
  }
  /****************************************************************************/
  /*                       Build a MAM and Partition it                       */
  else if (buildMAM && partition){
    fprintf(stderr,"************************************************************\n");
    fprintf(stderr,"BUILD A MAM & PARTITION IT\n");
    fprintf(stderr,"************************************************************\n");
    if (InGraph.Empty() || OutPart.Empty()){
      fprintf(stderr,"ERROR::Build MAM & Partition it\n");
      fprintf(stderr,"\t\t-igraph and -opart must be filled\n");
      usage();
      return -1;
    }
    Graph MAMGraph;
    fprintf(stderr,"=============== Build the MAM...\n");
    BuildMAM(InGraph,Motif,MAMGraph,OutMAM);
    fprintf(stderr,"...MAM Built ===============\n");
    fprintf(stderr,"=============== Partition the MAM...\n");
    Partition(MAMGraph,OutPart);
    fprintf(stderr,"...MAM Partitioned ===============\n");
  }
  /****************************************************************************/
  /*             Partition a MAM and Postprocess the partition                */
  else if (partition && postproc){
    fprintf(stderr,"************************************************************\n");
    fprintf(stderr,"PARTITION A MAM & POSTPROCESS THE PARTITION\n");
    fprintf(stderr,"************************************************************\n");
    if (InMAM.Empty() || InGraph.Empty() || OutPart.Empty()){
      fprintf(stderr,"ERROR::Partition a MAM & Postprocess the partition\n");
      fprintf(stderr,"\t\t-igraph or -isym and -imam and -opart must be filled\n");
      usage();
      return -1;
    }
    fprintf(stderr,"=============== Partition the MAM...\n");
    Graph MAMGraph(InMAM.CStr());
    if (OutPTmp.Empty()){
      OutPTmp = OutPart;
    }
    Partition(MAMGraph,OutPTmp);
    fprintf(stderr,"...MAM Partitioned ===============\n");
    fprintf(stderr,"=============== Postprocess the partition...\n");
    Graph SGraph;
    if (indir){
      BuildMAM(InGraph,"E2",SGraph,"");
    }
    else{
      SGraph = Graph(InGraph.CStr());
    }
    PostProc(SGraph,OutPTmp);
    fprintf(stderr,"...Partition Postprocessed ===============\n");
  }
  /****************************************************************************/
  /*                         Postprocess a partition                          */
  else if (postproc){
    fprintf(stderr,"************************************************************\n");
    fprintf(stderr,"POSTPROCESS A PARTITION\n");
    fprintf(stderr,"************************************************************\n");
    if (InPart.Empty() || InGraph.Empty() || OutPart.Empty()){
      fprintf(stderr,"ERROR::Postprocess a partition\n");
      fprintf(stderr,"\t\t-igraph or -isym and -ipart and -opart must be filled\n");
      usage();
      return -1;
    }
    fprintf(stderr,"=============== Postprocess the partition...\n");
    Graph SGraph;
    if (indir){
      BuildMAM(InGraph,"E2",SGraph,"");
    }
    else{
      SGraph = Graph(InGraph.CStr());
    }
    PostProc(SGraph,InPart);
    fprintf(stderr,"...Partition Postprocessed ===============\n");
  }
  /****************************************************************************/
  /*                          Partition a graph/MAM                           */
  else if (partition){
    fprintf(stderr,"************************************************************\n");
    fprintf(stderr,"PARTITION A DIRECTED GRAPH OR A MAM\n");
    fprintf(stderr,"************************************************************\n");
    if (InGraph.Empty() || OutPart.Empty()){
      fprintf(stderr,"ERROR::Partition a Graph or a MAM\n");
      fprintf(stderr,"\t\t-igraph or -isym and -opart must be filled\n");
      usage();
      return -1;
    }
    fprintf(stderr,"=============== Partition the MAM...\n");
    Graph AGraph;
    if (indir){
      BuildMAM(InGraph,"E2",AGraph,"");
    }
    else{
      AGraph = Graph(InGraph.CStr());
    }
    Partition(AGraph,OutPart);
    fprintf(stderr,"...MAM Partitioned ===============\n");
  }
  /****************************************************************************/
  /*                                 Build a MAM                              */
  else if (buildMAM){
    fprintf(stderr,"************************************************************\n");
    fprintf(stderr,"BUILD A MAM\n");
    fprintf(stderr,"************************************************************\n");
    if (InGraph.Empty() || OutMAM.Empty()){
      fprintf(stderr,"ERROR::Build a MAM\n");
      fprintf(stderr,"\t\t-igraph and -omam must be filled\n");
      usage();
      return -1;
    }
    Graph Dmy;
    fprintf(stderr,"=============== Build the MAM...\n");
    BuildMAM(InGraph,Motif,Dmy,OutMAM);
    fprintf(stderr,"...MAM Built ===============\n");
  }
  else{
    TExcept::Throw("Nonconsistent input parameters.");
  }
  return 0;
}
