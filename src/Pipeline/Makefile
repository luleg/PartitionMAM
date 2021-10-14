all: pipeline rmo

pipeline: benson louvain
	g++ -g3 -O0 -fopenmp -lrt -o pipeline pipeline.cpp louvain.o modularity.o quality.o graph_binary.o BensonGraphQuad.o BensonGraphTriangle.o BensonGraphEdge.o motifType.o ../UtilsSNAP/snap-core/Snap.o -I../UtilsSNAP/snap-core -I../UtilsSNAP/glib-core -I../UtilsSNAP/snap-adv -I../Louvain

benson:
	g++ -g3 -O0 -fopenmp -I../UtilsSNAP/snap-core/ -I../UtilsSNAP/glib-core -c ../UtilsSNAP/snap-adv/BensonGraphEdge.cpp ../UtilsSNAP/snap-adv/BensonGraphQuad.cpp ../UtilsSNAP/snap-adv/BensonGraphTriangle.cpp ../UtilsSNAP/snap-adv/motifType.cpp

louvain:
	g++ -g3 -O0 -fopenmp -I../UtilsSNAP/snap-core/ -I../UtilsSNAP/glib-core -c ../Louvain/graph_binary.cpp ../Louvain/quality.cpp ../Louvain/modularity.cpp ../Louvain/louvain.cpp

rmo:
	rm *.o
