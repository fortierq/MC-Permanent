#ifndef MH_MATCHING_H
#define MH_MATCHING_H

#include <vector>
#include <ogdf/basic/Graph.h>
#include <ogdf/graphalg/MinCostFlowReinelt.h>
#include <ogdf/energybased/FMMMLayout.h>
#include <ogdf/basic/GraphAttributes.h>
#include <string>
#include <iostream>
#include <time.h>

class MH_Matching
{
	bool accept(double acceptance);
	void increase();
public:
	double sum_degree;
	bool not_uniform;
	bool stats;
	ogdf::Graph G;
	ogdf::Graph state;
	std::vector<ogdf::node> idToNode_state;
	std::vector<ogdf::edge> idToEdge_state;
	std::vector<ogdf::Graph> visited; // in the order of the visit by the chain
	std::vector<ogdf::edge> idToEdge;
	std::vector<int> nVisit; // nVisit[i] is the number of visits of visited[i]
	double nBigMatching, nSmallMatching;
	MH_Matching(ogdf::Graph &G, int szMatching, int nVertices, bool stats, bool not_uniform); 
	void step();
	void init();
	int isVisited(); // returns the index of H if it was visited, -1 otherwise
	ogdf::edge contain(ogdf::Graph &S, ogdf::edge &e); // return the edge e in S or 0 if it is not in S
	void write_graph(ogdf::Graph &H, std::string file_name);
	double ratio();

private:
	int szMatching, nVertices;
	std::string file_name;
};

#endif
