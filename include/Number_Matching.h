#ifndef NUMBER_MATCHING_H
#define NUMBER_MATCHING_H

#include "MH_Matching.h"
#include <vector>
#include <ogdf/basic/Graph.h>
#include <ogdf/graphalg/MinCostFlowReinelt.h>
#include <ogdf/basic/GraphAttributes.h>
#include <string>
#include <iostream>
#include <time.h>
#include <fstream>
#include <iomanip>

class Number_Matching
{
private:
	ogdf::Graph G;
	int nVertices;
	int sz;
	std::vector<MH_Matching*> mhs;
	bool stats;
	bool use_uniform;

public:
	Number_Matching(ogdf::Graph &G, int nVertices, int sz, bool stats, bool use_uniform);
	double get_approx_matching(std::string file_name, int nIter, int step_write, double max_second);
	double get_ratio(int max);
	double get_all_ratio();
};

#endif