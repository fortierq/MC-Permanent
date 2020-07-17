#ifndef REJECTION_MATCHING_H
#define REJECTION_MATCHING_H

#include "MH_Matching.h"
#include <vector>
#include <ogdf/basic/Graph.h>
#include <ogdf/graphalg/MinCostFlowReinelt.h>
#include <ogdf/basic/GraphAttributes.h>
#include <string>
#include <iostream>
#include "binomial.h"
#include <algorithm>
#include <fstream>
#include <iomanip>

class Rejection_Matching
{
	std::vector<ogdf::node> idToNode;
public:
	ogdf::Graph &G;
	Rejection_Matching(ogdf::Graph &G);
	bool rdmSample(int szMatching); // return true iff a random sample of szMatching edges is a matching 
	double get_approx_matching(std::string file_name, int nSteps, int step_write, int szMatching, double max_second);
	double fact(int n);
	double get_ratio_accepted(int nSteps, int szMatching, double max_second);
};
#endif