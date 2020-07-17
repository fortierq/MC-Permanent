#ifndef IMPORTANCE_MATCHING_H
#define IMPORTANCE_MATCHING_H

#include "armadillo"
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <time.h>

class Importance_Matching
{
	arma::mat &A;
	double step(bool stats);
	bool stats;

public:
	std::vector<std::string> idToPerm;
	std::vector<int> idToFreqPerm;
	double ratio_accepted;
	Importance_Matching(arma::mat &A, bool stats);
	double get_approx(std::string file_name, int nSteps, int step_write, double max_second);
};

#endif