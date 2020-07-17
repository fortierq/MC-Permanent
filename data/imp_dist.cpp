#include "MH_Matching.h"
#include "Number_Matching.h"
#include <sstream>
#include "Rejection_Matching.h"
#include "Importance_Matching.h"
#include <iomanip>

using namespace arma;
using namespace std;
using namespace ogdf;

template<typename type> string toString(type x)
{
	std::ostringstream os;
	os <<fixed<<setprecision(10)<<x;
	return os.str();
}

void graph_to_mat(Graph &G, mat &M, int sz)
{
	M.set_size(sz, sz);
	M.fill(0);
	edge e;
	forall_edges(e, G)
	{
		if(e->source()->index() < sz)
			M(e->source()->index(), e->target()->index() - sz) = 1;
		else
			M(e->target()->index() - sz, e->source()->index()) = 1;
	}
}

void mat_to_graph(Graph &G, mat &M, int sz)
{
	vector<node> idToNode;

	for(int i = 0; i < 2*sz; i++)
		idToNode.push_back(G.newNode());

	for(int i = 0; i < sz; i++)
		for(int j = 0; j < sz; j++)
			if(M(i,j) == 1)
				G.newEdge(idToNode[i], idToNode[j+sz]);
}

void rand_one(mat &M, int sz, double density)
{
	M.set_size(sz, sz);
	M.fill(0.);
	for(int i = 0; i < sz; i++)
		for(int j = 0; j < sz; j++)
		{
			double u = ((double)rand()) / RAND_MAX;
			if (u < density)
				M(i,j) = 1.;
		}
}

int main()
{
	srand(time(NULL));
	int nNodes = 5;
	int szMatching = nNodes;

	double max_second = 20;
	int nSteps = 50000000;
	int step_write = 3000;/*nSteps / 500;*/
	
	//A.print("A:");
	//graph_to_mat(G, A, nNodes);
	//rand_one(A, 5);

	string dist_imp = "Dist_Imp = c(";

	mat A;
	Graph G;
	rand_one(A, nNodes, 0.9);
	Importance_Matching imp(A);
	imp.get_approx("imp.r", nSteps, step_write, max_second, true);
	
	for(int i = 0; i < imp.idToFreqPerm.size(); i++)
	{
		if(i == imp.idToFreqPerm.size() - 1)
			dist_imp += toString(imp.idToFreqPerm[i]) + ")\n";
		else
			dist_imp += (toString(imp.idToFreqPerm[i]) + ",");
	}
	
	ofstream outfile("distrib_imp.r", std::ios::out | std::ios::binary); 
	outfile<<dist_imp;
	outfile.close();
}