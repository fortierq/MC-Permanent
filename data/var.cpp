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
	int nNodes = 7;
	int szMatching = nNodes;

	double max_second = .1;
	int nSteps = 50000000;
	int step_write = 100;/*nSteps / 500;*/
	

	for(int i = 4; i <= 5; i+=1)
	{
		double d = ((double)i)/(double)10;
		mat A;
		Graph G;
		rand_one(A, nNodes, d);
		Importance_Matching imp(A, false);
		imp.get_approx("var_imp" + toString(i) + ".r", nSteps, step_write, max_second);
		mat_to_graph(G, A, nNodes);

		Rejection_Matching rej(G);
		cout<<"Rejection: "<<rej.get_approx_matching("var_rej" + toString(i) + ".r", nSteps, 
			step_write, szMatching, max_second)<<endl;

		Number_Matching nMatching(G, nNodes, szMatching, false, false);
		cout<<"MH: "<<nMatching.get_approx_matching("var_mh" + toString(i) + ".r", nSteps, step_write, max_second)<<endl;
	}
}
