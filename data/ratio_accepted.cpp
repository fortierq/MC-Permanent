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
	int nNodes = 10;
	int szMatching = nNodes;

	double max_second = 1;
	int nSteps = 50000000;
	int step_write = 3000;/*nSteps / 500;*/
	
	//A.print("A:");
	//graph_to_mat(G, A, nNodes);
	//rand_one(A, 5);

	string ratio_imp = "Rat_Imp = c(0";
	string ratio_rej = "Rat_Rej = c(0";
	string density = "D = c(0";

	double nTests = 10;
	for(double i = 3; i <= nTests; i++)
	{
		mat A;
		Graph G;
		double d = ((double)i)/(double)nTests;
		std::ostringstream os3;
		os3 <<fixed<<setprecision(1)<<d;
		density += string("," + os3.str());
		rand_one(A, nNodes, d);
		Importance_Matching imp(A);
		imp.get_approx("imp.r", nSteps, step_write, max_second);
		std::ostringstream os;
		os <<fixed<<setprecision(10)<<imp.ratio_accepted;
		ratio_imp += string(",") + os.str();
		mat_to_graph(G, A, nNodes);
		Rejection_Matching rej(G);
		std::ostringstream os2;
		os2 <<fixed<<setprecision(10)<< rej.get_ratio_accepted(nSteps, szMatching, max_second);
		ratio_rej += string(",") + os2.str();
	}
	ofstream outfile("ratio_accepted.r", std::ios::out | std::ios::binary); 
	outfile<<ratio_imp<<")\n"<<ratio_rej<<")\n"<<density<<")\n";
	outfile.close();
}