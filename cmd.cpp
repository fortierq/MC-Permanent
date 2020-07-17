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
	while(true)
	{
		srand(time(NULL));

		int nNodes = 7; // size matrix
		cout<<"Size of the matrix?"<<endl;
		cin>>nNodes;

		int szMatching = nNodes;

		double max_second = .1; // algorithms are run for max_second seconds
		cout<<"Maximum running time? (seconds)"<<endl;
		cin>>max_second;

		int nSteps = 50000000; // algorithms do at most nSteps steps

		int step_write = 100; // write every 100 steps
		cout<<"Write every what number of steps?"<<endl;
		cin>>step_write;

		double d = 0.7; // density
		cout<<"Density?"<<endl;
		cin>>d;

		mat A;
		Graph G;
		rand_one(A, nNodes, d); // we fill A randomly with ones

		cout<<endl<<"Approximation of the permanent of "<<endl;
		A.print("A:");
		cout<<endl;

		Importance_Matching imp(A, false);
		// run importance sampling and write result in imp.r
		cout<<"Importance: "<<imp.get_approx("imp.r", nSteps, step_write, max_second)<<endl; 

		mat_to_graph(G, A, nNodes); // construct the graph associated to A

		Rejection_Matching rej(G);
		cout<<"Rejection: "<<rej.get_approx_matching("rej.r", nSteps, 
			step_write, szMatching, max_second)<<endl;

		Number_Matching nMatching(G, nNodes, szMatching, false, false);
		cout<<"MH: "<<nMatching.get_approx_matching("mh.r", nSteps, step_write, max_second)<<endl;
		system("pause");
	}
}
