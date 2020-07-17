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
	os <<setprecision(100)<< (int)x;
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
	Graph G;
	int nNodes = 10;
	int szMatching = nNodes;
	/*vector<node> idToNode;
	for(int i = 0; i < 2*nNodes; i++)
		idToNode.push_back(G.newNode());

	for(int i = 0; i < nNodes; i++)
		for(int j = nNodes; j < 2*nNodes; j++)
			G.newEdge(idToNode[i], idToNode[j]);*/


	double max_second = 0.5;
	int nSteps = 50000000;
	int step_write = 1000;/*nSteps / 500;*/
	mat A;
	rand_one(A, nNodes, 0.7);
	A.print("A:");
	//graph_to_mat(G, A, nNodes);
	//rand_one(A, 5);

	Importance_Matching imp(A, false);
	cout<<"Importance: "<<imp.get_approx("imp.r", nSteps, step_write, max_second)<<endl;

	mat_to_graph(G, A, nNodes);
	Rejection_Matching rej(G);
	cout<<"Rejection: "<<rej.get_approx_matching("rej.r", nSteps, step_write, szMatching, max_second)<<endl;

	Number_Matching nMatching(G, nNodes, szMatching, false, false);
	cout<<"MH: "<<nMatching.get_approx_matching("MH.r", nSteps, step_write, max_second)<<endl;

	nSteps *= 10;
	step_write *= 5;
	max_second *= 5;
	Importance_Matching exact_imp(A, false);
	cout<<exact_imp.get_approx("exact_imp.r", nSteps, step_write, max_second)<<endl;
	Number_Matching exact_mh(G, nNodes, szMatching, false, false);
	cout<<exact_mh.get_approx_matching("exact_MH.r", nSteps, step_write, max_second)<<endl;

	/*MH_Matching mh(G, szMatching, nNodes, false);
	mh.init();
	
	for(int i = 0; i < 10; i++)
	{
		mh.step();
		if(i % 100 == 0)
			cout<<mh.ratio()<<endl;
		string file = "state";
		file += toString(i);
		file += ".gml";
		mh.write_graph(mh.state, file);

		string s = "gml2pic.exe ";
		s += file;
		system(s.c_str());
	}*/
	//cout<<endl<<endl;
	//int sum = 0;
	//for(int i = 0; i < mh.nVisit.size(); i++)
	//	sum += mh.nVisit[i];
	/*for(int i = 0; i < mh.nVisit.size(); i++)
		cout<<(double)mh.nVisit[i]/sum<<endl;*/
}