#include "Rejection_Matching.h"

using namespace std;
using namespace ogdf;

template<typename type> string toString(type x)
{
	std::ostringstream os;
	os <<fixed<<setprecision(0)<< x;
	return os.str();
}


Rejection_Matching::Rejection_Matching(Graph &G):
	G(G)
{
	node v;
	forall_nodes(v, G)
		idToNode.push_back(v);
}

bool Rejection_Matching::rdmSample(int szMatching)
{
	vector<int> idUsed;
	for(int i = 0; i < szMatching; i++)
	{
		node u = idToNode[i];
		node v_rdm = 0;
		do{
			int rdm = rand() % szMatching;
			v_rdm = idToNode[rdm + szMatching];
		}while(find(idUsed.begin(), idUsed.end(), v_rdm->index()) != idUsed.end());
		edge e = 0;
		forall_adj_edges(e, u)
			if(e->opposite(u) == v_rdm)
				goto found;
		return false;
found:
		idUsed.push_back(v_rdm->index());
	}
	return true;
}

double Rejection_Matching::fact(int n)
{
	double ret = 1;
	for(int i = 2; i <= n; i++)
		ret *= i;
	return ret;
}

double Rejection_Matching::get_ratio_accepted(int nSteps, int szMatching, double max_second)
{
	time_t t1 = clock();
	double nMatching = 0.;
	for(int i = 0; i < nSteps; i++)
	{
		if(max_second > 0)
		{
			time_t t3 = clock();
			if((t3 - t1)/1000. >= max_second)
			{
				nSteps = i;
				break;
			}
		}
		if(rdmSample(szMatching))
			nMatching++;
	}
	
	return ((double)nMatching / nSteps);
}
double Rejection_Matching::get_approx_matching(string file_name, int nSteps, int step_write, 
											   int szMatching, double max_second)
{
	time_t t1 = clock();
	double nMatching = 0.;
	string data = "M_Rej = c(0";

	for(int i = 1; i < nSteps; i++)
	{
		if(max_second > 0)
		{
			time_t t3 = clock();
			if((t3 - t1)/1000. >= max_second)
			{
				nSteps = i;
				break;
			}
		}
		if(rdmSample(szMatching))
			nMatching++;
		if(i % step_write == 0)
		{
			data += ",";
			data += toString(((double)nMatching / i)* fact(szMatching));
		}
	}
	data += ")\n";
	time_t t2 = clock();

	ofstream outfile(file_name.c_str(), std::ios::out | std::ios::binary); 
	outfile<<"# " + toString(nSteps) + " iterations, " + toString(szMatching) + " vertices\n";
	outfile<<"# time: "<<toString(t2-t1)<<"\n";
	outfile<<data;

	//cout<<nMatching<<" "<<nSteps<<endl;
	double ratio_acc = ((double)nMatching / nSteps);
	
	double approx = ratio_acc * fact(szMatching);//binomial(szMatching*szMatching, szMatching);
	return approx;
}