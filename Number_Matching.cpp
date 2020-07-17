#include "Number_Matching.h"

using namespace std;
using namespace ogdf;

Number_Matching::Number_Matching(Graph &G, int nVertices, int sz, bool stats, bool use_uniform):
	G(G), nVertices(nVertices), sz(sz), stats(stats), use_uniform(use_uniform)
{
	for(int k = 2; k <= sz; k++)
	{
		mhs.push_back(new MH_Matching(G, k, nVertices, stats, use_uniform));
		mhs[mhs.size() - 1]->init();
	}
}

double Number_Matching::get_ratio(int max)
{
	double ret = 1.;
	for(int i = 0; i < max; i++)
		ret *= mhs[i]->ratio();
	return ret;
}

double Number_Matching::get_all_ratio()
{
	return get_ratio(mhs.size());
}

template<typename type> string toString(type x)
{
	std::ostringstream os;
	os <<fixed<<setprecision(0)<< x;
	return os.str();
}

double Number_Matching::get_approx_matching(string file_name, int nSteps, int step_write, double max_second)
{
	time_t t1 = clock();
	int approx = G.numberOfEdges();

	vector<string> data;
	for(int i = 0; i < mhs.size(); i++)
		data.push_back("M" + toString(i+2) + "_MH = c(" + toString(approx));

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
		
		for(int j = 0; j < mhs.size(); j++)
			mhs[j]->step();
		if(i % step_write == 0)
		{
			for(int i = 0; i < mhs.size(); i++)
			{
				data[i] += ",";
				data[i] += toString(approx * get_ratio(i + 1));
			}
		}
		/*if(i % 100000 == 0)
			cout<<approx * get_all_ratio()<<endl;		*/
	}
	for(int i = 0; i < mhs.size(); i++)
		data[i] += ")\n";
	
	time_t t2 = clock();
	ofstream outfile(file_name.c_str(), std::ios::out | std::ios::binary); 
	outfile<<"# " + toString(nSteps) + " iterations, " + toString(nVertices) + " vertices\n";
	outfile<<"# time: "<<toString(t2-t1)<<"\n";
	string mean = string("EM = c(" + toString(approx));
	
	for(int i = 0; i < mhs.size(); i++)
	{
		outfile<<data[i];
		mean += ", mean(M" + toString(i+2) + ")";
	}
	mean += ")\n";
	outfile<<mean;
	outfile.close();

	return approx * get_all_ratio();
}
