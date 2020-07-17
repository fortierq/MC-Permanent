#include "Importance_Matching.h"

using namespace arma;
using namespace std;

template<typename type> string toString(type x)
{
	std::ostringstream os;
	os <<fixed<<setprecision(0)<< x;
	return os.str();
}

Importance_Matching::Importance_Matching(mat &A, bool stats):
	A(A), ratio_accepted(0.), stats(stats)
{
}


double Importance_Matching::get_approx(string file_name, int nSteps, int step_write, double max_second)
{
	double nAccepted = 0.;
	time_t t1 = clock();
	double approx = 0.;
	string data = string("M_Imp = c(" + toString(approx));

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
		double s = step(stats);
		if(s != 0) // if s == 0, we reject
		{
			nAccepted++;
			approx += 1/s;
		}
		if(i % step_write == 1) // i != 0
		{
			data += ("," + toString(approx/i));
		}
	}
	ratio_accepted = nAccepted / nSteps;
	approx /= nSteps;
	time_t t2 = clock();
	data += ")\n";

	ofstream outfile(file_name.c_str(), std::ios::out | std::ios::binary); 
	outfile<<"# " + toString(nSteps) + " iterations, " + toString(A.n_rows) + " vertices\n";
	outfile<<"# time: "<<toString(t2-t1)<<"\n";
	outfile<<data;

	//string mean = string(" = c(" + toString(approx));
	return approx;
}

double Importance_Matching::step(bool stats)
{
	mat tmp(A);
	double ret = 1.;
	string perm = "";

	for(int n = 0; n < tmp.n_cols; n++) // we choose in column number n
	{
		int row_chosen = 0;

		colvec d = sum(tmp,1);
		colvec p(d.n_elem);
		p.fill(0.);
		double sum = 0.;

		int rowWithOne = -1; // if != -1, this is the number where there is a row sum equal to 1

		for(int i = 0; i < d.n_elem; i++)
		{
			if(d(i) == 1 && (tmp(i, n) == 1))
			{
				if(rowWithOne != -1) // if there are two row sum equal to 1
					return 0.;
				rowWithOne = i;
			}
		}
		if(rowWithOne != -1)
		{
			row_chosen = rowWithOne;
			goto chosen;
		}

		for(int i = 0; i < d.n_elem; i++)
		{
			if(tmp(i, n) == 0)
				p(i) = 0.;
			else p(i) = 1. / (d(i) - 1);
			sum += p(i);
		}

		if(sum == 0.) // we can't chose any other row
			return 0.;

		p /= sum; 
		
		// we choose randomly
		double u = ((double)rand()) / RAND_MAX;
		//cout<<u<<endl;
		double sum_p = 0.;
		for(row_chosen = 0; row_chosen< d.n_elem; row_chosen++)
		{
			sum_p += p(row_chosen);
			if(u <= sum_p)
				break;
		}
		if(row_chosen == d.n_elem) // avoid rounding error
		{
			for(row_chosen = d.n_elem - 1; row_chosen >= 0 ; row_chosen--)
			{
				if(p(row_chosen) != 0.)
					break;
			}
		}
		ret *= p(row_chosen);

chosen:
		if(stats)
			perm += toString(row_chosen);
		
		for(int i = 0; i < tmp.n_rows; i++)
			tmp(i, n) = 0;
		for(int i = 0; i < tmp.n_cols; i++)
			tmp(row_chosen, i) = 0;
	}
	if(stats)
	{
		int id = 0;
		for(; id < idToPerm.size(); id++)
			if(idToPerm[id] == perm)
				break;
		if(id == idToPerm.size())
		{
			idToPerm.push_back(perm);
			idToFreqPerm.push_back(1);
		}
		else
			idToFreqPerm[id]++;
	}

	return ret;
}
