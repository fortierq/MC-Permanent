#include "MH_Matching.h"

using namespace std;
using namespace ogdf;

MH_Matching::MH_Matching(Graph &G, int szMatching, int nVertices, bool stats, bool not_uniform):
szMatching(szMatching), nVertices(nVertices), G(G), stats(stats), nBigMatching(0), not_uniform(not_uniform),
nSmallMatching(0.), sum_degree(0.)
{
	idToEdge.resize(G.numberOfEdges());
	edge e;
	forall_edges(e, G)
		idToEdge[e->index()] = e;
}

edge MH_Matching::contain(ogdf::Graph &S, ogdf::edge &e)
{
	edge e1;
	forall_edges(e1, S)
	{
		if( (e->source()->index() == e1->source()->index() && e->target()->index() == e1->target()->index())
			|| (e->source()->index() == e1->target()->index() && e->target()->index() == e1->source()->index()) )
			return e1;
	}
	return 0;
}

int MH_Matching::isVisited()
{
	for(int i = 0; i < visited.size(); i++)
	{
		if(visited[i].numberOfEdges() != state.numberOfEdges())
			continue;
		edge e;
		forall_edges(e, visited[i])
			if(contain(state, e) == 0) //if(idToEdge_state[e->index()] == 0)
				goto next;
		return i;
		next: continue;
	}
	return -1;
}

double MH_Matching::ratio()
{
	if(nSmallMatching == 0.)
		return 1.;
	return nBigMatching / nSmallMatching;
}

void MH_Matching::increase()
{
	if(!not_uniform)
	{
		if(state.numberOfEdges() == szMatching)
			nBigMatching++;
		else
			nSmallMatching++;
	}
	else
	{
		if(state.numberOfEdges() == szMatching)
			nBigMatching += sum_degree;
		else
			nSmallMatching += sum_degree;
	}
}

void MH_Matching::init()
{
	sum_degree = 0.;
	visited.clear();

	int nEdgesG = G.numberOfEdges();
	idToNode_state.resize(G.numberOfNodes(), 0);
	idToEdge_state.resize(nEdgesG, 0);

	node v;
	forall_nodes(v, G)
		idToNode_state[v->index()] = state.newNode(v->index());

	// Find first matching
	MinCostFlowReinelt maxFlow;
	node s = G.newNode();
	node t = G.newNode();
	
	forall_nodes(v, G)
	{
		if(v == s || v == t) continue;
		if(v->index() < nVertices) 
			G.newEdge(s, v);
		else
			G.newEdge(v, t);
	}
	
	EdgeArray<int> lowerBound(G, 0), upperBound(G, 1), cost(G, 0), flow(G, 0);
	NodeArray<int> supply(G, 0);
	
	supply[s] = szMatching;
	supply[t] = -szMatching;

	/*cout<<maxFlow.checkProblem(G, lowerBound, upperBound, supply)<<endl;*/
	maxFlow.call(G, lowerBound, upperBound, cost, supply, flow);

	write_graph(G, "G.gml");
	G.delNode(s);
	G.delNode(t);
	for(int i = 0; i < G.numberOfEdges(); i++)
	{
		//cout<<flow[i]<<endl;
		if(flow[i] == 1)
		{
			edge flowEdge = idToEdge[i];
			idToEdge_state[i] = state.newEdge(idToNode_state[flowEdge->source()->index()],
				idToNode_state[flowEdge->target()->index()], i);
			sum_degree += (flowEdge->source()->degree() + flowEdge->target()->degree());
		}
	}

	if(stats)
	{
		visited.push_back(state);
		nVisit.clear();
		nVisit.push_back(1);
	}
	
	increase();
}

void MH_Matching::write_graph(Graph &H, string file_name)
{
	GraphAttributes GA(H, GraphAttributes::nodeGraphics | GraphAttributes::edgeGraphics |
		GraphAttributes::nodeLabel | GraphAttributes::nodeColor | 
		GraphAttributes::edgeColor | GraphAttributes::edgeStyle | 
		GraphAttributes::nodeStyle | GraphAttributes::nodeTemplate);
	GA.directed(false);
	node v;
	forall_nodes(v, H)
	{
		GA.width(v) = GA.height(v) = 40.0;
		int index = v->index();
		if(index < nVertices)
		{
			GA.colorNode(v) = "#0000FF";
			GA.y(v) = 0 + index * 100;
			GA.x(v) = 0;
		}
		else
		{
			if(index < 2*nVertices)
			{
				GA.colorNode(v) = "#FF0000";
				GA.y(v) = 0 + (index - nVertices)* 100;
				GA.x(v) = 500;
			}
			else
				GA.colorNode(v) = "#000000";
		}
	}
	GA.writeGML(file_name.c_str());
}

bool MH_Matching::accept(double acceptance)
{
	double u = rand()/RAND_MAX;
	if(u <= acceptance)
		return true;
	return false;
}

void MH_Matching::step()
{
	double p = (double)rand()/RAND_MAX;
	
	if(p <= 0.5)
		return;

	double randEdge = rand() % idToEdge.size();
	edge e = idToEdge[randEdge];
	edge eInState = idToEdge_state[e->index()];

	if(eInState != 0)
	{
		if(state.numberOfEdges() == szMatching)
		{
			// even with not_uniform with accept with proba 1
			idToEdge_state[eInState->index()] = 0;
			state.delEdge(eInState);
			sum_degree -= (e->source()->degree() + e->target()->degree());
		}
	}
	else
	{
		if(state.numberOfEdges() == szMatching - 1)
		{
			node u = idToNode_state[e->source()->index()], v = idToNode_state[e->target()->index()];
			edge e_u = 0, e_v = 0;
			
			forall_adj_edges(e_u, u) { break; } 
			forall_adj_edges(e_v, v) { break; }
			
			if(e_u == 0 && e_v != 0)
			{
				if(not_uniform)
				{
					double lost_degree = e_v->source()->degree() + e_v->target()->degree();
					double add_degree = u->degree() + v->degree();
					
					double proposal = sum_degree + add_degree - lost_degree;
					if(!accept(sum_degree/proposal)) return;
					sum_degree = proposal;
				}
				idToEdge_state[e_v->index()] = 0;
				state.delEdge(e_v);
				edge newEdge = state.newEdge(u, v, e->index());
				idToEdge_state[newEdge->index()] = newEdge;
			}
			else if(e_u != 0 && e_v == 0)
			{
				if(not_uniform)
				{
					double lost_degree = e_u->source()->degree() + e_u->target()->degree();
					double add_degree = u->degree() + v->degree();
					double proposal = sum_degree + add_degree - lost_degree;
					if(!accept(sum_degree/proposal)) return;
					sum_degree = proposal;
				}
				idToEdge_state[e_u->index()] = 0;
				state.delEdge(e_u);
				edge newEdge = state.newEdge(u, v, e->index());
				idToEdge_state[newEdge->index()] = newEdge;
			}
			else if(e_u == 0 && e_v == 0)
			{
				if(not_uniform)
				{
					double add_degree = u->degree() + v->degree();
					double proposal = sum_degree + add_degree;
					if(!accept(sum_degree/proposal)) return;
					sum_degree = proposal;
				}
				edge newEdge = state.newEdge(u, v, e->index());
				idToEdge_state[newEdge->index()] = newEdge;
			}
		}
	}

	increase();
	if(stats)
	{
		int vis = isVisited();
		if(vis == -1)
		{
			visited.push_back(Graph(state));
			nVisit.push_back(1);
		}
		else
			nVisit[vis]++;
	}
}