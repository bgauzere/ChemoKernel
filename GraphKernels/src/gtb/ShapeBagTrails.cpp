/*-----------------------------------------------------------------------
 
 File        : ShapeBagTrails.cpp
 
 Description : Compute Bag of Trails with heuristics for shape graph
 
 Copyright  : Francois-Xavier Dupé - http://www.greyc.ensicaen.fr/~fdupe/
 
 This software is governed by the CeCILL  license under French law and
 abiding by the rules of distribution of free software.  You can  use,
 modify and/ or redistribute the software under the terms of the CeCILL
 license as circulated by CEA, CNRS and INRIA at the following URL
 "http://www.cecill.info".
 
 As a counterpart to the access to the source code and  rights to copy,
 modify and redistribute granted by the license, users are provided only
 with a limited warranty  and the software's author,  the holder of the
 economic rights,  and the successive licensors  have only  limited
 liability.
 
 In this respect, the user's attention is drawn to the risks associated
 with loading,  using,  modifying and/or developing or reproducing the
 software by the user in light of its specific status of free software,
 that may mean  that it is complicated to manipulate,  and  that  also
 therefore means  that it is reserved for developers  and  experienced
 professionals having in-depth computer knowledge. Users are therefore
 encouraged to load and test the software's suitability as regards their
 requirements in conditions enabling the security of their systems and/or
 data to be ensured and,  more generally, to use and operate it in the
 same conditions as regards security.
 
 The fact that you are presently reading this means that you have had
 knowledge of the CeCILL license and that you accept its terms.
 
 ---------------------------------------------------------------------*/

#include "ShapeBagTrails.h"
#include "VertexCovering.h"
#include "pandore.h"
#include <set>

using namespace std;
using namespace pandore;
using namespace cimg_library;

double ShapeBagTrails::computeSP ( const CImg<int> &histo, int partie, int cut, int debut, int fin )
{
	double mean = 0;
	double tmp2 = 0;
	double totalh = 0;
	
	if(partie == 1)
    {
		for(int i = debut+1 ; i < cut+1; i++)
        {
			mean += histo(i-1) * (double)(i)/histo.width();
			totalh += histo(i-1);
        }
		
		mean/=totalh;
		if(totalh == 0) mean = 0;
		//printf(" partie 1 = %g\n",totalh);
		//printf("Mean part 1 = %g\n",mean);
		for(int i = debut +1 ; i < cut+1; i++)
        {
			//printf("mean%.4g\n",mean);
			double sigma = (histo(i-1)*(double)(i)/histo.width()) - mean;
			//printf("sigma%.4g\n",sigma);
			tmp2 += sigma * sigma * histo(i-1);
			//printf("%.4g\n",tmp2);
        }
    }
	else if(partie == 2)
    {
		for(int i = cut+1 ; i < fin+1; i++)
        {
			mean += histo(i-1) * (double)(i)/histo.width();
			totalh += histo(i-1);
        }
        
		mean/=totalh;
		if(totalh == 0) mean = 0;
		//printf(" partie 2 = %g\n",totalh);
		//printf("Mean part 2 = %g\n",mean);
        
		for(int i = cut + 1 ; i < fin+1; i++)
        {
			//printf("mean%.4g\n",mean);
			double sigma = histo(i-1) *(double)(i)/histo.width() - mean;
			//printf("sigma %.4g\n",sigma);
			tmp2 += sigma * sigma * histo(i-1);
			//printf("%.4g\n",tmp2);
        }
    }
	//printf("SIGMA PARTIE %d %g\n",partie,tmp2);
	return tmp2 ;
}

void ShapeBagTrails::separate ( const CImg<int> &histo, int debut, int fin, double * stock, int nbcut )
{
	double max = -1e99;
	int maxt = 0;
	double d = 0;
	double sigmac = 0;
	double mean = 0;
	int totalh = 0;
	int nbcut_new = nbcut - 1; 
	
	for(int i = debut+1; i < fin+1; i++)
    {
		mean += histo(i-1) * (double)(i)/histo.width();
		totalh += histo(i-1);
    }
	mean /= totalh;
	if(totalh == 0) mean = 0;
	//printf("MEAN : %g\n",mean);
	
	for(int i = debut+1; i < fin+1; i++)
    {
		double tmp = ((histo(i-1)*(double)(i)/histo.width()) - mean);
		sigmac += tmp*tmp*histo(i-1);
    }
	//printf("SIGMA SQUARE : %g\n",sigmac);
    
	double sigma1[10];
	double sigma2[10];
	for(int k = debut; k <= fin;k++)
    { 
		//printf("\nSETTING t in: %d\n",k);
		sigma1[k] = roundl(computeSP(histo,1,k,debut,fin));
		sigma2[k] = roundl(computeSP(histo,2,k,debut,fin));
		
		double tmp3 = sigma2[k]+sigma1[k];
		//printf("Variance part %d = %g\n",1,sigma1[k]);
		//printf("Variance part %d = %g\n",2,sigma2[k]);
		//printf("Total Variance = %g\n",tmp3);
		
		d = sigmac - tmp3;
		//printf("Distance = %g\n",d);
        
		if(d > max )
        {
			max = d;
			maxt = k;
			//printf("New perfect cut: %d\n",maxt);
        }
    }
	
	//printf("Perfect cut: %d\n",maxt);
	//printf("Stocking Sigma\n");
	stock[debut] = sigma1[maxt];
	stock[maxt]  = sigma2[maxt];
	
	
	//printf("Adding : %d\n",maxt);
	
	// Find the 1-Column sigma
	for(int i= 1 ;i <= histo.width() ;i++)
    {
		if(i == histo.width() && stock[i] != -1e99)
			stock[i] = 0;
    }       
	
	//Reading cut stock
	// for(int i = 0; i < 10; i++)
	//   {
	//     if( stock[i] != -1e99 )
	// 	printf("%d sigma: %g\n",i,stock[i]);
	//   }
	
	//Get back sigma from board to see the highest
	double maxcut = -1e99;
	int cut_key = debut;
	
	for(int i= 0 ;i < histo.width();i++)
    {
		if(stock[i] > maxcut) 
		{
			maxcut = stock[i];
			cut_key = i;
		}
    }
	//printf("left cut: %d sigma: %g\n",cut_key,stock[cut_key]);
	
	int cut_limite = histo.width();
	
	//Get the next cut
	int i = cut_key+1;
	
	while (( i <= cut_limite ) && ( stock[i] == -1e99 ))
		++i;
	
	cut_limite = i;
	
	if ( cut_limite > histo.width() )
		cut_limite = histo.width();
	
	//printf("right cut: %d sigma: %g\n",cut_limite,stock[cut_limite]);
	//printf("cut(s) left %d\n",nbcut_new);
	
	if ( nbcut_new > 0 )
		separate(histo,cut_key,cut_limite,stock,nbcut_new);    
}

deque<trail> ShapeBagTrails::getGraphTrails ( Collection *col, int i, int maxsize )
{
	Graph3d *grp = 0;
	
	if ( useTree )
		grp = col->GETPOBJECT("tree",Graph3d);
	else
		grp = col->GETPOBJECT("graph",Graph3d);
    
	deque<trail> res;
	
	if ( (*grp)[i] == 0 ) return res;
	
	trail p;
	set<int> edges;
	deque<GEdge *> stack;
	bool stop = false;
	bool backtrack = false;
	
	p.t.push_back(i);
	int node = i;
	
	// The smallest trail...
	res.push_back(p);
	
	while ( ! stop )
    {
		GEdge *edge = 0;
		
		if ( backtrack )
		{
			edge = stack.back();
			stack.pop_back();
			backtrack = false;
		}
		else
			edge = (*grp)[node]->Neighbours();
		
		for ( ; edge != 0; edge = edge->Next() )
		{
			if ( edges.find(edge->Item()) != edges.end() )
				continue;
			
			if ( strictTree && useTree && (edge->Node() != (int)(*grp)[edge->Node()]->Item()) )
				continue;
			
			p.t.push_back(-edge->Item());
			p.t.push_back(edge->Node());
			
			//cout << node << " -- " << -edge->Item() << " -- " << edge->Node() << endl;
			
			if ( (int)p.t.size() < 2*(maxsize+1) )
				res.push_back(p);
			else
			{
				p.t.pop_back();
				p.t.pop_back();
				continue;
			}
			
			edges.insert(edge->Item());
			stack.push_back(edge->Next());
			
			node = edge->Node();
			break;
		}
		
		if (( edge == 0 ) && ( node == i ))
			stop = true;
		else if ( edge == 0 )
		{
			//cout << "Backtrack" << endl;
			
			p.t.pop_back();
			edges.erase(-p.t.back());
			p.t.pop_back();
			
			backtrack = true;
			node = p.t.back();
		}
    }
	
	return res;
}

deque<trail> ShapeBagTrails::filterTrail ( deque<trail> & bag, Collection * col, int redundancy )
{
	// 1- Compute the weights
	Graph3d * grp = col->GETPOBJECT("tree",Graph3d);
	CImg<double> pweight (bag.size(), 2, 1, 1, 0);
	
	for ( unsigned int i = 0; i < bag.size(); ++i )
    {
		pweight(i,0) = (double)i;
		pweight(i,1) = bag[i].weight;
    }
	
	// 2- Order the trails using their weight
	for ( int i = 1; i < pweight.width(); ++i )
		for ( int j = 0; j < i+1; ++j )
		{
			if ( pweight(i,1) > pweight(j,1) )
			{
				int itemp = (int)pweight(j,0);
				double wtemp = pweight(j,1);
				
				pweight(j,0) = pweight(i,0);
				pweight(j,1) = pweight(i,1);
				
				for ( int k = j+1; k < i+1; ++k )
				{
					int itemp2 = (int)pweight(k,0);
					double wtemp2 = pweight(k,1);
					
					pweight(k,0) = itemp;
					pweight(k,1) = wtemp;
					
					itemp = itemp2;
					wtemp = wtemp2;
				}
				break;
			}
		}
	
	// 3- Take trails until there is a covering and only of they cover a new vertex
	deque<trail> newbag;
	int nbbrch = col->GETVALUE("nb_branches",Long);
	CImg<int> cover (nbbrch+2,1,1,1,0);
	for ( unsigned int i = 0; i < bag.size(); ++i )
    {
		bool newcover = false;
		int index = (int)pweight(i,0);
		
		// Look if the trail covers a new edge without too much redundancy
		for ( unsigned int j = 0; j < bag[index].t.size()-1; j += 2 )
		{
			int n_edge = 0;
			GEdge * edge = (*grp)[bag[index].t[j]]->Neighbours();
			for ( ; edge != 0; edge = edge->Next() )
				if ( edge->Node() == bag[index].t[j+2] )
				{
					n_edge = edge->Item();
					break;
				}
			
			if ( cover(n_edge) < redundancy )
				newcover = true;
			else if ( cover(n_edge) >= redundancy )
			{
				newcover = false;
				break;
			}
		}
		
		// Do things if we need to insert the trail
		if ( newcover )
		{
			// Look if the same trail is in the bag (remove false redundancy)
			for ( unsigned int x = 0; x < newbag.size(); ++x )
				if ((( newbag[x].t.front() == bag[index].t.front() )  &&
					 ( newbag[x].t.back()  == bag[index].t.back()  )) ||
					(( newbag[x].t.front() == bag[index].t.back()  )  &&
					 ( newbag[x].t.back()  == bag[index].t.front() )))
				{
					newcover = false;
					break;
				}
		}
		
		if ( newcover )
		{ 
			newbag.push_back(bag[index]);
			
			// Update the cover table
			for ( unsigned int j = 0; j < bag[index].t.size()-1; j += 2 )
			{
				int n_edge = 0;
				GEdge * edge = (*grp)[bag[index].t[j]]->Neighbours();
				for ( ; edge != 0; edge = edge->Next() )
					if ( edge->Node() == bag[index].t[j+2] )
					{
						n_edge = edge->Item();
						break;
					}
				
				cover(n_edge) += 1;
			}
		}
    }
	
	return newbag;
}

deque<trail> ShapeBagTrails::percentTrail ( deque<trail> & bag, unsigned int nbp )
{
	// 1- Compute the weights
	CImg<double> pweight (bag.size(), 2, 1, 1, 0);
	
	for ( unsigned int i = 0; i < bag.size(); ++i )
    {
		pweight(i,0) = (double)i;
		pweight(i,1) = bag[i].weight;
    }
	
	// 2- Order the trails
	for ( int i = 1; i < pweight.width(); ++i )
		for ( int j = 0; j < i+1; ++j )
		{
			if ( pweight(i,1) > pweight(j,1) )
			{
				double wtemp = pweight(j,1);
				int itemp = (int)pweight(j,0);
				
				pweight(j,0) = pweight(i,0);
				pweight(j,1) = pweight(i,1);
				
				for ( int k = j+1; k < i+1; ++k )
				{
					double wtemp2 = pweight(k,1);
					int itemp2 = (int)pweight(k,0);
					
					pweight(k,0) = itemp;
					pweight(k,1) = wtemp;
					
					itemp = itemp2;
					wtemp = wtemp2;
				}
				break;
			}
		}
	
	// 3- Take only the best
	deque<trail> newbag;
	
	for ( unsigned int i = 0; (i < nbp) && (i < bag.size()); ++i )
		newbag.push_back(bag[(int)pweight(i,0)]);
	
	return newbag;
}

deque<trail> ShapeBagTrails::histogram ( const deque<trail> & ebag, const int nbcut )
{
	unsigned int size = ebag.size();
	deque<trail> bag = ebag; // Copy of the bag
	
	// 1 - Compute the histogram
	CImg<int> histo ( 10, 1, 1, 1, 0 );
	
	for ( unsigned int i = 0; i < size; ++i )
		histo((int)(bag[i].weight*10.0)) += 1;
	
	// 2 - Show the histogram
	// histo.display_graph("Histogram",3);
	
	// 3 - Cut the histogram into several parts
	double *tab_cut_sigma = new double [12];
    
	for(int i = 0; i < 12; i++)
		tab_cut_sigma[i] = -1e99;
	
	// Cut it
	separate(histo,0,10,tab_cut_sigma,nbcut);
    
	// Reading cut stock
	deque<int> cuts;
	cuts.push_back(0);
	for ( int i = 0; i < 10; i++ )
		if ( tab_cut_sigma[i] != -1e99 && i != 0 )
			cuts.push_back(i);
	
	delete[] tab_cut_sigma;
	
	// Partition the bag
	int label = 1;
	for ( int i = (int)cuts.size() - 1; i >= 0; --i )
    {
		for ( unsigned int j = 0; j < bag.size(); ++j )
			if (( (int)(bag[j].weight*10.0) >= cuts[i] ) && bag[j].l.empty() )
				bag[j].l.push_back(label);
		++label;
    }
	
	return bag;
}

ShapeBagTrails::ShapeBagTrails ( int size, float percent, double covering, bool useTree, bool strictTree, int nb_cuts )
: trailSize(size),percent(percent),covering(covering),useTree(useTree),strictTree(strictTree),nb_cuts(nb_cuts)
{
}

deque<trail> ShapeBagTrails::operator() ( Collection * col )
{
	Long idx = col->GETVALUE("btrails_idx", Long);
	
	if (idx != -1)
		return bags[idx];
	
	deque<trail> bag;
	Graph3d * grp   = 0;
	double totalw   = col->GETVALUE("sum_weight",Float);
	
	if ( useTree )
		grp = col->GETPOBJECT("tree",Graph3d);
	else
		grp = col->GETPOBJECT("graph",Graph3d);
	
	Float *imWeight = col->GETARRAY("weights",Float);
	
	// Trail are compose of the sequence node,edge,node,edge,node,...
	
	for ( int i = 0; i < grp->Size(); ++i )
		if ( (*grp)[i] != 0 )
		{
			if ( strictTree && useTree && (i != (int)(*grp)[i]->Item()) )
				continue;
			
			deque<trail> tmp = getGraphTrails(col,i,trailSize);
			bag.insert(bag.end(),tmp.begin(),tmp.end());
		}
	
	// Compute the weight of each path (and the sum weight)
	
	for ( unsigned int i = 0; i < bag.size(); ++i )
    {
		double weight = 0.0;
		for ( unsigned int j = 1; j < bag[i].t.size(); j += 2 )
			weight += imWeight[-bag[i].t[j]];
		bag[i].weight = weight/totalw;
    }
	
	// Filter the paths
	
	if ( covering > 0.0 )
    {
		// Filter the trail to an acceptable solution
		bag = filterTrail(bag,col,6);
		
		// Compute the optimal covering
		VertexCovering cover (col,bag);
		bag = cover.computeCovering(true,covering);
		
		// Complete the bags with the symmetric trails
		vector<trail> nbag;
		
		for ( unsigned int i = 0; i < bag.size(); ++i )
		{
			trail p;
			p.l = bag[i].l;
			p.weight = bag[i].weight;
			p.t.insert(p.t.end(),bag[i].t.rbegin(),bag[i].t.rend());
			nbag.push_back(p);
		}
		
		bag.insert(bag.end(),nbag.begin(),nbag.end());
		
		// If useTree then relabel the trails
		if ( useTree )
		{
			for ( unsigned j = 0; j < bag.size(); ++j )
				for ( unsigned i = 0; i < bag[j].t.size(); i += 2 )
					bag[j].t[i] = (*grp)[bag[j].t[i]]->Item();
		}
    }
	else if ( percent < 1.0 )
    {
		// If useTree then relabel the trails
		if ( useTree )
		{
			for ( unsigned j = 0; j < bag.size(); ++j )
				for ( unsigned i = 0; i < bag[j].t.size(); i += 2 )
					bag[j].t[i] = (*grp)[bag[j].t[i]]->Item();
		}
		
		unsigned int nbp = (unsigned int)((float)bag.size() * percent);
		if ( nbp == 0 ) nbp = 1;
		bag = percentTrail(bag,nbp);
    }
	else if ( useTree )
    {
		for ( unsigned j = 0; j < bag.size(); ++j )
			for ( unsigned i = 0; i < bag[j].t.size(); i += 2 )
				bag[j].t[i] = (*grp)[bag[j].t[i]]->Item();
    }
	
	// Labelled the paths
	if ( nb_cuts == 0 )
    {
		// This is the default rule
		for ( unsigned int j = 0; j < bag.size(); ++j )
			bag[j].l.push_back(1);
    }
	else
		bag = histogram(bag,nb_cuts);
	
	// Compute the sum of each part (given by labels)
	map<int,double> sumWes;
	for ( unsigned int i = 0; i < bag.size(); ++i )
		for ( unsigned int l = 0; l < bag[i].l.size(); ++l )
		{
			if ( sumWes.find(bag[i].l[l]) == sumWes.end() )
			{
				sumWes[bag[i].l[l]] = bag[i].weight;
			}
			else
			{
				sumWes[bag[i].l[l]] += bag[i].weight;
			}
		}
	
	for ( unsigned int i = 0; i < bag.size(); ++i )
		for ( unsigned int l = 0; l < bag[i].l.size(); ++l )
			bag[i].sumWeight = sumWes[bag[i].l[l]];
	
	col->SETVALUE("btrails_idx", Long, (Long) bags.size());
	bags.push_back(bag);
	
	return bag;
}
