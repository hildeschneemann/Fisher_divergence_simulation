#include "fisher.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <bitset>
#include <string>
#include <random>
using namespace std;

extern MTRand rnd;

// rec function: generates recombinant chromosome "res" from parental chromosomes "c1" and "c2"
// nbCo is the number of cross-overs


//---------------------------SCHEMATIC OVERVIEW -----------------------------------
//
//	position on chromosome: 0.0     0.2     0.4     0.6     0.8     1.0
//							--------------------------------------------
//		                     interval 1  | interval 2 |    interval 3
//							--------------------------------------------
//			   cross-overs:              |            |
//				 mutations:    X             X            X         X
//							--------------------------------------------
//  					c1:    1         |   1        |   1         0     > even intervals passed on ( interval 2)
//						c2:    0         |   0        |   0         1     > odd intervals passed on (interval 1 and 3)
// 	  			 	   res:    0         |   1        |   0         1


// freerec function is the simpler free recombination version in which all mutations segregate independently 



void rec(chr &res, chr &c1, chr &c2, double R, vector<vector<double> > &mutations, int nv)
{
	//cout << "doing rec" << mutations[0][2 + nv] << "\n";
	//cout << "before starting rec script: " << mutations[0][2 + nv] << "\n";
	vector<double> pos; //store positions of cross-overs
	double k;
	res.sel.clear(); //this is where the mutations passed on to the offspring will go
	pos.clear();
	
    int nbCo = int(poisdev(R)); //number of cross-overs

	//First we need to figure out the relative position of the cross-overs and the mutations
  	for (int j = 0; j < nbCo; j++) //generate a position for each cross-over
	{
		k = rnd.rand();
		pos.push_back(k);
	}
	pos.push_back(0); //store positions of cross-overs in pos vector


	if (pos.size() > 1) //if any cross-overs have occured, record which mutations lie in odd (1) and which lie in even (0) intervals
	{
		sort(pos.begin(), pos.end());
		if (!(mutations.empty()))
		{
			for (unsigned int i = 0; i < mutations.size(); i++) //go through mutations
			{
				if (pos.size() > 1)
				{
					for (unsigned int j = 1; j < pos.size(); j++) //go through cross-over positions
					{
						if(mutations[i][0] > pos[j - 1] && mutations[i][0] < pos[j]) //if a mutation lies between two cross-overs...
						{
							cout << "before assigning odd/even: " << mutations[i][3 + nv] << "\t" << (j % 2) << "\n";
							cout << "mutations size: " << mutations[i].size() << "\n";
							cout << "n: " << nv << "\n";
							
							mutations[i][1 + nv] = (j % 2); //record whether it lies in an odd or even interval 
							cout << "mutations size v2: " << mutations[i].size() << "\n";
							
							cout << "after assigning odd/even: " << mutations[i][2 + nv] << "\n";
						}
					}
				}
				if(mutations[i][0] > pos.back()) //if position of the mutation is greater than the last cross-over...
				{
					cout << "before assigning odd/even: " << mutations[i][1 + nv] << "\t" << (pos.size() % 2) << "\n";
					mutations[i][1 + nv] = (pos.size() % 2); //record whether it lies in an odd or even interval
					cout << "after assigning odd/even: " << mutations[i][1 + nv] << "\n";
				}
			}
		}
		
		if (!(c1.sel.empty()))
		{
			for (unsigned int i = 0; i < c1.sel.size(); i++) //include mutations from parent's chromsome 1 (c1) on even intervals between cross-overs
			{
				int index = c1.sel[i];
				if (mutations[index][1 + nv] == 0)
					res.sel.push_back(index); //add these mutations to offspring's recombined chromosome
			}

		}

		if (!(c2.sel.empty()))
		{
			for (unsigned int i = 0; i < c2.sel.size(); i++) //include mutations from parent's chromsome 2 (c2) on odd intervals between cross-overs 
			{
				int index = c2.sel[i];
				if (mutations[index][1 + nv] == 1)
					res.sel.push_back(index); //add these mutations to offspring's recombined chromosome
			}
		}
	}//end if any cross-overs
	
	if (pos.size() == 1) //if there are no cross-overs pick one chromosome to pass on at random
	{
		if(rnd.rand() > 0.5) //pick chromosome 1 (c1) to pass on
		{
			if (!(c1.sel.empty()))
			{
				for (unsigned int i = 0; i < c1.sel.size(); i++) //copy any mutations present on parental chromosome 1 (c1) to offspring (res)
				{
					int index = c1.sel[i];
					res.sel.push_back(index);
				}
			}
		}
		else //pick chromosome 2 (c2) to pass on
		{
			if (!(c2.sel.empty()))
            {
				for (unsigned int i = 0; i < c2.sel.size(); i++) //copy any mutations present on parental chromosome 2 (c2) to offspring (res)
				{
					int index = c2.sel[i];
					res.sel.push_back(index);
				}
			}
		}
	}
}



void freerec(chr &res, chr &c1, chr &c2)
{
	res.sel.clear();
	vector<int> hetero; //store which mutations are heterozygous in parent
	//include homozygous sites first
	set_intersection(c1.sel.begin(), c1.sel.end(), c2.sel.begin(), c2.sel.end(), back_inserter(res.sel)); //add homozygous mutations to offspring's chromosome
	set_symmetric_difference(c1.sel.begin(), c1.sel.end(), c2.sel.begin(), c2.sel.end(), back_inserter(hetero)); //store indices of heterozygous mutations
	if (!(hetero.empty()))
	{
		for (unsigned int i = 0; i < hetero.size(); i++) // loop through heterozygous mutations
		{
			if(rnd.rand() > 0.5) // heterozygous mutations are passed on with probability 0.5
				res.sel.push_back(hetero[i]);
		}
	}
}
