#include "fisher.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
//#include <boost/math/distributions.hpp>
#include <cmath>
#include <vector>
#include <algorithm>
#include <chrono>
using namespace std;

// Constants:
#define PIE 3.14159265359
#define EXITFLAG -1

// Hardcoded variables (need to be changed by hand and code recompiled):
#define ALPHA 0.5 // overall strength of selection
#define CYCLELENGTH 100000 // how many generations does it take for the moving optimum to complete a cycle (Ev=1)
#define UPDATELENGTH 100 // Report to the screen after this many generations

extern MTRand rnd;

/*****************************************************************************************************************************************
 *****************************************************************************************************************************************
 *****************************************************************************************************************************************
 ******** FUNCTION DEFINITIONS
 *****************************************************************************************************************************************
 *****************************************************************************************************************************************
 *****************************************************************************************************************************************/

/*************************************************************************
 ******** 	Return the position of the optimum for population "pop"
 ******** 	Ev=0: optimum fixed at 0 (the ancestral trait value)
 ******** 	Ev=-1: optimum fixed at 1 for population 1 and -1 for population 2
 ******** 	Ev=1: oscillating pattern over "CYCLELENGTH" generations
 ******** 	Ev=2: optimum fixed at 0 (ancestral state) for population 1, and fixed at 2 for population 2
 *************************************************************************/
double get_opt(int pop, int gen, int Ev, int cyclength, double kv)
{
	double tmp;
	
	if(Ev==0) // Fixed optimum at the origin
	{
		return(0.0);
	} 
	if(Ev==-1) // Optimum fixed at a unit length in opposite direction for population 1 and 2
	{
		return(pop == 1 ? 1.0 : -1.0);	    
	} 
	if(Ev==2) // Optimum fixed at 2 x unit length for population 1 and at ancestral state for population 2
	{    
		return(pop == 1 ? 2.0 : 0.0);
	}
	if(Ev==1)
	{    
		tmp = 0.5 + sin(2.0*PIE*gen/cyclength-PIE/2)*0.5;
		if(tmp>0.0)
			return(pow(tmp,2.0/kv));
		if(tmp<0.0)
			return(-pow(-tmp,2.0/kv));
		if(tmp==0.0)
			return(0.0);
	}
	exit(1);
}

/*****************************************************************************************************************************************
 *****************************************************************************************************************************************
 *****************************************************************************************************************************************
 ******** The main recursion
 *****************************************************************************************************************************************
 *****************************************************************************************************************************************
 *****************************************************************************************************************************************

Individual-based simulation of evolutionary divergence for two populations with gene flow between them under Fisher's geometrical model,
in a range of population genetic regimes

Input parameters:

** Population genetics parameters
Nv: size of the two populations combined; Total number of diploid individuals
mv: probability of migration
Lv: genome map length in Morgans (mean nb of cross-overs per meiosis); set to -1 for free recombination
Uv: mutation rate per diploid genome

** Fitness landscape
kv: curvature of fitness function
nv: number of phenotypic traits

** Distribution of mutation effects
Mv: Flag for mutational effect model: 1 "bottom up"; 0 "top down"
sv: The mean selection coefficient for random mutations in an optimal genotype
qv: mean phenotypic dominance effect (set to 1/2 in published work)
Fv: inflation of variance of distribution of phenotypic dominance (set to -1 for a uniform distribution)

** Demographic/environmental change
Ev: Flag for environmental change (i.e. position of the optimum): 
	-1  optimum fixed at 1 for population 1 and -1 for population 2
	0  (optimum is fixed to match ancestral state) 
	1  (oscillating over time) 
	2  optimum fixed at 0 (ancestral state) for population 1, and fixed at 2 for population 2

** Simulation parameters
Dv: Stop simulation after Dv fixations have occurred.

 *****************************************************************************************************************************************
 *****************************************************************************************************************************************
 *****************************************************************************************************************************************/
void recursion( int Nv, double mv, double Lv, double Uv,
                double kv, double sv, int nv, 
				double Fv, double qv, int Mv, 
				int Ev, int Dv, int repv)
{
	// declare variables:
	int i, j, k, loc, par1, par2, ind, nb, mut, indmut, index, migN1, migN2, migN2_lim, chrom;
	double w, d, sz2, mutsize, pos, Wmax;
	bool withrec;
	long long int maxgen;
   
	long long int gen = 0;
	int nb_seg = 0; //number of segregating mutations
	int nb_fix = 0;	//number of fixed variants

	
	double kd2 = kv / 2.0; // half of the fitness landscape curvature
	withrec = Lv == -1 ? false : true; // Flag for free recombination
	bool mono = true; //is the population monomorphic? skip to mutation step

	// compute various quantities
	int N_1 = Nv - 1;
	int twoN = 2*Nv; //total number of chromosomes
	int halfN = Nv/2; //number of individuals per population
	int halfN_1 = halfN - 1;
	double optval1 = get_opt(1,0,Ev,CYCLELENGTH,kv); //distance to optimum for population 1 at generation 0
	double optval2 = get_opt(2,0,Ev,CYCLELENGTH,kv); //distance to optimum for population 2 at generation 0
	
	// set up population of diploid invididuals
	chr * pop = new chr [twoN]; //store all genotypes in the population
	chr * temp = new chr [twoN]; // used to generate next generation
	chr * cp;

	// "Wtot" will hold the fitness of each individual:
	double * Wtot = new double [Nv];

	//set up lists to store information on mutations
	vector<double> mutline((3 + nv), 0.0); //store mutational effects of variants for each trait. elements are: position, mutational effects, allele count
	vector<vector <double> > mutations(0, mutline);
	vector<double> domline(nv, 0.0); 
	vector<vector <double> > delta(0, domline); //store dominance coefficients of variants for each trait
	vector<int> muttimes; //store generation at which variant arises
	vector<int> fixtimes; //store generation at which variant fixes
	vector<int> fixed; //store indices of fixed variants
	vector<int> lost; //store indices of lost variants
	vector<int> homo; //store indices of homozygous variants
	vector<int> hetero; // store indices of heterozygous variants
	vector<double>currentmutation(nv,0.0);
	vector<double>wildtype(nv,0.0);

	//reserve capacity for vectors that will grow
	delta.reserve(100 * nv);
	mutations.reserve(100 * nv);

	//record duration of simulation
	time_t start, finish;
	struct tm *ptr;
	start = time(0);
	
	
	//set up output file for fixed variants
	string fileName;
	stringstream nameF;
	nameF << "res_N" << Nv << "_m" << mv <<
					"_L" << Lv << "_U" << Uv << 
					"_k" << kv << "_s" << sv << "_n" << nv << "_F" << Fv << "_q" << qv << "_M" << Mv << 
					"_E" << Ev << "_rep" << repv << ".txt";
	nameF >> fileName;
	
	//set up output file for frequency of segregating variants
	string fileNamefreq;
    stringstream nameFreq;
	nameFreq << "freq_N" << Nv << "_m" << mv <<
					"_L" << Lv << "_U" << Uv << 
					"_k" << kv << "_s" << sv << "_n" << nv << "_F" << Fv << "_q" << qv << "_M" << Mv << 
					"_E" << Ev << "_rep" << repv << ".txt";
    nameFreq >> fileNamefreq;

	//set up output file with individual genotypes at the end of the simulation
	string fileNamepop;
	stringstream nameFpop;
	nameFpop << "pop_N" << Nv << "_m" << mv <<
					"_L" << Lv << "_U" << Uv << 
					"_k" << kv << "_s" << sv << "_n" << nv << "_F" << Fv << "_q" << qv << "_M" << Mv << 
					"_E" << Ev << "_rep" << repv << ".txt";
    nameFpop >> fileNamepop;

	
	//make header output file
	ofstream fout;
	ofstream fout2;
	ofstream fout3;
	
	fout.open(fileName);
        fout << "position\tmut_time\tfix_time\topt_dist\t";
        for (i = 0; i < nv; i++)
        {
                fout << "mut_effect" << i << "\t" << "dom" << i << "\t";
        }
	fout << endl;
	
	fout2.open(fileNamefreq);
	fout2 << "gen\tindex\tposition\tfreq_pop1\tfreq_pop2\tmut_time\tfix_time\topt_dist1\topt_dist2\t";
    for (i = 0; i < nv; i++)
    {
    	fout2 << "mut_effect" << i << "\t" << "dom" << i << "\t";
    }
    fout2 << endl;
	
	// set parameter for the distribution of mutation effect sizes		
	if(Mv) //Multivariate normal
		mutsize = get_mutsize_MVN(nv, kv, sv, ALPHA);		
	else //Exponential
		mutsize = get_mutsize_EXP(2.0, sv, ALPHA); // r_bar
	
	
	//   ---------------------------   START RECURSION    ---------------------------   
	// loop through generations:
	while(gen != EXITFLAG)
	{
		gen += 1;

		//   (1) ---------------------------   COMPUTE FITNESS    --------------------------- 
		if (mono == 0) //if the population is monomorphic, skip to mutation step. if not: compute fitnesses
		{
			Wmax = 0; //maximum fitness in the population
			optval1 = get_opt(1,((gen) % CYCLELENGTH),Ev,CYCLELENGTH,kv); //distance to optimum for population 1
			optval2 = get_opt(2,((gen) % CYCLELENGTH),Ev,CYCLELENGTH,kv); //distance to optimum for population 2
			
			for (j = 0; j < Nv; j++)  // for each individual
			{
				nb = 2 * j; //(nb = index of 1st chromosome individual j, nb+1 = index of 2nd chromosome individual j)
				sz2 = 0; //summed squared distance to optimum
				d = 0; //distance to optimum
				
				sort(pop[nb].sel.begin(), pop[nb].sel.end());
				sort(pop[nb + 1].sel.begin(), pop[nb + 1].sel.end());

				//remove fixed mutations from genotype
				for (i = 0; i < fixed.size(); i++)
				{
					if (find(pop[nb].sel.begin(), pop[nb].sel.end(), fixed[i]) == pop[nb].sel.end())
						cout << "'fixed' mut " << fixed[i] << " absent in chrom 1 ind " << j << endl;
					if (find(pop[nb + 1].sel.begin(), pop[nb + 1].sel.end(), fixed[i]) == pop[nb + 1].sel.end())
						cout << "'fixed' mut " << fixed[i] << " absent in chrom 2 ind " << j << endl;
					pop[nb].sel.erase(remove(pop[nb].sel.begin(), pop[nb].sel.end(), fixed[i]), pop[nb].sel.end());
					pop[nb + 1].sel.erase(remove(pop[nb + 1].sel.begin(), pop[nb + 1].sel.end(), fixed[i]), pop[nb + 1].sel.end());
				}

				//determine homozygous and heterozygous sites
				set_intersection(pop[nb].sel.begin(), pop[nb].sel.end(),pop[nb + 1].sel.begin(), pop[nb + 1].sel.end(), back_inserter(homo));
				set_symmetric_difference(pop[nb].sel.begin(), pop[nb].sel.end(),pop[nb + 1].sel.begin(), pop[nb + 1].sel.end(), back_inserter(hetero));
					
				//compute distance to optimum for each trait
				for (k = 0; k < nv; k++) //loop through n traits
				{
					d = wildtype[k];
					if (k == 0) //for 1st trait, consider current position of optimum
						if(j < halfN)
							d -= optval1;
						else if(j >= halfN)
							d -= optval2;
	
					for (i = 0; i < homo.size(); i++) //loop through homozygous sites
					{
						index = homo[i];
						d += mutations[index][k + 1];
					}
			
					for (i = 0; i < hetero.size(); i++) //loop through heterozygous sites
					{
						index = hetero[i];
						d += mutations[index][k + 1]*delta[index][k];
					}

					sz2 += d * d;
				}
				homo.clear();
				hetero.clear();
				
	
				// compute fitness:
				w = getfitness(sz2, kd2, ALPHA);
				Wtot[j] = w;
				if (Wmax < w)
					Wmax = w;
			} //end for loop going through indvidiuals
		
			// **** bookkeeping **** 
			//clear indices of fixed variants
			if (!(fixed.empty()))
			{
				for (i = 0; i < fixed.size(); i++)
					lost.push_back(fixed[i]);
				fixed.clear();
			}
			// **** end bookkeeping **** 
			
			//   (2) ---------------------------   SAMPLE NEXT GENERATION    --------------------------- 
			migN1 = int(binldev(mv, halfN)); // number of offspring in population 1 with parents from population 2
			migN2 = int(binldev(mv, halfN)); // number of offspring in population 1 with parents from population 2
			migN2_lim = migN2 + halfN;
			for (ind = 0; ind < Nv; ind++) //loop through offspring
			{ 
				// sampling first parent:
				do
				{
						par1 = int(rnd.randInt(halfN_1)); //pick random index for candidate parent 1
						
						//parents from offspring 0 to migN1 come from population 2 (migrant)
						//parents from offspring migN1 to halfN come from population 1 (non-migrant)
						//parents from offspring halfN to halfN+migN2 come from population 1 (migrant)
						//parents from offspring halfN+migN2 to N come from population 2 (non-migrant)
						
						if(ind < migN1 || ind > migN2_lim) // add halfN to candidate index if parent 1 comes from population 2
							par1 += halfN;
						
				} while (rnd.rand()> (Wtot[par1] / Wmax)); //accept candidate with probabiliy proportional to their fitness
				
				// recombination
				if (withrec)
								rec(temp[2*ind], pop[2*par1], pop[2*par1+1], Lv, mutations, nv);
				else
								freerec(temp[2*ind], pop[2*par1], pop[2*par1+1]);

				//sample second parent
				do
				{
					par2 = int(rnd.randInt(halfN_1)); //pick random index for candidate parent 2
					if(ind < migN1 || ind > migN2_lim)
							par2 += halfN;
				} while (rnd.rand()> (Wtot[par2] / Wmax)); //accept candidate with probabiliy proportional to their fitness
				
				// recombination
				if (withrec)
								rec(temp[2*ind+1], pop[2*par2], pop[2*par2+1], Lv, mutations, nv);
				else
								freerec(temp[2*ind+1], pop[2*par2], pop[2*par2+1]);
			}

		} //skip straight to mutation in case of monomorphic population (no need to recompute fitnesses)


		
		//   (3) ---------------------------   MUTATION    --------------------------- 
		mut = int(poisdev(2 * Uv * Nv)); // number of new mutations in entire population

		for (j = 0; j < mut; j++) //loop through new mutations to distribute over individuals
		{
			indmut = int(rnd.randInt(twoN - 1)); //pick random chromosome to receive mutation
			pos = rnd.rand(); //pick random position in the genome [0,1] ie infinite site model
			
			// info new mutation
			mutline[0] = pos;
			generatemutation(currentmutation,nv,mutsize,Mv,kv);
			for (k = 0; k < nv; k++) // generate mutational effect and dominance of new mutation for each trait
			{
				mutline[k + 1] = currentmutation[k];
				domline[k] = getdominance(qv, Fv);
			}
			mutline[1 + nv] = 0; //this is where the mutation's frequency in population 1 will go.
			mutline[2 + nv] = 0; //this is where the mutation's frequency in population 2 will go.
			
			if (lost.empty()) //generate new mutation index if none available to recycle
			{
				index = mutations.size();
				delta.push_back(domline);
				mutations.push_back(mutline); //add mutation to list of segregating mutations
				muttimes.push_back(gen); //store generation in which mutation arose
				fixtimes.push_back(0);
			}
			
			//to keep the keep the size of the mutations list managable, we re-use the position in the mutation list once a mutation is lost
			if (!(lost.empty())) //retrieve unused index if available
			{
				if (lost.size() > 1)
					sort(lost.begin(), lost.end());
				index = lost[0];
				lost.erase(remove(lost.begin(), lost.end(), index), lost.end()); //this index is now no longer available (until mutation is lost)
				mutations[index] = mutline; //add mutation to list of segregating mutations
				delta[index] = domline;
				muttimes[index] = gen; //store generation in which mutation arose
				fixtimes[index] = 0;
			}
			temp[indmut].sel.push_back(index); // add mutation to individual genotype
		}


		// **** bookkeeping **** 

		// update population:
		cp = pop;
		pop = temp;
		temp = cp;
		
		// **** end bookkeeping **** 
		
		
		
		
		
		//   (4) ---------------------------   COMPUTE ALLELE FREQUENCIES    --------------------------- 
		if ((mut > 0) || (mono==false)) //if no mutations are segregating (monomorphic population), no need to compute allele frequencies
		{
			mono = true;
			nb_seg = 0;

			//first set all allele counts to 0
			for (loc = 0; loc < mutations.size(); loc++)
			{
				mutations[loc][1 + nv] = 0; //allele count in population 1
				mutations[loc][2 + nv] = 0; //allele count in population 2
			}

			for (ind = 0; ind < Nv; ind++) //count up alleles of all individuals in population 1
			{	
				if (!(pop[ind].sel.empty()))
				{
					mono = false;
					for (loc = 0; loc < pop[ind].sel.size(); loc++)
					{	
						index = pop[ind].sel[loc];
						mutations[index][1 + nv] += 1;
					}
				}
			}

			for (ind = Nv; ind < twoN; ind++) //count up alleles of all individuals in population 2
			{	
				if (!(pop[ind].sel.empty()))
				{
					mono = false;
					for (loc = 0; loc < pop[ind].sel.size(); loc++)
					{	
						index = pop[ind].sel[loc];
						mutations[index][2 + nv] += 1;
					}
				}
			}
			
			// record fixed, lost and segregating variants and update wildtype
			for (loc = 0; loc < mutations.size(); loc++)
			{
				if ((mutations[loc][1 + nv] / Nv) == 1 && (mutations[loc][2 + nv] / Nv) == 1) // fixed variants
				{
					fixed.push_back(loc); //add recently fixed mutations to list of fixed mutations
					if( fixtimes[loc] == 0)
						fixtimes[loc] = gen; //record generation of fixation
					fout << mutations[loc][0] << "\t" << muttimes[loc] << "\t" << fixtimes[loc] << "\t" << optval1 << "\t" << optval2 << "\t"; //report fixation to output file
					for (k = 0; k < nv; k++)
					{
						wildtype[k] += mutations[loc][k + 1]; //update the wildtype
						fout << mutations[loc][k + 1] << "\t" << delta[loc][k] << "\t";
					}
					fout << endl;
					nb_fix += 1;
				}
				else if (mutations[loc][1 + nv] == 0 && mutations[loc][2 + nv] == 0) //lost variants
				{
					if (find(lost.begin(), lost.end(), loc) == lost.end()) //if locus was not yet on the lost list
						lost.push_back(loc);
				}
				else
					nb_seg += 1;
			}
		}
				
		
		if (gen % UPDATELENGTH == 0 || gen == EXITFLAG)    //display progress of simulation           
		{
			fprintf(stderr,"After this cycle, at gen %lld, we have fixed %d with %d segregating.\n",gen,nb_fix, nb_seg);
			for (loc = 0; loc < mutations.size(); loc++)
            {
				if((double(mutations[loc][1 + nv]) / double(Nv)) > 0 || (double(mutations[loc][2 + nv]) / double(Nv)) > 0) //no need to report lost mutants
				{
					fout2 << gen << "\t" << loc << "\t" << mutations[loc][0] << "\t" << (double(mutations[loc][1 + nv]) / double(Nv)) << "\t" << (double(mutations[loc][2 + nv]) / double(Nv)) << "\t" << muttimes[loc] << "\t" << fixtimes[loc] << "\t" << optval1 << "\t" << optval2 << "\t";
					for (k = 0; k < nv; k++)
					{
						fout2 << mutations[loc][k + 1] << "\t" << delta[loc][k] << "\t";
					}
					fout2 << endl;
				}
			}
		}		
		
		
		// **** bookkeeping **** 
		
		// stop simulation once D mutations fixed
		if (nb_fix >= Dv)
		{
			maxgen = gen;
			gen = EXITFLAG;
		}
			

		// **** end bookkeeping **** 
		
			
	} // end of recursion for loop
	//   ---------------------------   END RECURSION    --------------------------- 

	//record individuals from both populations
	fout3.open(fileNamepop);
	
	//make header genotype file
    fout3 << "ind\tpopulation\t";
    for (i = 0; i < mutations.size(); i++)
    {
        fout3 << "mut_" << i << "\t";
    }
    fout3 << endl;
	
	//record genotypes of individuals in population 1
	vector<int>genotype(mutations.size(),0); //declare vector to store genotype
	for (ind = 0; ind < halfN; ind++)
	{
		fill(genotype.begin(), genotype.end(), 0);
		fout3 << ind << "\t" << 1 << "\t"; //record individual and population number
		for(chrom = 0; chrom < 2; chrom++) //go through their two chromosomes
		{
			nb = 2*ind+chrom;
			if (!(pop[nb].sel.empty()))
			{
				for (loc = 0; loc < pop[nb].sel.size(); loc++)
				{
					index = pop[nb].sel[loc];
					genotype[index] += 1;
				}
			}
		}
		for(i = 0; i < mutations.size(); i++) //record genotype
			fout3 << genotype[i] << "\t";
		fout3 << endl;
	}

	for (ind = halfN; ind < Nv; ind++) //record individuals population 2
	{
		fill(genotype.begin(), genotype.end(), 0);
		fout3 << ind << "\t" << 2 << "\t"; //record individual and population number
		for(chrom = 0; chrom < 2; chrom++) //go through their two chromosomes
		{
			nb = 2*ind+chrom;
			if (!(pop[nb].sel.empty()))
			{
				for (loc = 0; loc < pop[nb].sel.size(); loc++)
				{
					index = pop[nb].sel[loc];
					genotype[index] += 1;
				}
			}
		}
		for(i = 0; i < mutations.size(); i++) //record genotype
			fout3 << genotype[i] << "\t";
		fout3 << endl;
	}

	//close output files
	fout.close();
    fout2.close();
	fout3.close();
	
	//record duration of simulation
	finish = time(0);
	int temps = int(difftime(finish, start));
	cout << maxgen << " generations have taken " << temps << " seconds\n";

	delete [] temp;
	delete [] Wtot;
	delete [] pop;

}
