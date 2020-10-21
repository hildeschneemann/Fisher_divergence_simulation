// main() function: reads parameter values in input file,
// and runs the simulation.

#include "fisher.h"
#include <iostream>
#include "MersenneTwister.h"
using namespace std;

// input and output files:
FILE * fichierE;

// random number generator (Mersenne Twister):
MTRand rnd;

int main(int argc, char * argv[])
{
	
	//set seed
	rnd.seed(123);
	
	// definitions of variables:
	int N, n, reps, M, E, D;
	double L, U, k, s, F, q;
	
	// opens input and output files:
	bool fin;
	ouvrirFichierE(argv[1]);
	fin = false;
	do
	{
		// reads parameter values;
		fin = lireFichier( N, L, U, k,  s,  n, F, q, M, E, D, reps);
		if (!fin)
		{
			// runs the simulation:
			for(int i=0; i<reps; i++)
				recursion( N, L, U, k, s, n, F, q, M, E, D, i+1);
		}
	} while (!fin);

	// closes files:
	fclose(fichierE);
	return 0 ;
}
