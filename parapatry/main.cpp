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
	rnd.seed();
	
	// definitions of variables:
		
	int N, n, reps, M, E, D;
	double L, U, k, s, F, q, m;
	
	// opens input and output files:

	bool fin;
	
	ouvrirFichierE(argv[1]);
	

	fin = false;

	do
	{
		// reads parameter values;
		fin = lireFichier( N, m, L, U, k, s, n, F, q, M, E, D, reps);
		if (!fin)
		{
			// runs the simulation:
			fprintf(stderr,"\tN=%d m=%f U=%f k=%f s=%f n=%d MVN=%d E=%d D=%d\n",N,m,U,k,s,n,M,E,D);
			for(int i=0; i<reps; i++)
				recursion( N, m, L, U, k,  s,  n, F, q, M, E, D, i+1);
		}
	} while (!fin);

	// closes files:
	fclose(fichierE);
	return 0 ;
}
