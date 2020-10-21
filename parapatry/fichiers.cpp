// Functions to open input and output files,
// read parameter values from input file and 
// write them in output file.

#include "fisher.h"
#include <iostream>
#include <fstream>
using namespace std;

extern FILE * fichierE;
extern FILE * fichierM;

// opens input file:
void ouvrirFichierE(char * param)    
{						 
	fichierE = fopen(param,"r");
}


// reads parameter values from input file,
// returns 1 if end of input file, else returns 0
bool lireFichier(int & Nr, double & mr, 
				double & Lr, double & Ur,
				double & kr, double & sr, int & nr, double & Fr, double & qr, int & Mv, int & Ev, int & Dv,
				int & repsr)
{					 
	int x;
	bool term;
	do {x = fgetc(fichierE);} while (!((x == '*') || (x == EOF)));
		// each parameter set must start with *
	if (x == EOF)
		term = true;
	else
	{
        fscanf(fichierE,"%d ",&Nr);
		fscanf(fichierE,"%lf ",&mr);
		fscanf(fichierE,"%lf ",&Lr);
		fscanf(fichierE,"%lf ",&Ur);	
		fscanf(fichierE,"%lf ",&kr);
		fscanf(fichierE,"%lf ",&sr);
		fscanf(fichierE,"%d ",&nr);
		fscanf(fichierE,"%lf ",&Fr);
		fscanf(fichierE,"%lf ",&qr);		
		fscanf(fichierE,"%d ",&Mv);
		fscanf(fichierE,"%d ",&Ev);
		fscanf(fichierE,"%d ",&Dv);
		fscanf(fichierE,"%d ",&repsr);
		
		term = false;
	}
	return term;
}
