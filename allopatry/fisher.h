// Header file: definitions of global variables, function prototypes
#ifndef FISHER_H
#define FISHER_H
#include <vector>
#include <iostream>
#include "MersenneTwister.h"
using namespace std;

// Global variables:
// "chr": represents a chromosome:
struct chr
{
	vector<int> sel; // selected loci (chain of 0 and 1)
};

// Function prototypes:
void ouvrirFichierE(char * param);
bool lireFichier(int & Nr, double & Lr, double & Ur,
				double & kr, double & sr, int & nr, double & Fr, double & qr, 
				int & Mr, int & Er, int & Dr, int & repsr);
void recursion(int Nv, double Lv, double Uv,
			double kv, double sv, int nv, double Fv, double qv, 
			int Mv, int Ev, int Dv, int repv);
double gammln(const double xx);
double poisdev(const double xm);
double gasdev();
double binldev(const double pp, const int n);
void rec(chr &res, chr &c1, chr &c2, double R, vector<vector<double> > &mutations, int n);
void freerec(chr &res, chr &c1, chr &c2);
double getfitness(double sz2, double kd2, double alpha);
double getdominance(double qv, double Fv);
double rexp(double mymean);
double factorial(double n);
double get_mutsize_MVN(int nv, double kv, double s_barv, double alpha);
double get_mutsize_EXP(double kv, double s_barv, double alpha);
void generatemutation(vector<double> & currentmutation, int nv, double mutsize, int Mv, double kv);
#endif