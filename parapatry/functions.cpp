#include "fisher.h"
#include <iostream>

using namespace std;

extern MTRand rnd;


/*****************************************************************************************************************************************
 *****************************************************************************************************************************************
 *****************************************************************************************************************************************
 ******** FUNCTION DEFINITIONS
 *****************************************************************************************************************************************
 *****************************************************************************************************************************************
 *****************************************************************************************************************************************/

/*************************************************************************
 ******** 	Isotropic fitness function:
 ******** 	sz2: summed squared distance from optimum
 ******** 	kd2: k/2 half of the curvature of the fitness function
 ******** 	alpha: tunes the overall strength of selection
 *************************************************************************/
double getfitness(double sz2, double kd2, double alpha)
{
        return(exp(-alpha * pow(sz2,kd2))); 
}

/*************************************************************************
 ******** 	Generate a random dominance coefficient for a new mutation:
 ******** 	Assumes a beta distribution with mean qv and variance qv(1-qv)Fv
 ******** 	Fv=-1 is a quicker way of implementing a uniform distribution
 *************************************************************************/
double getdominance(double qv, double Fv)
{
    if (Fv == 0) //without variation in dominance, dominance coefficient is equal to mean qv
    	return(qv);

    if (Fv == 1) //mutations are either completely dominant or completely recessive
    	return(rnd.rand() < qv ? 1 : 0);

    if (Fv == -1) //dominance coefficients are uniformly distributed
    	return(rnd.rand());
	if(Fv > 0 & Fv < 1) //draw dominance coefficient from beta distribution
	{
    	double alpha = qv * (1.0-Fv)/Fv;
        double beta = (1.0-Fv)*(1.0-qv)/Fv;
        //boost::math::beta_distribution<double> betadist(alpha, beta);
        //return(quantile(betadist, rnd.rand()));
    }
	exit(1);
}

/*************************************************************************
 ******** 	Exponentially distributed random number with mean "mymean"
 *************************************************************************/
double rexp(double mymean)
{
        return(-mymean*log(rnd.rand()));
}

/*************************************************************************
 ******** 	Factorials
 *************************************************************************/
double factorial(double n)
{
        return (n == 1.0 || n == 0.0) ? 1 : factorial(n - 1.0) * n;

}

/*************************************************************************
 ******** 	Calculate the standard deviation for the 
 ******** 	"bottom-up" mutation model (i.i.d normal on all traits)
 ******** 	This is set such that the mean selection coefficient of new 
 ******** 	mutations in an optimal genotype is s_bar.
 *************************************************************************/
double get_mutsize_MVN(int nv, double kv, double s_barv, double alpha)
{
	double sigv;
	sigv = s_barv/alpha;
	sigv *= exp(lgamma((double)nv/2.0) - lgamma(((double)nv+kv)/2.0));
	sigv = pow (sigv, 2.0/kv);
	sigv = sqrt(sigv/2.0);
	return(sigv);
}

/*************************************************************************
 ******** 	Calculate the parameter for the exponential distribution 
 ******** 	used in the "top-down" mutation model
 ******** 	This is set such that the mean selection coefficient of new
 ******** 	mutations in an optimal genotype is s_bar.
 *************************************************************************/
double get_mutsize_EXP(double kv, double s_barv, double alpha)
{
	double lamv;
	lamv = s_barv / alpha; 
	lamv /= factorial(kv);
	lamv = pow (lamv, 1.0/kv);
	return(lamv);
}

/*************************************************************************
 ******** 	Generate the n changes in trait values for a new mutation with
 ******** 	a random direction. Mv==1 for "bottom up", Mv=0 for "top down" 
 *************************************************************************/
void generatemutation(vector<double> & currentmutation, int nv, double mutsize, int Mv, double kv)
{
	if(Mv) //Multivariate normal
	{
		for(int i = 0; i<nv; i++)
			currentmutation[i] = mutsize * gasdev();
	} else { //Transformed exponential
		double r = 0.0;
		double curr_mutsize;
		for(int i = 0; i<nv; i++)
		{
			currentmutation[i] = gasdev();
			r += currentmutation[i]*currentmutation[i];
		}
		r = sqrt(r);
		// Do the transformation to fix the distribution of s for all values of k
		curr_mutsize = pow(rexp(mutsize),2.0/kv); 
		curr_mutsize /= r; 		
		
		for(int i = 0; i<nv; i++)
			currentmutation[i] *= curr_mutsize;
	}	
}