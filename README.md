# Fisher_divergence_simulation
This repository contains the scripts to perform individual-based simulations of divergence under Fisher's geometric model
############################################################################################

Individual-based simulation of evolutionary divergence under Fisher's geometrical model,
in a range of population genetic regimes

############################################################################################

For more details see: "The geometry and genetics of hybridization" , Appendix 2.
Schneemann H, Sanctis BD, Roze D, Bierne N, Welch JJ. 2020. The geometry and genetics of hybridization. Evolution (doi:10.1101/862235) 

For questions, email: hilde.schneemann@evobio.eu

############################################################################################

There are two versions of the code to simulate divergence in allopatry and parapatry respectively.
The allopatry version simulates a single population. 
The parapatry version simulates two populations with gene flow between them. 

############################################################################################
########                 COMPILATION                                                ########
############################################################################################


Before compilation, the boost library must be installed. (https://www.boost.org/)
To compile, make sure all *.cpp files and header files (*.h) are present in the same directory and run:

g++ *.cpp -I/usr/local/include/boost

############################################################################################
########                 RUNTIME                                                    ########
############################################################################################

To run:

./a.out inputfile.txt


############################################################################################
########                 INPUTFILES                                                 ########
############################################################################################

The inputfile needs to specify the following parameters (see example inputfile.txt)

Input parameters:

** Population genetics parameters
N: size of the population (post-bottleneck if implemented); Number of diploid individuals
m: probability of migration (NB only in case of parapatry)
L: genome map length in Morgans (mean nb of cross-overs per meiosis); set to -1 for free recombination
U: mutation rate per diploid genome

** Fitness landscape
k: curvature of fitness function
n: number of phenotypic traits

** Distribution of mutation effects
M: Flag for mutational effect model: 1 "bottom up"; 0 "top down"
s: The mean selection coefficient for random mutations in an optimal genotype
q: mean phenotypic dominance effect (set to 1/2 in published work)
F: inflation of variance of distribution of phenotypic dominance (set to -1 for a uniform distribution)

** Demographic/environmental change
E: Flag for environmental change (i.e. position of the optimum): 
	-1 (optimum is fixed to match ancestral state + initial bottleneck) 
	0  (optimum is fixed to match ancestral state) 
	1  (oscillating over time) 
	2  (fixed at 1, away from ancestral state)

** Simulation parameters
D:    Stop simulation after D fixations have occurred.
reps: Number of times simulation should be repeated



############################################################################################
########                 OUTPUTFILES                                                ########
############################################################################################

Output ALLOPATRY & PARAPATRY:

The allopatry version will produce one outputfile called "res_*.txt". 
This file lists the details of each fixation that has occured during the fixation, in the order in which they became fixed.

It has the following columns:
position: position of the fixed mutation in the genome (between 0 and 1; infinite-sites model)	
mut_time: generation in which the mutation arose
fix_time: generation in which the mutation became fixed
opt_dist: position of the optimum for trait 1

for each of the n traits
mut_effect: mutational effect size on this trait
dom: dominance coefficient for this trait


Output PARAPATRY:

For the parapatry version the "res_*.txt" file only lists mutations that have fixed in both populations. 
It also produces two additional outputfiles:
"freq_*.txt" lists the allele frequencies of segregating mutations in each population/deme every 100 generations,
(this interval can be changed in the recursion_parapatry.cpp script, under UPDATELENGTH).

Each row corresponds to one segregating mutation, and the columns report the following:
gen: generation of recording
index: index of this mutation in the mutations object
position: position of this mutation in the genome (between 0 and 1; infinite-sites model)	
freq_pop1: frequency of this mutation in population 1
freq_pop2: frequency of this mutation in population 2
mut_time: generation in which the mutation arose
fix_time: generation in which the mutation became fixed (should be 0 as segregating mutations are by definition not fixed)
opt_dist1: position of the optimum for trait 1 in population 1
opt_dist2: position of the optimum for trait 1 in population 2

for each of the n traits
mut_effect: mutational effect size on this trait
dom: dominance coefficient for this trait
 
 
 
The third file "pop_*.txt" contains the genotypes of each individual at the end of the simulation.

Here, each row corresponds to an individual. 
The first column indicates which population this individual belongs to. 
The subsequent columns correspond to the segregating mutations, with the number in the column name corresponding to the index of the mutation in the last recorded generation of the "freq_*.txt" file. 

A value of 0 indicates that this segregating mutation is absent in this individual.
A value of 1 indicates that this individual is heterozygous for this mutation.
A value of 2 indicates that this individual is homozygous for this mutation.


