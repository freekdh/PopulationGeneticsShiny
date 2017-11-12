#include <ostream>
#include "random.h"
#include <assert.h>
#include <fstream>
#include <Rcpp.h>
#include <bitset>
#include <vector>

// I am changing this file to see if git works..

// [[Rcpp::plugins(cpp11)]]

enum allele { A, rest };
enum Invaderallele { sex, nosex };

class Haploid{
public:
	Haploid(unsigned int p) {nMutations = p;}
	void addNmutations(unsigned int add) {
		nMutations += add;
	}
	unsigned int returnmutations() { 
		return nMutations;
	}

private:
	unsigned int nMutations;
};

class Diploid {
public:
	Diploid(const int &p, const int &q) { nMutations1 = p; nMutations2 = q; }

	Diploid(Diploid * const &p1, Diploid * const &p2) {
		int n = p1->returnmutations();
		nMutations1 = rnd::binomial(n, 0.5);

		n = p2->returnmutations();
		nMutations2 = rnd::binomial(n, 0.5);
	}

	void addNmutations(unsigned int chrom1, unsigned int chrom2) {
		nMutations1 += chrom1;
		nMutations2 += chrom2;
	}

	unsigned int returnmutations() {
		return nMutations1 + nMutations2;
	}

private:

	int nMutations1;
	int nMutations2;
};

class Individual {
public:
	Individual(bool p) { type = p; }
	bool ReturnType() { return type; }
private:
	bool type;
};

class Individual2 {
public:
	Individual2(const int &p, const int &q, const bool &i) { nMutations1 = p; nMutations2 = q; InvaderAllele = i; }

	Individual2(Individual2 * const &p1, Individual2 * const &p2) {
		int n = p1->returnmutations();
		nMutations1 = rnd::binomial(n, 0.5);

		n = p2->returnmutations();
		nMutations2 = rnd::binomial(n, 0.5);

		assert(p1->returnInvaderAllele() == p2->returnInvaderAllele());
		InvaderAllele = p1->returnInvaderAllele();
	}

	void addNmutations(unsigned int chrom1, unsigned int chrom2) {
		nMutations1 += chrom1;
		nMutations2 += chrom2;
	}

	void setInvaderAllele(bool p) { InvaderAllele = p; }
	bool returnInvaderAllele() { return InvaderAllele; }

	unsigned int returnmutations() {
		return nMutations1 + nMutations2;
	}

private:

	int nMutations1;
	int nMutations2;

	bool InvaderAllele;
};


// [[Rcpp::export]]
Rcpp::NumericVector MullerRatchetHaploid(const int &NGENERATIONS, const int &NINDIDIVUALS, const double &MRATE, const double &hs) {
	rnd::set_seed();
	assert(hs < 1);

	Haploid* parents[NINDIDIVUALS];
	Haploid* offspring[NINDIDIVUALS];

	for (int i = 0; i < NINDIDIVUALS; ++i) {
		parents[i] = new Haploid(0);
	}

	Rcpp::NumericVector data(NGENERATIONS);

	for (int i = 0; i < NGENERATIONS; ++i) {
		// mutations:
		for (int j = 0; j < NINDIDIVUALS; ++j) {
			parents[j]->addNmutations(rnd::poisson(MRATE));
		}

		// selection:
		rnd::discrete_distribution dist(NINDIDIVUALS);

		for (int j = 0; j < NINDIDIVUALS; ++j) {
		dist[j] = pow(1.0-hs, parents[j]->returnmutations());			
		}

		for (int j = 0; j < NINDIDIVUALS; ++j) {
			offspring[j] = parents[dist.sample()];
		}
		
		for (int j = 0; j < NINDIDIVUALS; ++j) {
			parents[j] = offspring[j];
		}

		// Collect data 
		
		int minimum = parents[0]->returnmutations();
		for (int j = 1; j < NINDIDIVUALS; ++j) {
			int temp = parents[j]->returnmutations();
			if (temp < minimum) minimum = temp;
		}
	
		data[i] = minimum;
	}

	return data;
}

// [[Rcpp::export]]
Rcpp::NumericVector MullerRatchetDiploid(const int &NGENERATIONS, const int &NINDIDIVUALS, const double &MRATE, const double &hs, const double &p_sex) {
	rnd::set_seed();
	assert(hs < 1);

	Diploid* parents[NINDIDIVUALS];
	Diploid* offspring[NINDIDIVUALS];

	for (int i = 0; i < NINDIDIVUALS; ++i) {
		parents[i] = new Diploid(0,0);
	}

	Rcpp::NumericVector data(NGENERATIONS);

	for (int i = 0; i < NGENERATIONS; ++i) {
		// mutations:
		for (int j = 0; j < NINDIDIVUALS; ++j) {
			parents[j]->addNmutations(rnd::poisson(MRATE), rnd::poisson(MRATE));
		}

		// selection:
		rnd::discrete_distribution dist(NINDIDIVUALS);

		for (int j = 0; j < NINDIDIVUALS; ++j) {
			dist[j] = pow(1.0 - hs, parents[j]->returnmutations());
		}

		for (int j = 0; j < NINDIDIVUALS; ++j) {
			offspring[j] = rnd::uniform() < p_sex ? new Diploid(parents[dist.sample()], parents[dist.sample()]) : new Diploid(*parents[dist.sample()]);
		}

		for (int j = 0; j < NINDIDIVUALS; ++j) {
			delete parents[j];
			parents[j] = offspring[j];

		}

		// Collect data 

		int minimum = parents[0]->returnmutations();
		for (int j = 1; j < NINDIDIVUALS; ++j) {
			int temp = parents[j]->returnmutations();
			if (temp < minimum) minimum = temp;
		}

		data[i] = minimum;
	}

	return data;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix DriftSimulation(const int &NGENERATIONS, const int &Nindividuals, const double&pA, const double&mutrate, const int &NREPS, const bool &ISDIPLOID) {
	rnd::set_seed();

	const int N_total = ISDIPLOID ? 2 * Nindividuals : Nindividuals;

	Individual *parents[N_total];	Individual *offspring[N_total];

	Rcpp::NumericMatrix data(NGENERATIONS, NREPS);

	for (int rep = 0; rep < NREPS; ++rep) {
		for (int i = 0; i < N_total; ++i) {
			parents[i] = rnd::uniform() < pA ? new Individual(A) : new Individual(rest);
		}

		// Safe initial population
		for (int i = 0; i < N_total; ++i) {
			if (static_cast<int>(parents[i]->ReturnType()) == A) ++data(0, rep);
		}
		data(0, rep) /= static_cast<double>(N_total);

		for (int j = 1; j < NGENERATIONS; ++j) {
			//create offspring
			for (int i = 0; i < N_total; ++i) {
				int p = rnd::integer(N_total);
				if (parents[p]->ReturnType() == A) {
					offspring[i] = (rnd::uniform() < mutrate) ? new Individual(rest) : new Individual(A);
					}
				offspring[i] = new Individual(parents[rnd::integer(N_total)]->ReturnType());
			}

			//offspring become parent
			for (int i = 0; i < N_total; ++i) {
				delete parents[i];
				parents[i] = offspring[i];
			}

			//safe files
			for (int i = 0; i < N_total; ++i) {
				if (static_cast<int>(parents[i]->ReturnType()) == A) ++data(j, rep);
			}

			data(j, rep) /= static_cast<double>(N_total);
		}
		for (int i = 0; i < N_total; ++i) {
			delete parents[i];
		}
	}
	return data;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix MullerRatchetInvader(int &NINDIDIVUALS, const double &MRATE, const double &hs, const int &LATENCY) {
	rnd::set_seed();
	assert(hs < 1);

	// class Individuls
	std::vector<Individual2*> parents(NINDIDIVUALS);
	std::vector<Individual2*> offspring(NINDIDIVUALS);
	std::vector<Individual2*>::iterator iter;

	for (iter = parents.begin(); iter != parents.end(); ++iter) {
		*iter = new Individual2(0, 0, nosex);
	}

	std::vector<int> data_min_mutation;
	std::vector<double> data_sexy;
	std::vector<double> data_avfitness;
	/*
	Individuals go through the following life stage events:
	- Mutations
	- Reproduction based on fertility selection
	*/

	// latency of introducing invader allele
	for (int i = 0; i < LATENCY; ++i) {
		// Mutation:
		for (int j = 0; j < NINDIDIVUALS; ++j) {
			parents[j]->addNmutations(rnd::poisson(MRATE), rnd::poisson(MRATE));
		}

		// Average fitness:
		rnd::discrete_distribution dist(NINDIDIVUALS);
		double AV_fitness = 0.0;
		for (int j = 0; j < NINDIDIVUALS; ++j) {
			dist[j] = pow(1.0 - hs, parents[j]->returnmutations());
			AV_fitness += dist[j];
		}
		AV_fitness /= static_cast<double>(NINDIDIVUALS);

		for (int j = 0; j < NINDIDIVUALS; ++j) {
			offspring[j] = new Individual2(*parents[dist.sample()]);
		}

		for (int j = 0; j < NINDIDIVUALS; ++j) {
			delete parents[j];
			parents[j] = offspring[j];
		}

		// Collect data 
		int minimum = parents[0]->returnmutations();
		for (int j = 1; j < NINDIDIVUALS; ++j) {
			int temp = parents[j]->returnmutations();
			if (temp < minimum) minimum = temp;
		}
		data_min_mutation.push_back(minimum);
		data_sexy.push_back(0.0);
		data_avfitness.push_back(AV_fitness);
	}

	// Introduce invader allele
	parents[rnd::integer(NINDIDIVUALS)]->setInvaderAllele(sex);
	std::vector<Individual2*> sexy;
	std::vector<Individual2*> nosexy;
	int Nsexy;
	int Nnosexy;

	do{
		// Splitpopulation and add mutations;
		sexy.clear(); nosexy.clear();
		for (iter = parents.begin(); iter != parents.end(); ++iter) {
			(*iter)->returnInvaderAllele() == sex ? sexy.push_back((*iter)) : nosexy.push_back((*iter));
			(*iter)->addNmutations(rnd::poisson(MRATE), rnd::poisson(MRATE));
		}

		// Sampling distribution:
		double AV_fitness = 0.0;
		rnd::discrete_distribution dist(NINDIDIVUALS);
		for (int j = 0; j < NINDIDIVUALS; ++j) {
			dist[j] = pow(1.0 - hs, parents[j]->returnmutations());
			AV_fitness += dist[j];
		}
		AV_fitness /= static_cast<double>(NINDIDIVUALS);

		// Create offspring:
		Nsexy = sexy.size();
		Nnosexy = nosexy.size();
		int parent1;
		for (int j = 0; j < NINDIDIVUALS; ++j) {
			parent1 = dist.sample();
			offspring[j] = parents[parent1]->returnInvaderAllele() == sex ? new Individual2(parents[parent1], sexy[rnd::integer(Nsexy)]) : new Individual2(*parents[parent1]);
			}

		for (iter = parents.begin(); iter != parents.end(); ++iter) {delete (*iter);}
		parents = offspring;

		// Collect data 
		int minimum = parents[0]->returnmutations();
		for (int j = 1; j < NINDIDIVUALS; ++j) {
			int temp = parents[j]->returnmutations();
			if (temp < minimum) minimum = temp;
		}

		data_min_mutation.push_back(minimum);
		data_sexy.push_back(static_cast<double>(Nsexy)/static_cast<double>(NINDIDIVUALS));
		data_avfitness.push_back(AV_fitness);
	} while (Nsexy != NINDIDIVUALS && Nnosexy != NINDIDIVUALS);

	const int NEXTRA = data_min_mutation.size(); // Extra generations after invader allele either fixated or not.

	// Continue simulation a bit more
	if (Nsexy == NINDIDIVUALS) {
		for (int i = 0; i < NEXTRA; ++i) {
			// Mutation:
			for (int j = 0; j < NINDIDIVUALS; ++j) {
				parents[j]->addNmutations(rnd::poisson(MRATE), rnd::poisson(MRATE));
			}

			// Selection:
			double AV_fitness = 0.0;
			rnd::discrete_distribution dist(NINDIDIVUALS);
			for (int j = 0; j < NINDIDIVUALS; ++j) {
				dist[j] = pow(1.0 - hs, parents[j]->returnmutations());
				AV_fitness += dist[j];
			}
			AV_fitness /= static_cast<double>(NINDIDIVUALS);

			int parent1;
			for (int j = 0; j < NINDIDIVUALS; ++j) {
				parent1 = dist.sample();
				offspring[j] = new Individual2(parents[parent1], sexy[rnd::integer(Nsexy)]);
			}

			for (int j = 0; j < NINDIDIVUALS; ++j) {
				delete parents[j];
				parents[j] = offspring[j];
			}

			// Collect data 
			int minimum = parents[0]->returnmutations();
			for (int j = 1; j < NINDIDIVUALS; ++j) {
				int temp = parents[j]->returnmutations();
				if (temp < minimum) minimum = temp;
			}
			data_min_mutation.push_back(minimum);
			data_sexy.push_back(1.0);
			data_avfitness.push_back(AV_fitness);
		}
	}
	else{ 
		for (int i = 0; i < NEXTRA; ++i) {
			// Mutation:
			for (int j = 0; j < NINDIDIVUALS; ++j) {
				parents[j]->addNmutations(rnd::poisson(MRATE), rnd::poisson(MRATE));
			}

			// Selection:
			double AV_fitness = 0.0;
			rnd::discrete_distribution dist(NINDIDIVUALS);
			for (int j = 0; j < NINDIDIVUALS; ++j) {
				dist[j] = pow(1.0 - hs, parents[j]->returnmutations());
				AV_fitness += dist[j];
			}
			AV_fitness /= static_cast<double>(NINDIDIVUALS);

			for (int j = 0; j < NINDIDIVUALS; ++j) {
				offspring[j] = new Individual2(*parents[dist.sample()]);
			}

			for (int j = 0; j < NINDIDIVUALS; ++j) {
				delete parents[j];
				parents[j] = offspring[j];
			}

			// Collect data 
			int minimum = parents[0]->returnmutations();
			for (int j = 1; j < NINDIDIVUALS; ++j) {
				int temp = parents[j]->returnmutations();
				if (temp < minimum) minimum = temp;
			}
			data_min_mutation.push_back(minimum);
			data_sexy.push_back(0.0);
			data_avfitness.push_back(AV_fitness);
		}
	}

	// Write to R
	const int N_total = data_min_mutation.size();
	Rcpp::NumericMatrix NmatrixR(N_total,3);
	
	for (unsigned int i = 0; i < data_min_mutation.size(); ++i) {
		NmatrixR(i, 0) = data_min_mutation[i];
		NmatrixR(i, 1) = data_sexy[i];
		NmatrixR(i, 2) = data_avfitness[i];
	}

	return NmatrixR;
}