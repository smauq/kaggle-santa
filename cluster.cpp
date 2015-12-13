#include <cstdio>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <ctime>
#include <cstring>
using namespace std;

#define MAX_WEIGHT 1000
#define SLEIGH_WEIGHT 10

#define N_GIFTS 2392
#define MAX_GIFTS N_GIFTS

#define MAX_POP 100

inline double get_weariness(int route[], int n_gifts, double dist[][MAX_GIFTS],
		     double weights[]){

    double weariness = 0;

    for(int i = 1; i < n_gifts; i++){
	weariness += dist[route[i]][route[i-1]];
    }

    return weariness;
}

inline int pop_fitness(int pop[][MAX_GIFTS], int n_pop, int n_gifts, double dist[][MAX_GIFTS], 
		double weights[], double fitness[]){

    int best = 0;

    for(int i = 0; i < n_pop; i++){
	fitness[i] = get_weariness(pop[i], n_gifts, dist, weights);

	if(fitness[i] < fitness[best]){
	    best = i;
	}
    }

    return best;

}

inline int tournament(double fitness[], int n_pop, int size, default_random_engine &generator){

    int best, idx;

    for(int i = 0; i < size; i++){
	idx = generator() % n_pop;

	if(i == 0 || fitness[idx] < fitness[best]){
	    best = idx;
	}
    }

    return best;
    
}

int tsp_genetic(double coord[][2], double weights[], vector<int> &gifts,
		int n_pop=100, int gen=100, double mutate=40){
    
    int n_gifts, i, j;
    n_gifts = gifts.size();

    // Random number generator
    default_random_engine generator;
    generator.seed(time(NULL));

    // Precompute distance matrix
    static double dist[MAX_GIFTS][MAX_GIFTS];
    
    for(i = 0; i < n_gifts; i++){
	for(j = i + 1; j < n_gifts; j++){
	    dist[i][j] = dist[j][i] = sqrt(pow(coord[gifts[i]][0] - coord[gifts[j]][0], 2) +
					   pow(coord[gifts[i]][1] - coord[gifts[j]][1], 2));
	}
    }

    // Initialize population
    static int pop[MAX_POP][MAX_GIFTS];

    for(i = 0; i < n_pop; i++){
	for(j = 0; j < n_gifts; j++){
	    pop[i][j] = j;
	}
	
	shuffle(*(pop+i), *(pop+i)+n_gifts, generator);
    }

    // Tournament size
    int n_champions = (int)round(n_pop / 14.0);
    
    // Generations
    static int new_pop[MAX_POP][MAX_GIFTS];
    int best, used[MAX_GIFTS], *parent[2], s, e, idx;
    double fitness[MAX_POP];

    for(int g = 0; g < gen; g++){
	best = pop_fitness(pop, n_pop, n_gifts, dist, weights, fitness);
	
	//printf("%lf\n", fitness[best]);

	// Elitism
	memcpy(new_pop[0], pop[best], sizeof(new_pop[0]));

	// Evolve a new population
	for(int k = 1; k < n_pop; k++){
	    // Crossover
	    for(i = 0; i < 2; i++){
		idx = tournament(fitness, n_pop, n_champions, generator);
		parent[i] = pop[idx];
	    }

	    s = generator() % n_gifts;
	    e = s + (generator() % (n_gifts - s)) + 1;

	    memset(used, 0, sizeof(used));
	    for(i = s; i < e; i++){
		used[parent[0][i]] = 1;
		new_pop[k][i] = parent[0][i];
	    }

	    idx = 0;
	    for(i = 0; i < s; i++){
		while(used[parent[1][idx]]){
		    idx++;	    
		}

		new_pop[k][i] = parent[1][idx];
		idx++;
	    }

	    for(i = e; i < n_gifts; i++){
		while(used[parent[1][idx]]){
		    idx++;	    
		}

		new_pop[k][i] = parent[1][idx];
		idx++;
	    }
	    
	    // Mutate
	    if((generator() % 100) <= mutate){
		i = generator() % n_gifts;
		j = generator() % n_gifts;

		swap(new_pop[k][i], new_pop[k][j]);
	    }
	}
	
	memcpy(pop, new_pop, sizeof(pop));
    }

    return 0;
}

int main(){

    int gifts[N_GIFTS], i;
    double coord[N_GIFTS][2], weights[N_GIFTS];

    // Input
    FILE *fp_gifts = fopen("3", "r");

    fscanf(fp_gifts, "%*s\n");
    for(i = 0; i < N_GIFTS; i++){
	//fscanf(fp_gifts, "%d,%lf,%lf,%lf", gifts+i, *(coord+i),
	//       *(coord+i)+1, weights+i);
	fscanf(fp_gifts, "%d %le %le", gifts+i, *(coord+i), *(coord+i)+1);
    }
    
    vector<int> labels;
    for(i = 0; i < N_GIFTS; i++){
	labels.push_back(i);
    }
    
    tsp_genetic(coord, weights, labels, 100, 300);
}
