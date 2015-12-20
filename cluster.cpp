#include <cstdio>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <ctime>
#include <cstring>
using namespace std;

const int EPS = 1e-6;

const int MAX_WEIGHT = 1000;
const int SLEIGH_WEIGHT = 10;

const int N_GIFTS = 100000;
const int MAX_GIFTS = 67;

const double NORTH_POLE[2] = {0.5 * M_PI, 0};
const double EARTH_RADIUS = 6372.8;

const int N_POP = 100;
const int MUTATE = 60;
const int TOURNAMENT = 14;

const double TOLERANCE = 1e-6;
const int CONVERGE = 1;
const int ATTEMPTS = 1;

struct Cluster {
    vector<int> gifts;
    double weight, cost;
};

default_random_engine generator;

inline double haversine(const double c1[], const double c2[]){

    double a, b, slon, slat;

    slat = sin((c2[0]-c1[0])/2);
    slon = sin((c2[1]-c1[1])/2);
    a = slat*slat+cos(c1[0])*cos(c2[0])*slon*slon;
    b = 2 * asin(sqrt(a)) * EARTH_RADIUS;

    return b;

}

inline double get_weariness(int route[], int n_gifts, double dist[][MAX_GIFTS],
			    double north[], double weights[]){
    
    int i;
    double weariness, weight;

    weight = SLEIGH_WEIGHT;
    for(i = 0; i < n_gifts; i++){
	weight += weights[i];
    }

    weariness = weight * north[route[0]];
    weight -= weights[route[0]];

    for(i = 1; i < n_gifts; i++){
	weariness += dist[route[i]][route[i-1]] * weight;
	weight -= weights[route[i]];
    }

    weariness += north[route[n_gifts-1]] * SLEIGH_WEIGHT;

    return weariness;
}

inline int pop_fitness(int pop[][MAX_GIFTS], int n_gifts, double dist[][MAX_GIFTS], 
		       double north[], double weights[], double fitness[]){

    int best = 0;

    for(int i = 0; i < N_POP; i++){
	fitness[i] = get_weariness(pop[i], n_gifts, dist, north, weights);

	if(fitness[i] < fitness[best]){
	    best = i;
	}
    }

    return best;

}

inline int tournament(double fitness[], default_random_engine &generator){

    int best, idx;

    for(int i = 0; i < TOURNAMENT; i++){
	idx = generator() % N_POP;

	if(i == 0 || fitness[idx] < fitness[best]){
	    best = idx;
	}
    }

    return best;
    
}

inline double tsp_optimized2(double coord[][2], double north[], double weights[], 
			     vector<int> &gifts){

    int a, b;
    double dist, cost_a, cost_b, cost_sleigh;

    a = gifts[0];
    b = gifts[1];

    dist = haversine(coord[a], coord[b]);
	
    cost_sleigh = SLEIGH_WEIGHT*(north[a]+dist+north[b]);
    cost_a = (weights[a]+weights[b])*north[a]+weights[b]*dist;
    cost_b = (weights[b]+weights[a])*north[b]+weights[a]*dist;

    if(cost_a < cost_b){
	return cost_sleigh + cost_a;
    } else {
	swap(gifts[0], gifts[1]);
	return cost_sleigh + cost_b;
    }

}

double tsp_genetic(double coord[][2], double all_north[], double all_weights[], 
		   vector<int> &gifts){
    
    int i, j, n_gifts = gifts.size();

    // Special optimization case for 2 clusters
    if(n_gifts == 2){
	return tsp_optimized2(coord, all_north, all_weights, gifts);
    }

    // Precompute and select relevant data
    static double dist[MAX_GIFTS][MAX_GIFTS];
    
    for(i = 0; i < n_gifts; i++){
	for(j = i + 1; j < n_gifts; j++){
	    dist[i][j] = dist[j][i] = haversine(coord[gifts[i]], coord[gifts[j]]);
	}
    }

    double north[MAX_GIFTS], weights[MAX_GIFTS];

    for(i = 0; i < n_gifts; i++){
	north[i] = all_north[gifts[i]];
	weights[i] = all_weights[gifts[i]];
    }

    // Initialize population
    static int pop[N_POP][MAX_GIFTS];

    for(i = 0; i < N_POP; i++){
	for(j = 0; j < n_gifts; j++){
	    pop[i][j] = j;
	}
	
	//shuffle(*(pop+i), *(pop+i)+n_gifts, generator);
    }
    
    // Generations
    static int new_pop[N_POP][MAX_GIFTS];
    int best, used[MAX_GIFTS], *parent[2], s, e, idx, converge;
    double fitness[N_POP], p_best;
    
    best = pop_fitness(pop, n_gifts, dist, north, weights, fitness);
    
    converge = 0;
    while(true){
	// Elitism
	memcpy(new_pop[0], pop[best], sizeof(new_pop[0]));

	// Evolve a new population
	for(int k = 1; k < N_POP; k++){
	    // Crossover
	    for(i = 0; i < 2; i++){
		idx = tournament(fitness, generator);
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
	    if((generator() % 100) <= MUTATE){
		i = generator() % n_gifts;
		j = generator() % n_gifts;

		swap(new_pop[k][i], new_pop[k][j]);
	    }
	}
	
	memcpy(pop, new_pop, sizeof(pop));

	// Convergance
	p_best = fitness[best];
	best = pop_fitness(pop, n_gifts, dist, north, weights, fitness);
	
	//printf("%lf\n", fitness[best]);
	if((1 - fitness[best] / p_best) < TOLERANCE){
	    converge++;
	    if(converge == CONVERGE){
		break;
	    }
	} else {
	    converge = 0;
	}
    }

    // Reorder
    for(i = 0; i < n_gifts; i++){
	pop[best][i] = gifts[pop[best][i]];
    }

    for(i = 0; i < n_gifts; i++){
	gifts[i] = pop[best][i];
    }

    return fitness[best];
}

int main(){

    int i, j;
    double coord[N_GIFTS][2], weights[N_GIFTS], north[N_GIFTS];

    generator.seed(time(NULL));

    // Input
    FILE *fp_gifts = fopen("gifts.csv", "r");

    fscanf(fp_gifts, "%*s\n");
    for(i = 0; i < N_GIFTS; i++){
	fscanf(fp_gifts, "%*d,%lf,%lf,%lf", *(coord+i), *(coord+i)+1, weights+i);
	coord[i][0] *= M_PI/180.0;
	coord[i][1] *= M_PI/180.0;
    }

    // Initialize
    int used[N_GIFTS];
    double total_cost;

    total_cost = 0;
    for(i = 0; i < N_GIFTS; i++){
	north[i] = haversine(NORTH_POLE, coord[i]);
	total_cost += north[i] * (2 * SLEIGH_WEIGHT + weights[i]);
	used[i] = 0;
    }

    printf("%12.0lf\n", total_cost);

    // Cluster
    vector<Cluster> clusters;
    int idx[N_GIFTS], count, region[N_GIFTS];

    for(i = 0; i < N_GIFTS; i++){
	idx[i] = i;
	
	region[i] = 0;
	if(coord[i][0] < -52.0/180.0*M_PI){
	    region[i] += 2;
	}

	if(coord[i][1] < -105.0/180.0*M_PI){
	    region[i] += 2;
	} else {
	    region[i] += 1;
	}
    }

    while(true){
	sort(idx, idx+N_GIFTS, [&used, &coord, &region](int a, int b){
		if(used[a]){
		    return false;
		} else if(used[b]){
		    return true;
		}

		if(region[a] < region[b]){
		    return true;
		} else if(region[a] > region[b]){
		    return false;
		} else {
		    return coord[a][1] < coord[b][1];
		}
	    });

	for(count = 0; count < MAX_GIFTS; count++){
	    if(used[idx[count]]){
		break;
	    }
	}

	sort(idx, idx+count, [&coord](int a, int b){
		return coord[a][0] > coord[b][0];
	    });
	
	Cluster cluster;
	cluster.weight = 0;

	for(i = 0; i < count; i++){
	    if(cluster.weight + weights[idx[i]] < MAX_WEIGHT){
		cluster.gifts.push_back(idx[i]);
		cluster.weight += weights[idx[i]];
		used[idx[i]] = 1;
	    }
	}

	if(cluster.gifts.size() > 0){
	    clusters.push_back(cluster);
	}

	if(count != MAX_GIFTS){
	    break;
	}
    }

    // Optimize using TSP
    double min_cost, cost;
    vector<int> route, min_route;

    total_cost = 0;
    for(i = 0; i < clusters.size(); i++){
	for(j = 0; j < ATTEMPTS; j++){
	    route = clusters[i].gifts;

	    cost = tsp_genetic(coord, north, weights, route);
	    
	    if(j == 0 || cost < min_cost){
		min_cost = cost;
		min_route = route;
	    }
	}

	total_cost += min_cost;
	clusters[i].gifts = min_route;
    }

    printf("%12.0lf\n", total_cost);

    // Output result
    FILE *fp_route = fopen("submission.csv", "w");
    fprintf(fp_route, "TripId,GiftId\n");
    for(i = 0; i < clusters.size(); i++){
	for(j = 0; j < clusters[i].gifts.size(); j++){
	    fprintf(fp_route, "%d,%d\n", i, clusters[i].gifts[j] + 1);
	}
    }

    return 0;
    
 }

