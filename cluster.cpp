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

const int N_GIFTS = 1000;
const int MAX_GIFTS = MAX_WEIGHT;

const double NORTH_POLE[2] = {0.5 * M_PI, 0};
const double EARTH_RADIUS = 6372.8;

const int N_POP = 100;
const int MUTATE = 40;
const int TOURNAMENT = 7;
const double TOLERANCE = 0.01;
const int CONVERGE = 7;

const int NEIGHBOURS = 10;

struct Cluster {
    vector<int> gifts;
    int target;
    double weight, cost, merge_cost;
    double center[2];
};

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

    // Random number generator
    default_random_engine generator;
    generator.seed(time(NULL));

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
	
	shuffle(*(pop+i), *(pop+i)+n_gifts, generator);
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

    int i, j, k, m, e;
    double coord[N_GIFTS][2], weights[N_GIFTS], north[N_GIFTS];

    // Input
    FILE *fp_gifts = fopen("gifts.csv", "r");

    fscanf(fp_gifts, "%*s\n");
    for(i = 0; i < N_GIFTS; i++){
	fscanf(fp_gifts, "%*d,%lf,%lf,%lf", *(coord+i), *(coord+i)+1, weights+i);
	coord[i][0] *= M_PI/180.0;
	coord[i][1] *= M_PI/180.0;
    }

    for(i = 0; i < N_GIFTS; i++){
	north[i] = haversine(NORTH_POLE, coord[i]);
    }

    // Initialize clusters
    Cluster clusters[N_GIFTS], new_clusters[N_GIFTS];
    int n_clusters = N_GIFTS;
    
    for(i = 0; i < N_GIFTS; i++){
	clusters[i].gifts.push_back(i);
	clusters[i].weight = weights[i];
	clusters[i].cost = north[i] * (2 * SLEIGH_WEIGHT + weights[i]);
    }
    
    // Merge
    double dist[N_GIFTS], cost;
    int idx[N_GIFTS], used[N_GIFTS], a, b, edge_count;
    vector<int> merge(MAX_GIFTS);

    for(e = 0; e < 25; e++){
	// Print current total
	cost = 0;
	for(i = 0; i < n_clusters; i++){
	    cost += clusters[i].cost;
	}

	printf("Epoch %2d: %12.0lf\n", e, cost);

	// Cluster centers
	for(i = 0; i < n_clusters; i++){
	    clusters[i].center[0] = 0;
	    clusters[i].center[1] = 0;

	    for(j = 0; j < clusters[i].gifts.size(); j++){
		clusters[i].center[0] += coord[clusters[i].gifts[j]][0];
		clusters[i].center[1] += coord[clusters[i].gifts[j]][1];
	    }

	    clusters[i].center[0] /= clusters[i].gifts.size();
	    clusters[i].center[1] /= clusters[i].gifts.size();
	}

	edge_count = 0;

	for(i = 0; i < n_clusters; i++){
	    // Distance to neighbours
	    for(j = 0; j < n_clusters; j++){
		if(i != j){
		    dist[j] = haversine(clusters[i].center, clusters[j].center);
		} else {
		    dist[j] = 2 * M_PI * EARTH_RADIUS;
		}

		idx[j] = j;
	    }

	    // Find the closest neighbours
	    nth_element(idx, idx+NEIGHBOURS, idx+n_clusters, [&dist](int a, int b) 
			{return dist[a] < dist[b];});

	    // Merge cost
	    clusters[i].merge_cost = 0;

	    for(m = 0; m < NEIGHBOURS; m++){
		j = idx[m];

		if(clusters[i].weight + clusters[j].weight > MAX_WEIGHT){
		    continue;
		}

		merge = clusters[i].gifts;
		merge.insert(merge.begin(), clusters[j].gifts.begin(), 
			     clusters[j].gifts.end());

		cost = tsp_genetic(coord, north, weights, merge) - 
		    (clusters[i].cost + clusters[j].cost);
		
		if(cost < clusters[i].merge_cost){
		    clusters[i].merge_cost = cost;
		    clusters[i].target = j;
		}
	    }

	    if(abs(clusters[i].merge_cost) > EPS){
		edge_count++;
	    }
	}

	// Merge
	if(edge_count == 0){
	    continue;
	}

	for(i = 0; i < n_clusters; i++){
	    idx[i] = i;
	}

	sort(idx, idx+n_clusters, [&clusters](int a, int b)
	     {return clusters[a].merge_cost < clusters[b].merge_cost;});

	memset(used, 0, sizeof(used));

	m = 0;
	for(i = 0; i < max(edge_count / 4, 1); i++){
	    a = idx[i];
	    b = clusters[a].target;

	    if(used[a] || used[b]){
		continue;
	    }
	    
	    new_clusters[m].gifts = clusters[a].gifts;
	    new_clusters[m].gifts.insert(new_clusters[m].gifts.end(), clusters[b].gifts.begin(), 
					 clusters[b].gifts.end());

	    new_clusters[m].weight = clusters[a].weight + clusters[b].weight;
	    new_clusters[m].cost = clusters[a].cost + clusters[b].cost + clusters[a].merge_cost;

	    m++;
	    used[a] = used[b] = 1;
	}

	for(i = 0; i < n_clusters; i++){
	    if(!used[i]){
		new_clusters[m] = clusters[i];
		m++;
	    }
	}

	for(i = 0; i < m; i++){
	    clusters[i] = new_clusters[i];
	}

	n_clusters = m;
    }

    return 0;
    
 }

