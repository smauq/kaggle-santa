#include <cstdio>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <ctime>
#include <cstring>
using namespace std;

const double EPS = 1e-6;

const int MAX_WEIGHT = 1000;
const int SLEIGH_WEIGHT = 10;

const int N_GIFTS = 100000;
const int MAX_GIFTS = 67;

const double NORTH_POLE[2] = {0.5 * M_PI, 0};
const double EARTH_RADIUS = 6371.0;

const int N_POP = 100;
const int MUTATE = 60;
const int TOURNAMENT = 14;

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

inline double get_weariness(vector<int> &gifts, double coord[][2],
			    double north[], double weights[]){
    
    int i, n_gifts = gifts.size();
    double weariness, weight;

    weight = SLEIGH_WEIGHT;
    for(i = 0; i < n_gifts; i++){
	weight += weights[gifts[i]];
    }

    weariness = weight * north[gifts[0]];
    weight -= weights[gifts[0]];

    for(i = 1; i < n_gifts; i++){
	weariness += haversine(coord[gifts[i]], coord[gifts[i-1]]) * weight;
	weight -= weights[gifts[i]];
    }

    weariness += north[gifts[n_gifts - 1]] * SLEIGH_WEIGHT;

    return weariness;
}

inline double get_weariness(int gifts[], int n_gifts, double dist[][MAX_GIFTS],
			    double north[], double weights[]){
    
    int i;
    double weariness, weight;

    weight = SLEIGH_WEIGHT;
    for(i = 0; i < n_gifts; i++){
	weight += weights[i];
    }

    weariness = weight * north[gifts[0]];
    weight -= weights[gifts[0]];

    for(i = 1; i < n_gifts; i++){
	weariness += dist[gifts[i]][gifts[i-1]] * weight;
	weight -= weights[gifts[i]];
    }

    weariness += north[gifts[n_gifts-1]] * SLEIGH_WEIGHT;

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
		   vector<int> &gifts, int optimize=10){
    
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
    }
    
    // Generations
    static int new_pop[N_POP][MAX_GIFTS];
    int best, used[MAX_GIFTS], *parent[2], s, e, idx, converge;
    double fitness[N_POP], p_best;
    
    best = pop_fitness(pop, n_gifts, dist, north, weights, fitness);
    
    for(converge = 0; converge < optimize; ){
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
	
	if((1 - fitness[best] / p_best) < EPS){
	    converge++;
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

    int i, j, k, progbar;
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

    printf("STR\t\t%12.0lf\n", total_cost);

    // Cluster
    fprintf(stderr, "CLT ");

    vector<Cluster> clusters;
    int idx[N_GIFTS], count, region[N_GIFTS], total_used;

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

    progbar = 1;
    total_cost = 0;
    total_used = 0;
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
		total_used++;
	    }
	}

	if(cluster.gifts.size() > 0){
	    cluster.cost = get_weariness(cluster.gifts, coord, north, weights);
	    total_cost += cluster.cost;
	    clusters.push_back(cluster);
	}

	if(total_used / (N_GIFTS/10.0) >= progbar){
	    fprintf(stderr, "|");
	    progbar++;
	}

	if(count != MAX_GIFTS){
	    break;
	}
    }

    fprintf(stderr, "\t%12.0lf\n", total_cost);

    // Optimize using TSP
    fprintf(stderr, "TSP ");

    double min_cost, cost;
    vector<int> route, min_route;

    progbar = 1;
    total_cost = 0;
    for(i = 0; i < clusters.size(); i++){
	for(j = 0; j < ATTEMPTS; j++){
	    route = clusters[i].gifts;

	    cost = tsp_genetic(coord, north, weights, route, 100);
	    
	    if(j == 0 || cost < min_cost){
		min_cost = cost;
		min_route = route;
	    }
	}

	total_cost += min_cost;
	clusters[i].gifts = min_route;
	clusters[i].cost = min_cost;

	if((i+1) / (clusters.size()/10.0) >= progbar){
	    fprintf(stderr, "|");
	    progbar++;
	}
    }

    fprintf(stderr, "\t%12.0lf\n", total_cost);

    // Split
    fprintf(stderr, "SPL ");

    double cost_a, cost_b, min_cost_a, min_cost_b;
    vector<int> route_a, route_b, min_route_a, min_route_b;

    progbar = 1;
    for(i = 0; i < clusters.size(); i++){
	min_cost = clusters[i].cost;
	for(j = 1; j < clusters[i].gifts.size(); j++){
	    route_a.clear();
	    for(k = 0; k < j; k++){
		route_a.push_back(clusters[i].gifts[k]);
	    }

	    route_b.clear();
	    for(k = j; k < clusters[i].gifts.size(); k++){
		route_b.push_back(clusters[i].gifts[k]);
	    }

	    cost_a = tsp_genetic(coord, north, weights, route_a, 10);
	    cost_b = tsp_genetic(coord, north, weights, route_b, 10);
	    
	    if(cost_a + cost_b < min_cost){
		min_cost = cost_a + cost_b;
		min_route_a = route_a;
		min_route_b = route_b;
		min_cost_a = cost_a;
		min_cost_b = cost_b;
	    }
	}

	if(abs(min_cost - clusters[i].cost) > EPS){
	    total_cost += min_cost - clusters[i].cost;
	    //printf("%4d %6.0lf\n", i, min_cost - clusters[i].cost);
	    
	    Cluster cluster;
	    cluster.gifts = min_route_a;
	    cluster.cost = min_cost_a;
	    clusters.push_back(cluster);

	    clusters[i].gifts = min_route_b;
	    clusters[i].cost = min_cost_b;
	}

	if((i+1) / (clusters.size()/10.0) >= progbar){
	    fprintf(stderr, "|");
	    progbar++;
	}
    }

    fprintf(stderr, "\t%12.0lf\n", total_cost);

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

