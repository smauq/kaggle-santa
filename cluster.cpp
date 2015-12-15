#include <cstdio>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <ctime>
#include <cstring>
using namespace std;

const int MAX_WEIGHT = 1000;
const int SLEIGH_WEIGHT = 10;

const int N_GIFTS = 10000;
const int MAX_GIFTS = MAX_WEIGHT;

const double NORTH_POLE[2] = {0.5 * M_PI, 0};
const double EARTH_RADIUS = 6372.8;

const int N_POP = 100;
const int MUTATE = 40;
const int TOURNAMENT = 7;
const double TOLERANCE = 0.01;

const int NEIGHBOURS = 5;

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

double tsp_genetic(double coord[][2], double all_north[], double all_weights[], int gifts[], 
		   int n_gifts, int route[]=NULL){
    
    int i, j;
    
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
    int best, used[MAX_GIFTS], *parent[2], s, e, idx;
    double fitness[N_POP], p_best;
    
    best = pop_fitness(pop, n_gifts, dist, north, weights, fitness);

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
	    break;
	}
    }

    if(route != NULL){
	for(i = 0; i < n_gifts; i++){
	    route[i] = gifts[pop[best][i]];
	}
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
    vector<int> clusters[N_GIFTS];
    int n_clusters = N_GIFTS;
    double cl_weights[N_GIFTS], cl_cost[N_GIFTS];
    
    for(i = 0; i < N_GIFTS; i++){
	clusters[i].push_back(i);
	cl_weights[i] = weights[i];
	cl_cost[i] = north[i] * (2 * SLEIGH_WEIGHT + weights[i]);
    }
    
    // Merge
    double cl_center[N_GIFTS][2], dist[N_GIFTS], cost;
    int idx[N_GIFTS*NEIGHBOURS], merge[MAX_GIFTS], size_i, size_j, used[N_GIFTS], edge_count, a, b;
    vector<double> merge_cost;
    vector<int> gift_a, gift_b;
        
    for(e = 0; e < 10; e++){
	// Print current total
	cost = 0;
	for(i = 0; i < n_clusters; i++){
	    cost += cl_cost[i];
	}

	printf("Epoch %2d: %12.0lf\n", e, cost);

	// Cluster centers
	for(i = 0; i < n_clusters; i++){
	    cl_center[i][0] = 0;
	    cl_center[i][1] = 0;

	    for(j = 0; j < clusters[i].size(); j++){
		cl_center[i][0] += coord[clusters[i][j]][0];
		cl_center[i][1] += coord[clusters[i][j]][1];
	    }

	    cl_center[i][0] /= clusters[i].size();
	    cl_center[i][1] /= clusters[i].size();
	}

	merge_cost.clear();
	gift_a.clear();
	gift_b.clear();

	for(i = 0; i < n_clusters; i++){
	    // Distance to neighbours
	    for(j = 0; j < n_clusters; j++){
		if(i != j){
		    dist[j] = haversine(cl_center[i], cl_center[j]);
		} else {
		    dist[j] = 2 * M_PI * EARTH_RADIUS;
		}

		idx[j] = j;
	    }

	    // Find the closest neighbours
	    nth_element(idx, idx+NEIGHBOURS, idx+n_clusters, [&dist](int a, int b) 
			{return dist[a] < dist[b];});

	    // Merge cost
	    size_i = clusters[i].size();
	    for(k = 0; k < size_i; k++){
		merge[k] = clusters[i][k];
	    }

	    for(m = 0; m < NEIGHBOURS; m++){
		j = idx[m];

		if(cl_weights[i] + cl_weights[j] > MAX_WEIGHT){
		    continue;
		}

		size_j = clusters[j].size();
		for(k = 0; k < size_j; k++){
		    merge[k + size_i] = clusters[j][k];
		}

		cost = tsp_genetic(coord, north, weights, merge, size_i+size_j);
		
		cost -= cl_cost[i] + cl_cost[j];

		if(cost < 0){
		    merge_cost.push_back(cost);
		    gift_a.push_back(i);
		    gift_b.push_back(j);
		}
	    }
	}

	// Merge
	edge_count = merge_cost.size();
	for(i = 0; i < edge_count; i++){
	    idx[i] = i;
	}

	sort(idx, idx+edge_count, [&merge_cost](int a, int b)
	     {return merge_cost[a] < merge_cost[b];});

	memset(used, 0, sizeof(used));

	vector<int> new_clusters[N_GIFTS];
	int new_n;
	double new_weights[N_GIFTS], new_cost[N_GIFTS];

	new_n = 0;
	for(i = 0; i < edge_count; i++){
	    a = gift_a[idx[i]];
	    b = gift_b[idx[i]];

	    if(used[a] || used[b]){
		continue;
	    }
	    
	    new_clusters[new_n] = clusters[a];
	    new_clusters[new_n].insert(new_clusters[new_n].end(), clusters[b].begin(), 
				       clusters[b].end());

	    new_weights[new_n] = cl_weights[a] + cl_weights[b];
	    new_cost[new_n] = cl_cost[a] + cl_cost[b] + merge_cost[idx[i]];
	    new_n++;

	    used[a] = used[b] = 1;
	}

	for(i = 0; i < n_clusters; i++){
	    if(!used[i]){
		new_clusters[new_n] = clusters[i];
		new_weights[new_n] = cl_weights[i];
		new_cost[new_n] = cl_cost[i];
		new_n++;
	    }
	}

	for(i = 0; i < new_n; i++){
	    clusters[i] = new_clusters[i];
	    cl_weights[i] = new_weights[i];
	    cl_cost[i] = new_cost[i];
	}

	n_clusters = new_n;
    }
    
    // Print result
    /*FILE *fp_submission = fopen("submission.csv", "w");
    fprintf(fp_submission, "GiftId,TripId\n");

    int route[MAX_GIFTS], gifts[MAX_GIFTS];
    
    cost = 0;
    for(i = 0; i < n_clusters; i++){
        tsp_genetic(coord, north, weights, &clusters[i][0], clusters[i].size(), 
			  100, 100, route);
    
	for(j = 0; j < clusters[i].size(); j++){
	    fprintf(fp_submission, "%d,%d\n", route[j]+1, i);
	}

	cost += cl_cost[i];
    }

    printf("%lf\n", cost);*/
    
 }

