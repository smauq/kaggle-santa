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
const int MAX_GIFTS = 100;

const double NORTH_POLE[2] = {0.5 * M_PI, 0};
const double EARTH_RADIUS = 6371.0;

const int N_POP = 100;
const int MUTATE = 60;
const int TOURNAMENT = 14;

const int ATTEMPTS = 1;
const int NEIGHBOURS = 3;

struct Cluster {
    vector<int> gifts;
    double weight, cost, center[2];
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

inline double get_weariness(vector<int> &gifts, double dist[][MAX_GIFTS],
			    double north[], double weights[]){
    
    int i, n_gifts = gifts.size();
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

inline int pop_fitness(vector<int> pop[], double dist[][MAX_GIFTS], double north[], 
		       double weights[], double fitness[]){

    int best = 0;

    for(int i = 0; i < N_POP; i++){
	fitness[i] = get_weariness(pop[i], dist, north, weights);

	if(fitness[i] < fitness[best]){
	    best = i;
	}
    }

    return best;

}

inline int tournament(double fitness[]){

    int best, idx;

    for(int i = 0; i < TOURNAMENT; i++){
	idx = generator() % N_POP;

	if(i == 0 || fitness[idx] < fitness[best]){
	    best = idx;
	}
    }

    return best;
    
}

inline double tsp_optimized2(vector<int> &gifts, double coord[][2], double north[], 
			     double weights[]){

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

double tsp_genetic(vector<int> &gifts, double coord[][2], double all_north[], 
		   double all_weights[], int optimize=10){
    
    int i, j, n_gifts = gifts.size();

    // Special optimization case for 2 gifts
    if(n_gifts == 2){
	return tsp_optimized2(gifts, coord, all_north, all_weights);
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
    static vector<int> pop[N_POP], new_pop[N_POP];

    for(i = 0; i < N_POP; i++){
	pop[i].resize(n_gifts);
	for(j = 0; j < n_gifts; j++){
	    pop[i][j] = j;
	}

	new_pop[i].resize(n_gifts);
    }
    
    // Generations
    vector<int> parent[2];
    int best, used[MAX_GIFTS], s, e, idx, converge;
    double fitness[N_POP], p_best;
    
    best = pop_fitness(pop, dist, north, weights, fitness);
    
    for(converge = 0; converge < optimize; ){
	// Elitism
	new_pop[0] = pop[best];

	// Evolve a new population
	for(int k = 1; k < N_POP; k++){
	    // Crossover
	    for(i = 0; i < 2; i++){
		idx = tournament(fitness);
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
	
	for(i = 0; i < N_POP; i++){
	    pop[i] = new_pop[i];
	}

	// Convergance
	p_best = fitness[best];
	best = pop_fitness(pop, dist, north, weights, fitness);
	
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

inline double get_weight(vector<int> &gifts, double weights[]){

    double weight = 0;

    for(int i = 0; i < gifts.size(); i++){
	weight += weights[gifts[i]];
    }

    return weight;

}

int main(){

    int i, j, k, m, p, h, progbar;
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

	for(count = 0; count < 67; count++){
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

	if(count != 67){
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

	    cost = tsp_genetic(route, coord, north, weights, 100);
	    
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

    double cost_a, cost_b;
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

	    cost_a = tsp_genetic(route_a, coord, north, weights, 10);
	    cost_b = tsp_genetic(route_b, coord, north, weights, 10);
	    
	    if(cost_a + cost_b < min_cost){
		min_cost = cost_a + cost_b;
		min_route_a = route_a;
		min_route_b = route_b;
	    }
	}

	if(min_cost < clusters[i].cost){
	    total_cost += min_cost - clusters[i].cost;
	    
	    Cluster cluster_a;
	    cluster_a.gifts = min_route_a;
	    cluster_a.cost = get_weariness(cluster_a.gifts, coord, north, weights);
	    cluster_a.weight = get_weight(cluster_a.gifts, weights);

	    clusters.push_back(cluster_a);

	    Cluster cluster_b;
	    cluster_b.gifts = min_route_b;
	    cluster_b.cost = get_weariness(cluster_b.gifts, coord, north, weights);
	    cluster_b.weight = get_weight(cluster_b.gifts, weights);
	    clusters.push_back(cluster_b);

	    clusters.erase(clusters.begin() + i);
	    i--;
	}

	if((i+1) / (clusters.size()/10.0) >= progbar){
	    fprintf(stderr, "|");
	    progbar++;
	}
    }

    fprintf(stderr, "\t%12.0lf\n", total_cost);
    
    // Swap
    for(i = 0; i < clusters.size(); i++){
	clusters[i].center[0] = 0;
	clusters[i].center[1] = 0;

	for(j = 0; j < clusters[i].gifts.size(); j++){
	    clusters[i].center[0] += coord[clusters[i].gifts[j]][0];
	    clusters[i].center[1] += coord[clusters[i].gifts[j]][1];
	}

	clusters[i].center[0] /= clusters[i].gifts.size();
	clusters[i].center[1] /= clusters[i].gifts.size();
    }

    int min_p, min_k;
    double dist[N_GIFTS];

    for(h = 0; h < 10; h++){
	fprintf(stderr, "S%02d ", h);

	vector<int> target_a, target_b, pos_a, pos_b;
	vector<double> swap_cost;
	
	progbar = 1;
	for(i = 0; i < clusters.size(); i++){
	    if(clusters[i].gifts.size() <= 1){
		printf("Warning!");
		continue;
	    }

	    // Distance to neighbours
	    for(j = 0; j < clusters.size(); j++){
		if(i != j){
		    dist[j] = haversine(clusters[i].center, clusters[j].center);
		} else {
		    dist[j] = 2 * M_PI * EARTH_RADIUS;
		}
		
		idx[j] = j;
	    }
	    
	    // Find the closest neighbours
	    nth_element(idx, idx+NEIGHBOURS, idx+clusters.size(), [&dist](int a, int b) 
			{return dist[a] < dist[b];});

	    for(m = 0; m < NEIGHBOURS; m++){
		j = idx[m];
		
		min_cost = clusters[i].cost + clusters[j].cost;

		for(k = 0; k < clusters[i].gifts.size(); k++){
		    if(clusters[j].weight + weights[clusters[i].gifts[k]] > MAX_WEIGHT){
			continue;
		    }

		    route_a = clusters[i].gifts;
		    route_a.erase(route_a.begin() + k);
		    
		    cost_a = get_weariness(route_a, coord, north, weights);
		    
		    for(p = 0; p < clusters[j].gifts.size(); p++){
			route_b = clusters[j].gifts;
			route_b.insert(route_b.begin() + p, clusters[i].gifts[k]);
			
			cost_b = get_weariness(route_b, coord, north, weights);
			
			if(cost_a + cost_b < min_cost){
			    min_cost = cost_a + cost_b;
			    min_k = k;
			    min_p = p;
			}
		    }
		}
	    
		if(min_cost < clusters[i].cost + clusters[j].cost){
		    swap_cost.push_back(min_cost - (clusters[i].cost + clusters[j].cost));
		    target_a.push_back(i);
		    target_b.push_back(j);
		    pos_a.push_back(min_k);
		    pos_b.push_back(min_p);
		}
	    }

	    if((i+1) / (clusters.size()/10.0) >= progbar){
		fprintf(stderr, "|");
		progbar++;
	    }
	}
	
	for(i = 0; i < swap_cost.size(); i++){
	    idx[i] = i;
	}
	
	sort(idx, idx+swap_cost.size(), [&swap_cost](int a, int b){
		return swap_cost[a] < swap_cost[b];
	    });

	for(i = 0; i < clusters.size(); i++){
	    used[i] = 0;
	}
	
	for(m = 0; m < max((int)swap_cost.size()>>2, 1); m++){
	    i = target_a[idx[m]];
	    j = target_b[idx[m]];
	    
	    if(used[i] || used[j]){
		continue;
	    }
	    
	    //printf("%lu %lu %d %d\n", clusters[i].gifts.size(), clusters[j].gifts.size(), 
	    //	   pos_a[idx[m]], pos_b[idx[m]]);
	    
	    clusters[j].gifts.insert(clusters[j].gifts.begin() + pos_b[idx[m]],
				     clusters[i].gifts[pos_a[idx[m]]]);

	    clusters[i].gifts.erase(clusters[i].gifts.begin() + pos_a[idx[m]]);
	    
	    clusters[i].weight = get_weight(clusters[i].gifts, weights);
	    clusters[j].weight = get_weight(clusters[j].gifts, weights);
	    
	    clusters[i].cost = tsp_genetic(clusters[i].gifts, coord, north, weights, 100);
	    clusters[j].cost = tsp_genetic(clusters[j].gifts, coord, north, weights, 100);
	    
	    used[i] = used[j] = 1;
	}

	total_cost = 0;
	for(i = 0; i < clusters.size(); i++){
	    total_cost += clusters[i].cost;
	}
	
	fprintf(stderr, "\t%12.0lf\n", total_cost);
    }

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

