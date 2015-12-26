#include <cstdio>
#include <vector>
#include <map>
#include <algorithm>
#include <random>
#include <cmath>
#include <ctime>
#include <cstring>
#include <sys/stat.h>
using namespace std;

const bool RESTART = false;

const int MAX_WEIGHT = 1000;
const int SLEIGH_WEIGHT = 10;

const int N_GIFTS = 100000;
const int MAX_GIFTS = 200;

const double NORTH_POLE[2] = {0.5 * M_PI, 0};
const double EARTH_RADIUS = 6371.0;

const int N_POP = 100;
const int MUTATE = 60;
const int TOURNAMENT = 14;

const int MAGIC_GIFTS = 67;
const double MAGIC_LAT = -52.0 / 180.0 * M_PI;
const double MAGIC_LON = -105.0 / 180.0 * M_PI;

const int TSP_CONVERGE = 1000;
const int TSP_ATTEMPTS = 10;

const int SPLIT_CONVERGE = 10;
const int SPLIT_ATTEMPTS = 5;

const int NEIGHBOURS = 256;
const int MOVE_JITTER = 3;
const int SWAP_JITTER = 1;
const int USE_SWAP = true;
const int SWAP_CONVERGE = 100;
const int SWAP_ATTEMPTS = 3;

const double EPS = 1e-6;

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

    int i, best = 0;

    for(i = 0; i < N_POP; i++){
	fitness[i] = get_weariness(pop[i], dist, north, weights);

	if(fitness[i] < fitness[best]){
	    best = i;
	}
    }

    return best;

}

inline int tournament(double fitness[]){

    int i, best, idx;

    for(i = 0; i < TOURNAMENT; i++){
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
		   double all_weights[], int converge){
    
    int i, j, k, n_gifts = gifts.size();

    // Special optimization cases
    if(n_gifts == 0){
	return 0;
    } else if(n_gifts == 1){
	return all_north[gifts[0]] * (2 * SLEIGH_WEIGHT + all_weights[gifts[0]]);
    } else if(n_gifts == 2){
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
    int best, used[MAX_GIFTS], s, e, idx, c;
    double fitness[N_POP], p_best;
    
    best = pop_fitness(pop, dist, north, weights, fitness);
    
    c = 0;
    while(c < converge){
	// Elitism
	new_pop[0] = pop[best];

	// Evolve a new population
	for(k = 1; k < N_POP; k++){
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
	    c++;
	} else {
	    c = 0;
	}
    }

    // Reorder
    for(i = 0; i < n_gifts; i++){
	pop[best][i] = gifts[pop[best][i]];
    }

    gifts = pop[best];

    return fitness[best];
}

double tsp_genetic(vector<int> &gifts, double coord[][2], double north[], 
		   double weights[], int converge, int attempts){

    int i;
    double cost, min_cost;
    vector<int> route, min_route;

    for(i = 0; i < attempts; i++){
	route = gifts;
	
	cost = tsp_genetic(route, coord, north, weights, converge);
	    
	if(i == 0 || cost < min_cost){
	    min_cost = cost;
	    min_route = route;
	}
    }

    gifts = min_route;

    return min_cost;

}

inline double get_weight(vector<int> &gifts, double weights[]){

    int i;
    double weight = 0;

    for(i = 0; i < gifts.size(); i++){
	weight += weights[gifts[i]];
    }

    return weight;

}

void basic_routing(vector<Cluster> &clusters, double coord[][2], double north[],
		   double weights[], double &total_cost){

    fprintf(stderr, "\tCLT ");

    int i, j, count, total_used, progbar;
    int idx[N_GIFTS], region[N_GIFTS], used[N_GIFTS];

    for(i = 0; i < N_GIFTS; i++){
	idx[i] = i;
	used[i] = 0;
	
	region[i] = 0;
	if(coord[i][0] < MAGIC_LAT){
	    region[i] += 2;
	}

	if(coord[i][1] < MAGIC_LON){
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

	for(count = 0; count < MAGIC_GIFTS; count++){
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

	if(count != MAGIC_GIFTS){
	    break;
	}
    }

    fprintf(stderr, "\t%12.0lf\n", total_cost);

}

void split(vector<Cluster> &clusters, double coord[][2], double north[],
	   double weights[], double &total_cost){

    fprintf(stderr, "\tSPL ");

    int i, j, k, progbar;
    double cost_a, cost_b, min_cost;
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

	    cost_a = tsp_genetic(route_a, coord, north, weights, 
				 SPLIT_CONVERGE, SPLIT_ATTEMPTS);

	    cost_b = tsp_genetic(route_b, coord, north, weights,
				 SPLIT_CONVERGE, SPLIT_ATTEMPTS);
	    
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

}

int swap_gifts(vector<Cluster> &clusters, double coord[][2], double north[],
		  double weights[], int neighbours[][NEIGHBOURS], double &total_cost){

    fprintf(stderr, "\tSWP ");

    int i, j, k, p, m, t, progbar, target;
    int used[N_GIFTS], gift_cluster[N_GIFTS], gift_position[N_GIFTS];
    double dist[N_GIFTS], cost, cost_a, cost_b, cost_c;
    vector<int> route_a, route_b, route_c, target_a, target_b, pos_a, pos_b, op, nn, idx;
    vector<double> swap_cost;

    for(i = 0; i < clusters.size(); i++){
	for(j = 0; j < clusters[i].gifts.size(); j++){
	    gift_cluster[clusters[i].gifts[j]] = i;
	    gift_position[clusters[i].gifts[j]] = j;
	}
    }
    
    progbar = 1;
    for(i = 0; i < clusters.size(); i++){
	for(k = 0; k < clusters[i].gifts.size(); k++){
	    // Swap
	    route_a = clusters[i].gifts;
	    route_a.erase(route_a.begin() + k);
		
	    cost_a = get_weariness(route_a, coord, north, weights);

	    for(m = 0; m < NEIGHBOURS; m++){
		target = neighbours[clusters[i].gifts[k]][m];
		
		j = gift_cluster[target];
		p = gift_position[target];

		if(i == j){
		    continue;
		}
	    
		// Move
		if(clusters[j].weight + weights[clusters[i].gifts[k]] <= MAX_WEIGHT){
		    for(t = max(p - MOVE_JITTER, 0); t <= p + MOVE_JITTER + 1 && t <= clusters[j].gifts.size(); t++){
			route_b = clusters[j].gifts;
			route_b.insert(route_b.begin() + t, clusters[i].gifts[k]);
			
			cost_b = get_weariness(route_b, coord, north, weights);
			cost = cost_a + cost_b - (clusters[i].cost + clusters[j].cost);
			
			if(cost < 0){
			    swap_cost.push_back(cost);
			    target_a.push_back(i);
			    target_b.push_back(j);
			    pos_a.push_back(k);
			    pos_b.push_back(t);
			    nn.push_back(m);
			    op.push_back(0);
			}
		    }
		}
		
		// Swap
		if(USE_SWAP){
		    for(t = max(p - SWAP_JITTER, 0); t <= p + SWAP_JITTER && t < clusters[j].gifts.size(); t++){
			if(clusters[i].weight - weights[clusters[i].gifts[k]] + weights[clusters[j].gifts[t]] > MAX_WEIGHT ||
			   clusters[j].weight - weights[clusters[j].gifts[t]] + weights[clusters[i].gifts[k]] > MAX_WEIGHT){
			    continue;
			}
			
			route_c = clusters[i].gifts;
			route_b = clusters[j].gifts;
			
			swap(route_c[k], route_b[t]);
			
			cost_c = get_weariness(route_c, coord, north, weights);
			cost_b = get_weariness(route_b, coord, north, weights);
			
			cost = cost_c + cost_b - (clusters[i].cost + clusters[j].cost);
			
			if(cost < 0){
			    swap_cost.push_back(cost);
			    target_a.push_back(i);
			    target_b.push_back(j);
			    pos_a.push_back(k);
			    pos_b.push_back(t);
			    nn.push_back(m);
			    op.push_back(1);
			}
		    }
		}
	    }
	}
	
	if((i+1) / (clusters.size()/10.0) >= progbar){
	    fprintf(stderr, "|");
	    progbar++;
	}
    }
    
    for(i = 0; i < swap_cost.size(); i++){
	idx.push_back(i);
    }
    
    sort(idx.begin(), idx.end(), [&swap_cost, &op](int a, int b){
	    return swap_cost[a] < swap_cost[b];
	});
    
    for(i = 0; i < clusters.size(); i++){
	used[i] = 0;
    }

    int nn_hist[NEIGHBOURS], op_hist[2];
    double decrease;

    op_hist[0] = op_hist[1] = 0;
    for(i = 0; i < NEIGHBOURS; i++){
	nn_hist[i] = 0;
    }
    
    decrease = 0;
    for(m = 0; m < swap_cost.size(); m++){
	k = idx[m];
	i = target_a[k];
	j = target_b[k];

	if(used[i] || used[j]){
	    continue;
	}
	
	decrease -= clusters[i].cost + clusters[j].cost;
	
	if(op[k] == 0){
	    clusters[j].gifts.insert(clusters[j].gifts.begin() + pos_b[k],
				     clusters[i].gifts[pos_a[k]]);
	    
	    clusters[i].gifts.erase(clusters[i].gifts.begin() + pos_a[k]);
	} else {
	    swap(clusters[i].gifts[pos_a[k]], clusters[j].gifts[pos_b[k]]);
	}
	
	op_hist[op[k]]++;

	clusters[i].weight = get_weight(clusters[i].gifts, weights);
	clusters[j].weight = get_weight(clusters[j].gifts, weights);
	
	clusters[i].cost = tsp_genetic(clusters[i].gifts, coord, north, weights,
				       SWAP_CONVERGE, SWAP_ATTEMPTS);
	
	clusters[j].cost = tsp_genetic(clusters[j].gifts, coord, north, weights,
				       SWAP_CONVERGE, SWAP_ATTEMPTS);
	
	decrease += clusters[i].cost + clusters[j].cost; 
	nn_hist[nn[k]]++;
	
	used[i] = used[j] = 1;
    }
    
    total_cost = 0;
    for(i = 0; i < clusters.size(); i++){
	if(clusters[i].gifts.size() == 0){
	    clusters.erase(clusters.begin() + i);
	}
	
	total_cost += clusters[i].cost;
    }
    
    fprintf(stderr, "\t%12.0lf %10.0lf |", total_cost, decrease);

    fprintf(stderr, " OP: %d / %d | N:", op_hist[0], op_hist[1]);
    for(i = 0; i < NEIGHBOURS; i++){
	fprintf(stderr, " %d", nn_hist[i]);
    }
    fprintf(stderr, "\n");

    return swap_cost.size();

}

void output_result(vector<Cluster> &clusters, char fname[]){

    int i, j;
    FILE *fp_result = fopen(fname, "w");
    
    fprintf(fp_result, "TripId,GiftId\n");
    for(i = 0; i < clusters.size(); i++){
	for(j = 0; j < clusters[i].gifts.size(); j++){
	    fprintf(fp_result, "%d,%d\n", i, clusters[i].gifts[j] + 1);
	}
    }

    fclose(fp_result);

}

int main(){

    int i, j, k, progbar;
    double coord[N_GIFTS][2], weights[N_GIFTS], north[N_GIFTS];
    vector<Cluster> clusters;

    generator.seed(time(NULL));

    // Input
    char fname[64];
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

    mkdir("submission", 0775);

    // Initialize
    double total_cost;

    if(RESTART){
	total_cost = 0;
	for(i = 0; i < N_GIFTS; i++){
	    total_cost += north[i] * (2 * SLEIGH_WEIGHT + weights[i]);
	}
    } else {
	int gift, r, prev_r;
	Cluster cluster;
	FILE *fp_best = fopen("basic.csv", "r");

	fscanf(fp_best, "%*s\n");
	for(i = 0; i < N_GIFTS; i++){
	    fscanf(fp_best, "%d,%d", &r, &gift);
	    gift--;

	    if(i != 0 && r != prev_r || i == N_GIFTS - 1){
		if(i == N_GIFTS - 1){
		    cluster.gifts.push_back(gift);
		}

		clusters.push_back(cluster);
		cluster.gifts.clear();
	    }

	    cluster.gifts.push_back(gift);

	    prev_r = r;
	}

	total_cost = 0;
	for(i = 0; i < clusters.size(); i++){
	    clusters[i].cost = get_weariness(clusters[i].gifts, coord, north, weights);
	    clusters[i].weight = get_weight(clusters[i].gifts, weights);
	    total_cost += clusters[i].cost;
	}
    }

    fprintf(stderr, "\tSTR\t\t%12.0lf\n", total_cost);

    if(RESTART){
	// Initial routing
	basic_routing(clusters, coord, north, weights, total_cost);

	// Optimize using TSP
	fprintf(stderr, "\tTSP ");
	
	progbar = 1;
	total_cost = 0;
	for(i = 0; i < clusters.size(); i++){
	    clusters[i].cost = tsp_genetic(clusters[i].gifts, coord, north, weights,
					   TSP_CONVERGE, TSP_ATTEMPTS);

	    total_cost += clusters[i].cost;

	    if((i+1) / (clusters.size()/10.0) >= progbar){
		fprintf(stderr, "|");
		progbar++;
	    }
	}

	fprintf(stderr, "\t%12.0lf\n", total_cost);
	
	// Split unoptimal clusters
	split(clusters, coord, north, weights, total_cost);

	sprintf(fname, "basic.csv");
	output_result(clusters, fname);
    }

    // Swap gifts among the routes
    static int neighbours[N_GIFTS][NEIGHBOURS];
    
    sprintf(fname, "nn%d", NEIGHBOURS);
    FILE *fp_nn = fopen(fname, "r");
    
    if(fp_nn == NULL){
	fprintf(stderr, "\tNN  ");

	int idx[N_GIFTS];
	double dist[N_GIFTS];

	progbar = 1;
	fp_nn = fopen(fname, "w");
	for(i = 0; i < N_GIFTS; i++){
	    for(j = 0; j < N_GIFTS; j++){
		if(i != j){
		    dist[j] = haversine(coord[i], coord[j]);
		} else {
		    dist[j] = 2 * M_PI * EARTH_RADIUS;
		}
		
		idx[j] = j;
	    }

	    sort(idx, idx+N_GIFTS, [&dist](int a, int b){
		    return dist[a] < dist[b];
		});

	    for(j = 0; j < NEIGHBOURS; j++){
		neighbours[i][j] = idx[j];
		fprintf(fp_nn, "%d ", idx[j]);
	    }
	    fprintf(fp_nn, "\n");
	    
	    if((i+1) / (N_GIFTS/10.0) >= progbar){
		fprintf(stderr, "|");
		progbar++;
	    }
	}
	fprintf(stderr, "\n");
    } else {
	for(i = 0; i < N_GIFTS; i++){
	    for(j = 0; j < NEIGHBOURS; j++){
		fscanf(fp_nn, "%d", &neighbours[i][j]);
	    }
	}
    }

    for(i = 0; ; i++){
	fprintf(stderr, "%d", i);
	if(swap_gifts(clusters, coord, north, weights, neighbours, total_cost) == 0){
	    break;
	}
	
	// Output result
	sprintf(fname, "submission/%03d.csv", i);
	output_result(clusters, fname);
    }
	
    return 0;
    
 }

