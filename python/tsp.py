import numpy as np
from random import random, randint, shuffle

def get_weariness(route, dist, weights):
    weariness = 0
    for i in range(0,len(route)):
        weariness += dist[route[i]][route[i-1]]

    return weariness

def pop_fitness(pop, dist, weights):
    fitness = []
    
    for i,trip in enumerate(pop):
        fitness.append(get_weariness(trip, dist, weights))
    
    return np.array(fitness)

def tournament(fitness, size):
    idx = []

    for i in range(0,size):
        idx.append(randint(0, len(fitness)-1))
        
    return idx[np.argmin(fitness[idx])]


def tsp_genetic(coord, weights, labels, n_pop=100, gen=100, mutation=0.4):
    n_gifts = coord.shape[0]

    # Precompute distance matrix
    dist = [[0]*n_gifts for i in range(0,n_gifts)]
    for i in range(0,n_gifts):
        for j in range(i+1,n_gifts):
            dist[i][j] = dist[j][i] = np.linalg.norm(coord[i] - coord[j])
    
    # Initialize population
    pop = []
    for i in range(0,n_pop):
        trip = range(0,n_gifts)
        shuffle(trip)
        pop.append(trip)

    # Tournament size
    n_champions = int(n_pop / 14 + 0.5)

    # Generations
    for g in range(0,gen):
        fitness = pop_fitness(pop, dist, weights)
        best = np.argmin(fitness)
        
        # Evolve a new population
        new_pop = [pop[best]]
        for k in range(0,n_pop-1):
            # Crossover
            parents = (pop[tournament(fitness, n_champions)],
                       pop[tournament(fitness, n_champions)])

            s = randint(0, n_gifts-1)
            e = randint(s, n_gifts-1) + 1

            child = [0]*n_gifts
            present = [0]*n_gifts
            for i in range(s,e):
                present[parents[0][i]] = 1
                child[i] = parents[0][i]

            idx = 0
            for i in range(0,s) + range(e,n_gifts):
                while present[parents[1][idx]]:
                    idx += 1

                child[i] = parents[1][idx]
                idx += 1

            # Mutate
            if random() <= mutation:
                i = randint(0,n_gifts-1)
                j = randint(0,n_gifts-1)

                child[i], child[j] = child[j], child[i]

            new_pop.append(child)

        pop = new_pop

    fitness = pop_fitness(pop, dist, weights)
    best = np.argmax(fitness)
        
    return (fitness[best], pop[best])
