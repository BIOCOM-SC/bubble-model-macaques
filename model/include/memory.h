#ifndef BM_MEMORY_H
#define BM_MEMORY_H

#include "lesions.h"
#include "clusters.h"
#include "lesions.h"
#include "results.h"
#include "simulation.h"
#include "random.h"


// Random-number state of one thread, on its own cache line (no false sharing between threads)
typedef struct t_rng_slot {
    _Alignas(64) s_random_context ctx;
} t_rng_slot;


typedef struct memory_pool {
    t_results ***results;                   // x2 per thread
    t_lesions ***lesions;                   // x3 per thread
    t_cluster ***clusters;                  // x2 per thread
    t_cluster_summary **cluster_summary;    // x1 per thread
    int **cluster_labels;                   // x1 per thread
    double **cumsum;                        // x1 per thread
    int ***UF;                              // x2 per thread
    t_limit_simu **limit_simu;              // x1 per thread
    t_rng_slot *rng;                        // x1 per thread: stream of the running simulation
    s_random_context *rng_sim;              // initial stream of each simulation of this rank
    int N_rng_sim;                          // allocated length of rng_sim
} t_memory_pool;


extern t_memory_pool memory_pool;

int initialize_memory_pool(void);
void free_memory_pool(void);


void free_lung(void);  // TODO put this in corresponding module?

void request_memory_lesions(t_lesions **lesions_1, t_lesions **lesions_2, t_lesions **lesions_3);
void request_memory_results(t_results **results_1, t_results **results_2);
void request_memory_UF(int **arr_1, int **arr_2);
void request_memory_clusters(t_cluster **cluster_1, t_cluster **cluster_2);
void request_memory_cluster_summary(t_cluster_summary **cluster_summary);
void request_memory_cluster_labels(int **cluster_labels);
void request_memory_limit_simu(t_limit_simu **limit_simu);
void request_memory_cumsum(double **cumsum);

int initialize_rng_streams(uint64_t seed);
void use_rng_stream(int ii);
double bm_uniform(void);

#endif
