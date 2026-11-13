#include "main.h"
#include "memory.h"
#include "lung.h"


int initialize_memory_pool(void) 
{
    int size_OMP = omp_get_max_threads();

    memory_pool.results = malloc(sizeof(t_results**) * size_OMP);
    memory_pool.lesions = malloc(sizeof(t_results**) * size_OMP);
    memory_pool.clusters = malloc(sizeof(t_cluster**) * size_OMP);
    memory_pool.cluster_summary = malloc(sizeof(t_cluster_summary*) * size_OMP);
    memory_pool.cluster_labels = malloc(sizeof(int*) * size_OMP);
    memory_pool.cumsum = malloc(sizeof(double*) * size_OMP);
    memory_pool.UF = malloc(sizeof(int**) * size_OMP);
    memory_pool.limit_simu = malloc(sizeof(t_limit_simu*) * size_OMP);

    if (!memory_pool.results || !memory_pool.results || \
        !memory_pool.clusters || !memory_pool.cluster_summary || \
        !memory_pool.cluster_labels || !memory_pool.cumsum || \
        !memory_pool.UF || !memory_pool.limit_simu) 
        return ERROR;

    for (int ii=0; ii<size_OMP; ii++) {
        memory_pool.results[ii] = malloc(sizeof(t_results *) * 2);
        if (!memory_pool.results[ii]) return ERROR;
        memory_pool.results[ii][0] = malloc(sizeof(t_results));
        memory_pool.results[ii][1] = malloc(sizeof(t_results));
        if (!memory_pool.results[ii][0] || !memory_pool.results[ii][1]) return ERROR;

        memory_pool.lesions[ii] = malloc(sizeof(t_results *) * 3);
        if (!memory_pool.lesions[ii]) return ERROR;
        memory_pool.lesions[ii][0] = malloc(sizeof(t_lesions));
        memory_pool.lesions[ii][1] = malloc(sizeof(t_lesions));
        memory_pool.lesions[ii][2] = malloc(sizeof(t_lesions));
        if (!memory_pool.lesions[ii][0] || !memory_pool.lesions[ii][1] || !memory_pool.lesions[ii][2]) return ERROR;

        memory_pool.clusters[ii] = malloc(sizeof(t_cluster *) * 2);
        if (!memory_pool.clusters[ii]) return ERROR;
        memory_pool.clusters[ii][0] = malloc(sizeof(t_cluster));
        memory_pool.clusters[ii][1] = malloc(sizeof(t_cluster));
        if (!memory_pool.clusters[ii][0] || !memory_pool.clusters[ii][1]) return ERROR;

        memory_pool.cluster_summary[ii] = malloc(sizeof(t_cluster_summary) * MAX_LESIONS);
        if (!memory_pool.cluster_summary[ii]) return ERROR;

        memory_pool.cluster_labels[ii] = malloc(sizeof(int) * MAX_LESIONS);
        if (!memory_pool.cluster_labels[ii]) return ERROR;

        memory_pool.cumsum[ii] = malloc(sizeof(double) * lung.N_terminals);
        if (!memory_pool.cumsum[ii]) return ERROR;

        memory_pool.UF[ii] = malloc(sizeof(int*) * 2); 
        if (!memory_pool.UF[ii]) return ERROR;
        memory_pool.UF[ii][0] = malloc(sizeof(int) * MAX_LESIONS);
        memory_pool.UF[ii][1] = malloc(sizeof(int) * MAX_LESIONS);
        if (!memory_pool.UF[ii][0] || !memory_pool.UF[ii][1]) return ERROR;

        if(!memory_pool.limit_simu) return ERROR;
        memory_pool.limit_simu[ii] = malloc(sizeof(t_limit_simu));
        if (!memory_pool.limit_simu[ii]) return ERROR;
    }

    memory_pool.rng = aligned_alloc(_Alignof(t_rng_slot), sizeof(t_rng_slot) * size_OMP);
    if (!memory_pool.rng) return ERROR;
    memory_pool.rng_sim = NULL;
    memory_pool.N_rng_sim = 0;
    return OK;
}


void free_memory_pool(void)
{
    int size_OMP = omp_get_max_threads();
    for (int ii=0; ii<size_OMP; ii++) {
        free(memory_pool.results[ii][0]);
        free(memory_pool.results[ii][1]);
        free(memory_pool.results[ii]);
        
        free(memory_pool.lesions[ii][0]);
        free(memory_pool.lesions[ii][1]);
        free(memory_pool.lesions[ii][2]);
        free(memory_pool.lesions[ii]);

        free(memory_pool.clusters[ii][0]);
        free(memory_pool.clusters[ii][1]);
        free(memory_pool.clusters[ii]);

        free(memory_pool.cluster_summary[ii]);

        free(memory_pool.cluster_labels[ii]);

        free(memory_pool.cumsum[ii]);

        free(memory_pool.UF[ii][0]);
        free(memory_pool.UF[ii][1]);
        free(memory_pool.UF[ii]);

        free(memory_pool.limit_simu[ii]);
    }

    free(memory_pool.results);
    free(memory_pool.lesions);
    free(memory_pool.clusters);
    free(memory_pool.cluster_summary);
    free(memory_pool.cluster_labels);
    free(memory_pool.cumsum);
    free(memory_pool.UF);
    free(memory_pool.limit_simu);
    free(memory_pool.rng);
    free(memory_pool.rng_sim);
}





void free_lung(void)
{
    for (int ii=0; ii<lung.N_lobes; ii++) {
        free(lung.lobes[ii].faces);
    }
    free(lung.lobes);
    free(lung.terminals);
    free(lung.distance);
}





void request_memory_lesions(t_lesions **lesions_1, t_lesions **lesions_2, t_lesions **lesions_3)
{
    int OMP_id = omp_get_thread_num();
    if (lesions_1) *lesions_1 = memory_pool.lesions[OMP_id][0];
    if (lesions_2) *lesions_2 = memory_pool.lesions[OMP_id][1];
    if (lesions_3) *lesions_3 = memory_pool.lesions[OMP_id][2];
}


void request_memory_results(t_results **results_1, t_results **results_2)
{
    int OMP_id = omp_get_thread_num();
    if (results_1) *results_1 = memory_pool.results[OMP_id][0];
    if (results_2) *results_2 = memory_pool.results[OMP_id][1];
}


void request_memory_UF(int **arr_1, int **arr_2)
{
    int OMP_id = omp_get_thread_num();
    if (arr_1) *arr_1 = memory_pool.UF[OMP_id][0];
    if (arr_2) *arr_2 = memory_pool.UF[OMP_id][1];
}


void request_memory_clusters(t_cluster **cluster_1, t_cluster **cluster_2)
{
    int OMP_id = omp_get_thread_num();
    if (cluster_1) *cluster_1 = memory_pool.clusters[OMP_id][0];
    if (cluster_2) *cluster_2 = memory_pool.clusters[OMP_id][1];
}


void request_memory_cluster_summary(t_cluster_summary **cluster_summary)
{
    int OMP_id = omp_get_thread_num();
    if (cluster_summary) *cluster_summary = memory_pool.cluster_summary[OMP_id];
}


void request_memory_cluster_labels(int **cluster_labels)
{
    int OMP_id = omp_get_thread_num();
    if (cluster_labels) *cluster_labels = memory_pool.cluster_labels[OMP_id];
}


void request_memory_limit_simu(t_limit_simu **limit_simu)
{
    int OMP_id = omp_get_thread_num();
    if (limit_simu) *limit_simu = memory_pool.limit_simu[OMP_id];
}


void request_memory_cumsum(double **cumsum)
{
    int OMP_id = omp_get_thread_num();
    if (cumsum) *cumsum = memory_pool.cumsum[OMP_id];
}


// RANDOM NUMBERS (random.h, xoshiro256**): one stream per simulation. Simulation j of the
// evaluation (global index j = start_MPI + ii, 0 <= j < N_tot) starts from the seed's state
// advanced by j jumps of 2^128 draws, i.e. context j of random_initialize_threads(seed, 1, N_tot):
// non-overlapping streams, and results depend only on (seed, N_tot, params), not on the number of
// ranks/threads or on which thread runs which simulation.
int initialize_rng_streams(uint64_t seed)
{
    int N = simu.N_MPI > 0 ? simu.N_MPI : 1;
    if (N > memory_pool.N_rng_sim) {
        s_random_context *aux = realloc(memory_pool.rng_sim, sizeof(s_random_context) * N);
        if (!aux) return ERROR;
        memory_pool.rng_sim = aux;
        memory_pool.N_rng_sim = N;
    }

    s_random_context ctx = random_initialize(seed);
    for (int jj=0; jj<simu.start_MPI; jj++) XOSHIRO256_jump(&ctx);
    for (int ii=0; ii<simu.N_MPI; ii++) {
        memory_pool.rng_sim[ii] = ctx;
        XOSHIRO256_jump(&ctx);
    }

    // Default stream of each thread (draws outside a simulation stream, e.g. -video): 2^192 ahead
    // of the simulation streams, then 2^128 apart.
    ctx = random_initialize(seed);
    XOSHIRO256_long_jump(&ctx);
    for (int tt=0; tt<omp_get_max_threads(); tt++) {
        memory_pool.rng[tt].ctx = ctx;
        XOSHIRO256_jump(&ctx);
    }
    return OK;
}


void use_rng_stream(int ii)
{   // The calling thread runs simulation ii (local index) of this rank: copy its initial stream
    memory_pool.rng[omp_get_thread_num()].ctx = memory_pool.rng_sim[ii];
}


double bm_uniform(void)
{   // Uniform in [0, 1), from the stream of the calling thread
    return random_uniform_double(&memory_pool.rng[omp_get_thread_num()].ctx);
}



