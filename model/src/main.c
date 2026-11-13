// Implementation of the Bubble Model on the Macaque's virtual lung.
// Started on 11/23, Fernando Muñoz
#include "main.h"
#include "memory.h"
#include "simulation.h"
#include "exp_err.h"
#include "lung.h"
#include "params.h"
#include "aux.h"
#include "plot.h"
#include <mpi.h>


time_t start_time;
t_exec_config simu;
t_config_lung lung;
t_params_bm params;
t_exp_results exp_results;


t_memory_pool memory_pool;
char server_dir[PATH_MAX];   // -server <dir> (absolute path)


void initialize_MPI(int argc, char **argv, int *rank_MPI, int *size_MPI) 
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, rank_MPI);
    MPI_Comm_size(MPI_COMM_WORLD, size_MPI);
}


void read_args(int argc, char **argv)
{   
    simu.N_tot = 100;
    simu.is_human = 0; 
    simu.dseg = 1;
    simu.plot = 1;
    if (argc > 1) {
        for (int ii=1; ii<argc; ii++) {
            if (strcmp(argv[ii],"-N") == 0 ) {
              simu.N_tot = atoi(argv[ii+1]);
            //} else if ( strcmp(argv[ii], "-dseg") == 0 ) {
            //  dseg = 1;
            } else if ( strcmp(argv[ii], "-dlobe") == 0 ) {
                simu.dseg = 0;
                simu.dlobe = 1;
            } else if ( strcmp(argv[ii],"-debug_seq") == 0 ) {
                simu.debug_seq = 1;
            } else if ( strcmp(argv[ii], "-debug_par") == 0 ) {
                simu.debug_par = 1;
            } else if ( strcmp(argv[ii], "-pin") == 0 ) {
                simu.inner_parallel = 1;
            } else if ( strcmp(argv[ii], "-video") == 0 ) {
                simu.video = 1;
            } else if ( strcmp(argv[ii], "-write_sample") == 0 ) {
                simu.write_sample = atoi(argv[ii+1]);
            } else if ( strcmp(argv[ii], "-append") == 0 ) {
                simu.append = 1;
            } else if ( strcmp(argv[ii], "-human") == 0 ) {
                simu.is_human = 1;
            } else if (strcmp(argv[ii], "-noplot") == 0 ) {
                simu.plot = 0;
            } else if (strcmp(argv[ii], "-seed") == 0 ) {
                simu.seed = strtoull(argv[ii+1], NULL, 10);
                simu.seed_given = 1;
            } else if (strcmp(argv[ii], "-server") == 0 ) {
                if (ii+1 >= argc || !realpath(argv[ii+1], server_dir))
                    display_and_exit("ERROR!! -server <dir>: directory not found\n");
                simu.server = 1;
            }
        }
    }
}


void choose_seed(int rank_MPI)
{   // Without -seed: a seed from the clock and the process id (rank 0), the same on every rank
    if (!simu.seed_given && rank_MPI == 0) {
        struct timespec ts;
        clock_gettime(CLOCK_REALTIME, &ts);
        simu.seed = ((uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec) ^ ((uint64_t)getpid() << 40);
    }
    MPI_Bcast(&simu.seed, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
}


void write_sample_and_exit(int rank_MPI)
{
    t_params_bm *samples;
    int N_samples;
    if (rank_MPI == 0) {
      N_samples = read_samples(&samples);
      if (simu.write_sample > N_samples) {
        printf("ERROR!! Only read %d samples, trying to write number %d...\n",  N_samples, simu.write_sample);
      } else {
        write_params(samples[simu.write_sample - 1]);
        printf("-> Succesfully wrote sample %d into %s\n", simu.write_sample, FILE_PARAMS);
      }
    }
    MPI_Finalize();
    exit(0);
}


void display_headers(int rank_MPI)
{
    if (rank_MPI == 0) {
        if (simu.is_human == 1) {
            printf("-> Chosen HUMAN lung, ");
        } else {
            printf("-> Chosen MACAQUE lung, ");
        }

        if (simu.dseg == 1) {
            if (simu.is_human == 1) {
                printf("limiting lesions to SECONDARY LUNG LOBULES\n");
            } else {
                printf("limiting lesions to SEGMENTS\n");
            }
        } else if (simu.dlobe == 1 ) {
            printf("limiting lesions to LOBES\n");
        }

        printf("      SEED: %llu\n", (unsigned long long)simu.seed);
        printf("      PARAMS: v0:%.3f, fvimm:%.3f, r0:%.3f, frimm:%.3f, beta:%.3f, nata:%.3f, n:%.3f, a:%.3f N0:%.3f, tc:%.3f, rmin:%.3f  \n", \
            params.v0, params.fvimm, params.r0, params.frimm, params.b, params.nata, params.n, params.a, params.n0, params.tc, params.rmin);
    }
}


void display_error_message(t_limit_simu *limit_simu)
{
    puts("ERROR!! There has been an error in the simulations:");
    if (limit_simu->max_nc != 0) puts("        Reached maximum NC");
    if (limit_simu->max_nm != 0) puts("        Reached maximum NM");
    if (limit_simu->max_vfc != 0) puts("        Reached maximum VFC");
    if (limit_simu->max_time != 0) puts("        Reached maximum simulation time");
}


void blank_file_errors_if_needed(int rank_MPI) 
{
    if (rank_MPI == 0) {
      if (simu.append == 0) {
        printf("WARNING!! Blanked file %s\n", FILE_ERROR);
        FILE *f_error = fopen(FILE_ERROR, "w");
        fclose(f_error);
      }
    }
}


void schedule_work_MPI(int rank_MPI, int size_MPI)
{
    int workload_MPI = simu.N_tot / size_MPI;
    int remainder_MPI = simu.N_tot % size_MPI;
    int start_MPI = rank_MPI * workload_MPI;
    int end_MPI = start_MPI + workload_MPI;
    if (rank_MPI == size_MPI - 1) {
        end_MPI += remainder_MPI; // Adjust for remainder
    }
    simu.N_MPI = end_MPI - start_MPI;
    simu.start_MPI = start_MPI;
}


void check_omp_cancel(int rank_MPI)
{
    if(rank_MPI == 0 && omp_get_cancellation()) {
      puts("OMP cancellations are enabled! Unsure if this is better or worse performance.");
      puts("Untested.");
      // puts("Enable by launching with OMP_CANCELLATION=true");
      // puts("Enabling cancellation and rerunning program...\n\n");
      // putenv("OMP_CANCELLATION=true");
      // execv(argv[0], argv);
    }
}


int initialize_main_pointers(int size_MPI, t_results **results_tot, t_results **syncMPI_results, 
                             t_results **arrayMPI_results, t_limit_simu **limit_simu_tot,
                             t_limit_simu **syncMPI_limit_simu, t_limit_simu **arrayMPI_limit_simu)
{
    *results_tot = malloc(sizeof(t_results));
    if (!results_tot) return ERROR;
    initialize_results(*results_tot);

    *syncMPI_results = malloc(sizeof(t_results));
    if (!syncMPI_results) return ERROR;
    initialize_results(*syncMPI_results);

    *arrayMPI_results = malloc(sizeof(t_results) * size_MPI);
    if (!arrayMPI_results) return ERROR;
    for (int ii=0; ii<size_MPI; ii++) {
        initialize_results(&(*arrayMPI_results)[ii]);
    }

    *limit_simu_tot = malloc(sizeof(t_limit_simu));
    if (!limit_simu_tot) return ERROR;
    initialize_limit_simu(*limit_simu_tot);

    *syncMPI_limit_simu= malloc(sizeof(t_limit_simu));
    if (!syncMPI_limit_simu) return ERROR;
    initialize_limit_simu(*syncMPI_limit_simu);

    *arrayMPI_limit_simu = malloc(sizeof(t_limit_simu) * size_MPI);
    for (int ii=0; ii<size_MPI; ii++) {
        initialize_limit_simu(&(*arrayMPI_limit_simu)[ii]);
    }
    return OK;
}


typedef struct t_main_buffers {
    t_results *results_tot, *syncMPI_results, *arrayMPI_results;
    t_limit_simu *limit_simu_tot, *syncMPI_limit_simu, *arrayMPI_limit_simu;
} t_main_buffers;


// One evaluation: N_tot simulations of ./config/params_bm.txt in the current directory, results in
// ./results. The sequence of the one-shot program, shared by the one-shot path and the server.
int evaluation(int rank_MPI, int size_MPI, t_main_buffers *b)
{
    int out;
    blank_file_errors_if_needed(rank_MPI);

    out = read_params_and_update_lung();
    if (out == ERROR) return ERROR;
    out = read_experimental_results();
    if (out == ERROR) return ERROR;

    display_headers(rank_MPI);

    schedule_work_MPI(rank_MPI, size_MPI);

    out = initialize_rng_streams(simu.seed);
    if (out == ERROR) return ERROR;
    initialize_results(b->syncMPI_results);         // accumulated on rank 0: reset every evaluation
    initialize_limit_simu(b->syncMPI_limit_simu);

    int global_out = 0;

    // Timer
    struct timespec t_start, t_end;
    clock_gettime(CLOCK_MONOTONIC, &t_start);
    start_time = t_start.tv_sec;

    // Run simulations depending on the choices
    MPI_Barrier(MPI_COMM_WORLD);
    if ( simu.debug_par == 1 && size_MPI == 1) {  
        printf("-> Now running %d simulations... ", simu.N_tot);
        printf("      MODE: debug parallel\n");
        global_out = DEBUG_simulations_parallel(b->results_tot, b->limit_simu_tot);
    } 
    else if ( simu.debug_seq == 1 && size_MPI == 1) {  
        printf("-> Now running %d simulations... ", simu.N_tot);
        printf("      MODE: debug sequential\n");
        global_out = DEBUG_simulations_sequential(b->results_tot, b->limit_simu_tot);
    }
    else if ( simu.video == 1 && size_MPI == 1) {
        printf("      MODE: video (tmax is changed to 250)\n");
        //params.tmax = 250;
        global_out = make_video();
    }
    else {  
        if (rank_MPI == 0) {
            printf("-> Now running %d simulations... ", simu.N_tot);
            printf("  MODE: Parallel on %d processes, each with %d threads\n", size_MPI, omp_get_max_threads());
        }
        out = simulations_parallel(b->results_tot, b->limit_simu_tot);

        MPI_Gather(b->results_tot, sizeof(t_results), MPI_BYTE, b->arrayMPI_results, sizeof(t_results), MPI_BYTE, 0, MPI_COMM_WORLD);
        MPI_Gather(b->limit_simu_tot, sizeof(t_limit_simu), MPI_BYTE, b->arrayMPI_limit_simu, sizeof(t_limit_simu), MPI_BYTE, 0, MPI_COMM_WORLD);
        MPI_Reduce(&out, &global_out, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    }

    if (rank_MPI == 0) {
        if ( 0 == global_out ) {
            for (int jj = 0; jj < size_MPI; jj++) {
                sum_results(&b->arrayMPI_results[jj], b->syncMPI_results);
            }
            //average and write results
            average_results(b->syncMPI_results, simu.N_tot);
            write_results(b->syncMPI_results);
            calculate_and_write_errors(b->syncMPI_results);
            if (simu.plot) plot_results(b->syncMPI_results);
        } else {
            for (int jj = 0; jj < size_MPI; jj++) {
                merge_limit_simu(&b->arrayMPI_limit_simu[jj], b->syncMPI_limit_simu);
            }
            display_error_message(b->syncMPI_limit_simu);
            write_null_errors();
        }
    }

    if (rank_MPI == 0) {
        clock_gettime(CLOCK_MONOTONIC, &t_end);
        double elapsed_time = time_difference_seconds(t_end, t_start);
        printf("-> TIMING: %f s\n", elapsed_time);
    }
    return OK;
}


// SERVER MODE (-server <dir>): set up once, then one evaluation per request. Protocol (files, atomic
// rename): the client writes <dir>/req.<id>.tmp and renames it to <dir>/req.<id>, with lines
// "workdir <absolute path>", "n <N>", "seed <uint64>", "domain lobe|segment" (or the single line
// "stop"); rank 0 reads and deletes it, all ranks run the evaluation in workdir, and rank 0 answers
// with <dir>/done.<id> ("ok", "error" for an invalid request, "stopped"). Removing <dir> also stops
// the server. A fatal error on any rank aborts every rank (MPI_Abort): the client sees the server die.
#define REQ_PATH 4096  // request paths (same on macOS and Linux; matches the %4095 format below)
typedef struct t_request {
    int stop;
    int N_tot;
    int dlobe;
    int has_seed;
    uint64_t seed;
    char id[128];
    char workdir[REQ_PATH];
} t_request;


void write_response(const char *id, const char *status)
{
    char tmp[PATH_MAX + 160], fin[PATH_MAX + 160];
    snprintf(tmp, sizeof(tmp), "%s/done.%s.tmp", server_dir, id);
    snprintf(fin, sizeof(fin), "%s/done.%s", server_dir, id);
    FILE *f = fopen(tmp, "w");
    if (!f) return;  // server directory removed
    fprintf(f, "%s\n", status);
    fclose(f);
    rename(tmp, fin);
}


int read_request(t_request *req)
{   // Rank 0: wait for the next request (poll every 2 ms), read it and delete it. ERROR if invalid.
    struct timespec pause = {0, 2000000};
    char name[128 + 8];
    memset(req, 0, sizeof(t_request));
    while (1) {
        DIR *dir = opendir(server_dir);
        if (!dir) {  // server directory removed: stop
            req->stop = 1;
            return OK;
        }
        name[0] = '\0';
        struct dirent *entry;
        while ((entry = readdir(dir)) != NULL) {
            size_t len = strlen(entry->d_name);
            if (strncmp(entry->d_name, "req.", 4) == 0 && len > 4 && len < sizeof(name) &&
                !(len > 8 && strcmp(entry->d_name + len - 4, ".tmp") == 0)) {
                snprintf(name, sizeof(name), "%s", entry->d_name);
                break;
            }
        }
        closedir(dir);
        if (name[0] != '\0') break;
        nanosleep(&pause, NULL);
    }

    snprintf(req->id, sizeof(req->id), "%s", name + 4);
    char path[PATH_MAX + 160], key[32], value[REQ_PATH];
    snprintf(path, sizeof(path), "%s/%s", server_dir, name);
    FILE *f = fopen(path, "r");
    if (!f) return ERROR;
    req->dlobe = 1;
    while (fscanf(f, "%31s", key) == 1) {
        if (strcmp(key, "stop") == 0) { req->stop = 1; continue; }
        if (fscanf(f, " %4095[^\n]", value) != 1) break;
        if (strcmp(key, "workdir") == 0) snprintf(req->workdir, sizeof(req->workdir), "%s", value);
        else if (strcmp(key, "n") == 0) req->N_tot = atoi(value);
        else if (strcmp(key, "seed") == 0) { req->seed = strtoull(value, NULL, 10); req->has_seed = 1; }
        else if (strcmp(key, "domain") == 0) req->dlobe = (strcmp(value, "segment") != 0);
    }
    fclose(f);
    unlink(path);
    if (req->stop) return OK;
    if (req->workdir[0] == '\0' || req->N_tot <= 0 || !req->has_seed) return ERROR;
    return OK;
}


void server(int rank_MPI, int size_MPI, t_main_buffers *b)
{
    if (rank_MPI == 0) {
        printf("-> SERVER: %d processes x %d threads, waiting for requests in %s\n",
               size_MPI, omp_get_max_threads(), server_dir);
        fflush(stdout);
    }
    t_request req;
    while (1) {
        if (rank_MPI == 0) {
            while (1) {
                int valid = read_request(&req);
                if (req.stop) break;
                if (valid == OK && chdir(req.workdir) == 0) break;
                printf("ERROR!! SERVER: invalid request %s (workdir '%s', N %d, seed given %d)\n",
                       req.id, req.workdir, req.N_tot, req.has_seed);
                fflush(stdout);
                write_response(req.id, "error");
            }
        }
        MPI_Bcast(&req, sizeof(t_request), MPI_BYTE, 0, MPI_COMM_WORLD);
        if (req.stop) break;
        if (rank_MPI != 0 && chdir(req.workdir) != 0) {
            printf("ERROR!! SERVER: rank %d cannot enter %s\n", rank_MPI, req.workdir);
            fflush(stdout);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        simu.N_tot = req.N_tot;
        simu.seed = req.seed;
        simu.dlobe = req.dlobe;
        simu.dseg = !req.dlobe;
        if (rank_MPI == 0) printf("\n### SERVER request %s: %s, N = %d\n", req.id, req.workdir, req.N_tot);

        if (evaluation(rank_MPI, size_MPI, b) == ERROR) {
            printf("ERROR!! SERVER: evaluation failed on rank %d (request %s)\n", rank_MPI, req.id);
            fflush(stdout);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if (rank_MPI == 0) write_response(req.id, "ok");
        fflush(stdout);
    }
    if (rank_MPI == 0) {
        if (req.id[0] != '\0') write_response(req.id, "stopped");
        printf("-> SERVER: stopped\n");
    }
}


int main(int argc, char **argv) 
{
    int out;   

    int rank_MPI, size_MPI;
    initialize_MPI(argc, argv, &rank_MPI, &size_MPI);

    check_omp_cancel(rank_MPI);

    read_args(argc, argv);
    if (simu.write_sample != 0) write_sample_and_exit(rank_MPI);  // Special case

    // Set-up, once (also in server mode): lung, memory pool, result buffers
    out = malloc_and_read_lung();
    if (out == ERROR) return ERROR;
    out = initialize_memory_pool();
    if (out == ERROR) return ERROR;
    t_main_buffers b = {0};
    out = initialize_main_pointers(size_MPI, &b.results_tot, &b.syncMPI_results, &b.arrayMPI_results, \
                        &b.limit_simu_tot, &b.syncMPI_limit_simu, &b.arrayMPI_limit_simu);
    if (out == ERROR) return ERROR;

    if (simu.server) {
        server(rank_MPI, size_MPI, &b);
    } else {
        choose_seed(rank_MPI);
        out = evaluation(rank_MPI, size_MPI, &b);
        if (out == ERROR) return ERROR;
    }

    free_memory_pool();
    free_lung();
    if (rank_MPI == 0) {
        printf("-> FINISHED THE PROGRAM!!!\n");
    }
    free(b.results_tot); free(b.syncMPI_results); free(b.arrayMPI_results);
    free(b.limit_simu_tot); free(b.syncMPI_limit_simu); free(b.arrayMPI_limit_simu);

    MPI_Finalize();
    return OK;
}

