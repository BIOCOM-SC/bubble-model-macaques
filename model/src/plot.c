
#include <math.h>
#include "main.h"
#include "results.h"
#include "params.h"
#include "fusion.h"
#include "reinfection.h"
#include "aux.h"
#include "memory.h"
#include "exp_err.h"
#include "gplot.h"


enum plot_variable {
    NC,
    NM,
    VFC,
    NMD,
    NMI,
    MAC
};


static void add_experimental(s_gplot *interface, enum plot_variable var)
{
    double exp_days[3] = {21, 53, 77};
    double exp[3][3];
    for (int ii=0; ii<3; ii++) {
        exp[0][ii] = exp_days[ii];
        switch (var) {
            case NC: 
                exp[1][ii] = exp_results.nc[ii];
                exp[2][ii] = exp_results.nc_std[ii];
                break;
            case NM:
                exp[1][ii] = exp_results.nm[ii];
                exp[2][ii] = exp_results.nm_std[ii];
                break;
            case VFC:
                exp[1][ii] = exp_results.vfc[ii];
                exp[2][ii] = exp_results.vfc_std[ii];
                break;
            case NMD:
                exp[1][ii] = exp_results.nmd[ii];
                exp[2][ii] = exp_results.nmd_std[ii];
                break;
            case NMI:
                exp[1][ii] = exp_results.nmi[ii];
                exp[2][ii] = exp_results.nmi_std[ii];
                break;
            case MAC:
                exp[1][ii] = exp_results.mac[ii];
                exp[2][ii] = exp_results.mac_std[ii];
                break;
        }
    }
    // Errorbars:
    gplot_draw_errorbars_2d(interface, 3, exp[0], exp[1], exp[2], .extra = "lc 8 lw 2 lt -1");
    // Points:
    gplot_draw_2d(interface, 3, exp[0], exp[1], .with = GPLOT_POINTS, .extra = "pt 5 lc 8 lw 4");
}


static void add_array_errorbars(s_gplot *interface, double *mean, double *std, int N, int skip, const char *config)
{
    for (int ii=0; ii<N; ii++) {
        if (ii % skip == 0) {
            double x = (double) ii;
            gplot_draw_errorbars_2d(interface, 1, &x, &mean[ii], &std[ii], .extra = config);
        }
    }
}


void plot_results(t_results *results)
{
    int N = (int)params.tmax;
    double days[N];
    for (int ii=0; ii<N; ii++) days[ii] = ii;  // index ii (as before; was an off-by-one loop writing days[N])

    char buff[1024];
    // gplot opens the 2x3 multiplot itself, so the overall title is a screen label (removed after the first panel)
    snprintf(buff, 1024, "set label 1 'v0:%.3f, fvimm:%.3f, r0:%.3f, frimm:%.3f, ß:%.3f, ρ:%.3f, n:%.3f, a:%.3f, N0:%.3f, tc:%.3f, rmin:%.3f' at screen 0.5,0.985 center font ',24'", \
          params.v0, params.fvimm, params.r0, params.frimm, params.b, params.nata, params.n, params.a, params.n0, params.tc, params.rmin);

    s_gplot *interface = gplot_start(GPLOT_PNG_2D, "./results/main_results.png",
                                     .size = {1920, 1080}, .font_size = 18, .layout = {2, 3});
    gplot_config(interface, buff, "set tmargin 4", "set grid", "set style line 1 lc rgb '#c8c8c8' lt 1 lw 1");
    
    // NC
    gplot_config(interface, "set xlabel 'time (days)'", "set ylabel 'Number'", 
                   "set title '# of consolidations (NC)'", "set yrange [0:30]");
    gplot_draw_2d(interface, N, days, results->nc.t1, .extra = "lc 6 lw 3");
    add_array_errorbars(interface, results->nc.t1, results->nc.t2, N, 5, "lc 3 lt -1");
    add_experimental(interface, NC);

    // NM
    gplot_next_subplot(interface, "unset label 1", "set xlabel 'time (days)'", "set ylabel 'Number'", 
                   "set title '# of micronodules (NM)'", "set yrange [0:150]");
    gplot_draw_2d(interface, N, days, results->nm.t1, .extra = "lc 6 lw 3");
    add_array_errorbars(interface, results->nm.t1, results->nm.t2, N, 5, "lc 3 lt -1");
    add_experimental(interface, NM);

    // VFC
    gplot_next_subplot(interface, "set xlabel 'time (days)'", "set ylabel 'Lung volume fraction'",
                    "set title 'Vol. fract. occupied by consol.'", "set yrange [0:0.14]");
    gplot_draw_2d(interface, N, days, results->vfc.t1, .extra = "lc 6 lw 3");
    add_array_errorbars(interface, results->vfc.t1, results->vfc.t2, N, 5, "lc 3 lt -1");
    add_experimental(interface, VFC);

    // NMD
    gplot_next_subplot(interface, "set xlabel 'time (days)'", "set ylabel 'Number'",
                "set title '# of daughter micronodules (NMD)'", "set yrange [0:100]");
    gplot_draw_2d(interface, N, days, results->nmd.t1, .extra = "lc 6 lw 3");
    add_array_errorbars(interface, results->nmd.t1, results->nmd.t2, N, 5, "lc 3 lt -1");
    add_experimental(interface, NMD);

    // NMI
    gplot_next_subplot(interface, "set xlabel 'time (days)'", "set ylabel 'Number'",
                    "set title '# of isolated micronodules (NMI)'", "set yrange [0:100]");
    gplot_draw_2d(interface, N, days, results->nmi.t1, .extra = "lc 6 lw 3");
    add_array_errorbars(interface, results->nmi.t1, results->nmi.t2, N, 5, "lc 3 lt -1");
    add_experimental(interface, NMI);

    // MAC
    gplot_next_subplot(interface, "set xlabel 'time (days)'", "set ylabel 'Main axis (mm)'",
                    "set title 'Main axis'", "set yrange [0:22]");
    gplot_draw_2d(interface, N, days, results->mac.t1, .extra = "lc 6 lw 3", .title = "MAC");
    gplot_draw_2d(interface, N, days, results->ma.t1, .extra = "lc 4 lw 3", .title = "MA(C+M)");
    add_array_errorbars(interface, results->mac.t1, results->mac.t2, N, 5, "lc 3 lt -1");
    add_experimental(interface, MAC);

    gplot_end(interface);


    /* Reinfections */
    interface = gplot_start(GPLOT_PNG_2D, "./results/reinfections.png", .size = {1920, 1080}, .font_size = 18);
    gplot_config(interface, "set yrange [0:180]", "set grid", "set style line 1 lc rgb '#c8c8c8' lt 1 lw 1");
    gplot_draw_2d(interface, N, days, results->nintra.t1, .extra = "lc 4 lw 3", .title = "intra");
    gplot_draw_2d(interface, N, days, results->nextra.t1, .extra = "lc 6 lw 3", .title = "extra");
    gplot_end(interface);
}


static void add_lung(s_gplot *interface)
{
    char dir_name[300];
    if (simu.is_human == 1) {
        snprintf(dir_name, 100, "%s/lobes", DIR_CONFIG_HUMAN);
    } else {
        if (simu.dseg == 1) {
            snprintf(dir_name, 100, "%s/segments", DIR_CONFIG_MACAQUE);
        } else {
            snprintf(dir_name, 100, "%s/lobes", DIR_CONFIG_MACAQUE);
        }
    }
  
    DIR *dir;
    dir = opendir(dir_name);
    if (!dir) display_and_exit("ERROR!! Unable to open directory (add_lung_as_data_blocks)");

    struct dirent *entry;
    entry = readdir(dir);
    int num_segments = 0;
    while (entry) {
        if (entry->d_type == DT_REG) {
            if (ends_with_txt(entry->d_name)) {
                char buff[256];
                snprintf(buff, 256, "DATA_%d", num_segments);
                char fname[256];
                snprintf(fname, 256, "%s/%s", dir_name, entry->d_name);
                gplot_draw_file(interface, fname, "skip 1", .with = GPLOT_POLYGONS, .extra = "fs transparent solid 0.1 fc 'salmon'");
                // add_datablock_from_file(interface, fname, 1, buff);
                num_segments++;
            }
        }
        entry = readdir(dir);
    }
    closedir(dir);
}


void add_lesions(s_gplot *interface, t_lesions *lesions)
{
    for (int ii=0; ii<lesions->N; ii++) {
        if (lesions->arr[ii].r > params.rmin) {
            gplot_draw_sphere_3d(interface, lesions->arr[ii].pos.x, lesions->arr[ii].pos.y, 
                                 lesions->arr[ii].pos.z, lesions->arr[ii].r,
                                 .with = GPLOT_POLYGONS, .extra = "fc 'black'");
        }
    }
}


void save_photogram(s_gplot *interface, t_lesions *lesions, int day)
{   
    char buffer[1024];
    static double rotation_z = 80;
    double scale = 1.4;

    rotation_z = fmod((rotation_z + 0.3), 360);
    snprintf(buffer, 1024, "set view 70, %f, %f, 1", rotation_z, scale);
    gplot_config(interface, buffer);
    snprintf(buffer, 1024, "set title 'Day %d' font ',48' offset 0,-7", day);
    gplot_config(interface, buffer);
    
    add_lung(interface);
    add_lesions(interface, lesions);
  
    for (int ii=0; ii<2; ii++) {
        rotation_z = fmod(rotation_z + 0.3, 360);
    }
}


static void set_config_lung(s_gplot *interface)
{
    if (simu.is_human == 1) {
        gplot_config(interface, "set xrange [-200:200]", "set yrange [0:350]", "set zrange [-275:25]");
    } else {
        gplot_config(interface, "set xrange [-85:85]", "set yrange [-105:65]", "set zrange [-110:60]");
    }
}


int make_video(void)
{
    s_gplot *interface = gplot_start(GPLOT_VIDEO_3D, "./results/video.mp4", .size = {1920, 1080}, .font_size = 18,
                                     .video_parallel_threads = omp_get_max_threads());
    gplot_config(interface, "unset tics", "unset border", "set xyplane relative 0");
    gplot_config(interface, "set parametric", "set urange [-pi/2:pi/2]", "set vrange [0:2*pi]",
                   "set samples 50", "set pm3d depthorder");
    set_config_lung(interface);

    t_lesions *lesions, *reinfected_lesions, *fusioned_lesions;
    request_memory_lesions(&lesions, &reinfected_lesions, &fusioned_lesions);
    initialize_lesions(lesions); 
    initialize_lesions(reinfected_lesions);
    initialize_lesions(fusioned_lesions);

    // Simulation
    int immune_state = 0;
    int n_frames = 0;
    double t = 0;
    double timm = normal_distribution(params.tc, params.tc * params.ftcdev);
    while (t <= params.tmax) {
        if ( (t >= timm) && (immune_state == 0) ) {
            activate_immune(lesions);
            immune_state = 1;
        }
        lesion_growth(lesions);
        int N_intra_it, N_extra_it;
        int N_reinf = endogenous_reinfection(lesions, reinfected_lesions, immune_state, &N_intra_it, &N_extra_it);
        if ( N_reinf == ERROR ) return ERROR;
        fusion(lesions, fusioned_lesions);
        synchronise_lesions(lesions, reinfected_lesions, fusioned_lesions);

        if ( (int)t < (int)(t+params.dt) || (int)(t+0.33) < (int)(t+params.dt+0.33) || (int)(t+0.66) < (int)(t+params.dt+0.66)) {
            // gplot keeps commands across frames, so the setup above need not be repeated.
            // gplot_start already opened frame 1; open a new frame only before each later one
            if (n_frames > 0) gplot_next_frame(interface);
            save_photogram(interface, lesions, (int)t);
            n_frames++;
        }
        t += params.dt;
    }
    puts("-> Rendering video frames...");
    gplot_end(interface);
    return OK;
}

