/**
Here we analyze the dump files
 */
#include <sys/stat.h>
#include "grid/octree.h"
#include "utils.h"
#include "profile6.h"

#define gCONST 9.81 // Gravitational constant
#define TREF 273.15 // Kelvin
#define RADIUS(s) (sqrt(sq(x - 500.) + sq(z - 1024.)) * ((y >= s) && (y <= (s + 1))))

struct sOutput
{
    double t0_slice;
    double dt_slices;
    char main_dir[200];
    char b_slicedir[250];
    char u_slicedir[250];
    char TV_slicedir[250];
    char vali_dir[250];
    char Rprof_dir[250];
};

struct sOutput out = {.dt_slices = 4.,
                      .t0_slice = 14400.,
                      .main_dir = "/scratch/ydai2/basilisk/D2M10/Bcooling/BV4"};

void sim_dir_create()
{
    sprintf(out.b_slicedir, "%s/b_slice/", out.main_dir);
    sprintf(out.u_slicedir, "%s/u_slice/", out.main_dir);
    sprintf(out.TV_slicedir, "%s/TV_slice/", out.main_dir);
    sprintf(out.vali_dir, "%s/VALI_sim/", out.main_dir);
    sprintf(out.vali_dir, "%s/Radial_prof/", out.main_dir);
    if (pid() == 0)
    {
        struct stat st = {0};
        if (stat(out.b_slicedir, &st) == -1)
        {
            mkdir(out.b_slicedir, 0777);
        }
        if (stat(out.u_slicedir, &st) == -1)
        {
            mkdir(out.u_slicedir, 0777);
        }
        if (stat(out.TV_slicedir, &st) == -1)
        {
            mkdir(out.TV_slicedir, 0777);
        }
        if (stat(out.vali_dir, &st) == -1)
        {
            mkdir(out.vali_dir, 0777);
        }
        if (stat(out.vali_dir, &st) == -1)
        {
            mkdir(out.vali_dir, 0777);
        }
    }
}

scalar b[];
scalar TV[];
vector u[];
// this will give 11 integer output with last one
double save_int = 100.;
int sec = 0;

int main(int argc, char *argv[])
{
    if (argc > 1)
        sec = atoi(argv[1]);
    printf("restore the dump file\n");
    sim_dir_create();
    // double t0 = out.t0_slice;
    for (double tt = out.t0_slice + sec * save_int; tt < out.t0_slice + (sec + 1) * save_int; tt+=out.dt_slices){
    // double tt = t0;
    char name[250];
    if (tt == out.t0_slice)
        sprintf(name, "%s/dump-%05d", out.main_dir, (int)tt);
    else    
        sprintf(name, "%s/dumpfields-%05d", out.main_dir, (int)tt);
    if (restore(file = name))
    {
        printf("this file at %g s is working\n", t);
    } else {
        return 1;
    }
    boundary(all);
    char fnameb[99];
    for (double zz = 0.; zz < 25.; zz++)
    {
        sprintf(fnameb, "profb%02g", zz);
        profile_equi({b}, RADIUS(zz), L0, fnameb);
    }
    }
}
