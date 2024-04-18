// This code test the goal function of different setting of fan

#include <sys/stat.h>
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "output_slices.h"

int maxlevel, minlevel;
double eps;
double TEND = 15000.;
double T_str = 12600.;
double T_end = 15000.;

#include "physics.h"
#include "fan.h"

char filedir[200] = "/scratch/ydai2/basilisk/D2M11/Bcooling/BV1/";
char W12dir[230];

int main()
{
  minlevel = 5;
  maxlevel = 11;
  a = av;
  L0 = 2048;
  X0 = Y0 = Z0 = 0;
  N = 256;
  foreach_dimension()
  {
    u.x.refine = refine_linear; // Momentum conserved     }
  }

  fan.prolongation = fraction_refine; // Fan is a volume fraction
  p.refine = p.prolongation = refine_linear;
  b.gradient = minmod2; // Flux limiter
  Boundary_C();
  leaf_BC();
  Pr = unityf;
  mu = Km;
  rot.fan = true;    // Yes we want a fan
  rot.rotate = true; // If we want it to rotate
  rot.start = T_str; // start time of rotor forcing
  rot.stop = T_end;  // stop time of rotor forcing
  rot.phit = 2 * M_PI / 288.;
  rot.ST_phit = 0.;  // Full rotation
  eps = .3;
  run();
}

void sim_dir_create(){
    sprintf(W12dir, "%sW12/", filedir);
    if (pid() == 0)
    {
        struct stat st = {0};
        if (stat(W12dir, &st) == -1)
        {
            mkdir(W12dir, 0777);
        }
    }
}

event init(t = 0)
{
  if (restore(file = "dump-12000")){
    N = 256;
    init_rotor();
    sim_dir_create();
  }
  else {
    init_physics();
    sim_dir_create();
    if (rot.fan)
    {
        init_rotor();
    }
    scalar n[];
    foreach ()
    {
        n[] = 0.;
        if (y <= Zh * 1.5)
        {
            n[] = noise();
            if (t >= rot.start && t < rot.stop && (sq(x - rot.x0) + sq(y - rot.y0) + sq(z - rot.z0) <= sq(rot.R * 1.5)))
            {
                n[] = noise();
            }
        }
    }
    while (adapt_wavelet((scalar *){u, b, n}, (double[]){eps, eps, eps, 0.35 * 9.81 / 273, 0.001}, maxlevel, minlevel).nf)
    {
        foreach ()
        {
            b[] = STRAT(y);
            u.x[] = WIND(y);
        }
    }
    foreach(){
        TV[] = 274.15;
        e120[] = Emin;
    }
  }
}

event adapt(i++)
{
    scalar n[];
    foreach(){
        n[] = 0.;
        if (y <= Zh * 1.5){
            n[] = noise();
            if (t >= rot.start && t < rot.stop && (sq(x - rot.x0) + sq(y - rot.y0) + sq(z - rot.z0) <= sq(rot.R * 1.5)))
            {
                n[] = noise();
            }
        }
    }
    adapt_wavelet((scalar *){u, b, n}, (double[]){eps, eps, eps, .35 * 9.81 / 273, 0.001}, maxlevel, minlevel);
}

double save_int = 600.;
event profiles(t += 1.)
{
    // ---- W12 output group (11 * 11): i = y (0.5..10.5), j = z (990.5..1100.5), x = 528 ----
    int W12_res = 11;
    // air temperature
    char nameW12_B[300];
    coord sliceW12 = {0., 1., 1.};
    sliceW12.x = (500. + 28.) / L0;
    snprintf(nameW12_B, 300, "%sW12_Bt=%02dx=528", W12dir, (int)(t / save_int));
    FILE *fpW12_B = fopen(nameW12_B, "a");
    output_W12(list = (scalar *){b}, fp = fpW12_B, n = W12_res, linear = true, plane = sliceW12);
    fclose(fpW12_B);

    // leaf temperature
    char nameW12_TV[300];
    snprintf(nameW12_TV, 300, "%sW12_TVt=%02dx=528", W12dir, (int)(t / save_int));
    FILE *fpW12_TV = fopen(nameW12_TV, "a");
    output_W12(list = (scalar *){TV}, fp = fpW12_TV, n = W12_res, linear = false, plane = sliceW12);
    fclose(fpW12_TV);

    char nameW12_U[300];
    snprintf(nameW12_U, 300, "%sW12_Ut=%02dx=528", W12dir, (int)(t / save_int));
    FILE *fpW12_U = fopen(nameW12_U, "a");
    output_W12(list = (scalar *){u.x}, fp = fpW12_U, n = W12_res, linear = true, plane = sliceW12);
    fclose(fpW12_U);

    char nameW12_V[300];
    snprintf(nameW12_V, 300, "%sW12_Vt=%02dx=528", W12dir, (int)(t / save_int));
    FILE *fpW12_V = fopen(nameW12_V, "a");
    output_W12(list = (scalar *){u.z}, fp = fpW12_V, n = W12_res, linear = true, plane = sliceW12);
    fclose(fpW12_V);

    char nameW12_W[300];
    snprintf(nameW12_W, 300, "%sW12_Wt=%02dx=528", W12dir, (int)(t / save_int));
    FILE *fpW12_W = fopen(nameW12_W, "a");
    output_W12(list = (scalar *){u.y}, fp = fpW12_W, n = W12_res, linear = true, plane = sliceW12);
    fclose(fpW12_W);




    coord sliceW12_e = {0., 1., 1.};
    sliceW12_e.x = (500. + 32.) / L0;
    char nameW12_B_e[300];
    snprintf(nameW12_B_e, 300, "%sW12_Bt=%02dx=532", W12dir, (int)(t / save_int));
    FILE *fpW12_B_e = fopen(nameW12_B_e, "a");
    output_W12(list = (scalar *){b}, fp = fpW12_B_e, n = W12_res, linear = true, plane = sliceW12_e);
    fclose(fpW12_B_e);

    // leaf temperature
    char nameW12_TV_e[300];
    snprintf(nameW12_TV_e, 300, "%sW12_TVt=%02dx=532", W12dir, (int)(t / save_int));
    FILE *fpW12_TV_e = fopen(nameW12_TV_e, "a");
    output_W12(list = (scalar *){TV}, fp = fpW12_TV_e, n = W12_res, linear = false, plane = sliceW12_e);
    fclose(fpW12_TV_e);

    char nameW12_U_e[300];
    snprintf(nameW12_U_e, 300, "%sW12_Ut=%02dx=532", W12dir, (int)(t / save_int));
    FILE *fpW12_U_e = fopen(nameW12_U_e, "a");
    output_W12(list = (scalar *){u.x}, fp = fpW12_U_e, n = W12_res, linear = true, plane = sliceW12_e);
    fclose(fpW12_U_e);

    char nameW12_V_e[300];
    snprintf(nameW12_V_e, 300, "%sW12_Vt=%02dx=532", W12dir, (int)(t / save_int));
    FILE *fpW12_V_e = fopen(nameW12_V_e, "a");
    output_W12(list = (scalar *){u.z}, fp = fpW12_V_e, n = W12_res, linear = true, plane = sliceW12_e);
    fclose(fpW12_V_e);

    char nameW12_W_e[300];
    snprintf(nameW12_W_e, 300, "%sW12_Wt=%02dx=532", W12dir, (int)(t / save_int));
    FILE *fpW12_W_e = fopen(nameW12_W_e, "a");
    output_W12(list = (scalar *){u.y}, fp = fpW12_W_e, n = W12_res, linear = true, plane = sliceW12_e);
    fclose(fpW12_W_e);
}

event dump_file1(t = 12000.; t += 288.; t <= TEND)
{
    char name[80];
    sprintf(name, "dump-%05d", (int)t);
    dump(file = name);
}

event progress(t = 12000.; t += 4.)
{
    char name[80];
    sprintf(name, "dumpfields-%05d", (int)t);
    dump(file = name, list = {b, TV, u});
}

event end(t = TEND)
{
}
