/**
Here we analyze the dump files
 */
#include <sys/stat.h>
#include "grid/octree.h"
#include "utils.h"
#include "output_slices.h"

#define gCONST 9.81 // Gravitational constant
#define TREF 273.15 // Kelvin

struct sOutput
{
    double t0_slice;
    double dt_slices;
    char main_dir[200];
    char b_slicedir[250];
    char u_slicedir[250];
    char TV_slicedir[250];
    char vali_dir[250];
};

struct sOutput out = {.dt_slices = 4.,
                      .t0_slice = 7800.,
                      .main_dir = "/projects/0/tdse0633/large_domain/DB_D2M11/Bcooling/BV0"};

void sim_dir_create()
{
    sprintf(out.b_slicedir, "%s/b_slice/", out.main_dir);
    sprintf(out.u_slicedir, "%s/u_slice/", out.main_dir);
    sprintf(out.TV_slicedir, "%s/TV_slice/", out.main_dir);
    sprintf(out.vali_dir, "%s/VALI_sim/", out.main_dir);
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
        sprintf(name, "%s/dump-%04d", out.main_dir, (int)tt);
    else    
        sprintf(name, "%s/dumpfields-%05d", out.main_dir, (int)tt);
    if (restore(file = name))
    {
        printf("this file at %g s is working\n", t);
    } else {
        return 1;
    }
    boundary(all);
    // ---- W12 output group (11 * 11): i = y (0.5..10.5), j = z (990.5..1100.5), x = 528.5 ----
    // int W12_res = 11;
    // // air temperature
    // char nameW12_B[300];
    // coord sliceW12 = {0., 1., 1.};
    // sliceW12.x = (500. + 32.5) / L0;
    // snprintf(nameW12_B, 300, "%sW12_Bt=%02dx=528", out.vali_dir, (int)((t - out.t0_slice) / save_int));
    // FILE *fpW12_B = fopen(nameW12_B, "a");
    // output_W12(list = (scalar *){b}, fp = fpW12_B, n = W12_res, linear = true, plane = sliceW12);
    // fclose(fpW12_B);

    // // leaf temperature
    // char nameW12_TV[300];
    // snprintf(nameW12_TV, 300, "%sW12_TVt=%02dx=528", out.vali_dir, (int)((t - out.t0_slice) / save_int));
    // FILE *fpW12_TV = fopen(nameW12_TV, "a");
    // output_W12(list = (scalar *){TV}, fp = fpW12_TV, n = W12_res, linear = false, plane = sliceW12);
    // fclose(fpW12_TV);

    // char nameW12_U[300];
    // snprintf(nameW12_U, 300, "%sW12_Ut=%02dx=528", out.vali_dir, (int)((t - out.t0_slice) / save_int));
    // FILE *fpW12_U = fopen(nameW12_U, "a");
    // output_W12(list = (scalar *){u.x}, fp = fpW12_U, n = W12_res, linear = true, plane = sliceW12);
    // fclose(fpW12_U);

    // char nameW12_V[300];
    // snprintf(nameW12_V, 300, "%sW12_Vt=%02dx=528", out.vali_dir, (int)((t - out.t0_slice) / save_int));
    // FILE *fpW12_V = fopen(nameW12_V, "a");
    // output_W12(list = (scalar *){u.z}, fp = fpW12_V, n = W12_res, linear = true, plane = sliceW12);
    // fclose(fpW12_V);

    // char nameW12_W[300];
    // snprintf(nameW12_W, 300, "%sW12_Wt=%02dx=528", out.vali_dir, (int)((t - out.t0_slice) / save_int));
    // FILE *fpW12_W = fopen(nameW12_W, "a");
    // output_W12(list = (scalar *){u.y}, fp = fpW12_W, n = W12_res, linear = true, plane = sliceW12);
    // fclose(fpW12_W);

    // ---  output b2 contour  ---
    int resb12 = 310;
    // coord sliceb2 = {1., 2., 1.};
    // char name_b2[300];
    // snprintf(name_b2, 300, "%sb2_t=%02d", out.vali_dir, (int)((t - out.t0_slice) / save_int));
    // FILE *fp_b2 = fopen(name_b2, "a");
    // b12output_slice(list = (scalar *){b}, fp = fp_b2, n = resb12, linear = true, plane = sliceb2);
    // fclose(fp_b2);

    coord sliceb3 = {1., 4., 1.};
    char name_b3[300];
    snprintf(name_b3, 300, "%sb4_t=%02d", out.vali_dir, (int)((t - out.t0_slice) / save_int));
    FILE *fp_b3 = fopen(name_b3, "a");
    b12output_slice(list = (scalar *){b}, fp = fp_b3, n = resb12, linear = true, plane = sliceb3);
    fclose(fp_b3);

    // // ---  horizontal slice  ---
    // int resL0 = L0;
    // coord slice = {1., 0., 1.};
    // for (double hh = 0.5; hh < 21; hh++){
    //     slice.y = hh / L0;
    //     int HH = hh * 10;

    //     char nameSlice_B[300];
    //     snprintf(nameSlice_B, 300, "%sB%03dt=%02d", out.b_slicedir, HH, (int)((t - out.t0_slice) / save_int));
    //     FILE *fp_SB = fopen(nameSlice_B, "a");
    //     output_slice(list = (scalar *){b}, fp = fp_SB, n = resL0, linear = true, plane = slice);
    //     fclose(fp_SB);

    //     char nameSlice_U[300];
    //     snprintf(nameSlice_U, 300, "%sU%03dt=%02d", out.u_slicedir, HH, (int)((t - out.t0_slice) / save_int));
    //     FILE *fp_SU = fopen(nameSlice_U, "a");
    //     output_slice(list = (scalar *){u.x}, fp = fp_SU, n = resL0, linear = true, plane = slice);
    //     fclose(fp_SU);

    //     char nameSlice_V[300];
    //     snprintf(nameSlice_V, 300, "%sV%03dt=%02d", out.u_slicedir, HH, (int)((t - out.t0_slice) / save_int));
    //     FILE *fp_SV = fopen(nameSlice_V, "a");
    //     output_slice(list = (scalar *){u.z}, fp = fp_SV, n = resL0, linear = true, plane = slice);
    //     fclose(fp_SV);

    //     char nameSlice_W[300];
    //     snprintf(nameSlice_W, 300, "%sW%03dt=%02d", out.u_slicedir, HH, (int)((t - out.t0_slice) / save_int));
    //     FILE *fp_SW = fopen(nameSlice_W, "a");
    //     output_slice(list = (scalar *){u.y}, fp = fp_SW, n = resL0, linear = true, plane = slice);
    //     fclose(fp_SW);
    // }

    // for (double hh = 0.5; hh < 3; hh++){
    //     slice.y = hh / L0;
    //     int HH = hh * 10;
    //     char nameSlice_TV[300];
    //     snprintf(nameSlice_TV, 300, "%sTV%03dt=%02d", out.TV_slicedir, HH, (int)((t - out.t0_slice) / save_int));
    //     FILE *fp_STV = fopen(nameSlice_TV, "a");
    //     output_slice(list = (scalar *){TV}, fp = fp_STV, n = resL0, linear = false, plane = slice);
    //     fclose(fp_STV);
    // }

    // // ---  centeral TX contour  xy plane with y_max = 100; z = L0 / 2; --- (x, y)
    // coord sliceV = {1., 1., 0.5};
    // char nameSlice_VB[300];
    // snprintf(nameSlice_VB, 300, "%sV_Bt=%02d", out.b_slicedir, (int)((t - out.t0_slice) / save_int));
    // FILE *fp_VB = fopen(nameSlice_VB, "a");
    // output_Vslice(list = (scalar *){b}, fp = fp_VB, n = resL0, linear = true, plane = sliceV);
    // fclose(fp_VB);

    // char nameSlice_VTV[300];
    // snprintf(nameSlice_VTV, 300, "%sV_TVt=%02d", out.TV_slicedir, (int)((t - out.t0_slice) / save_int));
    // FILE *fp_VTV = fopen(nameSlice_VTV, "a");
    // output_Vslice_TV(list = (scalar *){TV}, fp = fp_VTV, n = resL0, linear = false, plane = sliceV);
    // fclose(fp_VTV);

    // char nameSlice_VU[300];
    // snprintf(nameSlice_VU, 300, "%sV_Ut=%02d", out.u_slicedir, (int)((t - out.t0_slice) / save_int));
    // FILE *fp_VU = fopen(nameSlice_VU, "a");
    // output_Vslice(list = (scalar *){u.x}, fp = fp_VU, n = resL0, linear = true, plane = sliceV);
    // fclose(fp_VU);

    // char nameSlice_VV[300];
    // snprintf(nameSlice_VV, 300, "%sV_Vt=%02d", out.u_slicedir, (int)((t - out.t0_slice) / save_int));
    // FILE *fp_VV = fopen(nameSlice_VV, "a");
    // output_Vslice(list = (scalar *){u.z}, fp = fp_VV, n = resL0, linear = true, plane = sliceV);
    // fclose(fp_VV);

    // char nameSlice_VW[300];
    // snprintf(nameSlice_VW, 300, "%sV_Wt=%02d", out.u_slicedir, (int)((t - out.t0_slice) / save_int));
    // FILE *fp_VW = fopen(nameSlice_VW, "a");
    // output_Vslice(list = (scalar *){u.y}, fp = fp_VW, n = resL0, linear = true, plane = sliceV);
    // fclose(fp_VW);

    // // ---  centeral TZ contour  xy plane with y_max = 100; x = 499.; --- (y, z)
    // coord sliceZ = {0.5, 1., 1.};
    // sliceZ.x = 500. / L0;
    // char nameSlice_ZB[300];
    // snprintf(nameSlice_ZB, 300, "%sZ_Bt=%02d", out.b_slicedir, (int)((t - out.t0_slice) / save_int));
    // FILE *fp_ZB = fopen(nameSlice_ZB, "a");
    // output_Zslice(list = (scalar *){b}, fp = fp_ZB, n = resL0, linear = true, plane = sliceZ);
    // fclose(fp_ZB);

    // char nameSlice_ZTV[300];
    // snprintf(nameSlice_ZTV, 300, "%sZ_TVt=%02d", out.TV_slicedir, (int)((t - out.t0_slice) / save_int));
    // FILE *fp_ZTV = fopen(nameSlice_ZTV, "a");
    // output_Zslice_TV(list = (scalar *){TV}, fp = fp_ZTV, n = resL0, linear = false, plane = sliceZ);
    // fclose(fp_ZTV);

    // char nameSlice_ZU[300];
    // snprintf(nameSlice_ZU, 300, "%sZ_Ut=%02d", out.u_slicedir, (int)((t - out.t0_slice) / save_int));
    // FILE *fp_ZU = fopen(nameSlice_ZU, "a");
    // output_Zslice(list = (scalar *){u.x}, fp = fp_ZU, n = resL0, linear = true, plane = sliceZ);
    // fclose(fp_ZU);

    // char nameSlice_ZV[300];
    // snprintf(nameSlice_ZV, 300, "%sZ_Vt=%02d", out.u_slicedir, (int)((t - out.t0_slice) / save_int));
    // FILE *fp_ZV = fopen(nameSlice_ZV, "a");
    // output_Zslice(list = (scalar *){u.z}, fp = fp_ZV, n = resL0, linear = true, plane = sliceZ);
    // fclose(fp_ZV);

    // char nameSlice_ZW[300];
    // snprintf(nameSlice_ZW, 300, "%sZ_Wt=%02d", out.u_slicedir, (int)((t - out.t0_slice) / save_int));
    // FILE *fp_ZW = fopen(nameSlice_ZW, "a");
    // output_Zslice(list = (scalar *){u.y}, fp = fp_ZW, n = resL0, linear = true, plane = sliceZ);
    // fclose(fp_ZW);
    }
}
