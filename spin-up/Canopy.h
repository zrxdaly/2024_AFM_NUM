/**
## Canopy model
We implement canopy model based on Patton et al 2016; Boekee et al 2023. This parameterization should incorporate with large eddy simulation. This header file only solve energy balance for the leaf. 
To solve the flow part, the convective heat exchange (source term for temp) and wind drag should be including in other files (!!!)
*/
// #include "profile5b.h"
// #include "runge-kutta.h"                 // for now we use forward euler, later on maybe we will just use time integgrator from basilisk
// define where is the canopy
#define Zh 3. // canopy top height [m]   mind the float
scalar CL[]; // define the canopy fraction
#define CCan(x0, y0, z0) (((abs(x0-z0)) % 4 == 0) && (y0 < Zh) ? 1 : 0)
// #define CANOPY (y <= Zh)? 1:0              // create a mask for canopy layer
#define PAD(s) 2.0
// #define PAD(s) ((s <= Zh) ? 0.4 * sin(s / Zh * M_PI) : 0)
// #define PAD(s) (s <= Zh ? (s <= 0.7 ? (2. * (exp(s) * 0.49321 - 0.49321)) : (1.5 * (exp(-0.62366 * s) - exp(-3 * 0.62366)))) : 0)
// wind drag
#define Cd 0.20 // dimensionless drag coefficient (Shaw and Schumann 1992)

// default parameters
#define boltz 5.67E-8
#define Gconst 9.81
#define T_ref 273.15
#define Kd 0.024 // thermal conductivity of air [W m-1 K-1]
// other parameters
#define VF_s 0.1
#define VF_g VF_s
#define VF_l (1. - VF_s) // view factor of sky, ground, leaf
#define eps_s 0.8
#define eps_g 0.98
#define eps_l 0.96 // emissivity of sky, ground, leaf
// #define T_s 269.584
// #define T_g 275.961         // temperature of sky and surface in [kelvin]
#define T_s 268.584
#define T_g 272.961         // temperature of sky and surface in [kelvin]
#define Cp_l 2.0E6         // leaf heat capacity [J m^-3 K^-1]
#define Cp_a 1005.         // dry air heat capacity [J kg^-1 K^-1]
#define rho_a 1.27         // air density at near 10 C [kg m^-3]
#define dvis 1.718E-5      // dynamic viscosity of the air [N s m-2]
#define vis (dvis / rho_a) // kinematic viscosity of the air [m^2 s^-1]
// geometry of leaf
#define R_l 4E-2                  // radius of leaf surface [m]
#define L_l (2 * R_l)             // representative length scale [m]
#define d_l 2.0E-4                // thickness of the leaf  [m]
#define A_l (2. * M_PI * sq(R_l)) // surface area leaf      [m^2]
#define V_l (A_l / 2. * d_l)      // volume leaf            [m^3]
// leaf temperature
scalar TV[];                      // leaf temp [K]
scalar H[];                     // convective heat exchange [W m^-2]

void leaf_BC()
{
    TV.refine = refine_injection;
    TV.coarsen = refine_injection;
    H.refine = refine_injection;
    H.coarsen = refine_injection;
    CL.refine = refine_injection;
    CL.coarsen = refine_injection;
}

event leaf_flow(i++)
{
    // refine the canopy layer
    foreach ()
    {
        CL[] = CCan(x, y, z);
    }
    // radiation - Stefan–Boltzmann law
    // double Lwin, Lwout;
    boundary({CL});
    scalar Lwnet[];
    foreach(){
        Lwnet[] = 0.;
        if (CL[] > 0.)
        {
            double Lwin = 0.5 * VF_s * eps_s * boltz * pow(T_s, 4) +
                          0.5 * VF_g * eps_g * boltz * pow(T_g, 4) +
                          1. * VF_l * eps_l * boltz * pow(TV[], 4);
            double Lwout = eps_l * boltz * pow(TV[], 4);
            Lwnet[] = Lwin - Lwout;
        }
    }
    // Convective -- resistance of turbulent heat exchange -- Boekee et al 2023
    foreach ()
    {
        H[] = 0.;
        if (CL[] > 0.)
        {
            double T_a = b[] * T_ref / Gconst + T_ref;
            double gstar = Gconst * (TV[] - T_a) / T_a;

            // watch out for the sign here for gstar
            double M = sqrt(sq(u.x[]) + sq(u.y[]) + sq(u.z[]) + fabs(2 * L_l * gstar));
            double Re = fabs(M) * L_l / vis;
            double Nu = (Re > 2E4) ? 0.032 * pow(Re, 0.8) : 0.6 * pow(Re, 0.5);
            double rH = L_l / Nu / Kd * Cp_a * rho_a;

            H[] = Cp_a * rho_a / rH * (TV[] - T_a);
            // diagnose leaf temp with forward euler
            // for vegetation is negative, in the source term of b, it should be positive
            TV[] += dt * (Lwnet[] - H[]) * A_l / (Cp_l * V_l);
        }
    }
}