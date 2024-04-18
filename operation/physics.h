
#define CP 1005.	// C_p for air 
#define gCONST 9.81	// Gravitational constant
#define TREF 273.15     // Kelvin
#define karman 0.4      // von Karman constant 

#define roughY0u 0.1    // roughness wind length 
#define roughY0h 0.1     // roughness heat length

#define WIND(s) (s < 300? (max(0.20 * log(s / roughY0u), 0.)) : (s<1024 ? (1.6 + 0.01 * (s - 300)): 9.))   
#define QFLX -0.00024                                        // 0 (0.001 = 30wm2)
// #define BSURF (1.5 * b[] - 0.5 * b[0, 1])
#define BSURF ((karman * sqrt(sq(u.x[]) + sq(u.z[])))>0.01? b[]:(0.9 * b[]))
#define BCb (gCONST / TREF * 0.9)
#define GFLX (-Lambda * (BSURF - bd))
double Lambda = 0.0035, bd = 0.1; // Grass coupling
#define STRAT(s) (gCONST / TREF * (s <= 100 ? 9.5 : (9.5 + 0.003 * (s - 100))))

scalar b[];
scalar U[];
scalar e120[];
scalar * tracers = {b, e120};

#define CANOPY 1
#if CANOPY
#include "Canopy.h"
#endif

#if dimension == 3
#include "SGS_TKE.h"
#endif

face vector av[]; 

void Boundary_C(){
    b.nodump = false;

    u.r[bottom] = dirichlet(0.);
    u.n[bottom] = dirichlet(0.);
    u.t[bottom] = dirichlet(0.);

    // u.r[top] = dirichlet(WIND(y));
    u.r[top] = neumann(0.);
    u.n[top] = dirichlet(0.);
    u.t[top] = neumann(0.);

    periodic(left);

    // b[bottom] = BSURF;
    b[bottom] = dirichlet(BSURF);
    b[top] = dirichlet(STRAT(y));

    Evis[bottom] = dirichlet(0.); // Flux is explicitly calculated
    Evis[top] = dirichlet(0.);

    periodic(front);
}

void init_physics(){
    foreach ()
    {
        b[] = STRAT(y);
        u.x[] = WIND(y);
        u.y[] = 0.;
        u.z[] = 0.;
    }
}

event dampping(i += 1)
{
        double bf_L = L0/2.;
        double df = 0.15;
        double relaxtime = dt;
        foreach ()
        {
            double a = 0.;
            if (y > (L0 - bf_L))
            {
                a = df * sq(sin(M_PI / 2. * (y - L0 + bf_L) / bf_L));
            }
            b[] = b[] + a * (STRAT(y) - b[]) * relaxtime;
            u.x[] = u.x[] + a * (WIND(y) - u.x[]) * relaxtime;
            u.y[] = u.y[] - a * u.y[] * relaxtime;
            u.z[] = u.z[] - a * u.z[] * relaxtime;
        }
}


/* Gravity forcing */
event acceleration(i++){
	foreach_face(y){
		av.y[] = (b[] + b[0,-1])/2.;
	}
    foreach_face(x)
    {
        av.x[] = 4.4E-7;
    }
	#if CANOPY
        foreach ()
        {
            U[] = sqrt(sq(u.x[]) + sq(u.y[]) + sq(u.z[]));
        }
        foreach_face()
        {
            av.x[] = av.x[] - Cd * PAD(y) * (U[] + U[-1])/2. * (u.x[] + u.x[-1])/2. * (CL[] + CL[-1]) / 2.;
        }
	#endif
}

mgstats mgb;
/* Diffusion */
event tracer_diffusion(i++){
    scalar r[];
    foreach() {
        r[] = 0;
        if (y < Delta)
            r[] = (QFLX + GFLX)/Delta; // div needed as normalization
        #if CANOPY
        // if (y>Delta)
        // always don't forget about the unit
        // first (K m s-1) then (m2 s-3) then (m s-3)
        r[] = r[] + H[] / (Cp_a * rho_a) * (gCONST / TREF) * PAD(y) * CL[];
        #endif
	}
    double flx = 0, bt = 0;
    double fctr = rho_a * CP * TREF / gCONST;
    foreach_boundary(bottom, reduction(+:flx) reduction(+:bt)) {
        flx = flx + (QFLX + GFLX) * sq(Delta);
         bt = bt + BSURF * sq(Delta);
    }
    bt = bt/sq(L0);
    flx = flx/sq(L0);
    fprintf(stderr, "soil=%g %g %g %d\n", t, fctr*flx, bt * TREF / gCONST, i);
	mgb = diffusion(b, dt, Kh, r = r);
}

