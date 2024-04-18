#include "tracer.h"
#include "diffusion.h"

#define Emin 5.0E-5                    // limiting value for SGS TKE
#define kappa 0.4                      // von karman constant
#define ce1 0.19                       
#define ce2 0.51
#define Cm 0.12
#define Cn 0.76
#define eps1 1.0E-10                     // limiting value for temp gradient

face vector Km[], Kh[], Ke[];
(const) face vector Pr;
scalar Evis[]; // Cell Centered diffusivity

scalar e12p[];

/**
## Here i first diagnose Km using mixing length scale
## later on TKE
 * 
 */
event Km_tke(i++)
{
  // scalar lambda[];
  // stable or unstable condition !! this requries buoyancy field
  // how to transfer face gradient to cell center gradient and calculate the mixing length
  // - here we get face vector of buoyancy gradient
  // vector db[];
  // gradients({b}, {db});
  // scalar dbdz[];
  foreach ()
  {
    double lambda;
    double dbdz = (b[0, 1] - b[0, -1]) / (Delta * 2);
    // limiting the temp gradient value
    if (dbdz < 0 && dbdz >= -eps1)
      dbdz = -eps1;
    if (dbdz >= 0 && dbdz <= eps1)
      dbdz = eps1;

    if (dbdz <= 0)
    {
      lambda = Delta;
    }
    else
    {
      e120[] = (e120[] <= Emin) ? Emin : e120[];
      lambda = pow(pow(kappa * y, -1) + pow(Cn * e120[] / sqrt(fabs(dbdz)), -1), -1);
    }
    Evis[] = Cm * lambda * e120[];
    // cell center S^2/2  -the diagonal part
    /**
    // TKE prognostic -- shear production term sbshr
    $$\left(\frac{\partial \tilde{u}_j}{\partial x_i}+\frac{\partial \tilde{u}_i}{\partial x_j}\right) \frac{\partial \tilde{u}_i}{\partial x_j}$$
    */
    double tdef2 = 2 / sq(Delta) * (sq(u.x[1, 0, 0] - u.x[]) + sq(u.y[0, 1, 0] - u.y[]) + sq(u.z[0, 0, 1] - u.z[]));
    // other parts
    tdef2 += 0.25 * (sq((u.y[0, 0, 1] - u.y[-1, 0, 1]) / Delta + (u.x[0, 0, 1] - u.x[0, 0, 0]) / Delta) +
                    sq((u.y[0, 0, 0] - u.y[-1, 0, 0]) / Delta + (u.x[0, 0, 0] - u.x[0, 0, -1]) / Delta) +
                    sq((u.y[1, 0, 0] - u.y[0, 0, 0]) / Delta + (u.x[1, 0, 0] - u.x[1, 0, -1]) / Delta) +
                    sq((u.y[1, 0, 1] - u.y[0, 0, 1]) / Delta + (u.x[1, 0, 1] - u.x[1, 0, 0]) / Delta));
    tdef2 += 0.25 * (sq((u.x[0, 1, 0] - u.x[0, 0, 0]) / Delta + (u.z[0, 1, 0] - u.z[-1, 1, 0]) / Delta) +
                    sq((u.x[0, 0, 0] - u.x[0, -1, 0]) / Delta + (u.z[0, 0, 0] - u.z[-1, 0, 0]) / Delta) +
                    sq((u.x[1, 0, 0] - u.x[1, -1, 0]) / Delta + (u.z[1, 0, 0] - u.z[0, 0, 0]) / Delta) +
                    sq((u.x[1, 1, 0] - u.x[1, 0, 0]) / Delta + (u.z[1, 1, 0] - u.z[0, 1, 0]) / Delta));
    tdef2 += 0.25 * (sq((u.z[0, 0, 1] - u.z[0, 0, 0]) / Delta + (u.y[0, 0, 1] - u.y[0, -1, 1]) / Delta) +
                    sq((u.z[0, 0, 0] - u.z[0, 0, -1]) / Delta + (u.y[0, 0, 0] - u.y[0, -1, 0]) / Delta) +
                    sq((u.z[0, 1, 0] - u.z[0, 1, -1]) / Delta + (u.y[0, 1, 0] - u.y[0, 0, 0]) / Delta) +
                    sq((u.z[0, 1, 1] - u.z[0, 1, 0]) / Delta + (u.y[0, 1, 1] - u.y[0, 0, 1]) / Delta));

    // here we cancel out the e120 in the equations
    e12p[] = (Cm * lambda * tdef2 / 2 - 3 * Cm * lambda * dbdz / 2)      // shear and buoyancy
            - (ce1 + ce2 * lambda / Delta) * sq(e120[]) / (2 * lambda); // dissipation
    #if CANOPY
        e12p[] = e12p[] - 4. / 3. * Cd * PAD(y) * U[] * e120[] * CL[];
    #endif
  }
  boundary({Evis});
  
  foreach_face()
  {
    Km.x[] = (Evis[] + Evis[-1]) / 2;
    Kh.x[] = Km.x[] / Pr.x[];
    Ke.x[] = Km.x[] * 2.;
  }
  boundary({Km, Kh, Ke});
}

// diffusiona and advection terms happen here
mgstats mgb;
event tracer_diffusion(i++)
{
  mgb = diffusion(e120, dt, Ke, r = e12p);
}
