/*
    Phenomenological viscoelastic shell ("Model C" of Chahare, Ouzeri, Wilson,
    Bal et al., "Multiscale wrinkling and folding dynamics in epithelial
    shells", Supplementary Note 6) for the inflation -> hold -> deflation of an
    epithelial dome attached on a circular footprint.

    Unknowns
      posDHand  : X, Y, Z            shell mid-surface (subdivision surfaces)
      gloDHand  : P (+ 6 unused)     Lagrange multiplier of the enclosed-volume
                                     constraint (P/deltat = luminal pressure)
    Internal (Gauss-point) state, stored in paramStr->b_aux
      G^{IJ}     material (natural) metric, relaxes with viscosity fric2   Eq. (58)
      kbar_{IJ}  material curvature,        relaxes with viscosity fric3   Eq. (59)

    Energy per unit reference area (all moduli in tension units):
      w_mem(C,G)   neo-Hookean relative to G             (mu, lambda)       Eq. (55)
      w_pre(C)     neo-Hookean relative to the flat disc (mu_pre, lambda_pre)
                   = weak strain dependence of the pre-tension             Eq. (54)
      w_bend       h^2/4 C (k - kbar)(k - kbar), h = height_in / J_n       Eq. (56)
      sigma_a J    constant active pre-tension, sigma_a = gamma            Eq. (53)
    plus a viscous drag (fric) on the mid-surface velocity, a one-sided
    confinement below z = 0 (Kconf), and the volume constraint.

    Time discretisation: Rayleigh-type incremental minimisation; every term
    that has units of energy is multiplied by deltat, so the residual is
    (energy rate) * deltat^2 and the drag term carries no deltat.
*/

#ifndef AUX_VISCO_SHELL_H
#define AUX_VISCO_SHELL_H

#include <iostream>
#include <mpi.h>

#include "hl_DistributedMesh.h"
#include "hl_FillStructure.h"
#include "hl_ParamStructure.h"
#include "hl_DOFsHandler.h"
#include "hl_HiPerProblem.h"
#include "hl_Math.h"
#include "hl_Array.h"

using namespace std;
using namespace hiperlife;

using Teuchos::rcp;
using Teuchos::RCP;

struct MembParams
{
    enum RealParameters
    {
        deltat,
        // mechanics
        mu, lambda,             // passive (relaxing) membrane Lame moduli
        mu_pre, lambda_pre,     // non-relaxing pre-tension Lame moduli
        fact_elastic,           // multiplier on mu_pre, lambda_pre
        gamma,                  // active pre-tension sigma_a
        height_in,              // shell thickness (bending), thins as h/J
        fric,                   // drag coefficient
        fric2,                  // in-plane viscosity  eta_c  (metric relaxation)
        fric3,                  // bending viscosity   eta'_c (curvature relaxation)
        Kconf,                  // confinement z >= 0
        factor,                 // target volume / reference area (set every step)
        // optional increase of the drag near the substrate late in deflation
        fric_fact, width, thick,
        // relaxation integrator: 0 backward Euler (local Newton), 1 forward Euler
        forward_new
    };
    enum IntParameters
    {
        timestep,               // current step (set every step)
        time_target,            // G, kbar are reset to the reference while timestep < time_target
        vnstep, gap,            // inflation steps, hold steps
        control_fric            // drag increase starts at step vnstep+gap+control_fric
    };
    HL_PARAMETER_LIST DefaultValues
    {
        {"deltat", 0.0001},
        {"mu", 1.0}, {"lambda", 1.0},
        {"mu_pre", 1.5797}, {"lambda_pre", 1.2147},
        {"fact_elastic", 1.0},
        {"gamma", 1.0},
        {"height_in", 1.0},
        {"fric", 0.1}, {"fric2", 1.0}, {"fric3", 1.0},
        {"Kconf", 1.0},
        {"factor", 1.0},
        {"fric_fact", 1.0}, {"width", 0.0}, {"thick", 0.0},
        {"forward_new", 0.0}
    };
};

// ---------------------------------------------------------------------------
// Per-Gauss-point storage (indices into paramStr->a_aux / b_aux)
// a_aux: written by LS (mechanics) at the converged state, read by
//        LS_ReacDif (relaxation) and LS_ED (output projection)
// ---------------------------------------------------------------------------
enum AuxA
{
    A_METRIC   = 0,   // g_IJ             (4, row major)
    A_DWPDGP   = 4,   // dW_mem/dG        (4)  forward-Euler relaxation only
    A_IMETRICR = 8,   // reference g^IJ   (4)
    A_ICURVR   = 12,  // inverse reference curvature (4)
    A_CURVA    = 16,  // current curvature k_IJ (4)
    A_JACR     = 20,  // reference area element
    A_CONSTB   = 21,  // h^2/4 (bending prefactor)
    A_NORMAL   = 22,  // unit normal (3)          output
    A_THICK    = 25,  // current thickness h/J_n  output
    A_EA       = 26,  // membrane energy density  output
    A_EB       = 27,  // bending energy density   output
    A_PA       = 28,  // tension power            output
    A_TRSIG    = 29,  // trace of membrane stress output
    A_DEVSIG   = 30,  // norm of deviatoric stress output
    A_JAC      = 31,  // current area element
    N_AUXA     = 32
};
// b_aux: internal variables, written by LS_ReacDif, read by LS and LS_ED
enum AuxB
{
    B_GP    = 0,      // G^11, G^12, G^22     material metric (contravariant)
    B_GN    = 3,      // kbar_11, kbar_12, kbar_22  material curvature
    B_DA    = 6,      // in-plane viscous dissipation density   output
    B_DB    = 7,      // bending viscous dissipation density    output
    N_AUXB  = 8
};

void LS(hiperlife::FillStructure& fillStr);          // mechanics + volume constraint
void LS_ReacDif(hiperlife::FillStructure& fillStr);  // update of G, kbar (local, per Gauss point)
void LS_ED(hiperlife::FillStructure& fillStr);       // L2 projection of Gauss-point fields (output)

// one-sided cubic confinement keeping the shell above z = 0
inline double ConfinementPotential(hiperlife::Tensor::tensor<double,1>& x, double Kconf)
{
    return (x(2) < 0) ? -Kconf/3.0 * x(2)*x(2)*x(2) : 0.0;
}
inline double dConfinementPotential(hiperlife::Tensor::tensor<double,1>& x, double Kconf)
{
    return (x(2) < 0) ? -1*Kconf*x(2)*x(2) : 0.0;
}
inline double ddConfinementPotential(hiperlife::Tensor::tensor<double,1>& x, double Kconf)
{
    return (x(2) < 0) ? -2.0*Kconf*x(2) : 0.0;
}

#endif
