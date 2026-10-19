/*
    shell_model -- phenomenological viscoelastic shell ("Model C" of Chahare,
    Ouzeri, Wilson, Bal et al., "Multiscale wrinkling and folding dynamics in
    epithelial shells", Supplementary Note 6) for the inflation -> hold ->
    deflation of an epithelial dome on a circular footprint.
    Cleaned version of adam_shell: identical results (checked bit-for-bit).

    Unknowns
      posDHand : X, Y, Z            shell mid-surface (subdivision surfaces)
      gloDHand : P (+ 6 fixed)      Lagrange multiplier of the enclosed-volume
                                    constraint (P/deltat = luminal pressure)
    Internal (Gauss-point) state, paramStr->b_aux (11 per Gauss point)
      G^IJ      material (natural) metric, relaxes with viscosity fric2
      kbar_IJ   material curvature,        relaxes with viscosity fric3

    Energy per unit reference area:
      neo-Hookean membrane relative to G (mu, lambda)
      + weak neo-Hookean strain dependence of the pre-tension relative to the
        flat disc (moduli mu_pre, lambda_pre, times fact_elastic)
      + bending h^2/4 (k - kbar)(k - kbar) alpha,  h = height_in / J_n
      + constant active tension gamma (area term)
    plus drag on the mid-surface velocity (fric; fric*fric_fact near the base
    after the late-anchoring step fric_start + control_fric), a one-sided
    confinement below z = 0 (Kconf), an optional spring to the reference
    outside r = R - sp_gap (spring = 1, kspr) and the volume constraint.
    Every energy-like term is multiplied by deltat (incremental minimisation).

    Per-Gauss-point storage
      a_aux (58): written by LS at every assembly, read by LS_ReacDif / LS_ED
      b_aux (11): G, kbar (6) + output densities, written by LS_ReacDif
*/

#ifndef AUX_VISCO_SHELL_H
#define AUX_VISCO_SHELL_H

#include <iostream>
#include <vector>
#include <string>
#include <mpi.h>

#include "hl_DistributedMesh.h"
#include "hl_FillStructure.h"
#include "hl_DOFsHandler.h"
#include "hl_HiPerProblem.h"
#include "hl_ParamStructure.h"
#include "hl_Math.h"
#include "hl_Array.h"

using namespace std;
using namespace hiperlife;

constexpr int kNA = 58;   // a_aux per Gauss point
constexpr int kNB = 11;   // b_aux per Gauss point

struct MembParams
{
    enum RealParameters
    {
        deltat,
        mu, lambda,         // membrane moduli (relative to G)
        mu_pre, lambda_pre, // moduli of the strain dependence of the pre-tension
        fact_elastic,       // weight of the pre-tension strain dependence
        gamma,              // active tension (sigma_a)
        alpha,              // bending prefactor
        height_in,          // thickness
        fric,               // drag
        fric_fact,          // drag factor at the base after late anchoring
        width,              // sharpness of the base drag transition (tanh)
        thick,              // height of the base drag transition (z)
        fric2, fric3,       // viscosities of G and kbar
        forward_new,        // relaxation of G, kbar: 1 explicit, 0 implicit (local Newton)
        Kconf,              // confinement z >= 0
        R, sp_gap, kspr,    // optional spring outside r = R - sp_gap
        factor              // target volume / reference area (set every step)
    };
    enum IntParameters
    {
        timestep,
        time_target,        // G, kbar held at the flat reference while timestep < time_target
        fric_start,         // = vnstep + gap (start of the deflation)
        control_fric,       // late anchoring at timestep > fric_start + control_fric
        spring              // 1: spring on
    };
    HL_PARAMETER_LIST DefaultValues
    {
        {"deltat",0.0001}, {"mu",1.0}, {"lambda",1.0}, {"mu_pre",1.5797}, {"lambda_pre",1.2147}, {"fact_elastic",1.0}, {"gamma",1.0},
        {"alpha",1.0}, {"height_in",1.0}, {"fric",0.1}, {"fric_fact",1.0}, {"width",0.0},
        {"thick",0.01}, {"fric2",1.0}, {"fric3",0.0}, {"forward_new",0.0}, {"Kconf",0.0},
        {"R",1.0}, {"sp_gap",1.0}, {"kspr",0.0}, {"factor",0.0},
        {"timestep",0}, {"time_target",2}, {"fric_start",0}, {"control_fric",10000}, {"spring",0}
    };
};

void LS(hiperlife::FillStructure& fillStr);          // mechanics + volume constraint
void LS_ReacDif(hiperlife::FillStructure& fillStr);  // update of G, kbar (local, per Gauss point)
void LS_ED(hiperlife::FillStructure& fillStr);       // L2 projection of Gauss-point fields (output)

// one-sided cubic confinement keeping the shell above z = 0
double inline ConfinementPotential(hiperlife::Tensor::tensor<double, 1>& x, double Kconf)
{
    double pot{};
    if (x(2) < 0) pot = -Kconf / 3.0 * x(2) * x(2) * x(2);
    return pot;
}
double inline dConfinementPotential(hiperlife::Tensor::tensor<double, 1>& x, double Kconf)
{
    double dpot{};
    if (x(2) < 0) dpot = -1 * Kconf * x(2) * x(2);
    return dpot;
}
double inline ddConfinementPotential(hiperlife::Tensor::tensor<double, 1>& x, double Kconf)
{
    double ddpot{};
    if (x(2) < 0) ddpot = -2.0 * Kconf * x(2);
    return ddpot;
}

// ---- restart files of the Gauss-point state b_aux ---------------------------
// one binary file per rank, <prefix>.<rank>.gp; elements matched by the
// centroid of their nodes (tolerance tol; same number of MPI processes)
void GaussPointStateWrite(hiperlife::SmartPtr<hiperlife::DistributedMesh>& mesh,
    const std::vector<double>& b_aux, int gPts, const std::string& prefix);
bool GaussPointStateRead(hiperlife::SmartPtr<hiperlife::DistributedMesh>& mesh,
    std::vector<double>& b_aux, int gPts, const std::string& prefix, double tol);

#endif