/*
    Active-gel bilayer shell model of an epithelial dome (continuum active gel
    shell of Chahare, Ouzeri, Wilson, Bal et al., "Multiscale wrinkling and
    folding dynamics in epithelial shells", Supplementary Note 1).

    Kinematics (Cosserat-type shell)
      posDHand : X,Y,Z   mid-surface
                 H1..H3  director h (apical-basal axis, |h| ~ thickness)
                 LAM     Lagrange multiplier of the cell-volume constraint
                         J (h . n) = height_in   (thickness x area conserved)
      gloDHand : P       Lagrange multiplier of the enclosed volume (pressure)
    The apical / basal cortices are the offset surfaces x +/- h/2 (aap, bbn
    select the offsets) with metrics g+ / g-; the lateral cortex is sampled
    along three in-plane directions m1, m2, m3 (120 degrees apart), each with a
    2x2 metric CL built from the stretch along m_i and the director.

    Each cortical surface is a viscoelastic active gel:
      - neo-Hookean energy (Lame mu = lambda = young/(2(1+poisson))) relative
        to its material metric G (GP apical, GN basal, GL1..GL3 lateral),
      - G relaxes with viscosity fric2 (LS_ReacDif, per Gauss point),
      - active tension gamma_plus * rho_api (apical), gamma_minus * rho_bas
        (basal), gamma_l * f0/3 * rho_lat_i (lateral),
      - cortical density rho with turnover d rho/dt = k_p - k_d rho - rho tr(d)
        (LS_Rho, steady state rho = k_p/k_d).
    Plus: drag (fric), confinement z >= 0 (Kconf), enclosed-volume constraint,
    and, late in the deflation, a spring (kspr) that anchors the parts of the
    shell touching the substrate (z < 1.05 height_in) to where they are.

    Per-Gauss-point storage in paramStr:
      a_aux (88): mechanics -> relaxation / density / output   (written by LS)
      b_aux (18): material metrics G and dissipation          (LS_ReacDif)
      c_aux (9) : densities rho (LS_Rho) and anchor positions  (LS_node_end)
*/

#ifndef AUX_BILAYER_SHELL_H
#define AUX_BILAYER_SHELL_H

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

struct MembParams
{
    enum RealParameters
    {
        deltat,
        young, poisson,          // cortex elasticity: mu = lambda = young/(2(1+poisson))
        fric2,                   // cortex viscosity (relaxation of G)
        fact_elastic,            // (kept: multiplies a reference-elastic energy that is only reported)
        gamma_plus, gamma_minus, // apical / basal active tension
        gamma_l,                 // lateral active tension scale
        f0,                      // lateral weight (f0/3 per lateral direction)
        k_p, k_d,                // cortical turnover (polymerisation / depolymerisation)
        height_in,               // initial thickness
        aap, bbn,                // offset factors of the apical / basal surfaces (1 = +/- h/2)
        fric,                    // drag
        Kconf,                   // confinement z >= 0
        kspr, width,             // late anchoring spring and width of its z-transition
        factor,                  // target volume / reference area (set every step)
        forward_new              // relaxation integrator: 0 implicit, 1 explicit
    };
    enum IntParameters
    {
        timestep,                // current step (set every step)
        time_target,             // G kept at the reference while timestep < time_target
        fric_start,              // = vnstep + gap (start of the deflation)
        control_fric             // spring active for timestep > fric_start + control_fric
    };
    HL_PARAMETER_LIST DefaultValues
    {
        {"deltat", 0.0001},
        {"young", 1.0}, {"poisson", 0.25}, {"fric2", 1.0}, {"fact_elastic", 1.0},
        {"gamma_plus", 1.0}, {"gamma_minus", 1.0}, {"gamma_l", 1.0}, {"f0", 1.0},
        {"k_p", 1.0}, {"k_d", 1.0},
        {"height_in", 1.0}, {"aap", 1.0}, {"bbn", 1.0},
        {"fric", 0.1}, {"Kconf", 1.0}, {"kspr", 0.0}, {"width", 0.0},
        {"factor", 0.0}, {"forward_new", 0.0},
        {"timestep", 0}, {"time_target", 1}, {"fric_start", 0}, {"control_fric", 6000}
    };
};

void LS(hiperlife::FillStructure& fillStr);          // mechanics (Newton problem)
void LS_ReacDif(hiperlife::FillStructure& fillStr);  // relaxation of the material metrics (local)
void LS_Rho(hiperlife::FillStructure& fillStr);      // cortical density turnover (linear solve)
void LS_ED(hiperlife::FillStructure& fillStr);       // L2 projection of output fields
void LS_node_end(hiperlife::FillStructure& fillStr); // stores anchor positions for the late spring

// one-sided cubic confinement keeping the shell above z = 0
double inline ConfinementPotential(hiperlife::Tensor::tensor<double,1>& x, double Kconf)
{
    double pot{};
    if(x(2)<0)
        pot = -Kconf/3.0 * x(2)*x(2)*x(2);
    return pot;
}
double inline dConfinementPotential(hiperlife::Tensor::tensor<double,1>& x, double Kconf)
{
    double dpot{};
    if(x(2)<0)
        dpot =-1*Kconf*x(2)*x(2);
    return dpot;
}
double inline ddConfinementPotential(hiperlife::Tensor::tensor<double,1>& x, double Kconf)
{
    double ddpot{};
    if(x(2)<0)
        ddpot = -2.0*Kconf*x(2);
    return ddpot;
}

#endif
