/*
    *******************************************************************************
    Copyright (c) 2017-2021 Universitat Politècnica de Catalunya
    Authors: Daniel Santos-Olivan, Alejandro Torres-Sanchez and Guillermo Vilanova
    Contributors:
    *******************************************************************************
    This file is part of hiperlife - High Performance Library for Finite Elements
    Project homepage: https://git.lacan.upc.edu/HPLFEgroup/hiperlifelib.git
    Distributed under the MIT software license, see the accompanying
    file LICENSE or http://www.opensource.org/licenses/mit-license.php.
    *******************************************************************************
*/


#ifndef AUXPATFOR_H
#define AUXPATFOR_H


#include <iostream>
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

using Teuchos::rcp;
using Teuchos::RCP;

struct MembParams
{
    enum RealParameters
    {
        deltat,         // dparam[2]
        fric,           // dparam[3]
        young,          // dparam[8]
        poisson,        // dparam[9]
        thick,          // dparam[10]
        force,          // dparam[6], dparam[92]

        fric_fact,      // dparam[25]
        fric2,          // dparam[30]
        aap,            // dparam[32]
        gamma,          // dparam[21]

        R,              // dparam[33]
        nu,
        kappa,          // dparam[34]
        bbn,            // dparam[35]
        width,          // dparam[36]
        width_kap,

        tens_factor,    // dparam[39]
        height_in,      // dparam[75]
        Kconf,          // dparam[60]

        f0,             // dparam[81]
        g_ratio,        // dparam[82]
        bc_const,       // dparam[83]

        mu,             // dparam[84]
        lambda,         // dparam[85]
        alpha,          // dparam[67]

        sp_gap,         // dparam[90]

        a11,            // dparam[96]
        a12,            // dparam[97]
        a33,            // dparam[93]
        kap1,           // dparam[99]

        gamma_minus,    // dparam[18]
        gamma_plus,     // dparam[19]
        gamma_l,        // dparam[20]

        pn,             // dparam[69]
        fric3,          // dparam[11]

        apical,         // dparam[63]
        basal,          // dparam[64]

        spring,         // dparam[59]
        forward_new,    // dparam[54]

        Xmax,           // dparam[61]

        k_p,            // dparam[89]
        k_d,            // dparam[94]

        factor,
        vol_inc,
        tens1,

        kspr,           // dparam[88]
        fact_elastic,   // dparam[55]

        crypt,
        gamma_plus_ref,
        gamma_minus_ref,

        ang,
        ang_old
    };

    enum IntParameters
    {
        tens_start,     // dparam[22]

        fric_start,
        gap_step,

        control_fric,  // dparam[37]

        vn_gap,

        vnstep,
        gap,

        timeStep,
        timestep,

        case_sphere,   // dparam[7]

        choice,        // dparam[87]

        time_target    // dparam[98]
    };

    HL_PARAMETER_LIST DefaultValues
    {
        {"deltat",0.0001},
        {"fric",0.1},
        {"young",1.0},
        {"poisson",0.25},
        {"thick",0.01},
        {"force",0.0},

        {"fric_fact",1.0},
        {"fric2",1.0},
        {"aap",1.0},
        {"gamma",1.0},

        {"R",1.0},
        {"kappa",0.0},
        {"bbn",0.0},
        {"width",0.0},

        {"tens_factor",1.0},
        {"height_in",1.0},
        {"Kconf",1.0},

        {"f0",1.0},
        {"g_ratio",1.0},
        {"bc_const",1.0},

        {"mu",1.0},
        {"lambda",1.0},
        {"alpha",1.0},

        {"sp_gap",1.0},

        {"a11",1.0},
        {"a12",1.0},
        {"a33",1.0},
        {"kap1",1.0},

        {"gamma_minus",1.0},
        {"gamma_plus",1.0},
        {"gamma_l",1.0},

        {"pn",1.0},
        {"fric3",1.0},

        {"apical",0.0},
        {"basal",0.0},

        {"spring",0.0},
        {"forward_new",0.0},

        {"Xmax",1.0},

        {"k_p",0.8333},
        {"k_d",0.8333},

        {"factor",1.0},
        {"vol_inc",1.0},
        {"tens1",1.0},

        {"kspr",0.0},
        {"fact_elastic",1.0},

        {"crypt",1.0},
        {"gamma_plus_ref",1.0},
        {"gamma_minus_ref",1.0},

        {"ang",0.0},
        {"ang_old",0.0},

        {"tens_start",-10000},
        {"control_fric",10000},
        {"case_sphere",0},
        {"choice",100},
        {"time_target",1}
    };
};




void LS_ReacDif(hiperlife::FillStructure& fillStr);
void LS_ED(hiperlife::FillStructure& fillStr);
void LS_Rho(hiperlife::FillStructure& fillStr);

void LS(hiperlife::FillStructure& fillStr);

#endif


void RHS_Border(hiperlife::FillStructure& fillStr);




double inline ConfinementPotential(hiperlife::Tensor::tensor<double,1>& x, double Kconf)
{
    double pot{};
    if(x(2)<0)
    {
        pot = -Kconf/3.0 * x(2)*x(2)*x(2);
    }

    return pot;
}

double inline dConfinementPotential(hiperlife::Tensor::tensor<double,1>& x, double Kconf)
{
    using namespace hiperlife;
    using namespace hiperlife::Tensor;
    double dpot{};

    if(x(2)<0)
    {
        dpot =-1*Kconf*x(2)*x(2);
    }



    return dpot;
}

double  inline ddConfinementPotential(hiperlife::Tensor::tensor<double,1>& x, double Kconf)
{
    using namespace hiperlife;
    using namespace hiperlife::Tensor;
    double ddpot{};
    if(x(2)<0)
    {
        ddpot = -2.0*Kconf*x(2);
    }

    return ddpot;
}

double inline zz_max(hiperlife::Tensor::tensor<double,1>& Z_cord1,double size_n)
{
    using namespace hiperlife;
    using namespace hiperlife::Tensor;
    // find minimum

    double  Z_max = Z_cord1(0);
    // search num in inputArray from index 0 to elementCount-1
    for (int i = 0; i < size_n; i++)
    {
        if (Z_cord1(i) > Z_max)
        {
            Z_max = Z_cord1(i);
        }
    }

    return  Z_max;
}

double inline ConfinementPotential1(hiperlife::Tensor::tensor<double,1>& x, double Kconf)
{
    double pot{};
    pot = Kconf/2.0 * x(0)*x(0);
    return pot;
}

double inline dConfinementPotential1(hiperlife::Tensor::tensor<double,1>& x, double Kconf)
{
    using namespace hiperlife;
    using namespace hiperlife::Tensor;
    double dpot{};

    dpot =Kconf*x(0);

    return dpot;
}

double  inline ddConfinementPotential1(hiperlife::Tensor::tensor<double,1>& x, double Kconf)
{
    using namespace hiperlife;
    using namespace hiperlife::Tensor;
    double ddpot{};

    ddpot = Kconf;


    return ddpot;
}


// Define a cross product for 3D tensors
// cross product for 3D tensors
// cross product for 3D tensors

