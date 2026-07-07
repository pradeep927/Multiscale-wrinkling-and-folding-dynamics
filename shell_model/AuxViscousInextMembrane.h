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
        fric,
        young,
        poisson,
        thick,
        force,
        fric_fact,
        fric2,
        aap,
        R,
        nu,
        kappa,
        bbn,
        width,
        width_kap,
        tens_factor,
        height_in,
        Kconf,
        f0,
        g_ratio,
        bc_const,
        mu,
        lambda,
        sp_gap,
        a11,
        a12,
        a33,
        kap1,
        gamma_minus,
        gamma_plus,
        gamma_l,
        gamma,
        pn,
        fric3,
        apical,
        basal,
        spring,
        forward_new,
        Xmax,
        k_p,
        k_d,
        factor,
        vol_inc,
        tens1,
        kspr,
        fact_elastic
    };
    enum IntParameters
    {
        tens_start,
        gap_step,
        control_fric,
        vn_gap,
        timeStep,
        timestep,
        case_sphere,
        choice,
        time_target,
        vnstep,
        gap
    };
    HL_PARAMETER_LIST DefaultValues
  {
                        {"deltat", 0.0001},
                        {"fric", 0.1},
                        {"young", 1.0},

                        {"poisson", 0.25},
                        {"force", 0.0},
                        {"thick", 0.01},

                        {"tens_start", -10000},
                        {"fric_fact",1.0},
                        {"fric2",1.0},
                        {"aap",1.0},
                        {"R",10.0},

                        {"kappa",0.0},
                        {"bbn",0.0},
                        {"width",0.0},
                        {"tens_factor", 1.0},
                        {"height_in", 1.0},
                        {"Kconf", 1.0},
                        {"f0", 1.0},
                        {"g_ratio", 1.0},
                        {"bc_const", 1.0},
                        {"ang", 1.0},
                        {"ang_old", 1.0},
                        {"sp_gap", 1.0},
                        {"a11", 1.0},
                        {"a12", 1.0},
                        {"a33", 1.0},
                        {"kap1", 1.0},
                        {"crypt", 1.0},
                        {"gamma_minus", 1.0},
                        {"gamma_plus", 1.0},
                        {"gamma_l", 1.0},
                        {"factor", 1.0},
                        {"vol_inc", 1.0},
                        {"tens1", 1.0},
                       {"fact_elastic", 1.0}

                    };
};











void LS_ReacDif(hiperlife::FillStructure& fillStr);
void LS_ED(hiperlife::FillStructure & fillStr);

void LS(hiperlife::FillStructure& fillStr);

#endif


void RHS_Border(hiperlife::FillStructure & fillStr);



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

