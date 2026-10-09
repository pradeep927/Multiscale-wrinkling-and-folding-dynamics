// =============================================================================
// Phenomenological viscoelastic shell model of an epithelial dome
// (Model C, Supplementary Note 6 of Chahare, Ouzeri, Wilson, Bal et al.,
//  "Multiscale wrinkling and folding dynamics in epithelial shells").
//
// A flat disc of radius R, clamped at its edge (simply supported: position
// fixed, rotation free), is
//   1) inflated  : enclosed volume ramped linearly from 0 to v_target*R*(area)
//                  over vnstep steps,
//   2) held      : volume fixed for gap steps (the material metric G and the
//                  material curvature relax towards the inflated shape),
//   3) deflated  : volume ramped down linearly over v_rn_step steps (stopped
//                  at the fraction v_cont), with a smaller time step,
//   4) relaxed   : volume kept at the final value.
// The time step is adaptive (stepFactor), and is reduced to control_dt from
// control_steps steps before the deflation starts.
//
// Problems solved every step:
//   hiperProbl (posDHand + gloDHand): Newton-Raphson for the shell position and
//              the volume Lagrange multiplier (element filling LS)
//   viscoProbl (viscoDHand)        : local update of G and kbar at every Gauss
//              point (LS_ReacDif; no global system, only element loop)
//   ENProbl    (EDHand)            : L2 projection of Gauss-point fields to
//              nodes, only at output steps (LS_ED)
// =============================================================================

#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <Teuchos_RCP.hpp>
#include <math.h>

#include "hl_TypeDefs.h"
#include "hl_Geometry.h"
#include "hl_DistributedMesh.h"
#include "hl_FillStructure.h"
#include "hl_DOFsHandler.h"
#include "hl_HiPerProblem.h"
#include "hl_ConsistencyCheck.h"
#include <hl_NonlinearSolver_NewtonRaphson.h>
#include <hl_MeshLoader.h>
#include <hl_ConfigFile.h>
#include <hl_LinearSolver_Direct_MUMPS.h>
#include "hl_Parser.h"
#include "hl_ParamStructure.h"
#include "hl_Core.h"

#include "AuxViscousInextMembrane.h"

int main (int argc, char *argv[])
{
    using namespace std;
    using namespace hiperlife;
    using Teuchos::rcp;
    using Teuchos::RCP;
    using namespace hiperlife::Tensor;

    hiperlife::Init(argc, argv);
    const int myRank = hiperlife::MyRank();

    if (argc < 2)
    {
        cout << "Configuration file not provided!" << endl;
        hiperlife::Finalize();
        return 1;
    }
    ConfigFile config(argv[1]);

    // ---------------------------------------------------------------------
    // Parameters (see config_shell.cfg for the meaning of each one)
    // ---------------------------------------------------------------------
    // geometry
    double R{1.0}, height_in{0.05}, mesh_refine{0.0};
    string prefixMesh;
    int choice{10};
    // mechanics
    double mu{7.8372}, lambda{16.6482}, mu_pre{1.5797}, lambda_pre{1.2147}, fact_elastic{0.01};
    double gamma{5.0}, fric{0.01}, fric2{0.015}, fric3{0.02}, Kconf{0.0};
    double fric_fact{1}, width{80.0}, thick{0.0};
    int control_fric{10000};
    // protocol
    double v_target{0.667}, v_rn_step{1000.0}, v_cont{0.98}, control_dt{0.0001};
    int vnstep{500}, gap{2000}, control_steps{400}, time_target{1};
    // time stepping and numerics
    double deltat{}, deltatMax{}, totalTime{}, tSimu{}, stepFactor{}, forward_new{0.0};
    int totalSteps{}, nPrint{10}, ConsCheck{}, gPts{3}, MAXITER_NR{};
    double SOLTOL_NR{}, RESTOL_NR{};
    {
        config.readInto(R, "R");
        config.readInto(height_in, "height_in");
        config.readInto(mesh_refine, "mesh_refine");
        config.readInto(prefixMesh, "prefixMesh");
        config.readInto(choice, "choice");

        config.readInto(mu, "mu");
        config.readInto(lambda, "lambda");
        config.readInto(mu_pre, "mu_pre");
        config.readInto(lambda_pre, "lambda_pre");
        config.readInto(fact_elastic, "fact_elastic");
        config.readInto(gamma, "gamma");
        config.readInto(fric, "fric");
        config.readInto(fric2, "fric2");
        config.readInto(fric3, "fric3");
        config.readInto(Kconf, "Kconf");
        config.readInto(fric_fact, "fric_fact");
        config.readInto(width, "width");
        config.readInto(thick, "thick");
        config.readInto(control_fric, "control_fric");

        config.readInto(v_target, "v_target");
        config.readInto(vnstep, "vnstep");
        config.readInto(gap, "gap");
        config.readInto(v_rn_step, "v_rn_step");
        config.readInto(v_cont, "v_cont");
        config.readInto(control_steps, "control_steps");
        config.readInto(control_dt, "control_dt");
        config.readInto(time_target, "time_target");

        config.readInto(deltat, "deltat");
        config.readInto(deltatMax, "deltatMax");
        config.readInto(totalTime, "totalTime");
        config.readInto(totalSteps, "totalSteps");
        config.readInto(tSimu, "tSimu");
        config.readInto(stepFactor, "stepFactor");
        config.readInto(forward_new, "forward_new");
        config.readInto(nPrint, "nPrint");
        config.readInto(ConsCheck, "ConsCheck");
        config.readInto(gPts, "gPts");
        config.readInto(MAXITER_NR, "MAXITER_NR");
        config.readInto(SOLTOL_NR, "SOLTOL_NR");
        config.readInto(RESTOL_NR, "RESTOL_NR");
    }
    v_target=v_target*(R);   // target volume per unit footprint area ~ v_target * R

    SmartPtr<ParamStructure> paramStr = CreateParamStructure<MembParams>();
    paramStr->setRealParameter(MembParams::deltat, deltat);
    paramStr->setRealParameter(MembParams::mu, mu);
    paramStr->setRealParameter(MembParams::lambda, lambda);
    paramStr->setRealParameter(MembParams::mu_pre, mu_pre);
    paramStr->setRealParameter(MembParams::lambda_pre, lambda_pre);
    paramStr->setRealParameter(MembParams::fact_elastic, fact_elastic);
    paramStr->setRealParameter(MembParams::gamma, gamma);
    paramStr->setRealParameter(MembParams::height_in, height_in);
    paramStr->setRealParameter(MembParams::fric, fric);
    paramStr->setRealParameter(MembParams::fric2, fric2);
    paramStr->setRealParameter(MembParams::fric3, fric3);
    paramStr->setRealParameter(MembParams::Kconf, Kconf);
    paramStr->setRealParameter(MembParams::fric_fact, fric_fact);
    paramStr->setRealParameter(MembParams::width, width);
    paramStr->setRealParameter(MembParams::thick, thick);
    paramStr->setRealParameter(MembParams::forward_new, forward_new);
    paramStr->setIntParameter(MembParams::time_target, time_target);
    paramStr->setIntParameter(MembParams::vnstep, vnstep);
    paramStr->setIntParameter(MembParams::gap, gap);
    paramStr->setIntParameter(MembParams::control_fric, control_fric);

    if (myRank == 0)
    {
        cout << endl << "Parameters:" << endl;
        cout << "  footprint radius R               : " << R << endl;
        cout << "  thickness h (bending)            : " << height_in << "   R/h = " << R/height_in << endl;
        cout << "  membrane Lame mu, lambda         : " << mu << " / " << lambda << endl;
        cout << "  pre-tension Lame mu, lambda      : " << mu_pre*fact_elastic << " / " << lambda_pre*fact_elastic << endl;
        cout << "  active pre-tension gamma         : " << gamma << endl;
        cout << "  drag fric                        : " << fric << endl;
        cout << "  in-plane viscosity fric2         : " << fric2 << "   (tau ~ fric2/mu = " << fric2/mu << ")" << endl;
        cout << "  bending viscosity fric3          : " << fric3 << endl;
        cout << "  confinement Kconf                : " << Kconf << endl;
        cout << "  inflation / hold / deflation     : " << vnstep << " / " << gap << " / " << v_rn_step << " steps" << endl;
        cout << "  drag x" << fric_fact << " below z=" << thick << " from step " << vnstep+gap+control_fric << " (width " << width << ")" << endl;
        cout << "  relaxation integrator            : " << (forward_new > 0.1 ? "explicit" : "implicit") << endl;
        cout << "  consistency check                : " << ConsCheck << endl << endl;
    }

    // ---------------------------------------------------------------------
    // Mesh: unit-size disc (radius 0.34 in the mesh file) scaled to radius R
    // ---------------------------------------------------------------------
    RCP<DistributedMesh> posDisMesh, gloDisMesh;
    {
        RCP<MeshLoader> mesh = rcp(new MeshLoader);
        mesh->setMesh(ElemType::Triang, BasisFuncType::SubdivSurfs, 2);
        mesh->loadMesh(prefixMesh + ".vtk", MeshType::Sequential);
        mesh->transformFree([R](double x, double y)
        {
            x = R/0.34 * x;
            y = R/0.34 * y;
            return std::make_tuple(x, y);
        });

        if (mesh_refine < 1.0)
        {
            posDisMesh = rcp(new DistributedMesh);
            posDisMesh->setMesh(mesh);
            posDisMesh->setBalanceMesh(true);
            posDisMesh->Update();
        }
        else
        {
            // one level of uniform h-refinement
            RCP<DistributedMesh> coarse = rcp(new DistributedMesh);
            coarse->setMesh(mesh);
            coarse->setBalanceMesh(true);
            coarse->Update();
            posDisMesh = rcp(new DistributedMesh);
            posDisMesh->setHRefinement(1);
            posDisMesh->setMeshRelation(MeshRelation::hRefin, coarse);
            posDisMesh->setBalanceMesh(true);
            posDisMesh->Update();
        }
        posDisMesh->printFileLegacyVtk("pos_mesh_final_No_refine");

        gloDisMesh = rcp(new DistributedMesh);
        gloDisMesh->setMeshRelation(MeshRelation::GlobConstr, posDisMesh);
        gloDisMesh->setBalanceMesh(false);
        gloDisMesh->Update();
    }

    // ---------------------------------------------------------------------
    // DOF handlers
    // ---------------------------------------------------------------------
    RCP<DOFsHandler> posDHand, gloDHand, viscoDHand, EDHand;
    try
    {
        // shell position
        posDHand = rcp(new DOFsHandler(posDisMesh));
        posDHand->setNameTag("posDHand");
        posDHand->setDOFs({"X", "Y", "Z"});
        posDHand->setNodeAuxF({"Ux", "Uy", "Uz", "nx", "ny", "nz", "height",
                               "EA", "EB", "DA", "DB", "PA", "trace_sigma", "dev_sigma"});
        posDHand->Update();

        // only provides the element/Gauss-point loop for LS_ReacDif
        viscoDHand = rcp(new DOFsHandler(posDisMesh));
        viscoDHand->setNameTag("viscoDHand");
        viscoDHand->setDOFs({"dummy"});
        viscoDHand->Update();

        // nodal projection of Gauss-point output fields
        EDHand = rcp(new DOFsHandler(posDisMesh));
        EDHand->setNameTag("EDHand");
        EDHand->setDOFs({"nx","ny","nz","height","EA","EB","DA","DB","PA","trace_sigma","dev_sigma"});
        EDHand->Update();

        // P: volume Lagrange multiplier (FX..MZ are not used, kept fixed at 0)
        gloDHand = rcp(new DOFsHandler(gloDisMesh));
        gloDHand->setNameTag("gloDHand");
        gloDHand->setDOFs({"P", "FX", "FY", "FZ", "MX", "MY", "MZ"});
        gloDHand->Update();
    }
    catch (runtime_error &err)
    {
        cout << myRank << ": DOFsHandler could not be created " << err.what() << endl;
        hiperlife::Finalize();
        return 1;
    }

    // initial condition: flat disc; edge (crease) nodes clamped
    for (int i = 0; i < posDisMesh->loc_nPts(); i++)
    {
        posDHand->nodeDOFs->setValue("X", i, IndexType::Local, posDisMesh->nodeCoord(i, 0, IndexType::Local));
        posDHand->nodeDOFs->setValue("Y", i, IndexType::Local, posDisMesh->nodeCoord(i, 1, IndexType::Local));
        posDHand->nodeDOFs->setValue("Z", i, IndexType::Local, posDisMesh->nodeCoord(i, 2, IndexType::Local));
        if (posDisMesh->nodeCrease(i, IndexType::Local) > 0)
        {
            posDHand->setConstraint(0, i, IndexType::Local, 0.0);
            posDHand->setConstraint(1, i, IndexType::Local, 0.0);
            posDHand->setConstraint(2, i, IndexType::Local, 0.0);
        }
    }
    posDHand->nodeDOFs0->setValue(posDHand->nodeDOFs);
    posDHand->UpdateGhosts();

    gloDHand->setInitialCondition(0, 0.0);
    for (int k = 1; k <= 6; k++)
        gloDHand->setConstraint(k, 0.0);
    gloDHand->UpdateGhosts();

    // Gauss-point storage (see AuxA / AuxB in the header)
    paramStr->a_aux.resize(posDisMesh->loc_nElem()*gPts*N_AUXA);
    paramStr->b_aux.resize(posDisMesh->loc_nElem()*gPts*N_AUXB);

    // ---------------------------------------------------------------------
    // Problems
    // ---------------------------------------------------------------------
    SmartPtr<HiPerProblem> hiperProbl = Create<HiPerProblem>();
    SmartPtr<HiPerProblem> viscoProbl = Create<HiPerProblem>();
    SmartPtr<HiPerProblem> ENProbl    = Create<HiPerProblem>();
    try
    {
        hiperProbl->setParameterStructure(paramStr);
        hiperProbl->setConsistencyCheckDelta(1.E-8);
        hiperProbl->setConsistencyCheckTolerance(1.E-4);
        hiperProbl->setDOFsHandlers({posDHand,gloDHand});
        hiperProbl->setIntegration("Integ", {"posDHand","gloDHand"});
        hiperProbl->setCubatureGauss("Integ", gPts);
        if (ConsCheck == 1)
            hiperProbl->setElementFillings("Integ", ConsistencyCheck<LS>);
        else
            hiperProbl->setElementFillings("Integ", LS);
        // (E_tension1 is always 0; kept so that the columns of globalIntegrals.dat do not change)
        hiperProbl->setGlobalIntegrals({"Energy","Dissipation","Trd","volume","area","area_n","Epress","E1","E_tension1","P_tension1","E_elastic"});
        hiperProbl->Update();

        viscoProbl->setParameterStructure(paramStr);
        viscoProbl->setDOFsHandlers({viscoDHand});
        viscoProbl->setIntegration("Integ", {"viscoDHand"});
        viscoProbl->setCubatureGauss("Integ", gPts);
        viscoProbl->setElementFillings("Integ", LS_ReacDif);
        viscoProbl->setGlobalIntegrals({"Diss1","Diss2"});
        viscoProbl->Update();

        ENProbl->setParameterStructure(paramStr);
        ENProbl->setDOFsHandlers({EDHand});
        ENProbl->setIntegration("Integ", {"EDHand"});
        ENProbl->setCubatureGauss("Integ", gPts);
        ENProbl->setElementFillings("Integ", LS_ED);
        ENProbl->Update();
    }
    catch (runtime_error& err)
    {
        cout << myRank << ": HiPerProblem could not be created " << err.what() << endl;
        hiperlife::Finalize();
        return 1;
    }

    // solvers
    auto makeMUMPS = [](SmartPtr<HiPerProblem> pb)
    {
        RCP<MUMPSDirectLinearSolver> s = rcp(new MUMPSDirectLinearSolver());
        s->setHiPerProblem(pb);
        s->setMatrixType(MUMPSDirectLinearSolver::MatrixType::General);
        s->setAnalysisType(MUMPSDirectLinearSolver::AnalysisType::Parallel);
        s->setOrderingLibrary(MUMPSDirectLinearSolver::OrderingLibrary::Auto);
        s->setVerbosity(MUMPSDirectLinearSolver::Verbosity::None);
        s->setDefaultParameters();
        s->setWorkSpaceMemoryIncrease(1000);
        s->Update();
        return s;
    };
    RCP<MUMPSDirectLinearSolver> linSolver   = makeMUMPS(hiperProbl);
    RCP<MUMPSDirectLinearSolver> ENDirSolver = makeMUMPS(ENProbl);

    RCP<NewtonRaphsonNonlinearSolver> nonlinSolver = rcp(new NewtonRaphsonNonlinearSolver());
    nonlinSolver->setLinearSolver(linSolver);
    nonlinSolver->setMaxNumIterations(MAXITER_NR);
    nonlinSolver->setResTolerance(RESTOL_NR);
    nonlinSolver->setSolTolerance(SOLTOL_NR);
    nonlinSolver->setLineSearch(false);   // always plain Newton (no backtracking)
    nonlinSolver->setConvRelTolerance(false);
    nonlinSolver->setPrintIntermInfo(true);
    nonlinSolver->setPrintSummary(false);
    nonlinSolver->setResMaximum(1E4);
    nonlinSolver->setSolMaximum(1E4);
    nonlinSolver->setExitRelMaximum(1E4);
    nonlinSolver->Update();

    // initialise the Gauss-point storage: mechanics -> G = g_R^{-1}, kbar = k_R -> mechanics
    hiperProbl->FillGlobalIntegrals();
    viscoProbl->FillGlobalIntegrals();
    hiperProbl->FillGlobalIntegrals();

    // footprint area: choice = 10 -> circular footprint (3.141 R^2), else the mesh area
    double a0 = hiperProbl->globalIntegral("area_n");
    double base_area = (choice == 10) ? 3.141*R*R : a0;
    v_target=v_target/a0*base_area*0.999;
    double a_res = a0-base_area;
    double v0 = v_target*a0;        // fully inflated volume
    if (myRank == 0)
    {
        cout << std::scientific << std::setprecision(6);
        cout << "  Initial volume: " << hiperProbl->globalIntegral("volume") << endl;
        cout << "  Target volume : " << abs(v_target)*a0 << endl;
        cout << "  Area          : " << a0 << endl;
    }

    // globalIntegrals.dat (same columns as before; Zmax, stress, strain are dummy)
    ofstream gIntegFile;
    if (myRank == 0)
    {
        gIntegFile.open("globalIntegrals.dat");
        gIntegFile << "TS time deltat";
        for (auto g: hiperProbl->globIntegralNames())
            gIntegFile << " " << g;
        gIntegFile << " pressure Zmax stress strain";
        for (auto g: viscoProbl->globalIntegralNames())
            gIntegFile << " " << g;
        gIntegFile << endl;
    }

    int timeStep = 0;
    posDHand->printFileLegacyVtk("sol_dis." + to_string(timeStep), true);

    double factor = 0.0;            // enclosed volume / reference area at this step
    double v_r_step = 0.0;          // deflation steps done
    double press0 = 0.0;
    const double deltatMax1 = deltatMax;

    // ---------------------------------------------------------------------
    // Time loop
    // ---------------------------------------------------------------------
    while ((tSimu < totalTime) and (timeStep < totalSteps))
    {
        if (myRank == 0)
            cout << "TS: " << timeStep + 1 << " Time " << tSimu << " of " << totalTime << " with deltat=" << deltat << ": del_max: " << deltatMax << endl;

        // (1) inflation
        if (timeStep < vnstep)
        {
            factor = (timeStep+1.0)/vnstep*v_target;
            if (myRank == 0)
                cout << "inflation: fraction " << (timeStep+1.0)/vnstep << endl;
        }
        // small time step from control_steps steps before the deflation
        if (timeStep > vnstep+gap-control_steps)
        {
            deltat *= stepFactor;
            deltatMax = control_dt;
            if (myRank == 0)
                cout << "time-step reduction: dt " << deltat << " max dt " << deltatMax << endl;
        }
        // after the deflation: back to the normal time step
        if (timeStep > vnstep+gap+v_cont*v_rn_step)
        {
            deltat /= stepFactor;
            deltatMax = deltatMax1;
            if (myRank == 0)
                cout << "relaxation phase: dt " << deltat << " max dt " << deltatMax << endl;
        }
        // (3) deflation
        if (timeStep > vnstep+gap)
        {
            if (v_r_step < v_cont*v_rn_step)
            {
                factor = (1-(v_r_step+1.0)/v_rn_step)*v_target;
                v_r_step = v_r_step+1;
                if (myRank == 0)
                    cout << "deflation step " << v_r_step << ": volume fraction " << factor/v_target << endl;
            }
        }

        paramStr->setRealParameter(MembParams::factor, factor);
        paramStr->setIntParameter(MembParams::timestep, timeStep);

        if (deltat < 0.000001)
        {
            if (myRank == 0)
                cout << "stopping: time step too small" << endl;
            break;
        }

        // Newton-Raphson for the shell
        posDHand->nodeDOFs->setValue(posDHand->nodeDOFs0);
        gloDHand->nodeDOFs->setValue(gloDHand->nodeDOFs0);
        hiperProbl->UpdateGhosts();

        bool converged = nonlinSolver->solve();

        if (converged)
        {
            // Gauss-point data at the converged state, then relax G and kbar
            hiperProbl->FillGlobalIntegrals();
            viscoProbl->FillGlobalIntegrals();

            tSimu += deltat;
            timeStep += 1;
            press0 = gloDHand->nodeDOFs->getValue(0,0)/deltat;

            if (timeStep % nPrint == 0)
            {
                // nodal output fields
                ENDirSolver->solve();
                ENProbl->UpdateSolution();
                for (int i = 0; i < posDHand->mesh->loc_nPts(); i++)
                {
                    for (int d = 0; d < 3; d++)
                        posDHand->nodeAuxF->setValue(d, i, IndexType::Local,
                            posDHand->nodeDOFs->getValue(d,i,IndexType::Local) - posDisMesh->nodeCoord(i, d, IndexType::Local));
                    for (int k = 0; k < 11; k++)
                        posDHand->nodeAuxF->setValue(3+k, i, IndexType::Local, EDHand->nodeDOFs->getValue(k,i,IndexType::Local));
                }
            }

            posDHand->nodeDOFs0->setValue(posDHand->nodeDOFs);
            hiperProbl->UpdateGhosts();

            if (nonlinSolver->numberOfIterations() < MAXITER_NR and deltat < deltatMax)
                deltat /= stepFactor;
        }
        else
        {
            if (myRank == 0)
                cout << "Newton-Raphson did not converge: reducing the time step" << endl;
            deltat *= stepFactor;
            posDHand->nodeDOFs->setValue(posDHand->nodeDOFs0);
            gloDHand->nodeDOFs->setValue(gloDHand->nodeDOFs0);
        }

        // screen summary
        if (myRank == 0)
        {
            double v1 = hiperProbl->globalIntegral("volume");
            cout << std::scientific << std::setprecision(6);
            cout << "  volume " << v1 << "   volume fraction " << v1/v0 << endl;
            cout << "  area " << hiperProbl->globalIntegral("area") << "   reference area " << hiperProbl->globalIntegral("area_n")
                 << "   stretch " << (hiperProbl->globalIntegral("area")-a_res)/(hiperProbl->globalIntegral("area_n")-a_res) << endl;
            cout << "  pressure " << -gloDHand->nodeDOFs->getValue(0,0)/deltat
                 << "   excess area " << abs(hiperProbl->globalIntegral("area")-hiperProbl->globalIntegral("area_n"))/base_area << endl;
            cout << "  membrane energy " << hiperProbl->globalIntegral("Energy") << "   bending energy " << hiperProbl->globalIntegral("E1")
                 << "   pre-tension elastic energy " << hiperProbl->globalIntegral("E_elastic") << endl;
            cout << "  drag dissipation " << hiperProbl->globalIntegral("Dissipation")
                 << "   tension power " << hiperProbl->globalIntegral("P_tension1")
                 << "   viscous dissipation (in-plane / bending) " << viscoProbl->globalIntegral("Diss1") << " / " << viscoProbl->globalIntegral("Diss2") << endl;
        }

        // output
        if (timeStep % nPrint == 0)
        {
            posDHand->printFileLegacyVtk("sol_dis." + to_string(timeStep), true);
            if (myRank == 0)
            {
                gIntegFile << timeStep << " " << tSimu << " " << deltat;
                for (auto g: hiperProbl->globIntegrals())
                    gIntegFile << " " << g;
                gIntegFile << " " << press0 << " " << 1.0 << " " << 0.0 << " " << 0.0;
                for (auto g: viscoProbl->globIntegrals())
                    gIntegFile << " " << g;
                gIntegFile << endl;
            }
        }

        paramStr->setRealParameter(MembParams::deltat, deltat);
    }

    hiperlife::Finalize();
    return 0;
}
