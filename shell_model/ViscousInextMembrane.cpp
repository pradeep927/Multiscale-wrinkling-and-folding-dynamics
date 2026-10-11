/*
    shell_model -- viscoelastic shell of an epithelial dome (see
    AuxViscousInextMembrane.h for the model).

    Protocol (step counts, timestep = 0, 1, ...):
      inflation    timestep < vnstep          target volume ramped 0 -> v_target
      hold         vnstep .. vnstep+gap       G and kbar relax at fixed volume
      deflation    after vnstep+gap           volume ramped down over v_rn_step steps
                   (v_cont of it), dt reduced from vnstep+gap-control_steps
                   with maximum control_dt
      relaxation   after vnstep+gap+v_cont*v_rn_step (dt grows back)
    The clamped boundary is the mesh crease (outer rim) (+ optional zones,
    choice = 1000 / 2000). From fric_start + control_fric the nodes with
    z < thick in the current (start-of-step) configuration are clamped too.

    Output every nPrint steps: sol_dis.<step>.vtk, one line of globalIntegrals.dat.
    Restart files every hl_print steps: sol_pos.<step>.<rank>.hl,
    sol_gloCons.<step>.<rank>.hl, sol_gp.<step>.<rank>.gp (Gauss-point state:
    relaxed metric and curvature) and restart_info.<step>. Restart with
    restart = <step> and the same number of MPI processes.

    Cleaned version of adam_shell: the same results (bit-for-bit) with the
    same parameters.
*/
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <mpi.h>

#include "hl_Core.h"
#include "hl_TypeDefs.h"
#include "hl_DistributedMesh.h"
#include "hl_FillStructure.h"
#include "hl_DOFsHandler.h"
#include "hl_HiPerProblem.h"
#include "hl_ConsistencyCheck.h"
#include <hl_NonlinearSolver_NewtonRaphson.h>
#include <hl_MeshLoader.h>
#include <hl_ConfigFile.h>
#include <hl_LinearSolver_Direct_MUMPS.h>
#include "hl_ParamStructure.h"
#include "hl_Parser.h"
#include "AuxViscousInextMembrane.h"

int main(int argc, char** argv)
{
    using namespace std;
    using namespace hiperlife;

    hiperlife::Init(argc, argv);
    const int myRank = hiperlife::MyRank();
    if (argc < 2)
    {
        if (myRank == 0) cout << "Configuration file not provided!" << endl;
        hiperlife::Finalize();
        return 1;
    }
    ConfigFile config(argv[1]);
    auto rd = [&](auto& v, const char* name) { config.readInto(v, name); };

    // ---- mesh ----------------------------------------------------------------
    string prefixMesh;
    double R{1.0}, mesh_refine{0.0};
    int choice{100}, gPts{3};
    rd(prefixMesh, "prefixMesh"); rd(R, "R"); rd(mesh_refine, "mesh_refine");
    rd(choice, "choice"); rd(gPts, "gPts");

    // ---- solver ----------------------------------------------------------------
    int MAXITER_NR{15}, ConsCheck{0};
    double RESTOL_NR{5e-7}, SOLTOL_NR{5e-7};
    rd(MAXITER_NR, "MAXITER_NR"); rd(RESTOL_NR, "RESTOL_NR"); rd(SOLTOL_NR, "SOLTOL_NR");
    rd(ConsCheck, "ConsCheck");

    // ---- time stepping, output, restart ------------------------------------------
    double deltat0{0.0001}, deltatMax{0.0}, totalTime{0.0}, tSimu{0.0}, stepFactor{0.99};
    int totalSteps{0}, nPrint{10}, hl_print{1000}, restart{0};
    rd(deltat0, "deltat"); rd(deltatMax, "deltatMax"); rd(totalTime, "totalTime");
    rd(totalSteps, "totalSteps"); rd(tSimu, "tSimu"); rd(stepFactor, "stepFactor");
    rd(nPrint, "nPrint"); rd(hl_print, "hl_print"); rd(restart, "restart");

    // ---- protocol ----------------------------------------------------------------
    double v_target{0.667}, v_cont{0.98}, v_rn_step{1000.0}, control_dt{0.0001};
    int vnstep{500}, gap{2000}, control_steps{400}, time_target{1}, control_fric{10000};
    rd(v_target, "v_target"); rd(vnstep, "vnstep"); rd(gap, "gap");
    rd(v_rn_step, "v_rn_step"); rd(v_cont, "v_cont");
    rd(control_steps, "control_steps"); rd(control_dt, "control_dt");
    rd(time_target, "time_target"); rd(control_fric, "control_fric");

    // ---- mechanics -----------------------------------------------------------------
    double mu{7.8372}, lambda{16.6482}, gamma{1.0}, mu_pre{1.5797}, lambda_pre{1.2147},
           fact_elastic{0.01}, alpha{1.0}, height_in{0.05}, fric2{0.015}, fric3{0.02},
           forward_new{0.0}, fric{0.01}, Kconf{0.0}, fric_fact{1.0}, width{80.0}, thick{0.0},
           kspr{0.0}, sp_gap{0.058};
    int spring{0};
    rd(mu, "mu"); rd(lambda, "lambda"); rd(gamma, "gamma");
    rd(mu_pre, "mu_pre"); rd(lambda_pre, "lambda_pre"); rd(fact_elastic, "fact_elastic");
    rd(alpha, "alpha"); rd(height_in, "height_in");
    rd(fric2, "fric2"); rd(fric3, "fric3"); rd(forward_new, "forward_new");
    rd(fric, "fric"); rd(Kconf, "Kconf");
    rd(fric_fact, "fric_fact"); rd(width, "width"); rd(thick, "thick");
    rd(spring, "spring"); rd(kspr, "kspr"); rd(sp_gap, "sp_gap");

    // footprint of choice = 1000 (three-lobed, mesh units)
    const double r0_tri = 0.18, r_tri = 0.45, q_tri = 1.7;
    const int fric_start = vnstep + gap;

    // restart: state of the protocol at the saved step
    double vstep{0.0}, v_r_step{0.0};
    const double deltatMax1 = deltatMax;
    if (restart != 0)
    {
        std::ifstream rinfo("restart_info." + to_string(restart));
        if (!(rinfo >> tSimu >> deltat0 >> deltatMax >> vstep >> v_r_step))
        {
            if (myRank == 0) cout << "ERROR: restart_info." << restart << " not found" << endl;
            hiperlife::Finalize();
            return 1;
        }
    }

    SmartPtr<ParamStructure> paramStr = CreateParamStructure<MembParams>();
    paramStr->setRealParameter(MembParams::deltat, deltat0);
    paramStr->setRealParameter(MembParams::mu, mu);
    paramStr->setRealParameter(MembParams::lambda, lambda);
    paramStr->setRealParameter(MembParams::mu_pre, mu_pre);
    paramStr->setRealParameter(MembParams::lambda_pre, lambda_pre);
    paramStr->setRealParameter(MembParams::fact_elastic, fact_elastic);
    paramStr->setRealParameter(MembParams::gamma, gamma);
    paramStr->setRealParameter(MembParams::alpha, alpha);
    paramStr->setRealParameter(MembParams::height_in, height_in);
    paramStr->setRealParameter(MembParams::fric, fric);
    paramStr->setRealParameter(MembParams::fric_fact, fric_fact);
    paramStr->setRealParameter(MembParams::width, width);
    paramStr->setRealParameter(MembParams::thick, thick);
    paramStr->setRealParameter(MembParams::fric2, fric2);
    paramStr->setRealParameter(MembParams::fric3, fric3);
    paramStr->setRealParameter(MembParams::forward_new, forward_new);
    paramStr->setRealParameter(MembParams::Kconf, Kconf);
    paramStr->setRealParameter(MembParams::R, R);
    paramStr->setRealParameter(MembParams::sp_gap, sp_gap);
    paramStr->setRealParameter(MembParams::kspr, kspr);
    paramStr->setIntParameter(MembParams::timestep, restart);
    paramStr->setIntParameter(MembParams::time_target, time_target);
    paramStr->setIntParameter(MembParams::fric_start, fric_start);
    paramStr->setIntParameter(MembParams::control_fric, control_fric);
    paramStr->setIntParameter(MembParams::spring, spring);

    if (myRank == 0)
    {
        cout << endl << "Parameters:" << endl
             << "  R " << R << ", thickness " << height_in << " (R/h = " << R/height_in << ")" << endl
             << "  mu " << mu << ", lambda " << lambda << ", tension gamma " << gamma << endl
             << "  pre-tension moduli " << mu_pre << ", " << lambda_pre << " x fact_elastic " << fact_elastic << endl
             << "  viscosities fric2 " << fric2 << ", fric3 " << fric3 << ", drag fric " << fric << endl
             << "  inflation " << vnstep << " steps, hold " << gap << ", deflation " << v_rn_step
             << " (x " << v_cont << "), late anchoring from step " << fric_start + control_fric << endl
             << "  Kconf " << Kconf << ", spring " << spring << " (kspr " << kspr << ", sp_gap " << sp_gap << ")" << endl
             << "  consistency check " << ConsCheck << endl << endl;
    }

    // ======================================================================
    // Mesh and DOFs handlers
    // ======================================================================
    // posDHand AuxF: displacement, outputs of LS_ED, and "fixed" (1 = clamped
    // node, kept so that a restart re-applies the same constraints)
    const std::vector<string> posAux = {"Ux","Uy","Uz","nx","ny","nz","height","EA","EB","DA","DB","PA",
                                        "trace_sigma","dev_sigma","fixed"};
    const std::vector<string> edDOFs = {"nx","ny","nz","height","EA","EB","DA","DB","PA","trace_sigma","dev_sigma"};

    SmartPtr<DistributedMesh> posDisMesh, gloDisMesh;
    SmartPtr<DOFsHandler> posDHand, gloDHand, viscoDHand, EDHand;

    if (restart == 0)
    {
        SmartPtr<MeshLoader> mesh = Create<MeshLoader>();
        mesh->setMesh(ElemType::Triang, BasisFuncType::SubdivSurfs, 2);
        mesh->loadMesh(prefixMesh + ".vtk", MeshType::Sequential);
        // the mesh files are discs of radius 0.34
        mesh->transformFree([R](double x, double y)
        {
            x = R/0.34 * x;
            y = R/0.34 * y;
            return std::make_tuple(x, y);
        });
        if (mesh_refine < 1.0)
        {
            posDisMesh = Create<DistributedMesh>();
            posDisMesh->setMesh(mesh);
            posDisMesh->setBalanceMesh(true);
            posDisMesh->Update();
        }
        else
        {
            SmartPtr<DistributedMesh> coarse = Create<DistributedMesh>();
            coarse->setMesh(mesh);
            coarse->setBalanceMesh(true);
            coarse->Update();
            posDisMesh = Create<DistributedMesh>();
            posDisMesh->setHRefinement(1);
            posDisMesh->setMeshRelation(MeshRelation::hRefin, coarse);
            posDisMesh->setBalanceMesh(true);
            posDisMesh->Update();
        }
        gloDisMesh = Create<DistributedMesh>();
        gloDisMesh->setMeshRelation(MeshRelation::GlobConstr, posDisMesh);
        gloDisMesh->setBalanceMesh(false);
        gloDisMesh->Update();

        posDHand = Create<DOFsHandler>(posDisMesh);
        posDHand->setNameTag("posDHand");
        posDHand->setDOFs({"X", "Y", "Z"});
        posDHand->setNodeAuxF(posAux);
        posDHand->Update();
        for (int i = 0; i < posDisMesh->loc_nPts(); i++)
        {
            const double x = posDisMesh->nodeCoord(i, 0, IndexType::Local);
            const double y = posDisMesh->nodeCoord(i, 1, IndexType::Local);
            const double z = posDisMesh->nodeCoord(i, 2, IndexType::Local);
            posDHand->nodeDOFs->setValue("X", i, IndexType::Local, x);
            posDHand->nodeDOFs->setValue("Y", i, IndexType::Local, y);
            posDHand->nodeDOFs->setValue("Z", i, IndexType::Local, z);

            // clamped nodes: the rim (mesh crease) + optional zones
            bool fixed = posDisMesh->nodeCrease(i, IndexType::Local) > 0;
            if (choice == 1000)
            {
                if (y < r0_tri - r_tri + sp_gap) fixed = true;
                if (y - sqrt(3)*x - r0_tri - 2*r_tri + sp_gap > 0.0) fixed = true;
                if (y + sqrt(3)*x - r0_tri - 2*r_tri - sqrt(3)*q_tri + sp_gap > 0.0) fixed = true;
                if (y > sqrt(3)/2*q_tri + r0_tri + r_tri*0.5 - sp_gap &&
                    sqrt((x-0.5*q_tri)*(x-0.5*q_tri) + (y-r0_tri-sqrt(3)/2*q_tri)*(y-r0_tri-sqrt(3)/2*q_tri)) > r_tri - sp_gap) fixed = true;
                if (x > q_tri - sp_gap && y < r0_tri + 0.5*r_tri + sp_gap &&
                    sqrt((x-q_tri)*(x-q_tri) + (y-r0_tri)*(y-r0_tri)) > r_tri - sp_gap) fixed = true;
                if (x < 0.0 + sp_gap && y < r0_tri + 0.5*r_tri + sp_gap &&
                    sqrt(x*x + (y-r0_tri)*(y-r0_tri)) > r_tri - sp_gap) fixed = true;
            }
            if (choice == 2000 && sqrt(x*x + y*y) > R - sp_gap) fixed = true;
            posDHand->nodeAuxF->setValue("fixed", i, IndexType::Local, fixed ? 1.0 : 0.0);
        }
        posDHand->nodeDOFs0->setValue(posDHand->nodeDOFs);

        // pressure P (+ 6 rigid-body multipliers, kept fixed at 0)
        gloDHand = Create<DOFsHandler>(gloDisMesh);
        gloDHand->setNameTag("gloDHand");
        gloDHand->setDOFs({"P", "FX", "FY", "FZ", "MX", "MY", "MZ"});
        gloDHand->Update();
        gloDHand->setInitialCondition(0, 0.0);
    }
    else
    {
        if (myRank == 0) cout << "Starting from restart " << restart << " (t = " << tSimu << ", dt = " << deltat0 << ")" << endl;
        try
        {
            posDHand = Create<DOFsHandler>("posDHand");
            posDHand->setFilePrefix("sol_pos." + to_string(restart), OutputMode::Text);
            posDHand->setNumNodeAuxF(int(posAux.size()));
            posDHand->Update();
            posDisMesh = posDHand->mesh;
            // (a GlobConstr mesh cannot be built on a mesh read from file:
            //  gloDHand is read from its own file)
            gloDHand = Create<DOFsHandler>("gloDHand");
            gloDHand->setFilePrefix("sol_gloCons." + to_string(restart), OutputMode::Text);
            gloDHand->Update();
            gloDisMesh = gloDHand->mesh;
        }
        catch (runtime_error& err)
        {
            cout << myRank << ": restart failed: " << err.what() << endl;
            hiperlife::Finalize();
            return 1;
        }
        posDHand->nodeDOFs0->setValue(posDHand->nodeDOFs);
        gloDHand->nodeDOFs0->setValue(gloDHand->nodeDOFs);
    }
    for (int i = 0; i < posDisMesh->loc_nPts(); i++)
        if (posDHand->nodeAuxF->getValue("fixed", i, IndexType::Local) > 0.5)
            for (int d = 0; d < 3; d++)
                posDHand->setConstraint(d, i, IndexType::Local, 0.0);
    posDHand->UpdateGhosts();
    for (int k = 1; k < 7; k++)
        gloDHand->setConstraint(k, 0.0);
    gloDHand->UpdateGhosts();

    // carrier of the Gauss-point loop of LS_ReacDif (its DOF is never used)
    viscoDHand = Create<DOFsHandler>(posDisMesh);
    viscoDHand->setNameTag("viscoDHand");
    viscoDHand->setDOFs({"g"});
    viscoDHand->Update();
    viscoDHand->UpdateGhosts();

    // L2 projection of Gauss-point quantities (output)
    EDHand = Create<DOFsHandler>(posDisMesh);
    EDHand->setNameTag("EDHand");
    EDHand->setDOFs(edDOFs);
    EDHand->Update();
    EDHand->UpdateGhosts();

    paramStr->a_aux.assign(size_t(posDisMesh->loc_nElem())*gPts*kNA, 0.0);
    paramStr->b_aux.assign(size_t(posDisMesh->loc_nElem())*gPts*kNB, 0.0);
    if (restart != 0 && !GaussPointStateRead(posDisMesh, paramStr->b_aux, gPts, "sol_gp." + to_string(restart), 1.0e-6*R))
    {
        if (myRank == 0) cout << "ERROR: sol_gp." << restart << ".<rank>.gp missing or not matching the mesh" << endl;
        hiperlife::Finalize();
        return 1;
    }

    // ======================================================================
    // Problems and solvers
    // ======================================================================
    SmartPtr<HiPerProblem> hiperProbl = Create<HiPerProblem>();
    hiperProbl->setParameterStructure(paramStr);
    hiperProbl->setConsistencyCheckDelta(1.E-6);
    hiperProbl->setConsistencyCheckTolerance(1.E-4);
    hiperProbl->setDOFsHandlers({posDHand, gloDHand});
    hiperProbl->setIntegration("Integ", {"posDHand", "gloDHand"});
    hiperProbl->setCubatureGauss("Integ", gPts);
    if (ConsCheck == 1)
        hiperProbl->setElementFillings("Integ", ConsistencyCheck<LS>);
    else
        hiperProbl->setElementFillings("Integ", LS);
    hiperProbl->setGlobalIntegrals({"Energy","Dissipation","Trd","volume","area","area_n","Epress","E1",
                                    "E_tension1","P_tension1","E_elastic"});
    hiperProbl->Update();

    SmartPtr<HiPerProblem> viscoProbl = Create<HiPerProblem>();
    viscoProbl->setParameterStructure(paramStr);
    viscoProbl->setDOFsHandlers({viscoDHand});
    viscoProbl->setIntegration("Integ", {"viscoDHand"});
    viscoProbl->setCubatureGauss("Integ", gPts);
    viscoProbl->setElementFillings("Integ", LS_ReacDif);
    viscoProbl->setGlobalIntegrals({"Diss1","Diss2"});
    viscoProbl->Update();

    SmartPtr<HiPerProblem> ENProbl = Create<HiPerProblem>();
    ENProbl->setParameterStructure(paramStr);
    ENProbl->setDOFsHandlers({EDHand});
    ENProbl->setIntegration("Integ", {"EDHand"});
    ENProbl->setCubatureGauss("Integ", gPts);
    ENProbl->setElementFillings("Integ", LS_ED);
    ENProbl->Update();

    // (memIncrease: % extra MUMPS workspace)
    auto makeMumps = [](SmartPtr<HiPerProblem>& probl, int memIncrease)
    {
        SmartPtr<MUMPSDirectLinearSolver> s = Create<MUMPSDirectLinearSolver>();
        s->setHiPerProblem(probl);
        s->setMatrixType(MUMPSDirectLinearSolver::MatrixType::General);
        s->setAnalysisType(MUMPSDirectLinearSolver::AnalysisType::Parallel);
        s->setOrderingLibrary(MUMPSDirectLinearSolver::OrderingLibrary::Auto);
        s->setVerbosity(MUMPSDirectLinearSolver::Verbosity::None);
        s->setDefaultParameters();
        s->setWorkSpaceMemoryIncrease(memIncrease);
        s->Update();
        return s;
    };
    SmartPtr<MUMPSDirectLinearSolver> linSolver = makeMumps(hiperProbl, 1000);
    SmartPtr<MUMPSDirectLinearSolver> ENDirSolver = makeMumps(ENProbl, 100);   // small: the output projection only

    SmartPtr<NewtonRaphsonNonlinearSolver> nonlinSolver = Create<NewtonRaphsonNonlinearSolver>();
    nonlinSolver->setLinearSolver(linSolver);
    nonlinSolver->setMaxNumIterations(MAXITER_NR);
    nonlinSolver->setResTolerance(RESTOL_NR);
    nonlinSolver->setSolTolerance(SOLTOL_NR);
    nonlinSolver->setLineSearch(false);
    nonlinSolver->setConvRelTolerance(false);
    nonlinSolver->setPrintIntermInfo(true);
    nonlinSolver->setPrintSummary(false);
    nonlinSolver->setResMaximum(1E4);
    nonlinSolver->setSolMaximum(1E4);
    nonlinSolver->setExitRelMaximum(1E4);
    nonlinSolver->Update();

    // ED solve and copy of the outputs to posDHand's AuxF
    auto computeOutput = [&]()
    {
        ENDirSolver->solve();
        ENProbl->UpdateSolution();
        for (int i = 0; i < posDisMesh->loc_nPts(); i++)
        {
            for (int d = 0; d < 3; d++)
                posDHand->nodeAuxF->setValue(d, i, IndexType::Local,
                    posDHand->nodeDOFs->getValue(d, i, IndexType::Local) - posDisMesh->nodeCoord(i, d, IndexType::Local));
            for (int k = 0; k < int(edDOFs.size()); k++)
                posDHand->nodeAuxF->setValue(3 + k, i, IndexType::Local, EDHand->nodeDOFs->getValue(k, i, IndexType::Local));
        }
        hiperProbl->UpdateGhosts();
    };

    // initial assemblies: a_aux (LS), then G, kbar at the reference (LS_ReacDif)
    // (on a restart G, kbar come from sol_gp and only a_aux is rebuilt)
    hiperProbl->FillLinearSystem();
    if (restart == 0)
    {
        viscoProbl->FillLinearSystem();
        hiperProbl->FillLinearSystem();
    }

    // target volume: v_target * R * base area
    const double area_n = hiperProbl->globalIntegral("area_n");
    double base_area = area_n;
    if (choice == 1000) base_area = 3.141*r_tri*r_tri + 3*q_tri*r_tri + sqrt(3)/4*q_tri*q_tri;
    if (choice == 2000) base_area = 3.141*R*R;
    v_target = v_target*R;
    v_target = v_target/area_n*base_area*0.999;
    const double v0 = v_target*area_n;
    const double a_res = area_n - base_area;
    if (myRank == 0)
    {
        cout << std::scientific << std::setprecision(6);
        cout << "  Volume: " << hiperProbl->globalIntegral("volume") << endl;
        cout << "Target volume:  " << std::abs(v_target)*area_n << endl;
        cout << "  Area : phi: " << area_n << ": " << hiperProbl->globalIntegral("area") << endl;
    }

    ofstream gIntegFile;
    if (myRank == 0)
    {
        gIntegFile.open("globalIntegrals.dat", restart == 0 ? std::ios::out : std::ios::app);
        if (restart == 0)
        {
            gIntegFile << "TS time deltat";
            for (auto g : hiperProbl->globalIntegralNames()) gIntegFile << " " << g;
            gIntegFile << " pressure";
            for (auto g : viscoProbl->globalIntegralNames()) gIntegFile << " " << g;
            gIntegFile << endl;
        }
    }

    double& deltat  = paramStr->getRealParameter(MembParams::deltat);
    int& timestep   = paramStr->getIntParameter(MembParams::timestep);
    double& factor  = paramStr->getRealParameter(MembParams::factor);
    double press0 = 0.0;

    if (restart == 0)
    {
        computeOutput();
        posDHand->printFileLegacyVtk("sol_dis.0", true);
    }

    // ======================================================================
    // Time loop
    // ======================================================================
    while (tSimu < totalTime && timestep < totalSteps)
    {
        if (myRank == 0)
            cout << "TS: " << timestep + 1 << " Time " << tSimu << " of " << totalTime
                 << " with deltat=" << deltat << ": del_max: " << deltatMax << endl;

        // inflation
        if (timestep < vnstep)
        {
            factor = (timestep + 1.0)/vnstep*v_target;
            vstep = vstep + 1;
        }
        // smaller time step before and during the deflation
        if (timestep > vnstep + gap - control_steps)
        {
            deltat *= stepFactor;
            deltatMax = control_dt;
        }
        // final relaxation: the time step grows back
        if (timestep > vnstep + gap + v_cont*v_rn_step)
        {
            deltat /= stepFactor;
            deltatMax = deltatMax1;
        }
        // deflation
        if (timestep > vnstep + gap && v_r_step < v_cont*v_rn_step)
        {
            factor = (1 - (v_r_step + 1.0)/v_rn_step)*v_target;
            v_r_step = v_r_step + 1;
            if (myRank == 0) cout << "volume reduction step: " << v_r_step << " by factor: " << factor/v_target << endl;
        }
        if (deltat < 0.000001)
        {
            if (myRank == 0) cout << "stopping: dt too small" << endl;
            break;
        }
        // late anchoring: clamp the nodes close to the substrate
        if (timestep > fric_start + control_fric)
        {
            for (int i = 0; i < posDisMesh->loc_nPts(); i++)
                if (posDHand->nodeDOFs0->getValue(2, i, IndexType::Local) < thick)
                {
                    for (int d = 0; d < 3; d++)
                        posDHand->setConstraint(d, i, IndexType::Local, 0.0);
                    posDHand->nodeAuxF->setValue("fixed", i, IndexType::Local, 1.0);
                }
            posDHand->UpdateGhosts();
        }

        posDHand->nodeDOFs->setValue(posDHand->nodeDOFs0);
        gloDHand->nodeDOFs->setValue(gloDHand->nodeDOFs0);
        hiperProbl->UpdateGhosts();

        bool converged = nonlinSolver->solve();
        if (converged)
        {
            // a_aux at the converged state, then relaxation of G, kbar (b_aux)
            hiperProbl->FillLinearSystem();
            viscoProbl->FillLinearSystem();
            tSimu += deltat;
            timestep += 1;
            press0 = gloDHand->nodeDOFs->getValue(0, 0)/deltat;   // P/dt of this step
            posDHand->nodeDOFs0->setValue(posDHand->nodeDOFs);
            gloDHand->nodeDOFs0->setValue(gloDHand->nodeDOFs);
            hiperProbl->UpdateGhosts();

            // larger dt for the next step after an easy convergence
            if (nonlinSolver->numberOfIterations() < MAXITER_NR && deltat < deltatMax)
                deltat /= stepFactor;

            const bool doPrint = (timestep % nPrint == 0);
            const bool doRestart = (hl_print > 0 && timestep % hl_print == 0);
            if (doPrint || doRestart)
                computeOutput();

            if (doPrint)
            {
                posDHand->printFileLegacyVtk("sol_dis." + to_string(timestep), true);
                if (myRank == 0)
                {
                    gIntegFile << timestep << " " << tSimu << " " << deltat;
                    for (auto g : hiperProbl->globalIntegrals()) gIntegFile << " " << g;
                    gIntegFile << " " << press0;
                    for (auto g : viscoProbl->globalIntegrals()) gIntegFile << " " << g;
                    gIntegFile << endl;
                }
            }
            if (doRestart)
            {
                const string s = to_string(timestep);
                posDHand->printFile("sol_pos." + s, OutputMode::Text, true, tSimu);
                gloDHand->printFile("sol_gloCons." + s, OutputMode::Text, true, tSimu);
                GaussPointStateWrite(posDisMesh, paramStr->b_aux, gPts, "sol_gp." + s);
                if (myRank == 0)
                {
                    std::ofstream rinfo("restart_info." + s);
                    rinfo << std::setprecision(17) << tSimu << " " << deltat << " " << deltatMax << " "
                          << vstep << " " << v_r_step << endl;
                }
            }

        }
        else
        {
            if (myRank == 0) cout << "Newton-Raphson did not converge: dt reduced" << endl;
            deltat *= stepFactor;
            posDHand->nodeDOFs->setValue(posDHand->nodeDOFs0);
            gloDHand->nodeDOFs->setValue(gloDHand->nodeDOFs0);
            hiperProbl->UpdateGhosts();
        }

        if (myRank == 0)
        {
            const double v1 = hiperProbl->globalIntegral("volume");
            const double a1 = hiperProbl->globalIntegral("area");
            cout << std::scientific << std::setprecision(6);
            cout << "  VOLUME: " << v1 << "  volume fraction: " << v1/v0 << endl;
            cout << "  AREA: " << a1 << "  stretch: " << (a1 - a_res)/(area_n - a_res) << endl;
            cout << "  pressure: " << -gloDHand->nodeDOFs->getValue(0, 0)/deltat << endl;
            cout << "  strain energy: " << hiperProbl->globalIntegral("Energy")
                 << "  bending energy: " << hiperProbl->globalIntegral("E1") << endl;
            cout << "  drag dissipation: " << hiperProbl->globalIntegral("Dissipation")
                 << "  viscous dissipation: " << viscoProbl->globalIntegral("Diss1")
                 << "  curvature dissipation: " << viscoProbl->globalIntegral("Diss2") << endl;
        }
    }

    hiperlife::Finalize();
    return 0;
}
