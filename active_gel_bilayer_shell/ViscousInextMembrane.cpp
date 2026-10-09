// =============================================================================
// Active-gel bilayer shell model of an epithelial dome
// (continuum active gel shell, Supplementary Note 1 of Chahare, Ouzeri,
//  Wilson, Bal et al., "Multiscale wrinkling and folding dynamics in
//  epithelial shells"). See AuxViscousInextMembrane.h for the model.
//
// Protocol (same as the phenomenological shell_model):
//   1) inflation : enclosed volume ramped linearly 0 -> target over vnstep steps
//   2) hold      : gap steps (the cortices relax, turnover sets rho)
//   3) deflation : volume ramped down over v_rn_step steps (up to v_cont), with
//                  the time step limited to control_dt from control_steps
//                  steps before the deflation
//   4) relaxation at the final volume
//   From step vnstep+gap+control_fric on, the parts of the shell touching the
//   substrate are anchored by a spring (kspr) to their position at that step.
//
// Problems solved every converged step:
//   hiperProbl (posDHand + gloDHand): Newton-Raphson, element filling LS
//   viscoProbl (viscoDHand)         : relaxation of G at Gauss points (LS_ReacDif)
//   RhoProbl   (RhoDHand)           : density turnover, linear solve (LS_Rho)
//   ENProbl    (EDHand)             : output projection, only at print steps (LS_ED)
//   xyzProbl   (xyzDHand)           : stores the spring anchors once (LS_node_end)
// =============================================================================

#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
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
#include <hl_LocMongeParam.h>
#include "hl_Parser.h"
#include "hl_ParamStructure.h"
#include "hl_Core.h"

#include "AuxViscousInextMembrane.h"

int main (int argc, char** argv)
{
    using namespace std;
    using namespace hiperlife;
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
    // Parameters (see config_shell.cfg)
    // ---------------------------------------------------------------------
    // geometry and mesh
    double R{1.0}, height_in{0.05}, subdiv{1.0}, mesh_refine{0.0}, sp_gap{0.008};
    int choice{100};
    string prefixMesh;
    // cortex mechanics
    double young{25.0}, poisson{0.25}, fric2{0.015}, fact_elastic{0.01};
    double gamma_l{1.0}, gamma_sigma{2.0}, tens_factor{0.0}, f0_factor{1.0}, aap{1.0}, bbn{1.0};
    double k_p{0.8333}, k_d{0.8333};
    double fric{0.01}, Kconf{0.0}, kspr{0.0}, width{100.0};
    int control_fric{10000};
    double factor0{0.333};
    // protocol
    double v_target{0.73}, v_rn_step{1000.0}, v_cont{0.98}, control_dt{0.0001}, coeff{0.0};
    int vnstep{500}, gap{2000}, control_steps{400}, time_target{1};
    // time stepping and numerics
    double deltat1{0.0001}, deltatMax{}, totalTime{}, tSimu{}, stepFactor{}, forward_new{0.0};
    int totalSteps{}, nPrint{10}, ConsCheck{}, MAXITER_NR{};
    double gPts{}, SOLTOL_NR{}, RESTOL_NR{};
    {
        config.readInto(R, "R");
        config.readInto(height_in, "height_in");
        config.readInto(subdiv, "subdiv");
        config.readInto(mesh_refine, "mesh_refine");
        config.readInto(sp_gap, "sp_gap");
        config.readInto(choice, "choice");
        config.readInto(prefixMesh, "prefixMesh");

        config.readInto(young, "young");
        config.readInto(poisson, "poisson");
        config.readInto(fric2, "fric2");
        config.readInto(fact_elastic, "fact_elastic");
        config.readInto(gamma_l, "gamma_l");
        config.readInto(gamma_sigma, "gamma_sigma");
        config.readInto(tens_factor, "tens_factor");
        config.readInto(f0_factor, "f0_factor");
        config.readInto(aap, "aap");
        config.readInto(bbn, "bbn");
        config.readInto(k_p, "k_p");
        config.readInto(k_d, "k_d");
        config.readInto(fric, "fric");
        config.readInto(Kconf, "Kconf");
        config.readInto(kspr, "kspr");
        config.readInto(width, "width");
        config.readInto(control_fric, "control_fric");
        config.readInto(factor0, "factor");

        config.readInto(v_target, "v_target");
        config.readInto(vnstep, "vnstep");
        config.readInto(gap, "gap");
        config.readInto(v_rn_step, "v_rn_step");
        config.readInto(v_cont, "v_cont");
        config.readInto(coeff, "coeff");
        config.readInto(control_steps, "control_steps");
        config.readInto(control_dt, "control_dt");
        config.readInto(time_target, "time_target");

        config.readInto(deltat1, "deltat");
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
    v_target=v_target*(R);

    // active tensions: apical gamma+ and basal gamma- around gamma_l*gamma_sigma/2,
    // split by tens_factor; lateral weight f0 = 2 gamma_sigma f0_factor
    double gamma_del  = gamma_sigma*tens_factor;
    double gamma_plus = gamma_l*(gamma_sigma-gamma_del)/2;
    double gamma_minus= gamma_l*(gamma_sigma+gamma_del)/2;
    double f0         = 2*gamma_sigma*f0_factor;

    SmartPtr<ParamStructure> paramStr = CreateParamStructure<MembParams>();
    paramStr->setRealParameter(MembParams::young, young);
    paramStr->setRealParameter(MembParams::poisson, poisson);
    paramStr->setRealParameter(MembParams::fric2, fric2);
    paramStr->setRealParameter(MembParams::fact_elastic, fact_elastic);
    paramStr->setRealParameter(MembParams::gamma_plus, gamma_plus);
    paramStr->setRealParameter(MembParams::gamma_minus, gamma_minus);
    paramStr->setRealParameter(MembParams::gamma_l, gamma_l);
    paramStr->setRealParameter(MembParams::f0, f0);
    paramStr->setRealParameter(MembParams::k_p, k_p);
    paramStr->setRealParameter(MembParams::k_d, k_d);
    paramStr->setRealParameter(MembParams::height_in, height_in);
    paramStr->setRealParameter(MembParams::aap, aap);
    paramStr->setRealParameter(MembParams::bbn, bbn);
    paramStr->setRealParameter(MembParams::fric, fric);
    paramStr->setRealParameter(MembParams::Kconf, Kconf);
    paramStr->setRealParameter(MembParams::kspr, kspr);
    paramStr->setRealParameter(MembParams::width, width);
    paramStr->setRealParameter(MembParams::forward_new, forward_new);
    paramStr->setIntParameter(MembParams::time_target, time_target);
    paramStr->setIntParameter(MembParams::fric_start, vnstep + gap);
    paramStr->setIntParameter(MembParams::control_fric, control_fric);
    paramStr->setRealParameter(MembParams::deltat, deltat1);
    paramStr->setIntParameter(MembParams::timestep, 0);
    paramStr->setRealParameter(MembParams::factor, factor0);

    if (myRank == 0)
    {
        cout << endl << "Parameters:" << endl;
        cout << "  footprint radius R               : " << R << endl;
        cout << "  thickness h                      : " << height_in << "   R/h = " << R/height_in << endl;
        cout << "  cortex young, poisson            : " << young << " / " << poisson << "   (mu = lambda = " << 0.5*young/(1+poisson) << ")" << endl;
        cout << "  cortex viscosity fric2           : " << fric2 << "   (tau ~ 2 fric2 (1+nu)/E = " << 2*fric2*(1+poisson)/young << ")" << endl;
        cout << "  tensions apical / basal / lateral: " << gamma_plus << " / " << gamma_minus << " / " << gamma_l << " (f0 = " << f0 << ")" << endl;
        cout << "  turnover k_p, k_d                : " << k_p << " / " << k_d << "   (rho_ss = " << k_p/k_d << ", time " << 1/k_d << ")" << endl;
        cout << "  drag fric, confinement Kconf     : " << fric << " / " << Kconf << endl;
        cout << "  inflation / hold / deflation     : " << vnstep << " / " << gap << " / " << v_rn_step << " steps" << endl;
        cout << "  anchoring spring kspr " << kspr << " from step " << vnstep+gap+control_fric << " (width " << width << ")" << endl;
        cout << "  basis functions                  : " << (subdiv > 0.1 ? "subdivision" : "linear") << endl;
        cout << "  relaxation integrator            : " << (forward_new > 0.1 ? "explicit" : "implicit") << endl << endl;
    }

    // ---------------------------------------------------------------------
    // Mesh: disc of radius 0.25 in the mesh file, scaled to radius R
    // ---------------------------------------------------------------------
    SmartPtr<DistributedMesh> posDisMesh, gloDisMesh;
    {
        SmartPtr<MeshLoader> mesh = Create<MeshLoader>();
        if (subdiv > 0.1)
            mesh->setMesh(ElemType::Triang, BasisFuncType::SubdivSurfs, 2);
        else
            mesh->setMesh(ElemType::Triang, BasisFuncType::Lagrangian, 1);
        mesh->loadMesh(prefixMesh + ".vtk", MeshType::Sequential);
        mesh->transformFree([R](double x, double y)
        {
            x = R/0.25 * x;
            y = R/0.25 * y;
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
        posDisMesh->printFileLegacyVtk("pos_mesh_final_No_refine");

        gloDisMesh = Create<DistributedMesh>();
        gloDisMesh->setMeshRelation(MeshRelation::GlobConstr, posDisMesh);
        gloDisMesh->setBalanceMesh(false);
        gloDisMesh->Update();
    }

    // ---------------------------------------------------------------------
    // DOF handlers
    // ---------------------------------------------------------------------
    SmartPtr<DOFsHandler> posDHand, gloDHand, viscoDHand, EDHand, RhoDHand, xyzDHand;
    try
    {
        // element loop for LS_ReacDif; its AuxF holds the initial nodal normals
        viscoDHand = Create<DOFsHandler>(posDisMesh);
        viscoDHand->setNameTag("viscoDHand");
        viscoDHand->setDOFs({"dummy"});
        viscoDHand->setNodeAuxF({"e1x", "e1y", "e1z", "e2x", "e2y", "e2z", "nx", "ny", "nz"});
        viscoDHand->Update();

        // shell: position, director, thickness multiplier.
        // AuxF 7-9 (nxx,nyy,nzz) = reference nodal normal, READ by LS; the rest is output.
        posDHand = Create<DOFsHandler>(posDisMesh);
        posDHand->setNameTag("posDHand");
        posDHand->setDOFs({"X", "Y", "Z","H1", "H2", "H3","LAM"});
        posDHand->setNodeAuxF({"Ux", "Uy", "Uz","nx","ny","nz","height","nxx","nyy","nzz",
                               "rho_api","rho_bas","rho_lat1","rho_lat2","rho_lat3",
                               "EA","EB","EL","DA","DB","DL","PA","PB","PL","jacC"});
        posDHand->Update();

        EDHand = Create<DOFsHandler>(posDisMesh);
        EDHand->setNameTag("EDHand");
        EDHand->setDOFs({"nx","ny","nz","height","EA", "EB","EL","DA", "DB","DL","PA","PB","PL","jacC"});
        EDHand->Update();

        RhoDHand = Create<DOFsHandler>(posDisMesh);
        RhoDHand->setNameTag("RhoDHand");
        RhoDHand->setDOFs({"rho_api","rho_bas","rho_lat1","rho_lat2","rho_lat3"});
        RhoDHand->Update();

        // element loop for LS_node_end
        xyzDHand = Create<DOFsHandler>(posDisMesh);
        xyzDHand->setNameTag("xyzDHand");
        xyzDHand->setDOFs({"dummy"});
        xyzDHand->Update();

        gloDHand = Create<DOFsHandler>(gloDisMesh);
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

    LocMongeParam::computeLocal3DBasis(posDisMesh, viscoDHand->nodeAuxF);
    viscoDHand->UpdateGhosts();

    // initial condition: flat disc, director h = height_in * n, LAM = 0.001
    for (int i = 0; i < posDisMesh->loc_nPts(); i++)
    {
        double x = posDisMesh->nodeCoord(i, 0, IndexType::Local);
        double y = posDisMesh->nodeCoord(i, 1, IndexType::Local);
        double z = posDisMesh->nodeCoord(i, 2, IndexType::Local);
        posDHand->nodeDOFs->setValue("X", i, IndexType::Local, x);
        posDHand->nodeDOFs->setValue("Y", i, IndexType::Local, y);
        posDHand->nodeDOFs->setValue("Z", i, IndexType::Local, z);
        posDHand->nodeDOFs->setValue("LAM", i, IndexType::Local, 0.001);

        double nx1 = viscoDHand->nodeAuxF->getValue("nx",i,IndexType::Local);
        double ny1 = viscoDHand->nodeAuxF->getValue("ny",i,IndexType::Local);
        double nz1 = viscoDHand->nodeAuxF->getValue("nz",i,IndexType::Local);
        if(nz1<0.0)
            nz1=1.0;
        posDHand->nodeDOFs->setValue("H1", i, IndexType::Local, height_in*nx1);
        posDHand->nodeDOFs->setValue("H2", i, IndexType::Local, height_in*ny1);
        posDHand->nodeDOFs->setValue("H3", i, IndexType::Local, height_in*nz1);
        posDHand->nodeAuxF->setValue("nxx", i, IndexType::Local, nx1);
        posDHand->nodeAuxF->setValue("nyy", i, IndexType::Local, ny1);
        posDHand->nodeAuxF->setValue("nzz", i, IndexType::Local, nz1);

        // boundary conditions: edge (creases) clamped (all DOFs); outer ring r > R - sp_gap:
        // choice = 10 -> all DOFs fixed, otherwise only the position fixed
        bool edge  = posDisMesh->nodeCrease(i, IndexType::Local) > 0;
        bool outer = sqrt(x*x+y*y) > R-sp_gap;
        if (edge or (outer and choice == 10))
            for (int k = 0; k < 7; k++)
                posDHand->setConstraint(k, i, IndexType::Local, 0.0);
        else if (outer)
            for (int k = 0; k < 3; k++)
                posDHand->setConstraint(k, i, IndexType::Local, 0.0);
    }
    posDHand->nodeDOFs0->setValue(posDHand->nodeDOFs);
    posDHand->UpdateGhosts();

    // densities at steady state k_p/k_d
    for (int i = 0; i < posDisMesh->loc_nPts(); i++)
        for (int k = 0; k < 5; k++)
            RhoDHand->nodeDOFs->setValue(k, i, IndexType::Local, k_p/k_d);
    RhoDHand->nodeDOFs0->setValue(RhoDHand->nodeDOFs);
    RhoDHand->UpdateGhosts();

    gloDHand->setInitialCondition(0, 0.0);
    for (int k = 1; k <= 6; k++)
        gloDHand->setConstraint(k, 0.0);
    gloDHand->UpdateGhosts();

    // Gauss-point storage
    paramStr->a_aux.resize(posDisMesh->loc_nElem()*gPts*88);
    paramStr->b_aux.resize(posDisMesh->loc_nElem()*gPts*18);
    paramStr->c_aux.resize(posDisMesh->loc_nElem()*gPts*9);

    // ---------------------------------------------------------------------
    // Problems
    // ---------------------------------------------------------------------
    auto makeProblem = [&](SmartPtr<HiPerProblem> pb, SmartPtr<DOFsHandler> dh, string tag,
                           void (*fill)(hiperlife::FillStructure&), std::vector<string> integrals)
    {
        pb->setParameterStructure(paramStr);
        pb->setDOFsHandlers({dh});
        pb->setIntegration("Integ", {tag});
        pb->setCubatureGauss("Integ", gPts);
        pb->setElementFillings("Integ", fill);
        if (!integrals.empty())
            pb->setGlobalIntegrals(integrals);
        pb->Update();
    };
    SmartPtr<HiPerProblem> hiperProbl = Create<HiPerProblem>();
    SmartPtr<HiPerProblem> viscoProbl = Create<HiPerProblem>();
    SmartPtr<HiPerProblem> ENProbl    = Create<HiPerProblem>();
    SmartPtr<HiPerProblem> RhoProbl   = Create<HiPerProblem>();
    SmartPtr<HiPerProblem> xyzProbl   = Create<HiPerProblem>();
    try
    {
        hiperProbl->setParameterStructure(paramStr);
        hiperProbl->setConsistencyCheckDelta(1.E-6);
        hiperProbl->setConsistencyCheckTolerance(1.E-4);
        hiperProbl->setDOFsHandlers({posDHand,gloDHand});
        hiperProbl->setIntegration("Integ", {"posDHand","gloDHand"});
        hiperProbl->setCubatureGauss("Integ", gPts);
        if (ConsCheck == 1)
            hiperProbl->setElementFillings("Integ", ConsistencyCheck<LS>);
        else
            hiperProbl->setElementFillings("Integ", LS);
        // (E_tension1..3 hold the spring anchor x,y,z * area; kept for the columns of globalIntegrals.dat)
        hiperProbl->setGlobalIntegrals({"Energy","E1","E2","E3","Dissipation","Trd","volume","area","area_n","area_gn","area_gp","Lat_area","Epress","E_tension1","E_tension2","E_tension3","Ela_small","P_tension1","P_tension2","P_tension3"});
        hiperProbl->Update();

        makeProblem(viscoProbl, viscoDHand, "viscoDHand", LS_ReacDif, {"Diss1","Diss2","Diss3"});
        makeProblem(ENProbl,    EDHand,     "EDHand",     LS_ED,      {});
        makeProblem(RhoProbl,   RhoDHand,   "RhoDHand",   LS_Rho,     {});
        makeProblem(xyzProbl,   xyzDHand,   "xyzDHand",   LS_node_end,{});
    }
    catch (runtime_error& err)
    {
        cout << myRank << ": HiPerProblem could not be created " << err.what() << endl;
        hiperlife::Finalize();
        return 1;
    }

    auto makeMUMPS = [](SmartPtr<HiPerProblem> pb)
    {
        SmartPtr<MUMPSDirectLinearSolver> s = Create<MUMPSDirectLinearSolver>();
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
    SmartPtr<MUMPSDirectLinearSolver> linSolver    = makeMUMPS(hiperProbl);
    SmartPtr<MUMPSDirectLinearSolver> ENDirSolver  = makeMUMPS(ENProbl);
    SmartPtr<MUMPSDirectLinearSolver> RhoDirSolver = makeMUMPS(RhoProbl);

    SmartPtr<NewtonRaphsonNonlinearSolver> nonlinSolver = Create<NewtonRaphsonNonlinearSolver>();
    nonlinSolver->setLinearSolver(linSolver);
    nonlinSolver->setMaxNumIterations(MAXITER_NR);
    nonlinSolver->setResTolerance(RESTOL_NR);
    nonlinSolver->setSolTolerance(SOLTOL_NR);
    nonlinSolver->setLineSearch(false);            // always plain Newton
    nonlinSolver->setConvRelTolerance(false);
    nonlinSolver->setPrintIntermInfo(true);
    nonlinSolver->setPrintSummary(false);
    nonlinSolver->setResMaximum(1E4);
    nonlinSolver->setSolMaximum(1E4);
    nonlinSolver->setExitRelMaximum(1E4);
    nonlinSolver->Update();

    // initialise the Gauss-point storage (same order as the original code)
    hiperProbl->FillGlobalIntegrals();
    viscoProbl->FillGlobalIntegrals();
    RhoProbl->FillGlobalIntegrals();
    hiperProbl->FillGlobalIntegrals();
    xyzProbl->FillGlobalIntegrals();

    // circular footprint: base area 3.141 R^2
    double base_area = 3.141*R*R;
    double a0 = hiperProbl->globalIntegral("area_n");
    v_target = v_target/a0*base_area;
    double a_res = a0-base_area;
    double v0 = v_target*a0;
    if (myRank == 0)
    {
        cout << std::scientific << std::setprecision(6);
        cout << "  Initial volume: " << hiperProbl->globalIntegral("volume") << endl;
        cout << "  Target volume : " << abs(v_target)*a0 << endl;
        cout << "  Area          : " << a0 << "   apical " << hiperProbl->globalIntegral("area_gp") << "   basal " << hiperProbl->globalIntegral("area_gn") << endl;
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

    double& deltat   = paramStr->getRealParameter(MembParams::deltat);
    int&    timestep = paramStr->getIntParameter(MembParams::timestep);
    double& factor   = paramStr->getRealParameter(MembParams::factor);

    posDHand->printFileLegacyVtk("sol_dis." + to_string(timestep), true);

    double v_r_step = 0.0;
    double press0 = 0.0;
    const double deltatMax1 = deltatMax;

    // ---------------------------------------------------------------------
    // Time loop
    // ---------------------------------------------------------------------
    while ((tSimu < totalTime) and (timestep < totalSteps))
    {
        if (myRank == 0)
            cout << "TS: " << timestep + 1 << " Time " << tSimu << " of " << totalTime << " with deltat=" << deltat << ": del_max: " << deltatMax << endl;

        // (1) inflation
        if (timestep < vnstep)
        {
            factor = (timestep+1.0)/vnstep*v_target;
            if (myRank == 0)
                cout << "inflation: fraction " << (timestep+1.0)/vnstep << endl;
        }
        // small time step from control_steps steps before the deflation
        if (timestep > vnstep+gap-control_steps)
        {
            deltat *= stepFactor;
            deltatMax = control_dt;
        }
        // after the deflation: back to the normal time step
        if (timestep > vnstep+gap+v_cont*v_rn_step)
        {
            deltat /= stepFactor;
            deltatMax = deltatMax1;
        }
        // (3) deflation (coeff = 0: linear ramp)
        if (timestep > vnstep+gap)
        {
            if (v_r_step < v_cont*v_rn_step)
            {
                factor = (1-(v_r_step+1.0)/v_rn_step)*v_target;
                factor = factor*exp(coeff*deltat*v_r_step);
                v_r_step = v_r_step+1;
                if (myRank == 0)
                    cout << "deflation step " << v_r_step << ": volume fraction " << factor/v_target << endl;
            }
        }

        if (deltat < 0.000001)
        {
            if (myRank == 0)
                cout << "stopping: time step too small" << endl;
            break;
        }

        // store the spring anchors (positions at this step)
        if (timestep == vnstep + gap + control_fric)
        {
            if (myRank == 0)
                cout << "anchoring spring switched on" << endl;
            hiperProbl->FillGlobalIntegrals();
            xyzProbl->FillGlobalIntegrals();
            hiperProbl->FillGlobalIntegrals();
        }

        // Newton-Raphson for the shell
        posDHand->nodeDOFs->setValue(posDHand->nodeDOFs0);
        gloDHand->nodeDOFs->setValue(gloDHand->nodeDOFs0);
        hiperProbl->UpdateGhosts();
        RhoDHand->nodeDOFs0->setValue(RhoDHand->nodeDOFs);
        RhoDHand->UpdateGhosts();

        bool converged = nonlinSolver->solve();

        if (converged)
        {
            // Gauss-point data at the converged state, relax G, update densities
            hiperProbl->FillGlobalIntegrals();
            viscoProbl->FillGlobalIntegrals();
            RhoDirSolver->solve();
            RhoProbl->UpdateSolution();
            RhoProbl->FillGlobalIntegrals();        // densities -> c_aux

            tSimu += deltat;
            timestep += 1;
            press0 = gloDHand->nodeDOFs->getValue(0,0)/deltat;

            if (timestep % nPrint == 0)
            {
                // nodal output fields
                ENDirSolver->solve();
                ENProbl->UpdateSolution();
                for (int i = 0; i < posDHand->mesh->loc_nPts(); i++)
                {
                    // Ux,Uy,Uz: displacement of the last step
                    for (int d = 0; d < 3; d++)
                        posDHand->nodeAuxF->setValue(d, i, IndexType::Local,
                            posDHand->nodeDOFs->getValue(d,i,IndexType::Local) - posDHand->nodeDOFs0->getValue(d,i,IndexType::Local));
                    // nx,ny,nz,height
                    for (int k = 0; k < 4; k++)
                        posDHand->nodeAuxF->setValue(3+k, i, IndexType::Local, EDHand->nodeDOFs->getValue(k,i,IndexType::Local));
                    // rho_api, rho_bas, rho_lat1..3 (NB: rho_lat3 shows rho_lat2, as in the original code)
                    posDHand->nodeAuxF->setValue("rho_api",  i, IndexType::Local, RhoDHand->nodeDOFs->getValue(0,i,IndexType::Local));
                    posDHand->nodeAuxF->setValue("rho_bas",  i, IndexType::Local, RhoDHand->nodeDOFs->getValue(1,i,IndexType::Local));
                    posDHand->nodeAuxF->setValue("rho_lat1", i, IndexType::Local, RhoDHand->nodeDOFs->getValue(2,i,IndexType::Local));
                    posDHand->nodeAuxF->setValue("rho_lat2", i, IndexType::Local, RhoDHand->nodeDOFs->getValue(3,i,IndexType::Local));
                    posDHand->nodeAuxF->setValue("rho_lat3", i, IndexType::Local, RhoDHand->nodeDOFs->getValue(3,i,IndexType::Local));
                    // EA,EB,EL,DA,DB,DL,PA,PB,PL,jacC
                    for (int k = 4; k < 14; k++)
                        posDHand->nodeAuxF->setValue(11+k, i, IndexType::Local, EDHand->nodeDOFs->getValue(k,i,IndexType::Local));
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
            cout << "  area " << hiperProbl->globalIntegral("area") << "   apical " << hiperProbl->globalIntegral("area_gp")
                 << "   basal " << hiperProbl->globalIntegral("area_gn") << "   lateral " << hiperProbl->globalIntegral("Lat_area")
                 << "   stretch " << (hiperProbl->globalIntegral("area")-a_res)/(hiperProbl->globalIntegral("area_n")-a_res) << endl;
            cout << "  pressure " << -gloDHand->nodeDOFs->getValue(0,0)/deltat
                 << "   excess area " << abs(hiperProbl->globalIntegral("area")-hiperProbl->globalIntegral("area_n"))/base_area << endl;
            cout << "  elastic energy " << hiperProbl->globalIntegral("Energy") << "   (lateral " << hiperProbl->globalIntegral("E3") << ")" << endl;
            cout << "  tension power apical / basal / lateral " << hiperProbl->globalIntegral("P_tension1") << " / "
                 << hiperProbl->globalIntegral("P_tension2") << " / " << hiperProbl->globalIntegral("P_tension3") << endl;
            cout << "  viscous dissipation apical / basal / lateral " << viscoProbl->globalIntegral("Diss1") << " / "
                 << viscoProbl->globalIntegral("Diss2") << " / " << viscoProbl->globalIntegral("Diss3")
                 << "   drag " << hiperProbl->globalIntegral("Dissipation") << endl;
        }

        // output
        if (timestep % nPrint == 0)
        {
            posDHand->printFileLegacyVtk("sol_dis." + to_string(timestep), true);
            if (myRank == 0)
            {
                gIntegFile << timestep << " " << tSimu << " " << deltat;
                for (auto g: hiperProbl->globIntegrals())
                    gIntegFile << " " << g;
                gIntegFile << " " << press0 << " " << 1.0 << " " << 0.0 << " " << 0.0;
                for (auto g: viscoProbl->globIntegrals())
                    gIntegFile << " " << g;
                gIntegFile << endl;
            }
        }
    }

    hiperlife::Finalize();
    return 0;
}
