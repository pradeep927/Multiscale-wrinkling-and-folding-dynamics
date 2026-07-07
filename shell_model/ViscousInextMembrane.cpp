
// C++ headers
#include <iostream>
#include <mpi.h>

// Trilinos headers
#include <Teuchos_RCP.hpp>
#include <math.h>

// hiperlife headers
#include "hl_TypeDefs.h"
#include "hl_Geometry.h"
#include "hl_StructMeshGenerator.h"
#include "hl_DistributedMesh.h"
#include "hl_FillStructure.h"
#include "hl_DOFsHandler.h"
#include "hl_HiPerProblem.h"
#include "hl_ConsistencyCheck.h"
#include <hl_LoadBalance.h>
//#include <hl_LinearSolver_Direct_Amesos2.h>
#include <hl_NonlinearSolver_NewtonRaphson.h>
#include <hl_MeshLoader.h>
#include <hl_ConfigFile.h>
#include <hl_LinearSolver_Direct_MUMPS.h>
#include <hl_LocMongeParam.h>

// Header to auxiliary functions
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


    // **************************************************************//
    // *****                 INITIALIZATION                     *****//
    // **************************************************************//
    hiperlife::Init(argc, argv);
    const int myRank = hiperlife::MyRank();

    // Model parameters and solver options are read from a configuration file
    if (argc < 2)
    {
        cout << "Configuration file not provided!" << endl;
        MPI_Finalize();
        return 1;
    }
    // const char *config_filename = argv[1];
    ConfigFile config(argv[1]);


    //Cortex model parameters
    double  k_p{0.8333},k_d{0.8333}, factor{0.333},v_inc{0.0}, vol_start{0.0};
    {
        config.readInto(k_p, "k_p");
        config.readInto(k_d, "k_d");
        config.readInto(factor, "factor");
        config.readInto(v_inc, "v_inc");
        config.readInto(vol_start, "vol_start");


    }
    int aSteps{100},control_fric{10000}, time_target{1}, case_sphere{0};
    double  R{1.0},thick{},Ainc{},aap{1.0},kappa{0.0},apical{0.0},subdiv{1.0},basal{0.0},bbn{1.0},width{80.0},bound{0.0},bc_const{1.0},ang_old{0.0},mesh_refine{0.0},mom_steps{100.0},sp_gap{0.058},a11{0.51},a12{0.49},a33{0.8},kap1{1.0},fact_elastic{0.01};
    {
        config.readInto(R, "R");
        config.readInto(thick, "thick");
        config.readInto(aap, "aap");
        config.readInto(kappa, "kappa");
        config.readInto(aSteps, "aSteps");
        config.readInto(bbn,"bbn");
        config.readInto(width,"width");
        config.readInto(bound,"bound");
        config.readInto(bc_const,"bc_const");
        config.readInto(control_fric,"control_fric");
        config.readInto(ang_old,"ang_old");
        config.readInto(mesh_refine,"mesh_refine");
        config.readInto(mom_steps,"mom_steps");
        config.readInto(sp_gap,"sp_gap");
        config.readInto(fact_elastic,"fact_elastic");
        config.readInto(subdiv,"subdiv");
        config.readInto(case_sphere,"case_sphere");

        config.readInto(time_target,"time_target");

        config.readInto(a11,"a11");
        config.readInto(a12,"a12");
        config.readInto(a33,"a33");
        config.readInto(kap1, "kap1");

        config.readInto(apical,"apical");
        config.readInto(basal,"basal");

    }
    int spring{0};
    double mu{7.8372}, lambda{16.6482}, fric{0.01},Xmax{1.0},Ymax{1.0},forward_new{0.0},kspr{0.0},fric3{0.02},stretching{0.0}, young{25.0}, poisson{0.25},my_work{1.0},gamma_minus{1.0},gamma_sigma{2.0},gamma_del{0.0},floating{0.0},gamma_plus{1.0},force{0.0},Kconf{0.0},gamma_l{1.0},f0{0.9},g_ratio{0.8},v_target{0.667},vstep{0.0}, v_r_step{0.0},coeff{0.0},ang{0.0},tens_factor{0.0},out_choice{0.0};
    {
        config.readInto(fric, "fric");
        config.readInto(young, "young");
        config.readInto(lambda, "lambda");
        config.readInto(mu, "mu");


        config.readInto(poisson, "poisson");
        config.readInto(force, "force");
        config.readInto(Kconf, "Kconf");
        config.readInto(ang, "ang");
        config.readInto(gamma_del, "gamma_del");

        config.readInto(forward_new, "forward_new");

        config.readInto(tens_factor,"tens_factor");
        config.readInto(f0, "f0");
        config.readInto(g_ratio, "g_ratio");
        config.readInto(v_target, "v_target");
        config.readInto(vstep, "vstep");
        config.readInto(v_r_step, "v_r_step");
        config.readInto(coeff, "coeff");
        config.readInto(gamma_minus, "gamma_minus");
        config.readInto(gamma_plus, "gamma_plus");
        config.readInto(floating, "floating");
        config.readInto(Xmax, "Xmax");
        config.readInto(my_work, "my_work");
        config.readInto(Ymax, "Ymax");
        config.readInto(gamma_sigma, "gamma_sigma");
        config.readInto(fric3, "fric3");
        config.readInto(spring, "spring");


        config.readInto(gamma_l, "gamma_l");
        config.readInto(kspr, "kspr");


        config.readInto(stretching, "stretching");


    }

    v_target=v_target*(R);
    //Simulation time settings
    int nSave{},hl_print{1000},stretch_start{26000};
    double deltat{}, deltatMax{}, gamma{5.0},square{0.0},totalTime{}, tSimu{}, stepFactor{},stretch_xx{0.001},pn{1.0},f0_factor{1.0},fric_fact{1},height_in{0.05},fric2{0.015}, v_rn_step{1000.0}, vol_inc{1.0},control_dt{0.0001},v_cont{0.98},out_fact{1.2},thin_shell{0.0}, case_tension{0.0};
    int restart{}, totalSteps{}, nPrint{10},ConsCheck{},nstep{50} ,vnstep{500},tens_start{2},choice{100}, control_steps{400}, gap{2000}, test_jac{0},stretch_step{100};
    {
        config.readInto(deltat, "deltat");
        config.readInto(deltatMax, "deltatMax");
        config.readInto(totalTime, "totalTime");
        config.readInto(totalSteps, "totalSteps");
        config.readInto(tSimu, "tSimu");
        config.readInto(stepFactor, "stepFactor");
        config.readInto(restart, "restart");
        config.readInto(nPrint, "nPrint");
        config.readInto(ConsCheck, "ConsCheck");
        config.readInto(nstep, "nstep");
        config.readInto(vnstep, "vnstep");
        config.readInto(tens_start, "tens_start");
        config.readInto(height_in, "height_in");
        config.readInto(v_rn_step, "v_rn_step");
        config.readInto(v_cont, "v_cont");

        config.readInto(f0_factor, "f0_factor");
        config.readInto(case_tension, "case_tension");
        config.readInto(stretch_step, "stretch_step");
        config.readInto(stretch_xx, "stretch_xx");
        config.readInto(pn, "pn");
        config.readInto(stretch_start, "stretch_start");

        config.readInto(square, "square");

        config.readInto(out_fact, "out_fact");
        config.readInto(hl_print, "hl_print");

        config.readInto(thin_shell, "thin_shell");



        config.readInto(choice, "choice");
        config.readInto(fric2, "fric2");

        config.readInto(gap, "gap");
        config.readInto(test_jac, "test_jac");

        config.readInto(control_steps, "control_steps");
        config.readInto(control_dt, "control_dt");

    }






        gamma_del=gamma_sigma*tens_factor;
        gamma_plus= gamma_l*(gamma_sigma-gamma_del)/2;
        gamma_minus=gamma_l*(gamma_sigma+gamma_del)/2;
       gamma=0.5*(gamma_plus+gamma_minus);


        f0=2*gamma_sigma*f0_factor;


    //Numerical parameters
    double  gPts{};
    int MAXITER_AZ{}, MAXITER_NR{};
    double SOLTOL_NR{}, RESTOL_NR{}, SOLTOL_AZ{};
    string prefixMesh;
    {
        config.readInto(gPts, "gPts");
        config.readInto(MAXITER_NR, "MAXITER_NR");
        config.readInto(SOLTOL_NR, "SOLTOL_NR");
        config.readInto(RESTOL_NR, "RESTOL_NR");
        config.readInto(MAXITER_AZ, "MAXITER_AZ");
        config.readInto(SOLTOL_AZ, "SOLTOL_AZ");
        config.readInto(prefixMesh, "prefixMesh");
    }



        double    ppp=0.15;
        double qqq=1.7;
        double rrr=0.3;
        double r0_tri=0.18;//0.2, for trg_r2: q=2.1,r0=0.18,r=0.45, r_out=0.72; trg_r1: r0=0.2,q=1.7,r=0.5,r_out=0.8; trg_r11: r0=0.18,q=1.7,r3=0.45,r_out=0.72; tr_r1: r0=0.2,q=1.7,r=0.35, r_out=0.5;
        double r_tri=0.45;//0.35
        double q_tri=1.7;
 SmartPtr<ParamStructure> paramStr = CreateParamStructure<MembParams>();

paramStr->setRealParameter(MembParams::deltat, deltat);
paramStr->setRealParameter(MembParams::fric, fric);
paramStr->setRealParameter(MembParams::young, young);
paramStr->setRealParameter(MembParams::poisson, poisson);
paramStr->setRealParameter(MembParams::thick, thick);
paramStr->setRealParameter(MembParams::force, force);

paramStr->setIntParameter(MembParams::tens_start, tens_start);

paramStr->setRealParameter(MembParams::fric_fact, fric_fact);
paramStr->setRealParameter(MembParams::fric2, fric2);
paramStr->setRealParameter(MembParams::pn, pn);
paramStr->setRealParameter(MembParams::fric3, fric3);

paramStr->setRealParameter(MembParams::aap, aap);
paramStr->setRealParameter(MembParams::R, R);
paramStr->setRealParameter(MembParams::kappa, kappa);
paramStr->setRealParameter(MembParams::bbn, bbn);
paramStr->setRealParameter(MembParams::width, width);

paramStr->setIntParameter(MembParams::control_fric, control_fric);

paramStr->setRealParameter(MembParams::apical, apical);
paramStr->setRealParameter(MembParams::basal, basal);

paramStr->setIntParameter(MembParams::time_target, time_target);
paramStr->setRealParameter(MembParams::spring, spring);

paramStr->setRealParameter(MembParams::tens_factor, tens_factor);
paramStr->setRealParameter(MembParams::fact_elastic, fact_elastic);

paramStr->setIntParameter(MembParams::case_sphere, case_sphere);

paramStr->setRealParameter(MembParams::height_in, height_in);
paramStr->setRealParameter(MembParams::Kconf, Kconf);

paramStr->setRealParameter(MembParams::f0, f0);
paramStr->setRealParameter(MembParams::g_ratio, g_ratio);
paramStr->setRealParameter(MembParams::bc_const, bc_const);

paramStr->setRealParameter(MembParams::mu, mu);
paramStr->setRealParameter(MembParams::lambda, lambda);

paramStr->setRealParameter(MembParams::kspr, kspr);
paramStr->setRealParameter(MembParams::sp_gap, sp_gap);

paramStr->setRealParameter(MembParams::a11, a11);
paramStr->setRealParameter(MembParams::a12, a12);
paramStr->setRealParameter(MembParams::a33, a33);
paramStr->setRealParameter(MembParams::kap1, kap1);
    
paramStr->setRealParameter(MembParams::gamma_minus, gamma_minus);
paramStr->setRealParameter(MembParams::gamma_plus, gamma_plus);
paramStr->setRealParameter(MembParams::gamma_l, gamma_l);

paramStr->setRealParameter(MembParams::forward_new, forward_new);

paramStr->setIntParameter(MembParams::choice, choice);
 paramStr->setIntParameter(MembParams::vnstep, vnstep);
    paramStr->setIntParameter(MembParams::gap, gap);

    paramStr->setRealParameter(MembParams::gamma, gamma);

paramStr->setRealParameter(MembParams::Xmax, Xmax);

paramStr->setRealParameter(MembParams::k_p, k_p);
paramStr->setRealParameter(MembParams::k_d, k_d);
        //Output
        if (myRank == 0)
        {
            cout << endl << "Parameters: " << endl;
            cout << "Radius:   " << R << endl;
            cout << "thickness " <<height_in << endl;


            cout<<"gamma: "<<paramStr->getRealParameter(MembParams::gamma_minus)<<"g+: "<<paramStr->getRealParameter(MembParams::gamma_plus)<<endl;
            cout << endl;
            cout << "elastic parameters: lambda  " << lambda << "elastic parameters: mu  " << mu <<"gamma"<<gamma<<endl;

            cout << "friction:   " << fric << endl;
            cout << "friction 3: curvature viscoelasticity:   " << fric3  << "time: curvature viscoelasticity:   " << fric3/10 << endl;
            cout << "Total  volume reduction  step " << v_rn_step<< endl;
            cout << "tension: " << gamma<< endl;

            cout << "young:    " << young << endl;
            cout << "poisson: " << poisson<< endl;
            cout << "Applied force " << force<< endl;
            cout << "Total  force step " << nstep<< endl;
            cout << "Total  volume step " << vnstep<< endl;
            //cout << "Lateral parameter:f0 " << f0<< endl;
           cout << "spring radius gap " << sp_gap<< endl;
            cout << "spring constant " << kspr<< endl;
            cout << "a11 " << a11<<" a12: " << a12<< " a33: magnitude: " << a33<< " strength " <<width  <<endl;

            cout << endl;
            cout << "visco_dissipation constant:   " << fric2 << endl;
            cout << "basis function: 1 for loop, 0 for linear  " << subdiv << endl;
            cout << "mesh refinement 1:yes " <<mesh_refine<< endl;

            cout << "Consistency check:  " << ConsCheck << endl;
            cout << "Friction control step:  " << vnstep+gap+control_fric<<" width: "<<width<< endl;

            cout << "Confinement potential:  " << Kconf<< endl;


            cout << "Tension factor:  " << tens_factor<< "  :tension plus :"  <<   gamma_plus<<  " : tension minus: " <<   gamma_minus<<endl;
            cout << "Surface tension lateral:  " << gamma_l<< " :gamma_sigma:  " << gamma_sigma<<" :Surface tension delta:  " << gamma_del<<endl;

            cout << "Radius to thickness ratio:  " << R/height_in<< endl;

            cout << "tens_start:  " << tens_start<< endl;
            cout << "Elastic factor:  " << fact_elastic<< endl;
            cout << "Newton Method: 1.0 forward, 0 backward: " << forward_new<< endl;

            cout << "aap=1,bbn=1 normal mid:aap=2,bbn=0 up :aap=0,bbn=2 down::: " << aap << ": "<<bbn<< endl;
            cout << "xmax input:  " << Xmax<< "Ymax input:  " << Ymax<<endl;

            cout << "visco_elastic time:  " <<2*fric2*(1+ poisson)/young << endl;
            cout << endl;
        }


    //bool  ConsCheck = false;

    // Time related parameters
    string sol_prefixMesh = "sol";
    double vol_step=1;
    RCP<DistributedMesh> posDisMesh,gloDisMesh,tensDisMesh;
    RCP<DOFsHandler> posDHand, gloDHand,viscoDHand,EDHand ;

    if (restart == 0)
    {

        // **************************************************************//
        // *****                   MESH CREATION                    *****//
        // **************************************************************//

        // Load mesh with MeshLoader
        RCP <MeshLoader> mesh = rcp(new MeshLoader);
        mesh->setMesh(ElemType::Triang, BasisFuncType::SubdivSurfs, 2);
        mesh->loadMesh(prefixMesh + ".vtk", MeshType::Sequential);



        mesh->transformFree([R](double x, double y)
          {
              x = R/0.34 * x;//major_axis 0.25
              y = R/0.34 * y;//minor_axis 0.25
              return std::make_tuple(x, y);
          });



if(mesh_refine<1.0)
{
    // Distribute mesh for position
    posDisMesh = rcp(new DistributedMesh);
    posDisMesh->setMesh(mesh);
    posDisMesh->setBalanceMesh(true);
    posDisMesh->Update();
    posDisMesh->printFileLegacyVtk("pos_mesh_final_No_refine");

}
else
{

    // Distribute mesh for tension
    tensDisMesh = rcp(new DistributedMesh);
    tensDisMesh->setMesh(mesh);
    tensDisMesh->setBalanceMesh(true);
    tensDisMesh->Update();
    tensDisMesh->printFileLegacyVtk("pos_mesh_before_refine");


    // Distribute mesh for positions
    posDisMesh = rcp(new DistributedMesh);
    posDisMesh->setHRefinement(1);
    posDisMesh->setMeshRelation(MeshRelation::hRefin, tensDisMesh);
    posDisMesh->setBalanceMesh(true);
    posDisMesh->Update();
    posDisMesh->printFileLegacyVtk("pos_mesh_after_refined");

}

        // Distribute mesh global constraints
        gloDisMesh = rcp(new DistributedMesh);
        gloDisMesh->setMeshRelation(MeshRelation::GlobConstr, posDisMesh);
        gloDisMesh->setBalanceMesh(false);
        gloDisMesh->Update();

        gloDisMesh->printFileLegacyVtk("global_after");

        // **************************************************************//
        // *****               DOFHANDLER CREATION                  *****//
        // **************************************************************//

// CREATE VISCO DOF HAND
        viscoDHand = rcp(new DOFsHandler(posDisMesh));
        try
        {
            viscoDHand->setNameTag("viscoDHand");
            viscoDHand->setDOFs({"G11P", "G12P", "G22P", "G11N", "G12N", "G22N"});
            viscoDHand->setNodeAuxF({"e1x", "e1y", "e1z", "e2x", "e2y", "e2z", "nx", "ny", "nz"});
            viscoDHand->Update();
        }
        catch (runtime_error &err)
        {
            cout << myRank << ": VISCO DOFHandler could not be created " << err.what() << endl;
            MPI_Finalize();
            return 1;
        }

        for (int i = 0; i < posDisMesh->loc_nPts(); i++)
        {
            double x = posDisMesh->nodeCoord(i, 0, IndexType::Local);
            double y = posDisMesh->nodeCoord(i, 1, IndexType::Local);
            double z = posDisMesh->nodeCoord(i, 2, IndexType::Local);


            viscoDHand->nodeDOFs->setValue("G11P", i, IndexType::Local, 1.0);
            viscoDHand->nodeDOFs->setValue("G12P", i, IndexType::Local,0.0);
            viscoDHand->nodeDOFs->setValue("G22P", i, IndexType::Local, 1.0);
            viscoDHand->nodeDOFs->setValue("G11N", i, IndexType::Local,1.0);
            viscoDHand->nodeDOFs->setValue("G12N", i, IndexType::Local, 0.0);
            viscoDHand->nodeDOFs->setValue("G22N", i, IndexType::Local,1.0);

        }
        // InitializeNodalValues(posDHand, 3);
        viscoDHand->nodeDOFs0->setValue(viscoDHand->nodeDOFs);


        LocMongeParam::computeLocal3DBasis(posDisMesh, viscoDHand->nodeAuxF);// assign auxiliary
        viscoDHand->UpdateGhosts();
//POSITION DOF
        posDHand = rcp(new DOFsHandler(posDisMesh));
        try
        {
            posDHand->setNameTag("posDHand");
            posDHand->setDOFs({"X", "Y", "Z"});
            posDHand->setNodeAuxF({"Ux", "Uy", "Uz","nx","ny","nz","height","nxx","nyy","nzz", "EA","EB","DA","DB","PA","trace_sigma","dev_sigma"});
            posDHand->Update();
        }
        catch (runtime_error &err)
        {
            cout << myRank << ": DOFHandler could not be created " << err.what() << endl;
            MPI_Finalize();
            return 1;
        }
        //set initial condition;

        for (int i = 0; i < posDisMesh->loc_nPts(); i++)
        {
            double x = posDisMesh->nodeCoord(i, 0, IndexType::Local);
            double y = posDisMesh->nodeCoord(i, 1, IndexType::Local);
            double z = posDisMesh->nodeCoord(i, 2, IndexType::Local);

            posDHand->nodeDOFs->setValue("X", i, IndexType::Local, x);
            posDHand->nodeDOFs->setValue("Y", i, IndexType::Local, y);
            //posDHand->nodeDOFs->setValue("Z", i, IndexType::Local, z);

            double rad=sqrt(x*x+y*y);
            posDHand->nodeDOFs->setValue("Z", i, IndexType::Local, z + 0.00* (rad * rad - R*R));


            double nx1 = viscoDHand->nodeAuxF->getValue("nx",i,IndexType::Local);
            double ny1 = viscoDHand->nodeAuxF->getValue("ny",i,IndexType::Local);
            double nz1 = viscoDHand->nodeAuxF->getValue("nz",i,IndexType::Local);

            if(nz1<0.0)
                nz1=1.0;

            posDHand->nodeAuxF->setValue("nxx", i, IndexType::Local, nx1);
            posDHand->nodeAuxF->setValue("nyy", i, IndexType::Local, ny1);
            posDHand->nodeAuxF->setValue("nzz", i, IndexType::Local, nz1);



        }
        // InitializeNodalValues(posDHand, 3);
        posDHand->nodeDOFs0->setValue(posDHand->nodeDOFs);

        // find maximum z
        tensor<double, 1> X_cord1(posDisMesh->loc_nPts());

        for (int i = 0; i < posDisMesh->loc_nPts(); i++)
        {
            double x3 = posDHand->nodeDOFs->getValue(0, i, IndexType::Local);
            X_cord1(i) = x3;
        }


        //local maximum
        double X_max1 = zz_max(X_cord1, posDisMesh->loc_nPts());
        // 'global_max' will be the result on the root process.
        double global_max;
        //   MPI_Reduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

        // Perform the reduction operation MPI_MAX to find the global maximum.
        MPI_Reduce(&X_max1, &global_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        double X_max = global_max;
        if (myRank == 0)
            cout << "The maximum X  level is calaculated: " << X_max << endl;


        //Boundary conditions
        for(int i = 0; i < posDisMesh->loc_nPts(); i++)
        {
            double x = posDisMesh->nodeCoord(i, 0, IndexType::Local);
            double y = posDisMesh->nodeCoord(i, 1, IndexType::Local);
            double z = posDisMesh->nodeCoord(i, 2, IndexType::Local);

            int crease = posDisMesh->nodeCrease(i, hiperlife::IndexType::Local);
            if (crease > 0)
            {
               posDHand->setConstraint(0, i, hiperlife::IndexType::Local, 0.0);
               posDHand->setConstraint(1, i, hiperlife::IndexType::Local, 0.0);
               posDHand->setConstraint(2, i, hiperlife::IndexType::Local, 0.0);
            }


          
	   
	   
        
          // base_area= 3.141*r_tri*r_tri+3*(q_tri)*r_tri+sqrt(3)/4*q_tri*q_tri;
        //base_area=6.54584;
    
                if(choice==1000)
           {
                              if (y<r0_tri-r_tri+sp_gap)
                    {
                        posDHand->setConstraint(0, i, hiperlife::IndexType::Local, 0.0);
                        posDHand->setConstraint(1, i, hiperlife::IndexType::Local, 0.0);
                        posDHand->setConstraint(2, i, hiperlife::IndexType::Local, 0.0);
                    }

                    if ((y-sqrt(3)*x-r0_tri-2*r_tri+sp_gap)>0.0)
                    {
                        posDHand->setConstraint(0, i, hiperlife::IndexType::Local, 0.0);
                        posDHand->setConstraint(1, i, hiperlife::IndexType::Local, 0.0);
                        posDHand->setConstraint(2, i, hiperlife::IndexType::Local, 0.0);
                    }

                    if ((y+sqrt(3)*x-r0_tri-2*r_tri-(sqrt(3))*q_tri+sp_gap)>0.0) //done
                    {
                        posDHand->setConstraint(0, i, hiperlife::IndexType::Local, 0.0);
                        posDHand->setConstraint(1, i, hiperlife::IndexType::Local, 0.0);
                        posDHand->setConstraint(2, i, hiperlife::IndexType::Local, 0.0);

                    }


                    if(   (y>(sqrt(3)/2*q_tri+r0_tri+r_tri*0.5-sp_gap)) && (sqrt((x-0.5*q_tri)*(x-0.5*q_tri)+(y-r0_tri-sqrt(3)/2*q_tri)*(y-r0_tri-sqrt(3)/2*q_tri))>r_tri-sp_gap))

                    {
                        posDHand->setConstraint(0, i, hiperlife::IndexType::Local, 0.0);
                        posDHand->setConstraint(1, i, hiperlife::IndexType::Local, 0.0);
                        posDHand->setConstraint(2, i, hiperlife::IndexType::Local, 0.0);
                    }

                    if( (x>q_tri-sp_gap) && (y<r0_tri+0.5*r_tri+sp_gap) && (sqrt((x-q_tri)*(x-q_tri)+(y-r0_tri)*(y-r0_tri))>r_tri-sp_gap))
                    {
                        posDHand->setConstraint(0, i, hiperlife::IndexType::Local, 0.0);
                        posDHand->setConstraint(1, i, hiperlife::IndexType::Local, 0.0);
                        posDHand->setConstraint(2, i, hiperlife::IndexType::Local, 0.0);
                    }


                    if( (x<0.0+sp_gap) && (y<r0_tri+0.5*r_tri+sp_gap) && (sqrt((x-0.0)*(x-0.0)+(y-r0_tri)*(y-r0_tri))>r_tri-sp_gap))
                    {
                        posDHand->setConstraint(0, i, hiperlife::IndexType::Local, 0.0);
                        posDHand->setConstraint(1, i, hiperlife::IndexType::Local, 0.0);
                        posDHand->setConstraint(2, i, hiperlife::IndexType::Local, 0.0);

                    }
	   
	   
	   
	   }


           if (choice ==2000)
                {
                    if (sqrt(x*x+y*y) > R-sp_gap)
                    {
                        // boundary condition

                        posDHand->setConstraint(0, i, hiperlife::IndexType::Local, 0.0);
                        posDHand->setConstraint(1, i, hiperlife::IndexType::Local, 0.0);
                        posDHand->setConstraint(2, i, hiperlife::IndexType::Local, 0.0);


                    }
                }


            if (spring==0.5)
            {

                if (choice ==10)
                {
                    if (sqrt(x*x+y*y) > R-sp_gap)
                    {
                        // boundary condition

                        posDHand->setConstraint(0, i, hiperlife::IndexType::Local, 0.0);
                        posDHand->setConstraint(1, i, hiperlife::IndexType::Local, 0.0);
                        posDHand->setConstraint(2, i, hiperlife::IndexType::Local, 0.0);


                    }
                }

                if (choice ==11)
                {
                    // boundary condition

                    if (x< (0+sp_gap) || x>(Xmax-sp_gap) )
                    {
                        // boundary condition

                        posDHand->setConstraint(0, i, hiperlife::IndexType::Local, 0.0);
                        posDHand->setConstraint(1, i, hiperlife::IndexType::Local, 0.0);
                        posDHand->setConstraint(2, i, hiperlife::IndexType::Local, 0.0);
                        //posDHand->setConstraint(3, i, hiperlife::IndexType::Local, 0.0);
                        // posDHand->setConstraint(4, i, hiperlife::IndexType::Local, 0.0);
                        // posDHand->setConstraint(5, i, hiperlife::IndexType::Local, 0.0);
                        // posDHand->setConstraint(6, i, hiperlife::IndexType::Local, 0.0);

                    }
                }


                if (choice ==12)  // choice 1 clamp
                {
                    if (x< (0+sp_gap) || x>(Xmax-sp_gap) )
                    {
                        // boundary condition

                        posDHand->setConstraint(0, i, hiperlife::IndexType::Local, 0.0);
                        posDHand->setConstraint(1, i, hiperlife::IndexType::Local, 0.0);
                        posDHand->setConstraint(2, i, hiperlife::IndexType::Local, 0.0);
                        //posDHand->setConstraint(3, i, hiperlife::IndexType::Local, 0.0);
                        // posDHand->setConstraint(4, i, hiperlife::IndexType::Local, 0.0);
                        // posDHand->setConstraint(5, i, hiperlife::IndexType::Local, 0.0);
                        // posDHand->setConstraint(6, i, hiperlife::IndexType::Local, 0.0);

                    }

		        if (y< (0+sp_gap) || y>(1-sp_gap) )
                    {
                        // boundary condition

                        posDHand->setConstraint(0, i, hiperlife::IndexType::Local, 0.0);
                        posDHand->setConstraint(1, i, hiperlife::IndexType::Local, 0.0);
                        posDHand->setConstraint(2, i, hiperlife::IndexType::Local, 0.0);
                        //posDHand->setConstraint(3, i, hiperlife::IndexType::Local, 0.0);
                        // posDHand->setConstraint(4, i, hiperlife::IndexType::Local, 0.0);
                        // posDHand->setConstraint(5, i, hiperlife::IndexType::Local, 0.0);
                        // posDHand->setConstraint(6, i, hiperlife::IndexType::Local, 0.0);

                    }

                }


                    if (choice ==100 && crease > 0)
                    {
                        // boundary condition

                        posDHand->setConstraint(0, i, hiperlife::IndexType::Local, 0.0);
                        posDHand->setConstraint(1, i, hiperlife::IndexType::Local, 0.0);
                        posDHand->setConstraint(2, i, hiperlife::IndexType::Local, 0.0);
                    }

            }

        }
            /*posDHand->setConstraint(6,0.0);

               posDHand->setConstraint(3,0.0);
               posDHand->setConstraint(4, 0.0);
               posDHand->setConstraint(5, 0.0);*/


        posDHand->UpdateGhosts();
        posDHand->printFileLegacyVtk("checkinitialC", true);


        viscoDHand->printFileLegacyVtk("checkinitialV", true);

        EDHand = rcp(new DOFsHandler(posDisMesh));
        try
        {
            EDHand->setNameTag("EDHand");
            EDHand->setDOFs({"nx","ny","nz","height","EA","EB","DA","DB","PA","trace_sigma","dev_sigma"});
            EDHand->Update();
        }
        catch (runtime_error &err)
        {
            cout << myRank << ": DOFHandler could not be created " << err.what() << endl;
            MPI_Finalize();
            return 1;
        }

        EDHand->UpdateGhosts();
        for (int i = 0; i < posDisMesh->loc_nPts(); i++)
        {
            //
        }

        EDHand->nodeDOFs0->setValue(EDHand->nodeDOFs);

        EDHand->UpdateGhosts();
        EDHand->printFileLegacyVtk("checkinitialE", true);

//

        gloDHand = rcp(new DOFsHandler(gloDisMesh));
        try
        {

            gloDHand->setNameTag("gloDHand");
            gloDHand->setDOFs({"P", "FX", "FY", "FZ", "MX", "MY", "MZ"});
            gloDHand->Update();
        }
        catch (runtime_error &err)
        {
            cout << myRank << ": DOFHandler could not be created " << err.what() << endl;
            MPI_Finalize();
            return 1;
        }

        gloDHand->setInitialCondition(0, 0.0);
        // gloDHand->setConstraint(0,0.0);


       // gloDHand->setConstraint(0,0.0);

        gloDHand->setConstraint(1, 0.0);
        gloDHand->setConstraint(2, 0.0);
        gloDHand->setConstraint(3, 0.0);
        gloDHand->setConstraint(4, 0.0);
        gloDHand->setConstraint(5, 0.0);
        gloDHand->setConstraint(6, 0.0);
        gloDHand->UpdateGhosts();

        string fMesh = "sol_gloCons";
        gloDHand->printFile(fMesh, OutputMode::Text, true);


    }

    else
    {
        if (myRank == 0)
        cout<<"Starting from a restart" <<restart<<endl;
        //Cortex
        try
        {
            posDHand = rcp(new DOFsHandler("posDHand"));
            string fMesh = "sol_pos." + to_string(restart);
            posDHand->setFilePrefix(fMesh, OutputMode::Text);
            posDHand->setNumNodeAuxF(14);
            posDHand->Update();
        }
        catch(runtime_error err)
        {
            cout << myRank << ":  Position DOFsHandler could not be created. " << err.what() << endl;

            MPI_Finalize();
            return 1;
        }


//VISCO
        try
        {
            viscoDHand = rcp(new DOFsHandler("viscoDHand"));
            string fMesh = "sol_visco." + to_string(restart);
            viscoDHand->setFilePrefix(fMesh, OutputMode::Text);
            viscoDHand->setNumNodeAuxF(9);
            viscoDHand->Update();
        }
        catch(runtime_error err)
        {
            cout << myRank << ":  visco DOFsHandler could not be created. " << err.what() << endl;

            MPI_Finalize();
            return 1;
        }
//ENERGY HAND
        try
        {
            EDHand = rcp(new DOFsHandler("EDHand"));
            string fMesh = "sol_EN." + to_string(restart);
            EDHand->setFilePrefix(fMesh, OutputMode::Text);
            EDHand->Update();
        }
        catch(runtime_error err)
        {
            cout << myRank << ":  ENERGY DOFsHandler could not be created. " << err.what() << endl;

            MPI_Finalize();
            return 1;
        }



        //Global constraints
        try
        {
            gloDHand = rcp(new DOFsHandler("gloDHand"));
            string fMesh = "sol_gloCons." + to_string(restart);
            gloDHand->setFilePrefix(fMesh, OutputMode::Text);
            gloDHand->Update();
        }
        catch(runtime_error err)
        {
            cout << myRank << ":   gloCons DOFsHandler could not be created. " << err.what() << endl;

            MPI_Finalize();
            return 1;
        }


        //disMesh
        posDisMesh = posDHand->mesh;

    }

    
/*
    for (int i = 0; i < posDisMesh->loc_nPts(); i++) //FIXME: assign random initial?
    {
        auto x = posDisMesh->nodeCoords(i,IndexType::Local);
       // Z_cord(i)=x[2]+postHand->nodeAuxF->getValue("Uz",i,IndexType::Local);
        posDisMesh->_nodeData->setValue(0,i,IndexType::Local,x[0]+posDHand->nodeAuxF->getValue("Ux",i,IndexType::Local));
        posDisMesh->_nodeData->setValue(1,i,IndexType::Local,x[1]+posDHand->nodeAuxF->getValue("Uy",i,IndexType::Local));
        posDisMesh->_nodeData->setValue(2,i,IndexType::Local,x[2]+posDHand->nodeAuxF->getValue("Uz",i,IndexType::Local));
    }

    posDisMesh->printFileLegacyVtk("newMesh1");

    
return 0; */


    //Resize auxiliary fields
    paramStr->a_aux.resize(posDisMesh->loc_nElem()*gPts*58);//just a definition
    paramStr->b_aux.resize(posDisMesh->loc_nElem()*gPts*11);

    // **************************************************************//
    // *****               HIPERPROBLEM CREATION                *****//
    // **************************************************************//
    SmartPtr<HiPerProblem> hiperProbl = Create<HiPerProblem>();
    try
    {
        // Set paramStructure
        hiperProbl->setParameterStructure(paramStr);
        hiperProbl->setConsistencyCheckDelta(1.E-8);
        hiperProbl->setConsistencyCheckTolerance(1.E-4);

        // Set DOFHandler
        hiperProbl->setDOFsHandlers({posDHand,gloDHand});

        // Set Integration
        hiperProbl->setIntegration("Integ", {"posDHand","gloDHand"});
        hiperProbl->setCubatureGauss("Integ", gPts);
        //hiperProbl->setElementFillings("Integ", LS);
        if (ConsCheck == 1)
        {
            hiperProbl->setElementFillings("Integ", ConsistencyCheck<LS>);
        }
        else
        {
            hiperProbl->setElementFillings("Integ", LS);
        }
        hiperProbl->setGlobalIntegrals({"Energy","Dissipation","Trd","volume","area","area_n","Epress","E1","E_tension1","P_tension1","E_elastic"});

        // Update
        hiperProbl->Update();
    }
    catch (runtime_error& err)
    {
        cout << myRank << ": HiPerProblem could not be created " << err.what() << endl;
        MPI_Finalize();
        return 1;
    }


    // create visco problem

    SmartPtr<HiPerProblem> viscoProbl= Create<HiPerProblem>();
    viscoProbl->setParameterStructure(paramStr);
    viscoProbl->setDOFsHandlers({viscoDHand});
    viscoProbl->setIntegration("Integ", {"viscoDHand"});
    viscoProbl->setCubatureGauss("Integ", gPts);
    viscoProbl->setElementFillings("Integ", LS_ReacDif);
    viscoProbl->setGlobalIntegrals({"Diss1","Diss2"});
    viscoProbl->Update();

    // create ENERGY problem

    SmartPtr<HiPerProblem> ENProbl= Create<HiPerProblem>();
    ENProbl->setParameterStructure(paramStr);
    ENProbl->setDOFsHandlers({EDHand});
    ENProbl->setIntegration("Integ", {"EDHand"});
    ENProbl->setCubatureGauss("Integ", gPts);
    ENProbl->setElementFillings("Integ", LS_ED);
    ENProbl->Update();



    // **************************************************************//
    // *****                 SOLVERS' CREATION                  *****//
    // **************************************************************//


    RCP<MUMPSDirectLinearSolver> linSolver = rcp(new MUMPSDirectLinearSolver());
    linSolver->setHiPerProblem(hiperProbl);
    linSolver->setMatrixType(MUMPSDirectLinearSolver::MatrixType::General);
    linSolver->setAnalysisType(MUMPSDirectLinearSolver::AnalysisType::Parallel);
    linSolver->setOrderingLibrary(MUMPSDirectLinearSolver::OrderingLibrary::Auto);
    linSolver->setVerbosity(MUMPSDirectLinearSolver::Verbosity::None);
    linSolver->setDefaultParameters();
    linSolver->setWorkSpaceMemoryIncrease(1000);
    linSolver->Update();

    RCP<NewtonRaphsonNonlinearSolver> nonlinSolver = rcp(new NewtonRaphsonNonlinearSolver());
    nonlinSolver->setLinearSolver(linSolver);
    nonlinSolver->setMaxNumIterations(MAXITER_NR);
    nonlinSolver->setResTolerance(RESTOL_NR);
    nonlinSolver->setSolTolerance(SOLTOL_NR);
    nonlinSolver->setLineSearch(true);
    nonlinSolver->setConvRelTolerance(false);
    nonlinSolver->setPrintIntermInfo(true);
    nonlinSolver->setPrintSummary(false);
    nonlinSolver->setResMaximum(1E4);
    nonlinSolver->setSolMaximum(1E4);
    nonlinSolver->setExitRelMaximum(1E4);
    nonlinSolver->Update();

    RCP<MUMPSDirectLinearSolver>viscoDirSolver= rcp(new MUMPSDirectLinearSolver());
    viscoDirSolver->setHiPerProblem(viscoProbl);
    viscoDirSolver->setMatrixType(MUMPSDirectLinearSolver::MatrixType::General);
    viscoDirSolver->setAnalysisType(MUMPSDirectLinearSolver::AnalysisType::Parallel);
    viscoDirSolver->setOrderingLibrary(MUMPSDirectLinearSolver::OrderingLibrary::Auto);
    viscoDirSolver->setVerbosity(MUMPSDirectLinearSolver::Verbosity::None);
    viscoDirSolver->setDefaultParameters();
    viscoDirSolver->setWorkSpaceMemoryIncrease(1000);
    viscoDirSolver->Update();


    RCP<MUMPSDirectLinearSolver>ENDirSolver= rcp(new MUMPSDirectLinearSolver());
    ENDirSolver->setHiPerProblem(ENProbl);
    ENDirSolver->setMatrixType(MUMPSDirectLinearSolver::MatrixType::General);
    ENDirSolver->setAnalysisType(MUMPSDirectLinearSolver::AnalysisType::Parallel);
    ENDirSolver->setOrderingLibrary(MUMPSDirectLinearSolver::OrderingLibrary::Auto);
    ENDirSolver->setVerbosity(MUMPSDirectLinearSolver::Verbosity::None);
    ENDirSolver->setDefaultParameters();
    ENDirSolver->setWorkSpaceMemoryIncrease(1000);
    ENDirSolver->Update();


    hiperProbl->FillLinearSystem();



    viscoProbl->FillLinearSystem();




    hiperProbl->FillLinearSystem();


    ENProbl->FillLinearSystem();

    double a0{3.141};

  //  double  base_area= 3.141*R*R;
    double  base_area= hiperProbl->globalIntegral("area_n");

    if(choice==1000)
     {
	     base_area= 3.141*r_tri*r_tri+3*(q_tri)*r_tri+sqrt(3)/4*q_tri*q_tri;


    }

 if(choice==2000)
     {
             base_area= 3.141*R*R;


    }

if(choice==10)
     {
             base_area= 3.141*R*R;


    }

    v_target=v_target/hiperProbl->globalIntegral("area_n")*base_area*0.999;
    if (myRank == 0)
    {
        cout << " " << " Volume: " << hiperProbl->globalIntegral("volume") << endl;
        cout << std::scientific << std::setprecision(6);

        a0 = hiperProbl->globalIntegral("area_n");
        cout << "Target volume:  " << abs(v_target)*a0<< endl;

        cout << " " << " Area : phi: " << hiperProbl->globalIntegral("area_n")<< ": "<< hiperProbl->globalIntegral("area") << endl;
    }

double a_res=a0-base_area;



    // globalintegral
    //Open file to write global integrals and write headers
    ofstream gIntegFile;
    gIntegFile.open ("globalIntegrals.dat");
    gIntegFile << "TS time deltat";
    //print
    for (auto g:  hiperProbl->globIntegralNames())
        gIntegFile << " " << g;
        gIntegFile << " pressure" ;
        gIntegFile << " Zmax" ;
        gIntegFile << " stress" ;
        gIntegFile << " strain" ;
    for (auto g:  viscoProbl->globalIntegralNames())
        gIntegFile << " " << g;

    gIntegFile << endl;

    double fstep=0;
   // double  factor=1;
    double astep=0;
    double mstep=0;
    double press0=0.0;
    double stress=0.0;
    double strain=0.0;
    double Z_max=1.0;
    // Time loop
    int timeStep=restart;
    int stretch_ind=0;

//print initial

    string solName = "sol_dis." + to_string(timeStep);
    posDHand->printFileLegacyVtk(solName,true);
    solName = "sol_visco." + to_string(timeStep);
    viscoDHand->printFileLegacyVtk(solName,true);

    int fric_start = gap+vnstep;

     double force1=0.0;
    double v0=v_target*a0;

    double deltatMax1=deltatMax;
    double factor1=factor;
    //  factor=v_target;
    if (myRank == 0)
        cout<<"input factor for volume 0.333: "<< factor<<endl;

    while ((tSimu < totalTime) and (timeStep < totalSteps))
    {
        // Print info
        if (myRank == 0)
        {
            cout<< "TS: " << timeStep + 1  << " Time " << tSimu << " of " << totalTime << " with deltat=" << deltat <<": del_max: "<<deltatMax<< endl;
            cout << "Starting Newton-Raphson iteration for membrane evolution" << endl;
        }




        // reduce volume in decreaments
        if(timeStep<vnstep)
        {
            factor=(timeStep+1.0)/vnstep*v_target;
            force1=(timeStep+1.0)/vnstep*force;



            vstep=vstep+1;
            if (myRank == 0)
                cout<< "volume increased step: "<< vstep<<" fraction: "<<(timeStep+1.0)/vnstep<<endl;
            vol_inc=1;
        }

        if(timeStep>vnstep+gap-control_steps)
        {
            deltat *= stepFactor;
            deltatMax=control_dt;
            if (myRank == 0)
                cout << "time reduction step: " << timeStep-(vnstep+gap-control_steps)<<" dt: " << deltat<< " :max dt:  " << deltatMax <<endl;
        }

        if(timeStep>vnstep+gap+v_cont*v_rn_step)
        {
            deltat /= stepFactor;
            deltatMax=deltatMax1;

            if (myRank == 0)
                cout << "trelaxing pahse: " << timeStep<<" dt: " << deltat<< " :max dt:  " << deltatMax <<endl;
        }


        if(timeStep>vnstep+gap)
        {
            if (v_r_step<v_cont*v_rn_step)
            {
                factor=(1-(v_r_step+1.0)/v_rn_step)*v_target;
               // factor=factor*exp(coeff*deltat*v_r_step);
                v_r_step=v_r_step+1;
                if (myRank == 0)
                    cout << "volume reduction step: " << v_r_step <<"by factor: "<<factor/v_target<< endl;

            }



        }

        paramStr->setRealParameter(MembParams::factor,factor);
        paramStr->setRealParameter(MembParams::force,force1);


         if (myRank == 0)
           cout<< "factor: "<<  paramStr->getRealParameter(MembParams::factor)<<endl;

        paramStr->setRealParameter(MembParams::vol_inc,vol_inc);
        paramStr->setIntParameter(MembParams::timestep,timeStep);


        if (deltat<0.000001)
        {
            if (myRank == 0)
                cout << "stopping time dt: " << endl;
            return 0;
        }



        // Prepare for solver
        posDHand->nodeDOFs->setValue(posDHand->nodeDOFs0);
        gloDHand->nodeDOFs->setValue(gloDHand->nodeDOFs0);
        hiperProbl->UpdateGhosts();
        //Prepare for solver
        viscoDHand->nodeDOFs0->setValue(viscoDHand->nodeDOFs);
        viscoDHand->UpdateGhosts();

        EDHand->nodeDOFs0->setValue(EDHand->nodeDOFs);
        EDHand->UpdateGhosts();

        // Solve the problem
        bool converged = nonlinSolver->solve();

        // Print info
        if (myRank == 0)
            cout << "Finished Newton-Raphson iteration" << endl;

        // Check convergence
        if (converged)
        {

            hiperProbl->FillLinearSystem();



            viscoProbl->FillLinearSystem();

            if (myRank == 0)
                cout << "  Solving Energy problem: "<<endl;
            // ENERGY PROBLEM
            bool ENConverged = ENDirSolver->solve();

            if(ENConverged)
            {
                if (myRank == 0)
                    cout << "ENERGY Direct solver converged: " << endl;
            }
            ENProbl->UpdateSolution();
           // ENProbl->FillLinearSystem();

            tSimu += deltat;
            timeStep += 1;

            press0= gloDHand->nodeDOFs->getValue(0,0)/deltat;



            for (int i = 0; i< posDHand->mesh->loc_nPts();i++)
            {

                double x0 = posDHand->nodeDOFs0->getValue(0,i,IndexType::Local);
                double y0 = posDHand->nodeDOFs0->getValue(1,i,IndexType::Local);
                double z0 = posDHand->nodeDOFs0->getValue(2,i,IndexType::Local);

                double x = posDHand->nodeDOFs->getValue(0,i,IndexType::Local);
                double y = posDHand->nodeDOFs->getValue(1,i,IndexType::Local);
                double z = posDHand->nodeDOFs->getValue(2,i,IndexType::Local);

            double xr = posDisMesh->nodeCoord(i, 0, IndexType::Local);
            double yr = posDisMesh->nodeCoord(i, 1, IndexType::Local);
            double zr = posDisMesh->nodeCoord(i, 2, IndexType::Local);

                posDHand->nodeAuxF->setValue("Ux",i,IndexType::Local,(x-xr));
                posDHand->nodeAuxF->setValue("Uy",i,IndexType::Local,(y-yr));
                posDHand->nodeAuxF->setValue("Uz",i,IndexType::Local,(z-zr));
                //updated lagrangian
                /*   posDHand->mesh->_nodeData->setValue(0,i,IndexType::Local,x);
                   posDHand->mesh->_nodeData->setValue(1,i,IndexType::Local,y);
                   posDHand->mesh->_nodeData->setValue(2,i,IndexType::Local,z);*/


                double nxx = EDHand->nodeDOFs->getValue(0,i,IndexType::Local);
                double nyy = EDHand->nodeDOFs->getValue(1,i,IndexType::Local);
                double nzz = EDHand->nodeDOFs->getValue(2,i,IndexType::Local);
                double height1 = EDHand->nodeDOFs->getValue(3,i,IndexType::Local);

                double EA1 = EDHand->nodeDOFs->getValue(4,i,IndexType::Local);
                double EB1 = EDHand->nodeDOFs->getValue(5,i,IndexType::Local);

                double DA1 = EDHand->nodeDOFs->getValue(6,i,IndexType::Local);
                double DB1 = EDHand->nodeDOFs->getValue(7,i,IndexType::Local);
                double PA1 = EDHand->nodeDOFs->getValue(8,i,IndexType::Local);


                double trace_sigma1 = EDHand->nodeDOFs->getValue(9,i,IndexType::Local);
                double dev_sigma1 = EDHand->nodeDOFs->getValue(10,i,IndexType::Local);


                posDHand->nodeAuxF->setValue("nx",i,IndexType::Local,nxx);
                posDHand->nodeAuxF->setValue("ny",i,IndexType::Local,nyy);
                posDHand->nodeAuxF->setValue("nz",i,IndexType::Local,nzz);
                posDHand->nodeAuxF->setValue("height",i,IndexType::Local,height1);

                posDHand->nodeAuxF->setValue("EA",i,IndexType::Local,EA1);
                posDHand->nodeAuxF->setValue("EB",i,IndexType::Local,EB1);

                posDHand->nodeAuxF->setValue("DA",i,IndexType::Local,DA1);
                posDHand->nodeAuxF->setValue("DB",i,IndexType::Local,DB1);
                posDHand->nodeAuxF->setValue("PA",i,IndexType::Local,PA1);

                posDHand->nodeAuxF->setValue("trace_sigma",i,IndexType::Local,trace_sigma1);
                posDHand->nodeAuxF->setValue("dev_sigma",i,IndexType::Local,dev_sigma1);


            }

            posDHand->nodeDOFs0->setValue(posDHand->nodeDOFs);
           // viscoDHand->nodeDOFs0->setValue(viscoDHand->nodeDOFs);
            EDHand->nodeDOFs0->setValue(EDHand->nodeDOFs);

            hiperProbl->UpdateGhosts();
          //  viscoProbl->UpdateGhosts();
            ENProbl->UpdateGhosts();

            if (timeStep%nPrint==0)
            {
                string solName;
                solName = sol_prefixMesh + to_string(timeStep);
                if (myRank == 0)
                    cout << "Printing file " <<  solName << endl;
                // ppHand->printFileLegacyVtk(solName,true);

            }

            if (nonlinSolver->numberOfIterations() < MAXITER_NR and deltat < deltatMax)
            {
                deltat /= stepFactor;
                if (myRank == 0)
                    cout << "dividing dt "  << endl;
            }

        }


        else
        {
            if (myRank == 0)
                std::cout << "ERROR::EXIT from Newton-Raphson but linear solver converged" << endl;

            // Modify time-step
            deltat *= stepFactor;

            // Restore initial values
            posDHand->nodeDOFs->setValue(posDHand->nodeDOFs0);
            gloDHand->nodeDOFs->setValue(gloDHand->nodeDOFs0);
            viscoDHand->nodeDOFs->setValue(viscoDHand->nodeDOFs0);
            EDHand->nodeDOFs->setValue(EDHand->nodeDOFs0);
        }



        // Compute factor to reduce the volume
        double v1 = hiperProbl->globalIntegral("volume");
        double a1 = hiperProbl->globalIntegral("area");


        if (myRank == 0)
        {
            cout << std::scientific << std::setprecision(6);
            cout << " " << " VOLUME: " << hiperProbl->globalIntegral("volume") << endl;
            cout << " " << " AREA jac: " << hiperProbl->globalIntegral("area") << endl;
            cout << " " << " AREA_0 : jacR: " << hiperProbl->globalIntegral("area_n") << endl;
            cout << " " << " AREA_ jac_C*jacR: " << hiperProbl->globalIntegral("Trd") << endl;
            cout << " " << " stretch: " << (hiperProbl->globalIntegral("area")-a_res)/(hiperProbl->globalIntegral("area_n")-a_res) << endl;

            cout << " " << " volume fraction: " << v1/v0<< endl;
            // cout << " " << " volume needed to be reduced by: " << pow(Vf,timeStep)<< endl;
            cout << " " << " Pressure work: " << -1*hiperProbl->globalIntegral("Epress") << endl;

            cout << " " << " Total strain energy: " <<  hiperProbl->globalIntegral("Energy")  <<   "  Bending energy: "<< hiperProbl->globalIntegral("E1")<<endl;
            cout << " " << " Total energy: " <<  hiperProbl->globalIntegral("Energy")  + hiperProbl->globalIntegral("E1")<<endl;
            cout << " " << " Small Elastic energy: " <<  hiperProbl->globalIntegral("E_elastic") <<endl;

            double p0= gloDHand->nodeDOFs->getValue(0,0)/deltat;
            cout << " " << " pressure : " << -p0 << endl;
            cout << " " << "Excess area: " << abs(hiperProbl->globalIntegral("area")-hiperProbl->globalIntegral("area_n"))/base_area<< endl;


            cout << " " << " Frictional Dissipation: " << hiperProbl->globalIntegral("Dissipation") << endl;

            cout << " " << " Power Input: " << viscoProbl->globalIntegral("P_tension1") << endl;

            cout << " " << " Viscous Dissipation: " << viscoProbl->globalIntegral("Diss1") << endl;

            cout << " " << " curvature Dissipation: " << viscoProbl->globalIntegral("Diss2") << endl;


        }





        // Print solution
        if (timeStep%nPrint==0)
        {
            string solName = "sol_dis." + to_string(timeStep);
            posDHand->printFileLegacyVtk(solName,true);
            //solName = "sol_gloCons." + to_string(timeStep);
            //gloDHand->printFileLegacyVtk(solName,true);
            // solName = "sol_visco." + to_string(timeStep);
            // viscoDHand->printFileLegacyVtk(solName,true);

        }




        if (timeStep % nPrint == 0)
        {
            //In the time loop print global integrals
            //Write global integrals
            if (myRank == 0)
            {
                gIntegFile << timeStep << " " << tSimu << " " << deltat;
                for (auto g:  hiperProbl->globIntegrals())
                    gIntegFile << " " << g;
                gIntegFile << " " << press0;
                gIntegFile << " " << Z_max;
                gIntegFile << " " << stress;
                gIntegFile << " " << strain;
                for (auto g:  viscoProbl->globIntegrals())
                    gIntegFile << " " << g;
                gIntegFile << endl;
            }
        }

        if(timeStep>600)
        {
            if (timeStep%(hl_print) ==0)
            {

                string solName = "sol_pos." + to_string(timeStep);
                posDHand->printFile(solName, OutputMode::Text, true, tSimu);


                solName = "sol_gloCons." + to_string(timeStep);
                gloDHand->printFile(solName, OutputMode::Text, true, tSimu);

                solName = "sol_visco." + to_string(timeStep);
                viscoDHand->printFile(solName, OutputMode::Text, true, tSimu);

                solName = "sol_EN." + to_string(timeStep);
                EDHand->printFile(solName, OutputMode::Text, true, tSimu);


            }
        }





        // Update time-related quantities
        paramStr->setRealParameter(MembParams::deltat,deltat);

    }



    // **************************************************************//
    // *****                    FINALIZE                        *****//
    // **************************************************************//

    MPI_Finalize();
    return 0;
}
