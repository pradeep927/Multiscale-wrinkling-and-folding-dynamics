

/// hiperlife headers
#include "hl_FillStructure.h"
#include "hl_Geometry.h"
#include "hl_SurfLagrParam.h"
#include "hl_Tensor.h"
//#include "hl_LinearSolver_Direct_Amesos2.h"
#include <hl_LinearSolver_Direct_MUMPS.h>
#include <fstream>
#include "hl_Parser.h"

/// Header to auxiliary functions
#include "AuxViscousInextMembrane.h"

void LS(hiperlife::FillStructure& fillStr)
{

    using namespace hiperlife;
    using namespace hiperlife::Tensor;

    //------------------------------------------------------------------
    // [2] Input variables

    //[1.1] Positions
    auto &subFill = (fillStr)["posDHand"];
    int nDim = subFill.nDim;
    int pDim = subFill.pDim;
    int eNN = subFill.eNN;
    int numDOFs = subFill.numDOFs;
    int numAuxF = subFill.numAuxF;


    wrapper<double, 2> nborDOFs0(subFill.nborDOFs0.data(), eNN, numDOFs);
    wrapper<double, 2> nborDOFs(subFill.nborDOFs.data(), eNN, numDOFs);
    wrapper<double, 2> nborCoords(subFill.nborCoords.data(), eNN, nDim);
    wrapper<double,2> nborAuxF(subFill.nborAuxF.data(),eNN,numAuxF);

    wrapper<double, 1> bf(subFill.nborBFs(), eNN);
    wrapper<double, 2> Dbf(subFill.nborBFsGrads(), eNN, pDim);
    wrapper<double, 3> DDbf(subFill.nborBFsHess(), eNN, pDim, pDim);

    //[1.3] Global constraints
    auto &g_subFill = (fillStr)["gloDHand"];
    int g_numDOFs = g_subFill.numDOFs;

    wrapper<double, 1> gDOFs(g_subFill.nborDOFs.data(), g_numDOFs);
    double pressure = gDOFs(0);
    tensor<double, 1> F = gDOFs(range(1, 3));
    tensor<double, 1> L = gDOFs(range(4, 6));

    //[1.4] Parameters


 double deltat = fillStr.getRealParameter(MembParams::deltat);
 double fric = fillStr.getRealParameter(MembParams::fric);
 double factor = fillStr.getRealParameter(MembParams::factor);
 double young = fillStr.getRealParameter(MembParams::young)*deltat;



    double nu = fillStr.getRealParameter(MembParams::poisson);

 double thick = fillStr.getRealParameter(MembParams::thick);
 double Kconf = fillStr.getRealParameter(MembParams::Kconf) * deltat;

int vnstep=    fillStr.getIntParameter(MembParams::vnstep);
 int gap=   fillStr.getIntParameter(MembParams::gap);


    double gamma = fillStr.getRealParameter(MembParams::gamma)* deltat;
    int tens_start = fillStr.getIntParameter(MembParams::tens_start);
    int timestep = fillStr.getIntParameter(MembParams::timestep);//timestep
    int case_sphere = fillStr.getIntParameter(MembParams::case_sphere);

 double R    = fillStr.getRealParameter(MembParams::R);
 double width   = fillStr.getRealParameter(MembParams::width);
 
    double kappa = fillStr.getRealParameter(MembParams::kappa) * deltat;
    double aap    = fillStr.getRealParameter(MembParams::aap);
    double bbn    = fillStr.getRealParameter(MembParams::bbn);
    double fric3    = fillStr.getRealParameter(MembParams::fric3);

     double apical    = fillStr.getRealParameter(MembParams::apical);
    double basal    = fillStr.getRealParameter(MembParams::basal);

    double a11= fillStr.getRealParameter(MembParams::a11);
    double a12= fillStr.getRealParameter(MembParams::a12);
    double a33= fillStr.getRealParameter(MembParams::a33);
    double kap1= fillStr.getRealParameter(MembParams::kap1);

     double force=fillStr.getRealParameter(MembParams::force)*deltat;

    double fric_fact = fillStr.getRealParameter(MembParams::fric_fact);


    double fric2 = fillStr.getRealParameter(MembParams::fric2);
    double height_in = fillStr.getRealParameter(MembParams::height_in);
    int choice = fillStr.getIntParameter(MembParams::choice);
    double Xmax = fillStr.getRealParameter(MembParams::Xmax);
    int control_fric = fillStr.getIntParameter(MembParams::control_fric);


    double Lagrangian{};
    double f0 = fillStr.getRealParameter(MembParams::f0);
    double g_ratio = fillStr.getRealParameter(MembParams::g_ratio);
    double kspr= fillStr.getRealParameter(MembParams::kspr) * deltat;
    double sp_gap= fillStr.getRealParameter(MembParams::sp_gap) ;

    double tens_factor   = fillStr.getRealParameter(MembParams::tens_factor);
    int spring= fillStr.getIntParameter(MembParams::spring);
    int fric_start = fillStr.getIntParameter(MembParams::fric_start);


    double crypt=fillStr.getRealParameter(MembParams::crypt);
    double gamma_minus=fillStr.getRealParameter(MembParams::gamma_minus) * deltat;
    double gamma_plus=fillStr.getRealParameter(MembParams::gamma_plus) * deltat;
    double gamma_l=fillStr.getRealParameter(MembParams::gamma_l) * deltat;
    double fact_elastic=fillStr.getRealParameter(MembParams::fact_elastic);

    double pn=fillStr.getRealParameter(MembParams::pn) ;

    double gamma_plus_ref=fillStr.getRealParameter(MembParams::gamma_plus_ref)* deltat;
    double gamma_minus_ref=fillStr.getRealParameter(MembParams::gamma_minus_ref) * deltat;

    int numDOFs1=3;
    //OUTPUTS
    wrapper<double, 2> Bp(fillStr.Bk(0).data(), eNN, numDOFs);
    wrapper<double, 1> Bg(fillStr.Bk(1).data(), g_numDOFs);

    wrapper<double, 4> App(fillStr.Ak(0, 0).data(), eNN, numDOFs, eNN, numDOFs);
    wrapper<double, 3> Apg(fillStr.Ak(0, 1).data(), eNN, numDOFs, g_numDOFs);
    wrapper<double, 3> Agp(fillStr.Ak(1, 0).data(), g_numDOFs, eNN, numDOFs);


    //------------------------------------------------------------------
  // [2] Compute variables
    // [2.1] Compute geometry in the reference
    tensor<double, 1> xRv(nDim);
    tensor<double, 1> xRu(nDim);
    tensor<double, 1> xRuu(nDim);
    tensor<double, 1> xRuv(nDim);
    tensor<double, 1> xRvv(nDim);
    tensor<double, 1> normalR(nDim);
    tensor<double, 2> metricR(pDim, pDim);
    tensor<double, 2> curvatureR(pDim, pDim);
    tensor<double, 3> cristSymR(pDim, pDim, pDim);
    SurfLagrParam::ChristoffelSymbols(cristSymR, curvatureR, metricR, normalR, xRu, xRv, xRuu, xRvv, xRuv, eNN,
                                      nborCoords, Dbf, DDbf);
    tensor<double, 1> xR = bf * nborCoords;
    tensor<double, 2> imetricR = metricR.inv();
    double jacR = sqrt(metricR.det());

    tensor<double, 2> nborDOFs_xyz0(eNN, 3);
    tensor<double, 2> nborDOFs_xyz(eNN, 3);

     nborDOFs_xyz0 = nborDOFs0(all,range(0,2));//dof:({"offset","c1","c2","cDens"});
     nborDOFs_xyz = nborDOFs(all,range(0,2));//dof:({"offset","c1","c2","cDens"});

   // cout <<"numDOFs "<< nborDOFs_xyz<<endl;




    //[2.1] Previous time-step
    tensor<double, 1> xu_n(nDim);
    tensor<double, 1> xv_n(nDim);
    tensor<double, 1> normal_n(nDim);
    tensor<double, 2> metric_n(pDim, pDim);
    SurfLagrParam::Normal(normal_n, metric_n, xu_n, xv_n, eNN,  nborDOFs_xyz0, Dbf);

    tensor<double, 1> xuu_n(nDim);
    tensor<double, 1> xuv_n(nDim);
    tensor<double, 1> xvv_n(nDim);
    tensor<double, 2> curva_n(pDim, pDim);
    tensor<double, 3> christsym_n(pDim, pDim, pDim);
    SurfLagrParam::ChristoffelSymbols(christsym_n, curva_n, metric_n, normal_n, xu_n, xv_n, xuu_n, xvv_n, xuv_n, eNN,
                                      nborDOFs_xyz0, Dbf, DDbf);

    tensor<double, 1> x_n = bf *  nborDOFs_xyz0;
    tensor<double, 2> imetric_n = metric_n.inv();



    double jac_n = sqrt(metric_n.det());
    double xnormal_n = x_n * normal_n;

    // UNIFORMLY DISTRIBUTED FORCE
    // UNIFORMLY DISTRIBUTED FORCE
    tensor<double, 1> Coords = bf * nborCoords;


    double rad=sqrt(xR(0)*xR(0)+xR(1)*xR(1));




       // fric = fric * fric_fact + (fric * fric_fact - fric) * tanh(1.5*width * (1.0*thick - x_n(2)));
        //double fact=0.5+0.5* tanh(width * (1.4*height_in - x_n(2)));

 double arg = width * (1.05 * height_in - x_n(2)); //1.05
 double fact = 0.5 + 0.5 * std::tanh(arg);


        double kspr_in=kspr*fact;





    //[2.2] Current time-step
    tensor<double, 1> xu(nDim);
    tensor<double, 1> xv(nDim);
    tensor<double, 1> xuu(nDim);
    tensor<double, 1> xuv(nDim);
    tensor<double, 1> xvv(nDim);
    tensor<double, 1> normal(nDim);
    tensor<double, 2> metric(pDim, pDim);
    tensor<double, 2> curva(pDim, pDim);
    tensor<double, 3> christsym(pDim, pDim, pDim);
    SurfLagrParam::ChristoffelSymbols(christsym, curva, metric, normal, xu, xv, xuu, xvv, xuv, eNN,  nborDOFs_xyz, Dbf,
                                      DDbf);


    tensor<double, 1> x = bf *  nborDOFs_xyz;




    tensor<double, 2> imetric = metric.inv();
    double jac = sqrt(metric.det());
    double meancurva = product(imetric, curva, {{0, 0},{1, 1}});
    double xnormal = x * normal;
    double xRnormalR = xR * normalR;
    //3D unit tensor
    tensor<double, 2> Id(3, 3);
    Id(0, 0) = 1;
    Id(0, 1) = 0;
    Id(0, 2) = 0;
    Id(1, 0) = 0;
    Id(1, 1) = 1;
    Id(1, 2) = 0;
    Id(2, 0) = 0;
    Id(2, 1) = 0;
    Id(2, 2) = 1;

    // rate of deformation tensor
   // tensor<double, 2> rodt = metric - metric_n;
   // tensor<double, 2> rodt_CnCn = imetric_n * rodt * imetric_n;

    //[2.3] Other variables
  //  double tr_rodt = (jac - jac_n) / jac_n;


    //------------------------------------------------------------------

    //------------------------------------------------------------------
    // [3.1] First derivatives

    tensor<double, 4> d_metric(eNN, numDOFs1, pDim, pDim);
    SurfLagrParam::d_Metric(d_metric, Dbf, xu, xv);

    tensor<double, 4> d_metric_CC = product(product(d_metric, imetric, {{2, 0}}), imetric, {{2, 0}});
    tensor<double, 4> d_metric_CnCn = product(product(d_metric, imetric_n, {{2, 0}}), imetric_n, {{2, 0}});

    tensor<double, 4> d_imetric = -1.0 * d_metric_CC;

    tensor<double, 2> d_jac(eNN, numDOFs1);
    SurfLagrParam::d_Jac(d_jac, Dbf, xu, xv, normal);


    tensor<double, 3> d_normal(eNN, numDOFs1, nDim);
    SurfLagrParam::d_Normal(d_normal, Dbf, xu, xv, normal, jac, d_jac);

    tensor<double, 4> d_curva(eNN, numDOFs1, pDim, pDim);
    SurfLagrParam::d_Curva(d_curva, DDbf, xuu, xuv, xvv, normal, d_normal);

    tensor<double, 2> d_meancurva = product(imetric, d_curva, {{0, 2},{1, 3}}) + product(curva, d_imetric, {{0, 2},{1, 3}});



    //------------------------------------------------------------------
    // [4] second derivatives

    tensor<double, 6> dd_metric(eNN, numDOFs1, eNN, numDOFs1, pDim, pDim);
    SurfLagrParam::dd_Metric(dd_metric, Dbf);

    tensor<double, 6> dd_metric_CC = product(product(dd_metric, imetric, {{4, 0}}), imetric, {{4, 0}});

    tensor<double, 6> dd_imetric = -1.0 * dd_metric_CC -product(product(d_metric, d_imetric, {{2, 2}}).transpose({0, 1, 3, 4, 5, 2}),imetric, {{5, 0}}) -product(product(d_metric, imetric, {{2, 0}}), d_imetric, {{2, 2}}).transpose({0, 1, 3, 4, 2, 5});

    tensor<double, 4> dd_jac(eNN, numDOFs1, eNN, numDOFs1);
    SurfLagrParam::dd_Jac(dd_jac, Dbf, xu, xv, normal, d_normal);

    tensor<double, 5> dd_normal(eNN, numDOFs1, eNN, numDOFs1, nDim);
    SurfLagrParam::dd_Normal(dd_normal, Dbf, jac, d_jac, dd_jac, normal, d_normal);

    tensor<double, 6> dd_curva(eNN, numDOFs1, eNN, numDOFs1, pDim, pDim);
    SurfLagrParam::dd_Curva(dd_curva, DDbf, xuu, xuv, xvv, d_normal, dd_normal);


   // tensor<double, 4> dd_meancurva = product(imetric, dd_curva, {{0, 4},{1, 5}}) + product(d_imetric, d_curva, {{2, 2},{3, 3}}) +product(d_curva, d_imetric, {{2, 2},{3, 3}}) + product(curva, dd_imetric, {{0, 4},{1, 5}});



    // SHELL ANALYSIS

    //dphi0



    //cout<<"jacobian "<< jac<< " height: "<< height0 << " jacR: " <<jacR <<endl;
 // cout<<"print: "<< "height: "<<height0<<endl;

    double *auxiliary_b = &fillStr.paramStr->b_aux[(subFill.loc_elemID * subFill.cubaInfo.iPts + subFill.kPt) * 18];

    double G11PG = auxiliary_b[0];
    double G12PG = auxiliary_b[1];
    double G22PG = auxiliary_b[2];
    double G11NG = auxiliary_b[3];
    double G12NG = auxiliary_b[4];
    double G22NG = auxiliary_b[5];
    double GL11 = auxiliary_b[6];
    double GL22 = auxiliary_b[7];

    double GL33 = auxiliary_b[8];
    double GL44 = auxiliary_b[9];
    double GL55 = auxiliary_b[10];
    double GL66 = auxiliary_b[11];

    double GL12 = auxiliary_b[12];
    double GL34 = auxiliary_b[13];
    double GL56 = auxiliary_b[14];

  //density
double *auxiliary_c = &fillStr.paramStr->c_aux[(subFill.loc_elemID * subFill.cubaInfo.iPts + subFill.kPt) * 9];

double rho_api=auxiliary_c[0];
double rho_bas=auxiliary_c[1];
double rho_lat1=auxiliary_c[2];
double rho_lat2=auxiliary_c[3];
double rho_lat3=auxiliary_c[4];

 double xx_target=auxiliary_c[5];
double yy_target=auxiliary_c[6];
double zz_target=auxiliary_c[7];
double jac_target=auxiliary_c[8];

     tensor<double, 1> x_target(3);
  x_target(0)=xx_target;
  x_target(1)=yy_target;
  x_target(2)=zz_target;







    tensor<double, 2> GP(pDim, pDim);
    tensor<double, 2> GN(pDim, pDim);
    tensor<double, 2> GL1(pDim, pDim);
    tensor<double, 2> GL2(pDim, pDim);
    tensor<double, 2> GL3(pDim, pDim);

    GP(0, 0) = G11PG; // G ^AB // imetricR(0,0)
    GP(0, 1) = G12PG;
    GP(1, 0) = G12PG;
    GP(1, 1) = G22PG;
    GN(0, 0) = G11NG;
    GN(0, 1) = G12NG;
    GN(1, 0) = G12NG;
    GN(1, 1) = G22NG;

    GL1(0, 0) = GL11;
    GL1(0, 1) = GL12;
    GL1(1, 0) = GL12;
    GL1(1, 1) = GL22;

    GL2(0, 0) = GL33;
    GL2(0, 1) = GL34;
    GL2(1, 0) = GL34;
    GL2(1, 1) = GL44;


    GL3(0, 0) = GL55;
    GL3(0, 1) = GL56;
    GL3(1, 0) = GL56;
    GL3(1, 1) = GL66;


    tensor<double, 2> iGP = GP.inv(); //G_AB
    tensor<double, 2> iGN = GN.inv();

    tensor<double, 2> iGL1 = GL1.inv();
    tensor<double, 2> iGL2 = GL2.inv();
    tensor<double, 2> iGL3 = GL3.inv();


    tensor<double, 2> metric_GR(pDim, pDim); //gp ^A_b
    tensor<double, 2> metric_GPR(pDim, pDim); //gp ^A_b
    tensor<double, 2> metric_GNR(pDim, pDim);//gp ^A_b

    tensor<double, 2> metric_g(pDim, pDim); //gn ^a_b
    tensor<double, 2> metric_g0(pDim, pDim);//gn0 ^a_b
    tensor<double, 2> metric_gp(pDim, pDim); //gn ^a_b
    tensor<double, 2> metric_gp0(pDim, pDim);//gn0 ^a_b
    tensor<double, 2> metric_gn(pDim, pDim); //gn ^a_b
    tensor<double, 2> metric_gn0(pDim, pDim);//gn0 ^a_b

// PROCESS h
    tensor<double,2> H_nodes0(eNN, 3);//dof:({"offset","c1","c2","cDens"});
    tensor<double,2> H_nodes(eNN, 3);//dof:({"offset","c1","c2","cDens"});
    tensor<double,1> lam_nodes(eNN);//dof:({"offset","c1","c2","cDens"});
    tensor<double,2> NR_node(eNN, 3);//Auxdof: "nx", "ny", "nz"

    NR_node = nborAuxF(all,range(7,9));//Auxdof: "nx", "ny", "nz"

    H_nodes0 = nborDOFs0(all,range(3,5));//dof:({"offset","c1","c2","cDens"});
    H_nodes = nborDOFs(all,range(3,5));//dof:({"offset","c1","c2","cDens"});
    lam_nodes = nborDOFs(all,6);//dof:({"offset","c1","c2","cDens"});
    double lam=bf*lam_nodes;

    tensor<double, 1> H0(3);// 3 by1
    tensor<double, 1> H (3);// 3 by 1
    tensor<double, 1> HR (3);// 3 by 1
    tensor<double, 1> NR (3);// 3 by 1

    HR=height_in*normalR;
    H0 = bf * H_nodes0;// 3 by1
    H = bf * H_nodes;// 3 by 1
    NR = bf * NR_node;// 3 by 1

    tensor<double, 2> D_H0(nDim, pDim);// 3by 2
    tensor<double, 2> D_H(nDim, pDim);// 3 by 2
    tensor<double, 2> D_HR(nDim, pDim);// 3 by 2


    D_H0=product(Dbf, H_nodes0, {{0,0}}).transpose({1, 0});
    D_H=product(Dbf, H_nodes, {{0,0}}).transpose({1, 0});
    D_HR=product(Dbf, height_in*NR_node, {{0,0}}).transpose({1, 0});
    tensor<double, 1> Hu(nDim);
    tensor<double, 1> Hv(nDim);
    Hu=D_H(all,0);
    Hv=D_H(all,1);

   // cout<<"print: "<<NR-normalR   <<endl;
    // STRETCH CALCULATION  AND EVOLVE HEIGHT
    double jacC = jac/jacR; //
    double jacC_n = jac_n/jacR; //
    double height0 =product(H0, normal_n, {{0, 0}});//height_in/jac_n

    double mu=0.5*young/(1+nu) ;
   //
/*
    double theta1 = std::atan2(xR(1), xR(0));

   if (theta1 < 0.0)
    {
        theta1 += 2.0 * 3.14159;
    }

     mu = mu * (1.0 + 0.02 * std::cos(5.0 * theta1) / (2.0 * 3.14159));*/

 //

  double lambda=mu ;

    tensor<double, 2> DphiR(nDim, 2);
    DphiR(all, 0) = xRu(all);
    DphiR(all, 1) = xRv(all);
    tensor<double, 2> Dphi(nDim, 2);
    Dphi(all, 0) = xu(all);
    Dphi(all, 1) = xv(all);
    tensor<double, 2> Dphi0(nDim, 2);
    Dphi0(all, 0) = xu_n(all);
    Dphi0(all, 1) = xv_n(all);

 tensor<double, 2> Id2(pDim, pDim); //C_ab
    Id2(0,0)=1.0;
    Id2(0,1)=0.0;
    Id2(1,0)=0.0;
    Id2(1,1)=1.0;

    metric_GR = metricR; //Cp_ab
    metric_GPR = metricR + 0.5*aap* product(D_HR,DphiR,{{0,0}})+aap* 0.5* product(DphiR,D_HR,{{0,0}})+0.25*aap*aap*product(D_HR,D_HR,{{0,0}}) ;
    metric_GNR = metricR - 0.5* bbn*product(D_HR,DphiR,{{0,0}})- bbn*0.5* product(DphiR,D_HR,{{0,0}})+0.25*bbn*bbn*product(D_HR,D_HR,{{0,0}}) ; //Cp_ab

    metric_g = metric; //Cp_ab
    metric_g0 = metric_n; //Cp_ab

    metric_gp =metric + 0.5*aap* product(D_H,Dphi,{{0,0}})+ 0.5*aap* product(Dphi,D_H,{{0,0}})+0.25*aap*aap*product(D_H,D_H,{{0,0}}) ; //Cp_ab
    metric_gn = metric - 0.5* bbn*product(D_H,Dphi,{{0,0}})- 0.5*bbn* product(Dphi,D_H,{{0,0}})+0.25*bbn*bbn*product(D_H,D_H,{{0,0}}) ; //Cp_ab

    metric_gp0 = metric_n + 0.5* aap*product(D_H0,Dphi0,{{0,0}})+ 0.5*aap* product(Dphi0,D_H0,{{0,0}})+0.25*aap*aap*product(D_H0,D_H0,{{0,0}}) ; //Cp_ab
    metric_gn0 =  metric_n - 0.5*bbn* product(D_H0,Dphi0,{{0,0}})- 0.5*bbn* product(Dphi0,D_H0,{{0,0}})+0.25*bbn*bbn*product(D_H0,D_H0,{{0,0}}); //Cp_ab


    tensor<double,2> imetric_GR(pDim,pDim); //C+^ab
    tensor<double,2> imetric_GPR(pDim,pDim); //C+^ab
    tensor<double,2> imetric_GNR(pDim,pDim); //C+^ab

    tensor<double,2> imetric_g(pDim,pDim);
    tensor<double,2> imetric_g0(pDim,pDim);
    tensor<double,2> imetric_gp(pDim,pDim);
    tensor<double,2> imetric_gp0(pDim,pDim);
    tensor<double,2> imetric_gn(pDim,pDim);
    tensor<double,2> imetric_gn0(pDim,pDim);


    imetric_GR=metric_GR.inv();
    imetric_GPR=metric_GPR.inv();
    imetric_GNR=metric_GNR.inv();

    imetric_g=metric_g.inv();
    imetric_g0=metric_g0.inv();

    imetric_gp=metric_gp.inv();
    imetric_gp0=metric_gp0.inv();

    imetric_gn=metric_gn.inv();
    imetric_gn0=metric_gn0.inv();

    double jac_GR     = sqrt(metric_GR.det()); // J plus // sqrt(det gP)
    double jac_GPR       = sqrt(metric_GPR.det()); // J plus // sqrt(det gP)
    double jac_GNR       = sqrt(metric_GNR.det());


    double jac_g       = sqrt(metric_g.det()); // J plus // sqrt(det gP)
    double jac_g0       = sqrt(metric_g0.det()); // J plus // sqrt(det gP)

    double jac_gp       = sqrt(metric_gp.det()); // J plus // sqrt(det gP)
    double jac_gp0       = sqrt(metric_gp0.det()); // J plus // sqrt(det gP)


    double jac_gn       = sqrt(metric_gn.det()); // J plus // sqrt(det gP)
    double jac_gn0       = sqrt(metric_gn0.det()); // J plus // sqrt(det gP)


    //REAL EVOLVING
    tensor<double,2> Id33(3,3);
    Id33(0,0)=1.0;
    Id33(0,1)=0.0;
    Id33(0,2)=0.0;

    Id33(1,0)=0.0;
    Id33(1,1)=1.0;
    Id33(1,2)=0.0;

    Id33(2,0)=0.0;
    Id33(2,1)=0.0;
    Id33(2,2)=1.0;

    tensor<double,4> d_metric_g(eNN,numDOFs1,pDim,pDim);//X_I_alpa  ab
    tensor<double,4> d_imetric_g(eNN,numDOFs1,pDim,pDim);
    tensor<double,6> dd_metric_g(eNN,numDOFs1,eNN,numDOFs1,pDim,pDim);//X_J_beta X_I_alpa   ab

    tensor<double,4> d_metric_gn(eNN,numDOFs1,pDim,pDim);
    tensor<double,4> dh_metric_gn(eNN,numDOFs1,pDim,pDim); //h_I_alpa  ab

    tensor<double,6> dd_metric_gn(eNN,numDOFs1,eNN,numDOFs1,pDim,pDim);
    tensor<double,6> dh_d_metric_gn(eNN,numDOFs1,eNN,numDOFs1,pDim,pDim);
    tensor<double,6> d_dh_metric_gn(eNN,numDOFs1,eNN,numDOFs1,pDim,pDim);
    tensor<double,6> dh_dh_metric_gn(eNN,numDOFs1,eNN,numDOFs1,pDim,pDim);

    tensor<double,4> d_imetric_gn(eNN,numDOFs1,pDim,pDim);
    tensor<double,4> dh_imetric_gn(eNN,numDOFs1,pDim,pDim);

    tensor<double,4> dh_metric_g(eNN,numDOFs1,pDim,pDim);

    tensor<double,4> d_metric_gp(eNN,numDOFs1,pDim,pDim);
    tensor<double,4> dh_metric_gp(eNN,numDOFs1,pDim,pDim);


    tensor<double,6> dd_metric_gp(eNN,numDOFs1,eNN,numDOFs1,pDim,pDim);
    tensor<double,6> dh_d_metric_gp(eNN,numDOFs1,eNN,numDOFs1,pDim,pDim);
    tensor<double,6> d_dh_metric_gp(eNN,numDOFs1,eNN,numDOFs1,pDim,pDim);
    tensor<double,6> dh_dh_metric_gp(eNN,numDOFs1,eNN,numDOFs1,pDim,pDim);

    tensor<double,4> d_imetric_gp(eNN,numDOFs1,pDim,pDim);
    tensor<double,4> dh_imetric_gp(eNN,numDOFs1,pDim,pDim);

    dh_metric_g=0.0*d_metric;
    d_metric_g=d_metric;
    dd_metric_g=dd_metric;
    //metric_gp =metric + 0.5* product(D_H,Dphi,{{0,0}})+ 0.5* product(Dphi,D_H,{{0,0}})+0.25*product(D_H,D_H,{{0,0}}) ; //Cp_ab

    d_metric_gp= d_metric + 0.5*aap* outer(D_H,Dbf).transpose({2,0,1, 3})+ 0.5*aap* outer(D_H,Dbf).transpose({2,0,3, 1});//correct
    d_metric_gn = d_metric- 0.5*bbn* outer(D_H,Dbf).transpose({2,0,1, 3})- 0.5*bbn* outer(D_H,Dbf).transpose({2,0,3, 1}); //correct


    dh_metric_gp=0.5*aap* outer(Dbf,Dphi).transpose({0,2,3, 1})+ 0.5* aap*outer(Dbf,Dphi).transpose({0,2,1, 3})+0.25*aap*aap*outer(Dbf,D_H).transpose({0,2,1, 3})+0.25*aap*aap*outer(D_H,Dbf).transpose({2,0,1, 3});// correct
    dh_metric_gn=-0.5*bbn* outer(Dbf,Dphi).transpose({0,2,3, 1})- 0.5*bbn* outer(Dbf,Dphi).transpose({0,2,1, 3})+0.25*bbn*bbn*outer(Dbf,D_H).transpose({0,2,1, 3})+0.25*bbn*bbn*outer(D_H,Dbf).transpose({2,0,1, 3});// correct


    //metric + 0.5* product(D_H,Dphi,{{0,0}})+ 0.5* product(Dphi,D_H,{{0,0}})+0.25*product(D_H,D_H,{{0,0}})
    dd_metric_gp = dd_metric ; //correct
    dh_d_metric_gp=aap*0.5*outer(outer(Dbf,Id33),Dbf).transpose({0,2,4,3,1, 5})+0.5*aap*outer(outer(Dbf,Id33),Dbf).transpose({0,2,4,3,5, 1});//correct

    d_dh_metric_gp=0.5*aap*outer(Dbf,outer(Dbf,Id33)).transpose({2,5,0,4,1, 3})+aap*0.5*outer(Dbf,outer(Dbf,Id33)).transpose({2,5,0,4,3, 1});//correct
    dh_dh_metric_gp=0.25*aap*aap*outer(Dbf,outer(Dbf,Id33)).transpose({2,5,0,4,1, 3})+0.25*aap*aap*outer(outer(Dbf,Id33),Dbf).transpose({0,3,4,2,1, 5});//correct

    dd_metric_gn =  dd_metric; //Cp_ab
    dh_d_metric_gn=-0.5*bbn*outer(outer(Dbf,Id33),Dbf).transpose({0,2,4,3,1, 5})-0.5*bbn*outer(outer(Dbf,Id33),Dbf).transpose({0,2,4,3,5, 1});//correct

    d_dh_metric_gn=-0.5*bbn*outer(Dbf,outer(Dbf,Id33)).transpose({2,5,0,4,1, 3})-0.5*bbn*outer(Dbf,outer(Dbf,Id33)).transpose({2,5,0,4,3, 1});//correct
    dh_dh_metric_gn=0.25*bbn*bbn*outer(Dbf,outer(Dbf,Id33)).transpose({2,5,0,4,1, 3})+0.25*bbn*bbn*outer(outer(Dbf,Id33),Dbf).transpose({0,3,4,2,1, 5});//correct



    //




    tensor<double,4> d_metric_gp_CC  = product(product(imetric_gp,d_metric_gp,{{1,2}}).transpose({1,2,0,3}),imetric_gp,{{3,0}});
    d_imetric_gp = -1.0 * d_metric_gp_CC;

    tensor<double,4> d_metric_gn_CC  = product(product(imetric_gn,d_metric_gn,{{1,2}}).transpose({1,2,0,3}),imetric_gn,{{3,0}});
    d_imetric_gn = -1.0 * d_metric_gn_CC;

    tensor<double,4> d_metric_g_CC  = product(product(imetric_g,d_metric_g,{{1,2}}).transpose({1,2,0,3}),imetric_g,{{3,0}});
    d_imetric_g = -1.0 * d_metric_g_CC;

    tensor<double,4> dh_metric_gp_CC  = product(product(imetric_gp,dh_metric_gp,{{1,2}}).transpose({1,2,0,3}),imetric_gp,{{3,0}});
    dh_imetric_gp = -1.0 * dh_metric_gp_CC;

    tensor<double,4> dh_metric_gn_CC  = product(product(imetric_gn,dh_metric_gn,{{1,2}}).transpose({1,2,0,3}),imetric_gn,{{3,0}});
    dh_imetric_gn = -1.0 * dh_metric_gn_CC;

    //test
   // cout<<"print: "<<d_imetric_g-d_imetric<<endl;

// JACOBIAN AND DERIVATIVES
    //REAL EVOLVING

    tensor<double,2> d_jac_gp(eNN,numDOFs1);
    d_jac_gp=0.5*jac_gp* product(imetric_gp,d_metric_gp,{{0,2},{1,3}});// Good P by _gp, N by _gn, C by _g
    tensor<double,2> d_jac_gn(eNN,numDOFs1);
    d_jac_gn=0.5*jac_gn* product(imetric_gn,d_metric_gn,{{0,2},{1,3}});// Good
    tensor<double,2> dh_jac_gp(eNN,numDOFs1);

    dh_jac_gp=0.5*jac_gp* product(imetric_gp,dh_metric_gp,{{0,2},{1,3}});// Good P by _gp, N by _gn, C by _g
    tensor<double,2> dh_jac_gn(eNN,numDOFs1);
    dh_jac_gn=0.5*jac_gn* product(imetric_gn,dh_metric_gn,{{0,2},{1,3}});// Good

    //double jac

    //dd_jac
    tensor<double,4> dd_jac_gp(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,2> dd_jac_gp1(eNN,numDOFs1);
    dd_jac_gp1=product(imetric_gp,d_metric_gp,{{0,2},{1,3}});// dgP/dx:gP^-1
    dd_jac_gp=0.5* jac_gp*product(imetric_gp,dd_metric_gp,{{0,4},{1,5}})+ 0.5*jac_gp*product(d_imetric_gp,d_metric_gp,{{2,2},{3,3}})+0.5*outer(d_jac_gp,dd_jac_gp1);// good correct


    tensor<double,4> dd_jac_gn(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,2> dd_jac_gn1(eNN,numDOFs1);
    dd_jac_gn1=product(imetric_gn,d_metric_gn,{{0,2},{1,3}});// dgP/dx:gP^-1
    dd_jac_gn=0.5* jac_gn*product(imetric_gn,dd_metric_gn,{{0,4},{1,5}})+ 0.5*jac_gn*product(d_imetric_gn,d_metric_gn,{{2,2},{3,3}})+0.5*outer(d_jac_gn,dd_jac_gn1);// good correct

//dh_d_jac
    tensor<double,4> dh_d_jac_gp(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,2>dh_d_jac_gp1(eNN,numDOFs1);
    dh_d_jac_gp1=product(imetric_gp,d_metric_gp,{{0,2},{1,3}});// dgP/dx:gP^-1
    dh_d_jac_gp=0.5* jac_gp*product(imetric_gp,dh_d_metric_gp,{{0,4},{1,5}})+ 0.5*jac_gp*product(dh_imetric_gp,d_metric_gp,{{2,2},{3,3}})+0.5*outer(dh_jac_gp,dh_d_jac_gp1);// good correct

    tensor<double,4> dh_d_jac_gn(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,2>dh_d_jac_gn1(eNN,numDOFs1);
    dh_d_jac_gn1=product(imetric_gn,d_metric_gn,{{0,2},{1,3}});// dgP/dx:gP^-1
    dh_d_jac_gn=0.5* jac_gn*product(imetric_gn,dh_d_metric_gn,{{0,4},{1,5}})+ 0.5*jac_gn*product(dh_imetric_gn,d_metric_gn,{{2,2},{3,3}})+0.5*outer(dh_jac_gn,dh_d_jac_gn1);// good correct

//d_dh_jac
    tensor<double,4> d_dh_jac_gp(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,2>d_dh_jac_gp1(eNN,numDOFs1);
    d_dh_jac_gp1=product(imetric_gp,dh_metric_gp,{{0,2},{1,3}});// dgP/dx:gP^-1
    d_dh_jac_gp=0.5* jac_gp*product(imetric_gp,d_dh_metric_gp,{{0,4},{1,5}})+ 0.5*jac_gp*product(d_imetric_gp,dh_metric_gp,{{2,2},{3,3}})+0.5*outer(d_jac_gp,d_dh_jac_gp1);// good correct


    tensor<double,4> d_dh_jac_gn(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,2>d_dh_jac_gn1(eNN,numDOFs1);
    d_dh_jac_gn1=product(imetric_gn,dh_metric_gn,{{0,2},{1,3}});// dgP/dx:gP^-1
    d_dh_jac_gn=0.5* jac_gn*product(imetric_gn,d_dh_metric_gn,{{0,4},{1,5}})+ 0.5*jac_gn*product(d_imetric_gn,dh_metric_gn,{{2,2},{3,3}})+0.5*outer(d_jac_gn,d_dh_jac_gn1);// good correct

//dh_dh_jac

    tensor<double,4> dh_dh_jac_gp(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,2>dh_dh_jac_gp1(eNN,numDOFs1);
    dh_dh_jac_gp1=product(imetric_gp,dh_metric_gp,{{0,2},{1,3}});// dgP/dx:gP^-1
    dh_dh_jac_gp=0.5* jac_gp*product(imetric_gp,dh_dh_metric_gp,{{0,4},{1,5}})+ 0.5*jac_gp*product(dh_imetric_gp,dh_metric_gp,{{2,2},{3,3}})+0.5*outer(dh_jac_gp,dh_dh_jac_gp1);// good correct


    tensor<double,4> dh_dh_jac_gn(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,2>dh_dh_jac_gn1(eNN,numDOFs1);
    dh_dh_jac_gn1=product(imetric_gn,dh_metric_gn,{{0,2},{1,3}});// dgP/dx:gP^-1
    dh_dh_jac_gn=0.5* jac_gn*product(imetric_gn,dh_dh_metric_gn,{{0,4},{1,5}})+ 0.5*jac_gn*product(dh_imetric_gn,dh_metric_gn,{{2,2},{3,3}})+0.5*outer(dh_jac_gn,dh_dh_jac_gn1);// good correct

    //calculation of material stretch and lateral strain tensors
    tensor<double, 1> G1(nDim);
    G1 = xRu(all);
    tensor<double, 1> G2(nDim);
    G2 = xRv(all);


tensor<double,1> M1(pDim);
tensor<double,1> M2(pDim);
tensor<double,1> M3(pDim);


    tensor<double, 1> M11(3); //C_ab
    M11(0)=1.0;
    M11(1)=0.0;
    M11(2)=0.0;
    tensor<double, 1> M22(3); //C_ab
    M22(0)=-0.5;
    M22(1)=sqrt(3)/2;
    M22(2)=0.0;

    tensor<double, 1> M33(3); //C_ab
    M33(0)=-0.5;
    M33(1)=-sqrt(3)/2;
    M33(2)=0.0;

    double MG1_1=G1*M11;
    double MG2_1=G2*M11;
    double MG1_2=G1*M22;
    double MG2_2=G2*M22;
    double MG1_3=G1*M33;
    double MG2_3=G2*M33;

    tensor<double, 1> MF1(2); //C_ab
    MF1(0)=MG1_1;
    MF1(1)=MG2_1;

    tensor<double, 1> MF2(2); //C_ab
    MF2(0)=MG1_2;
    MF2(1)=MG2_2;
    tensor<double, 1> MF3(2); //C_ab
    MF3(0)=MG1_3;
    MF3(1)=MG2_3;


//  components in the basis of G

    M1=imetricR*MF1;
    M2=imetricR*MF2;
    M3=imetricR*MF3;



//new done
    tensor<double, 1> mm1(nDim); //C_ab
    tensor<double, 1> mm2(nDim); //C_ab
    tensor<double, 1> mm3(nDim); //C_ab

    mm1=M1(0)*xu+M1(1)*xv;
    mm2=M2(0)*xu+M2(1)*xv;
    mm3=M3(0)*xu+M3(1)*xv;


    tensor<double, 1> mm1_n(nDim); //C_ab
    tensor<double, 1> mm2_n(nDim); //C_ab
    tensor<double, 1> mm3_n(nDim); //C_ab
    mm1_n=M1(0)*xu_n+M1(1)*xv_n;
    mm2_n=M2(0)*xu_n+M2(1)*xv_n;
    mm3_n=M3(0)*xu_n+M3(1)*xv_n;

    tensor<double,3> d_mm1(eNN,numDOFs1,nDim);
    tensor<double,3> d_mm2(eNN,numDOFs1,nDim);
    tensor<double,3> d_mm3(eNN,numDOFs1,nDim);

    d_mm1=M1(0)*outer(Dbf(all,0),Id33)+M1(1)*outer(Dbf(all,1),Id33);
    d_mm2=M2(0)*outer(Dbf(all,0),Id33)+M2(1)*outer(Dbf(all,1),Id33);
    d_mm3=M3(0)*outer(Dbf(all,0),Id33)+M3(1)*outer(Dbf(all,1),Id33);


    double shear1=mm1*H/height_in;
    double shear2=mm2*H/height_in;
    double shear3=mm3*H/height_in;

    double shear1_n=mm1_n*H0/height_in;
    double shear2_n=mm2_n*H0/height_in;
    double shear3_n=mm3_n*H0/height_in;


    double CM1= product(M1,product(metric,M1,{{1,0}}),{{0,0}});
    double CM2= product(M2,product(metric,M2,{{1,0}}),{{0,0}});
    double CM3= product(M3,product(metric,M3,{{1,0}}),{{0,0}});

    double CM1_n= product(M1,product(metric_n,M1,{{1,0}}),{{0,0}});
    double CM2_n= product(M2,product(metric_n,M2,{{1,0}}),{{0,0}});
    double CM3_n= product(M3,product(metric_n,M3,{{1,0}}),{{0,0}});


    tensor<double,2> d_CM1(eNN,numDOFs1);
    tensor<double,2> d_CM2(eNN,numDOFs1);
    tensor<double,2> d_CM3(eNN,numDOFs1);
    d_CM1=product(M1,product(d_metric,M1,{{3,0}}),{{0,2}});
    d_CM2=product(M2,product(d_metric,M2,{{3,0}}),{{0,2}});
    d_CM3=product(M3,product(d_metric,M3,{{3,0}}),{{0,2}});

    tensor<double,4> dd_CM1(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dd_CM2(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dd_CM3(eNN,numDOFs1,eNN,numDOFs1);
    dd_CM1=product(M1,product(dd_metric,M1,{{5,0}}),{{0,4}});
    dd_CM2=product(M2,product(dd_metric,M2,{{5,0}}),{{0,4}});
    dd_CM3=product(M3,product(dd_metric,M3,{{5,0}}),{{0,4}});


    double CL22=product(H,H,{{0,0}})/(height_in*height_in);
    double CL22_n=product(H0,H0,{{0,0}})/(height_in*height_in);

    tensor<double,2> d_CL22(eNN,numDOFs1);
    d_CL22=0.0*outer(bf,H);
    tensor<double,2> dh_CL22(eNN,numDOFs1);
    dh_CL22=2.0*outer(bf,H)/(height_in*height_in);



    tensor<double,4> dd_CL22(eNN,numDOFs1,eNN,numDOFs1);
    dd_CL22=0.0*outer(outer(bf,bf),Id33).transpose({0,2,1,3});

    tensor<double,4> dh_d_CL22(eNN,numDOFs1,eNN,numDOFs1);
    dh_d_CL22=0.0*outer(outer(bf,bf),Id33).transpose({0,2,1,3});


    tensor<double,4> d_dh_CL22(eNN,numDOFs1,eNN,numDOFs1);
    d_dh_CL22=0.0*outer(outer(bf,bf),Id33).transpose({0,2,1,3});

    tensor<double,4> dh_dh_CL22(eNN,numDOFs1,eNN,numDOFs1);
    dh_dh_CL22=2.0*outer(outer(bf,bf),Id33).transpose({0,2,1,3})/(height_in*height_in);


//
    tensor<double,2> d_shear1(eNN,numDOFs1);
    d_shear1=d_mm1*H/height_in;
    tensor<double,2> d_shear2(eNN,numDOFs1);
    d_shear2=d_mm2*H/height_in;
    tensor<double,2> d_shear3(eNN,numDOFs1);
    d_shear3=d_mm3*H/height_in;

    tensor<double,2> dh_shear1(eNN,numDOFs1);
    dh_shear1=outer(bf,mm1)/height_in;
    tensor<double,2> dh_shear2(eNN,numDOFs1);
    dh_shear2=outer(bf,mm2)/height_in;
    tensor<double,2> dh_shear3(eNN,numDOFs1);
    dh_shear3=outer(bf,mm3)/height_in;

    tensor<double,4> dh_d_shear1(eNN,numDOFs1,eNN,numDOFs1);
    dh_d_shear1=outer(bf,d_mm1).transpose({0,3,1,2})/height_in;
    tensor<double,4> dh_d_shear2(eNN,numDOFs1,eNN,numDOFs1);
    dh_d_shear2=outer(bf,d_mm2).transpose({0,3,1,2})/height_in;
    tensor<double,4> dh_d_shear3(eNN,numDOFs1,eNN,numDOFs1);
    dh_d_shear3=outer(bf,d_mm3).transpose({0,3,1,2})/height_in;

    tensor<double,4> d_dh_shear1(eNN,numDOFs1,eNN,numDOFs1);
    d_dh_shear1=outer(bf,d_mm1).transpose({1,2,0,3})/height_in;
    tensor<double,4> d_dh_shear2(eNN,numDOFs1,eNN,numDOFs1);
    d_dh_shear2=outer(bf,d_mm2).transpose({1,2,0,3})/height_in;
    tensor<double,4> d_dh_shear3(eNN,numDOFs1,eNN,numDOFs1);
    d_dh_shear3=outer(bf,d_mm3).transpose({1,2,0,3})/height_in;

    tensor<double,2> CL1(pDim,pDim);
    CL1(0,0)=CM1;
    CL1(0,1)=shear1;
    CL1(1,0)=shear1;
    CL1(1,1)=CL22;

    double jacCL1=sqrt(CM1*CL22-shear1*shear1);
    double jacCL1_n=sqrt(CM1_n*CL22_n-shear1_n*shear1_n);

    tensor<double,2> iCL1=CL1.inv();


    tensor<double,2> CL2(pDim,pDim);
    CL2(0,0)=CM2;
    CL2(0,1)=shear2;
    CL2(1,0)=shear2;
    CL2(1,1)=CL22;

    double jacCL2=sqrt(CM2*CL22-shear2*shear2);
    double jacCL2_n=sqrt(CM2_n*CL22_n-shear2_n*shear2_n);


    tensor<double,2> iCL2=CL2.inv();

    tensor<double,2> CL3(pDim,pDim);
    CL3(0,0)=CM3;
    CL3(0,1)=shear3;
    CL3(1,0)=shear3;
    CL3(1,1)=CL22;

    double jacCL3=sqrt(CM3*CL22-shear3*shear3);
    double jacCL3_n=sqrt(CM3_n*CL22_n-shear3_n*shear3_n);

    tensor<double,2> iCL3=CL3.inv();

    tensor<double,2> d_jacC=d_jac/jacR;
    tensor<double,4> dd_jacC=dd_jac/jacR;
//till here updated

    tensor<double,2> d_jacCL1(eNN,numDOFs1);
    d_jacCL1= 0.5/(jacCL1)*(d_CM1*CL22-2*shear1*d_shear1);//CORRECT

    tensor<double,2> dh_jacCL1(eNN,numDOFs1);
    dh_jacCL1=0.5/(jacCL1)*( dh_CL22*CM1-2.0*shear1*dh_shear1);//correct

    tensor<double,4> dd_jacCL1(eNN,numDOFs1,eNN,numDOFs1);
    dd_jacCL1= 0.5/(jacCL1)*(dd_CM1*CL22-2*outer(d_shear1,d_shear1))-0.5/(jacCL1*jacCL1)*outer(d_jacCL1,d_CM1*CL22-2*shear1*d_shear1);//correct

    tensor<double,4> dh_d_jacCL1(eNN,numDOFs1,eNN,numDOFs1);
    dh_d_jacCL1=-0.5/(jacCL1*jacCL1)*outer(dh_jacCL1,d_CM1*CL22-2*shear1*d_shear1)+0.5/(jacCL1)*(outer(dh_CL22,d_CM1)-2*outer(dh_shear1,d_shear1)-2*shear1*dh_d_shear1);//correct
    tensor<double,4> d_dh_jacCL1(eNN,numDOFs1,eNN,numDOFs1);
    d_dh_jacCL1=-0.5/(jacCL1*jacCL1)*outer(d_jacCL1,(dh_CL22*CM1-2.0*shear1*dh_shear1))+0.5/(jacCL1)*(outer(d_CM1,dh_CL22)-2*outer(d_shear1,dh_shear1)-2*shear1*d_dh_shear1);//correct

    tensor<double,4> dh_dh_jacCL1(eNN,numDOFs1,eNN,numDOFs1);
    dh_dh_jacCL1=-0.5/(jacCL1*jacCL1)*outer(dh_jacCL1,(dh_CL22*CM1-2.0*shear1*dh_shear1))+0.5/(jacCL1)*(dh_dh_CL22*CM1-2.0*outer(dh_shear1,dh_shear1));//correct


    tensor<double,2> d_jacCL2(eNN,numDOFs1);
    d_jacCL2= 0.5/(jacCL2)*(d_CM2*CL22-2*shear2*d_shear2);//CORRECT

    tensor<double,2> dh_jacCL2(eNN,numDOFs1);
    dh_jacCL2=0.5/(jacCL2)*( dh_CL22*CM2-2.0*shear2*dh_shear2);//correct

    tensor<double,4> dd_jacCL2(eNN,numDOFs1,eNN,numDOFs1);
    dd_jacCL2= 0.5/(jacCL2)*(dd_CM2*CL22-2*outer(d_shear2,d_shear2))-0.5/(jacCL2*jacCL2)*outer(d_jacCL2,d_CM2*CL22-2*shear2*d_shear2);//correct

    tensor<double,4> dh_d_jacCL2(eNN,numDOFs1,eNN,numDOFs1);
    dh_d_jacCL2=-0.5/(jacCL2*jacCL2)*outer(dh_jacCL2,d_CM2*CL22-2*shear2*d_shear2)+0.5/(jacCL2)*(outer(dh_CL22,d_CM2)-2*outer(dh_shear2,d_shear2)-2*shear2*dh_d_shear2);//correct
    tensor<double,4> d_dh_jacCL2(eNN,numDOFs1,eNN,numDOFs1);
    d_dh_jacCL2=-0.5/(jacCL2*jacCL2)*outer(d_jacCL2,(dh_CL22*CM2-2.0*shear2*dh_shear2))+0.5/(jacCL2)*(outer(d_CM2,dh_CL22)-2*outer(d_shear2,dh_shear2)-2*shear2*d_dh_shear2);//correct

    tensor<double,4> dh_dh_jacCL2(eNN,numDOFs1,eNN,numDOFs1);
    dh_dh_jacCL2=-0.5/(jacCL2*jacCL2)*outer(dh_jacCL2,(dh_CL22*CM2-2.0*shear2*dh_shear2))+0.5/(jacCL2)*(dh_dh_CL22*CM2-2.0*outer(dh_shear2,dh_shear2));//correct


    tensor<double,2> d_jacCL3(eNN,numDOFs1);
    d_jacCL3= 0.5/(jacCL3)*(d_CM3*CL22-2*shear3*d_shear3);//CORRECT

    tensor<double,2> dh_jacCL3(eNN,numDOFs1);
    dh_jacCL3=0.5/(jacCL3)*( dh_CL22*CM3-2.0*shear3*dh_shear3);//correct

    tensor<double,4> dd_jacCL3(eNN,numDOFs1,eNN,numDOFs1);
    dd_jacCL3= 0.5/(jacCL3)*(dd_CM3*CL22-2*outer(d_shear3,d_shear3))-0.5/(jacCL3*jacCL3)*outer(d_jacCL3,d_CM3*CL22-2*shear3*d_shear3);//correct

    tensor<double,4> dh_d_jacCL3(eNN,numDOFs1,eNN,numDOFs1);
    dh_d_jacCL3=-0.5/(jacCL3*jacCL3)*outer(dh_jacCL3,d_CM3*CL22-2*shear3*d_shear3)+0.5/(jacCL3)*(outer(dh_CL22,d_CM3)-2*outer(dh_shear3,d_shear3)-2*shear3*dh_d_shear3);//correct
    tensor<double,4> d_dh_jacCL3(eNN,numDOFs1,eNN,numDOFs1);
    d_dh_jacCL3=-0.5/(jacCL3*jacCL3)*outer(d_jacCL3,(dh_CL22*CM3-2.0*shear3*dh_shear3))+0.5/(jacCL3)*(outer(d_CM3,dh_CL22)-2*outer(d_shear3,dh_shear3)-2*shear3*d_dh_shear3);//correct

    tensor<double,4> dh_dh_jacCL3(eNN,numDOFs1,eNN,numDOFs1);
    dh_dh_jacCL3=-0.5/(jacCL3*jacCL3)*outer(dh_jacCL3,(dh_CL22*CM3-2.0*shear3*dh_shear3))+0.5/(jacCL3)*(dh_dh_CL22*CM3-2.0*outer(dh_shear3,dh_shear3));//correct


//needed to be

    tensor<double,4> d_CL1(eNN,numDOFs1,pDim,pDim);
    d_CL1(all,all,0,0)=d_CM1;
    d_CL1(all,all,1,0)=d_shear1;
    d_CL1(all,all,0,1)=d_shear1;
    d_CL1(all,all,1,1)=d_jacC*0.0;

    tensor<double,4> dh_CL1(eNN,numDOFs1,pDim,pDim);
    dh_CL1(all,all,0,0)=d_jacC*0.0;
    dh_CL1(all,all,1,0)=dh_shear1;
    dh_CL1(all,all,0,1)=dh_shear1;
    dh_CL1(all,all,1,1)=dh_CL22;

    tensor<double,6> dd_CL1(eNN,numDOFs1,eNN,numDOFs1,pDim,pDim);
    dd_CL1(all,all,all,all,0,0)=dd_CM1;
    dd_CL1(all,all,all,all,0,1)=0.0*dd_jacC;
    dd_CL1(all,all,all,all,1,0)=0.0*dd_jacC;
    dd_CL1(all,all,all,all,1,1)=dd_jacC*0.0;

    tensor<double,6> dh_d_CL1(eNN,numDOFs1,eNN,numDOFs1,pDim,pDim);
    dh_d_CL1(all,all,all,all,0,0)=dd_jacC*0.0;
    dh_d_CL1(all,all,all,all,0,1)=dh_d_shear1;
    dh_d_CL1(all,all,all,all,1,0)=dh_d_shear1;
    dh_d_CL1(all,all,all,all,1,1)=dd_CL22*0.0;

    tensor<double,6> d_dh_CL1(eNN,numDOFs1,eNN,numDOFs1,pDim,pDim);
    d_dh_CL1(all,all,all,all,0,0)=dd_jacC*0.0;
    d_dh_CL1(all,all,all,all,0,1)=d_dh_shear1;
    d_dh_CL1(all,all,all,all,1,0)=d_dh_shear1;
    d_dh_CL1(all,all,all,all,1,1)=dd_CL22*0.0;

    tensor<double,6> dh_dh_CL1(eNN,numDOFs1,eNN,numDOFs1,pDim,pDim);
    dh_dh_CL1(all,all,all,all,0,0)=dd_jacC*0.0;
    dh_dh_CL1(all,all,all,all,0,1)=0.0*dd_jacC;
    dh_dh_CL1(all,all,all,all,1,0)=0.0*dd_jacC;
    dh_dh_CL1(all,all,all,all,1,1)=dh_dh_CL22;
// case2
    tensor<double,4> d_CL2(eNN,numDOFs1,pDim,pDim);
    d_CL2(all,all,0,0)=d_CM2;
    d_CL2(all,all,1,0)=d_shear2;
    d_CL2(all,all,0,1)=d_shear2;
    d_CL2(all,all,1,1)=d_jacC*0.0;

    tensor<double,4> dh_CL2(eNN,numDOFs1,pDim,pDim);
    dh_CL2(all,all,0,0)=d_jacC*0.0;
    dh_CL2(all,all,1,0)=dh_shear2;
    dh_CL2(all,all,0,1)=dh_shear2;
    dh_CL2(all,all,1,1)=dh_CL22;

    tensor<double,6> dd_CL2(eNN,numDOFs1,eNN,numDOFs1,pDim,pDim);
    dd_CL2(all,all,all,all,0,0)=dd_CM2;
    dd_CL2(all,all,all,all,0,1)=0.0*dd_jacC;
    dd_CL2(all,all,all,all,1,0)=0.0*dd_jacC;
    dd_CL2(all,all,all,all,1,1)=dd_jacC*0.0;

    tensor<double,6> dh_d_CL2(eNN,numDOFs1,eNN,numDOFs1,pDim,pDim);
    dh_d_CL2(all,all,all,all,0,0)=dd_jacC*0.0;
    dh_d_CL2(all,all,all,all,0,1)=dh_d_shear2;
    dh_d_CL2(all,all,all,all,1,0)=dh_d_shear2;
    dh_d_CL2(all,all,all,all,1,1)=dd_CL22*0.0;

    tensor<double,6> d_dh_CL2(eNN,numDOFs1,eNN,numDOFs1,pDim,pDim);
    d_dh_CL2(all,all,all,all,0,0)=dd_jacC*0.0;
    d_dh_CL2(all,all,all,all,0,1)=d_dh_shear2;
    d_dh_CL2(all,all,all,all,1,0)=d_dh_shear2;
    d_dh_CL2(all,all,all,all,1,1)=dd_CL22*0.0;

    tensor<double,6> dh_dh_CL2(eNN,numDOFs1,eNN,numDOFs1,pDim,pDim);
    dh_dh_CL2(all,all,all,all,0,0)=dd_jacC*0.0;
    dh_dh_CL2(all,all,all,all,0,1)=0.0*dd_jacC;
    dh_dh_CL2(all,all,all,all,1,0)=0.0*dd_jacC;
    dh_dh_CL2(all,all,all,all,1,1)=dh_dh_CL22;
    //case 3


    tensor<double,4> d_CL3(eNN,numDOFs1,pDim,pDim);
    d_CL3(all,all,0,0)=d_CM3;
    d_CL3(all,all,1,0)=d_shear3;
    d_CL3(all,all,0,1)=d_shear3;
    d_CL3(all,all,1,1)=d_jacC*0.0;

    tensor<double,4> dh_CL3(eNN,numDOFs1,pDim,pDim);
    dh_CL3(all,all,0,0)=d_jacC*0.0;
    dh_CL3(all,all,1,0)=dh_shear3;
    dh_CL3(all,all,0,1)=dh_shear3;
    dh_CL3(all,all,1,1)=dh_CL22;

    tensor<double,6> dd_CL3(eNN,numDOFs1,eNN,numDOFs1,pDim,pDim);
    dd_CL3(all,all,all,all,0,0)=dd_CM3;
    dd_CL3(all,all,all,all,0,1)=0.0*dd_jacC;
    dd_CL3(all,all,all,all,1,0)=0.0*dd_jacC;
    dd_CL3(all,all,all,all,1,1)=dd_jacC*0.0;

    tensor<double,6> dh_d_CL3(eNN,numDOFs1,eNN,numDOFs1,pDim,pDim);
    dh_d_CL3(all,all,all,all,0,0)=dd_jacC*0.0;
    dh_d_CL3(all,all,all,all,0,1)=dh_d_shear3;
    dh_d_CL3(all,all,all,all,1,0)=dh_d_shear3;
    dh_d_CL3(all,all,all,all,1,1)=dd_CL22*0.0;

    tensor<double,6> d_dh_CL3(eNN,numDOFs1,eNN,numDOFs1,pDim,pDim);
    d_dh_CL3(all,all,all,all,0,0)=dd_jacC*0.0;
    d_dh_CL3(all,all,all,all,0,1)=d_dh_shear3;
    d_dh_CL3(all,all,all,all,1,0)=d_dh_shear3;
    d_dh_CL3(all,all,all,all,1,1)=dd_CL22*0.0;

    tensor<double,6> dh_dh_CL3(eNN,numDOFs1,eNN,numDOFs1,pDim,pDim);
    dh_dh_CL3(all,all,all,all,0,0)=dd_jacC*0.0;
    dh_dh_CL3(all,all,all,all,0,1)=0.0*dd_jacC;
    dh_dh_CL3(all,all,all,all,1,0)=0.0*dd_jacC;
    dh_dh_CL3(all,all,all,all,1,1)=dh_dh_CL22;



    //end new good





    // jacobian of viscoelastic metric

    double jacGP       = (GP.det()); // J of GPP (G11P*G22P-G12P*G12P)
    double jacGN       = (GN.det()); // J of GPP (G11P*G22P-G12P*G12P)

    double jacGL1      =  (GL1.det()); // J of GPP (G11P*G22P-G12P*G12P)
    double jacGL2      =  (GL2.det()); // J of GPP (G11P*G22P-G12P*G12P)
    double jacGL3      =  (GL3.det()); // J of GPP (G11P*G22P-G12P*G12P)



//first and third invariants of visco elastic tensor
    double I1_gp=product(GP,metric_gp ,{{0,0},{1,1}}) ;
    double I1_gn=product(GN,metric_gn ,{{0,0},{1,1}}) ;

    double I1L1=product(GL1,CL1 ,{{0,0},{1,1}}) ;
    double I1L2=product(GL2,CL2 ,{{0,0},{1,1}}) ;
    double I1L3=product(GL3,CL3 ,{{0,0},{1,1}}) ;

//

    double I3_gp=jac_gp*jac_gp* jacGP;
    double I3_gn=jac_gn*jac_gn*jacGN;


    double I3L1=jacCL1*jacCL1*jacGL1;
    double I3L2=jacCL2*jacCL2*jacGL2;
    double I3L3=jacCL3*jacCL3*jacGL3;



    double I3_gp_SR=sqrt(jacGP)*jac_gp;
    double I3_gn_SR=sqrt(jacGN)*jac_gn;

    double I3L1_SR=sqrt(jacGL1)*jacCL1;
    double I3L2_SR=sqrt(jacGL2)*jacCL2;
    double I3L3_SR=sqrt(jacGL3)*jacCL3;

//d_I1, d_I3
    tensor<double,2> d_I1_gp(eNN,numDOFs1);
    tensor<double,2> d_I1_gn(eNN,numDOFs1);
    tensor<double,2> d_I3_gp(eNN,numDOFs1);
    tensor<double,2> d_I3_gn(eNN,numDOFs1);

    tensor<double,2> d_I3_gp_SR(eNN,numDOFs1);
    tensor<double,2> d_I3_gn_SR(eNN,numDOFs1);

    tensor<double,2> d_I1L1(eNN,numDOFs1);
    tensor<double,2> d_I1L2(eNN,numDOFs1);
    tensor<double,2> d_I1L3(eNN,numDOFs1);

    tensor<double,2> d_I3L1(eNN,numDOFs1);
    tensor<double,2> d_I3L2(eNN,numDOFs1);
    tensor<double,2> d_I3L3(eNN,numDOFs1);

    tensor<double,2> d_I3L1_SR(eNN,numDOFs1);
    tensor<double,2> d_I3L2_SR(eNN,numDOFs1);
    tensor<double,2> d_I3L3_SR(eNN,numDOFs1);


    d_I1_gp=product(GP,d_metric_gp ,{{0,2},{1,3}}) ;
    d_I1_gn=product(GN,d_metric_gn ,{{0,2},{1,3}}) ;

    d_I1L1=product(GL1,d_CL1 ,{{0,2},{1,3}}) ;
    d_I1L2=product(GL2,d_CL2 ,{{0,2},{1,3}}) ;
    d_I1L3=product(GL3,d_CL3 ,{{0,2},{1,3}}) ;



    d_I3_gp=2*jac_gp*d_jac_gp* jacGP;
    d_I3_gn=2*jac_gn*d_jac_gn* jacGN;

    d_I3L1=2*jacCL1*d_jacCL1* jacGL1;
    d_I3L2=2*jacCL2*d_jacCL2* jacGL2;
    d_I3L3=2*jacCL3*d_jacCL3* jacGL3;



    d_I3_gp_SR=sqrt(jacGP)*d_jac_gp;
    d_I3_gn_SR=sqrt(jacGN)*d_jac_gn;

    d_I3L1_SR=sqrt(jacGL1)*d_jacCL1;
    d_I3L2_SR=sqrt(jacGL2)*d_jacCL2;
    d_I3L3_SR=sqrt(jacGL3)*d_jacCL3;

    //dh_I1, dh_I3
    tensor<double,2> dh_I1_gp(eNN,numDOFs1);
    tensor<double,2> dh_I1_gn(eNN,numDOFs1);
    tensor<double,2> dh_I3_gp(eNN,numDOFs1);
    tensor<double,2> dh_I3_gn(eNN,numDOFs1);

    tensor<double,2> dh_I3_gp_SR(eNN,numDOFs1);
    tensor<double,2> dh_I3_gn_SR(eNN,numDOFs1);

    tensor<double,2> dh_I1L1(eNN,numDOFs1);
    tensor<double,2> dh_I1L2(eNN,numDOFs1);
    tensor<double,2> dh_I1L3(eNN,numDOFs1);

    tensor<double,2> dh_I3L1(eNN,numDOFs1);
    tensor<double,2> dh_I3L2(eNN,numDOFs1);
    tensor<double,2> dh_I3L3(eNN,numDOFs1);

    tensor<double,2> dh_I3L1_SR(eNN,numDOFs1);
    tensor<double,2> dh_I3L2_SR(eNN,numDOFs1);
    tensor<double,2> dh_I3L3_SR(eNN,numDOFs1);


    dh_I1_gp=product(GP,dh_metric_gp ,{{0,2},{1,3}}) ;
    dh_I1_gn=product(GN,dh_metric_gn ,{{0,2},{1,3}}) ;

    dh_I1L1=product(GL1,dh_CL1 ,{{0,2},{1,3}}) ;
    dh_I1L2=product(GL2,dh_CL2 ,{{0,2},{1,3}}) ;
    dh_I1L3=product(GL3,dh_CL3 ,{{0,2},{1,3}}) ;



    dh_I3_gp=2*jac_gp*dh_jac_gp* jacGP;
    dh_I3_gn=2*jac_gn*dh_jac_gn* jacGN;

    dh_I3L1=2*jacCL1*dh_jacCL1* jacGL1;
    dh_I3L2=2*jacCL2*dh_jacCL2* jacGL2;
    dh_I3L3=2*jacCL3*dh_jacCL3* jacGL3;



    dh_I3_gp_SR=sqrt(jacGP)*dh_jac_gp;
    dh_I3_gn_SR=sqrt(jacGN)*dh_jac_gn;

    dh_I3L1_SR=sqrt(jacGL1)*dh_jacCL1;
    dh_I3L2_SR=sqrt(jacGL2)*dh_jacCL2;
    dh_I3L3_SR=sqrt(jacGL3)*dh_jacCL3;



//dd_I1,dd_I3
    tensor<double,4> dd_I1_gp(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dd_I1_gn(eNN,numDOFs1,eNN,numDOFs1);

    tensor<double,4> dd_I1L1(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dd_I1L2(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dd_I1L3(eNN,numDOFs1,eNN,numDOFs1);




    tensor<double,4> dd_I3_gp(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dd_I3_gn(eNN,numDOFs1,eNN,numDOFs1);

    tensor<double,4> dd_I3L1(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dd_I3L2(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dd_I3L3(eNN,numDOFs1,eNN,numDOFs1);



    tensor<double,4> dd_I3_gp_SR(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dd_I3_gn_SR(eNN,numDOFs1,eNN,numDOFs1);

    tensor<double,4> dd_I3L1_SR(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dd_I3L2_SR(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dd_I3L3_SR(eNN,numDOFs1,eNN,numDOFs1);


    dd_I1_gp=product(GP,dd_metric_gp ,{{0,4},{1,5}}) ;
    dd_I1_gn=product(GN,dd_metric_gn ,{{0,4},{1,5}}) ;
    dd_I3_gp=(2*jac_gp*dd_jac_gp+2*outer(d_jac_gp,d_jac_gp))* jacGP;
    dd_I3_gn=(2*jac_gn*dd_jac_gn+2*outer(d_jac_gn,d_jac_gn))* jacGN;
    dd_I3_gp_SR=sqrt(jacGP)*dd_jac_gp;
    dd_I3_gn_SR=sqrt(jacGN)*dd_jac_gn;



    dd_I1L1=product(GL1,dd_CL1 ,{{0,4},{1,5}}) ;
    dd_I1L2=product(GL2,dd_CL2 ,{{0,4},{1,5}}) ;
    dd_I1L3=product(GL3,dd_CL3 ,{{0,4},{1,5}}) ;


    dd_I3L1=(2*jacCL1*dd_jacCL1+2*outer(d_jacCL1,d_jacCL1))* jacGL1;
    dd_I3L2=(2*jacCL2*dd_jacCL2+2*outer(d_jacCL2,d_jacCL2))* jacGL2;
    dd_I3L3=(2*jacCL3*dd_jacCL3+2*outer(d_jacCL3,d_jacCL3))* jacGL3;


    dd_I3L1_SR=sqrt(jacGL1)*dd_jacCL1;
    dd_I3L2_SR=sqrt(jacGL2)*dd_jacCL2;
    dd_I3L3_SR=sqrt(jacGL3)*dd_jacCL3;
//dh_d_I1,dh_d_I3

    tensor<double,4> dh_d_I1_gp(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dh_d_I1_gn(eNN,numDOFs1,eNN,numDOFs1);

    tensor<double,4> dh_d_I1L1(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dh_d_I1L2(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dh_d_I1L3(eNN,numDOFs1,eNN,numDOFs1);




    tensor<double,4> dh_d_I3_gp(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dh_d_I3_gn(eNN,numDOFs1,eNN,numDOFs1);

    tensor<double,4> dh_d_I3L1(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dh_d_I3L2(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dh_d_I3L3(eNN,numDOFs1,eNN,numDOFs1);



    tensor<double,4> dh_d_I3_gp_SR(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dh_d_I3_gn_SR(eNN,numDOFs1,eNN,numDOFs1);

    tensor<double,4> dh_d_I3L1_SR(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dh_d_I3L2_SR(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dh_d_I3L3_SR(eNN,numDOFs1,eNN,numDOFs1);



    dh_d_I1_gp=product(GP,dh_d_metric_gp ,{{0,4},{1,5}}) ;
    dh_d_I1_gn=product(GN,dh_d_metric_gn ,{{0,4},{1,5}}) ;
    dh_d_I3_gp=(2*jac_gp*dh_d_jac_gp+2*outer(dh_jac_gp,d_jac_gp))* jacGP;
    dh_d_I3_gn=(2*jac_gn*dh_d_jac_gn+2*outer(dh_jac_gn,d_jac_gn))* jacGN;
    dh_d_I3_gp_SR=sqrt(jacGP)*dh_d_jac_gp;
    dh_d_I3_gn_SR=sqrt(jacGN)*dh_d_jac_gn;



    dh_d_I1L1=product(GL1,dh_d_CL1 ,{{0,4},{1,5}}) ;
    dh_d_I1L2=product(GL2,dh_d_CL2 ,{{0,4},{1,5}}) ;
    dh_d_I1L3=product(GL3,dh_d_CL3 ,{{0,4},{1,5}}) ;


    dh_d_I3L1=(2*jacCL1*dh_d_jacCL1+2*outer(dh_jacCL1,d_jacCL1))* jacGL1;
    dh_d_I3L2=(2*jacCL2*dh_d_jacCL2+2*outer(dh_jacCL2,d_jacCL2))* jacGL2;
    dh_d_I3L3=(2*jacCL3*dh_d_jacCL3+2*outer(dh_jacCL3,d_jacCL3))* jacGL3;


    dh_d_I3L1_SR=sqrt(jacGL1)*dh_d_jacCL1;
    dh_d_I3L2_SR=sqrt(jacGL2)*dh_d_jacCL2;
    dh_d_I3L3_SR=sqrt(jacGL3)*dh_d_jacCL3;
//d_dh_I1,d_dh_I3

    tensor<double,4> d_dh_I1_gp(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> d_dh_I1_gn(eNN,numDOFs1,eNN,numDOFs1);

    tensor<double,4> d_dh_I1L1(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> d_dh_I1L2(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> d_dh_I1L3(eNN,numDOFs1,eNN,numDOFs1);

    tensor<double,4> d_dh_I3_gp(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> d_dh_I3_gn(eNN,numDOFs1,eNN,numDOFs1);

    tensor<double,4> d_dh_I3L1(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> d_dh_I3L2(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> d_dh_I3L3(eNN,numDOFs1,eNN,numDOFs1);

    tensor<double,4> d_dh_I3_gp_SR(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> d_dh_I3_gn_SR(eNN,numDOFs1,eNN,numDOFs1);

    tensor<double,4> d_dh_I3L1_SR(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> d_dh_I3L2_SR(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> d_dh_I3L3_SR(eNN,numDOFs1,eNN,numDOFs1);



    d_dh_I1_gp=product(GP,d_dh_metric_gp ,{{0,4},{1,5}}) ;
    d_dh_I1_gn=product(GN,d_dh_metric_gn ,{{0,4},{1,5}}) ;
    d_dh_I3_gp=(2*jac_gp*d_dh_jac_gp+2*outer(d_jac_gp,dh_jac_gp))* jacGP;
    d_dh_I3_gn=(2*jac_gn*d_dh_jac_gn+2*outer(d_jac_gn,dh_jac_gn))* jacGN;
    d_dh_I3_gp_SR=sqrt(jacGP)*d_dh_jac_gp;
    d_dh_I3_gn_SR=sqrt(jacGN)*d_dh_jac_gn;



    d_dh_I1L1=product(GL1,d_dh_CL1 ,{{0,4},{1,5}}) ;
    d_dh_I1L2=product(GL2,d_dh_CL2 ,{{0,4},{1,5}}) ;
    d_dh_I1L3=product(GL3,d_dh_CL3 ,{{0,4},{1,5}}) ;


    d_dh_I3L1=(2*jacCL1*d_dh_jacCL1+2*outer(d_jacCL1,dh_jacCL1))* jacGL1;
    d_dh_I3L2=(2*jacCL2*d_dh_jacCL2+2*outer(d_jacCL2,dh_jacCL2))* jacGL2;
    d_dh_I3L3=(2*jacCL3*d_dh_jacCL3+2*outer(d_jacCL3,dh_jacCL3))* jacGL3;


    d_dh_I3L1_SR=sqrt(jacGL1)*d_dh_jacCL1;
    d_dh_I3L2_SR=sqrt(jacGL2)*d_dh_jacCL2;
    d_dh_I3L3_SR=sqrt(jacGL3)*d_dh_jacCL3;

    //dh_dh_I1,dh_dh_I3

    tensor<double,4> dh_dh_I1_gp(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dh_dh_I1_gn(eNN,numDOFs1,eNN,numDOFs1);

    tensor<double,4> dh_dh_I1L1(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dh_dh_I1L2(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dh_dh_I1L3(eNN,numDOFs1,eNN,numDOFs1);




    tensor<double,4> dh_dh_I3_gp(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dh_dh_I3_gn(eNN,numDOFs1,eNN,numDOFs1);

    tensor<double,4> dh_dh_I3L1(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dh_dh_I3L2(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dh_dh_I3L3(eNN,numDOFs1,eNN,numDOFs1);



    tensor<double,4> dh_dh_I3_gp_SR(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dh_dh_I3_gn_SR(eNN,numDOFs1,eNN,numDOFs1);

    tensor<double,4> dh_dh_I3L1_SR(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dh_dh_I3L2_SR(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dh_dh_I3L3_SR(eNN,numDOFs1,eNN,numDOFs1);



    dh_dh_I1_gp=product(GP,dh_dh_metric_gp ,{{0,4},{1,5}}) ;
    dh_dh_I1_gn=product(GN,dh_dh_metric_gn ,{{0,4},{1,5}}) ;
    dh_dh_I3_gp=(2*jac_gp*dh_dh_jac_gp+2*outer(dh_jac_gp,dh_jac_gp))* jacGP;
    dh_dh_I3_gn=(2*jac_gn*dh_dh_jac_gn+2*outer(dh_jac_gn,dh_jac_gn))* jacGN;
    dh_dh_I3_gp_SR=sqrt(jacGP)*dh_dh_jac_gp;
    dh_dh_I3_gn_SR=sqrt(jacGN)*dh_dh_jac_gn;



    dh_dh_I1L1=product(GL1,dh_dh_CL1 ,{{0,4},{1,5}}) ;
    dh_dh_I1L2=product(GL2,dh_dh_CL2 ,{{0,4},{1,5}}) ;
    dh_dh_I1L3=product(GL3,dh_dh_CL3 ,{{0,4},{1,5}}) ;


    dh_dh_I3L1=(2*jacCL1*dh_dh_jacCL1+2*outer(dh_jacCL1,dh_jacCL1))* jacGL1;
    dh_dh_I3L2=(2*jacCL2*dh_dh_jacCL2+2*outer(dh_jacCL2,dh_jacCL2))* jacGL2;
    dh_dh_I3L3=(2*jacCL3*dh_dh_jacCL3+2*outer(dh_jacCL3,dh_jacCL3))* jacGL3;


    dh_dh_I3L1_SR=sqrt(jacGL1)*dh_dh_jacCL1;
    dh_dh_I3L2_SR=sqrt(jacGL2)*dh_dh_jacCL2;
    dh_dh_I3L3_SR=sqrt(jacGL3)*dh_dh_jacCL3;
//


 // cout<<"print: "<<kspr/deltat<< ": "<< <<endl;


 double jacPRR=jac_gp0*rho_api;
 double jacNRR=jac_gn0*rho_bas;

 double jacLRR1=jacCL1_n*jacR*rho_lat1;
 double jacLRR2=jacCL2_n*jacR*rho_lat2;
 double jacLRR3=jacCL3_n*jacR*rho_lat3;

    ////small_elasticity not evolve


    //  viscoelastic metric

    double jacGP_ref       = 1/(jac_GPR*jac_GPR); // J of GPP (G11P*G22P-G12P*G12P)
    double jacGN_ref       = 1/(jac_GNR*jac_GNR); // J of GPP (G11P*G22P-G12P*G12P)


//first and third invariants
    double I1_gp_ref=product(imetric_GPR,metric_gp ,{{0,0},{1,1}}) ;
    double I1_gn_ref=product(imetric_GNR,metric_gn ,{{0,0},{1,1}}) ;

    double I3_gp_ref=jac_gp*jac_gp* jacGP_ref;
    double I3_gn_ref=jac_gn*jac_gn*jacGN_ref;

    double I3_gp_SR_ref=sqrt(jacGP_ref)*jac_gp;
    double I3_gn_SR_ref=sqrt(jacGN_ref)*jac_gn;


//d_I1_ref

    tensor<double,2> d_I1_gp_ref(eNN,numDOFs1);
    tensor<double,2> d_I1_gn_ref(eNN,numDOFs1);
    tensor<double,2> d_I3_gp_ref(eNN,numDOFs1);
    tensor<double,2> d_I3_gn_ref(eNN,numDOFs1);

    tensor<double,2> d_I3_gp_SR_ref(eNN,numDOFs1);
    tensor<double,2> d_I3_gn_SR_ref(eNN,numDOFs1);



    d_I1_gp_ref=product(imetric_GPR,d_metric_gp ,{{0,2},{1,3}}) ;
    d_I1_gn_ref=product(imetric_GNR,d_metric_gn ,{{0,2},{1,3}}) ;



    d_I3_gp_ref=2*jac_gp*d_jac_gp* jacGP_ref;
    d_I3_gn_ref=2*jac_gp*d_jac_gn* jacGN_ref;


    d_I3_gp_SR_ref=sqrt(jacGP_ref)*d_jac_gp;
    d_I3_gn_SR_ref=sqrt(jacGN_ref)*d_jac_gn;

//dh_I1_ref

    tensor<double,2> dh_I1_gp_ref(eNN,numDOFs1);
    tensor<double,2> dh_I1_gn_ref(eNN,numDOFs1);
    tensor<double,2> dh_I3_gp_ref(eNN,numDOFs1);
    tensor<double,2> dh_I3_gn_ref(eNN,numDOFs1);

    tensor<double,2> dh_I3_gp_SR_ref(eNN,numDOFs1);
    tensor<double,2> dh_I3_gn_SR_ref(eNN,numDOFs1);



    dh_I1_gp_ref=product(imetric_GPR,dh_metric_gp ,{{0,2},{1,3}}) ;
    dh_I1_gn_ref=product(imetric_GNR,dh_metric_gn ,{{0,2},{1,3}}) ;



    dh_I3_gp_ref=2*jac_gp*dh_jac_gp* jacGP_ref;
    dh_I3_gn_ref=2*jac_gp*dh_jac_gn* jacGN_ref;


    dh_I3_gp_SR_ref=sqrt(jacGP_ref)*dh_jac_gp;
    dh_I3_gn_SR_ref=sqrt(jacGN_ref)*dh_jac_gn;



    //dd_I1_ref
    tensor<double,4> dd_I1_gp_ref(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dd_I1_gn_ref(eNN,numDOFs1,eNN,numDOFs1);

    tensor<double,4> dd_I3_gp_ref(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dd_I3_gn_ref(eNN,numDOFs1,eNN,numDOFs1);

    tensor<double,4> dd_I3_gp_SR_ref(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dd_I3_gn_SR_ref(eNN,numDOFs1,eNN,numDOFs1);

    dd_I1_gp_ref=product(imetric_GPR,dd_metric_gp ,{{0,4},{1,5}}) ;
    dd_I1_gn_ref=product(imetric_GNR,dd_metric_gn ,{{0,4},{1,5}}) ;
    dd_I3_gp_ref=(2*jac_gp*dd_jac_gp+2*outer(d_jac_gp,d_jac_gp))* jacGP_ref;
    dd_I3_gn_ref=(2*jac_gn*dd_jac_gn+2*outer(d_jac_gn,d_jac_gn))* jacGN_ref;
    dd_I3_gp_SR_ref=sqrt(jacGP_ref)*dd_jac_gp;
    dd_I3_gn_SR_ref=sqrt(jacGN_ref)*dd_jac_gn;

    //dh_d_I1_ref
    tensor<double,4> dh_d_I1_gp_ref(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dh_d_I1_gn_ref(eNN,numDOFs1,eNN,numDOFs1);

    tensor<double,4> dh_d_I3_gp_ref(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dh_d_I3_gn_ref(eNN,numDOFs1,eNN,numDOFs1);

    tensor<double,4> dh_d_I3_gp_SR_ref(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dh_d_I3_gn_SR_ref(eNN,numDOFs1,eNN,numDOFs1);

    dh_d_I1_gp_ref=product(imetric_GPR,dh_d_metric_gp ,{{0,4},{1,5}}) ;
    dh_d_I1_gn_ref=product(imetric_GNR,dh_d_metric_gn ,{{0,4},{1,5}}) ;
    dh_d_I3_gp_ref=(2*jac_gp*dh_d_jac_gp+2*outer(dh_jac_gp,d_jac_gp))* jacGP_ref;
    dh_d_I3_gn_ref=(2*jac_gn*dh_d_jac_gn+2*outer(dh_jac_gn,d_jac_gn))* jacGN_ref;
    dh_d_I3_gp_SR_ref=sqrt(jacGP_ref)*dh_d_jac_gp;
    dh_d_I3_gn_SR_ref=sqrt(jacGN_ref)*dh_d_jac_gn;

    //d_dh_I1_ref
    tensor<double,4> d_dh_I1_gp_ref(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> d_dh_I1_gn_ref(eNN,numDOFs1,eNN,numDOFs1);

    tensor<double,4> d_dh_I3_gp_ref(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> d_dh_I3_gn_ref(eNN,numDOFs1,eNN,numDOFs1);

    tensor<double,4> d_dh_I3_gp_SR_ref(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> d_dh_I3_gn_SR_ref(eNN,numDOFs1,eNN,numDOFs1);

    d_dh_I1_gp_ref=product(imetric_GPR,d_dh_metric_gp ,{{0,4},{1,5}}) ;
    d_dh_I1_gn_ref=product(imetric_GNR,d_dh_metric_gn ,{{0,4},{1,5}}) ;
    d_dh_I3_gp_ref=(2*jac_gp*d_dh_jac_gp+2*outer(d_jac_gp,dh_jac_gp))* jacGP_ref;
    d_dh_I3_gn_ref=(2*jac_gn*d_dh_jac_gn+2*outer(d_jac_gn,dh_jac_gn))* jacGN_ref;
    d_dh_I3_gp_SR_ref=sqrt(jacGP_ref)*d_dh_jac_gp;
    d_dh_I3_gn_SR_ref=sqrt(jacGN_ref)*d_dh_jac_gn;

    //dh_dh_I1_ref
    tensor<double,4> dh_dh_I1_gp_ref(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dh_dh_I1_gn_ref(eNN,numDOFs1,eNN,numDOFs1);

    tensor<double,4> dh_dh_I3_gp_ref(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dh_dh_I3_gn_ref(eNN,numDOFs1,eNN,numDOFs1);

    tensor<double,4> dh_dh_I3_gp_SR_ref(eNN,numDOFs1,eNN,numDOFs1);
    tensor<double,4> dh_dh_I3_gn_SR_ref(eNN,numDOFs1,eNN,numDOFs1);

    dh_dh_I1_gp_ref=product(imetric_GPR,dh_dh_metric_gp ,{{0,4},{1,5}}) ;
    dh_dh_I1_gn_ref=product(imetric_GNR,dh_dh_metric_gn ,{{0,4},{1,5}}) ;
    dh_dh_I3_gp_ref=(2*jac_gp*dh_dh_jac_gp+2*outer(dh_jac_gp,dh_jac_gp))* jacGP_ref;
    dh_dh_I3_gn_ref=(2*jac_gn*dh_dh_jac_gn+2*outer(dh_jac_gn,dh_jac_gn))* jacGN_ref;
    dh_dh_I3_gp_SR_ref=sqrt(jacGP_ref)*dh_dh_jac_gp;
    dh_dh_I3_gn_SR_ref=sqrt(jacGN_ref)*dh_dh_jac_gn;


    //
    double mu_ela=mu*fact_elastic ;
    double lambda_ela=lambda*fact_elastic;

    double E_small =(0.5*lambda_ela*log(I3_gp_SR_ref)*log(I3_gp_SR_ref)-mu_ela*log(I3_gp_SR_ref)+0.5*mu_ela*(I1_gp_ref-2) )*jacPRR+ (0.5*lambda_ela*log(I3_gn_SR_ref)*log(I3_gn_SR_ref)-mu_ela*log(I3_gn_SR_ref)+0.5*mu_ela*(I1_gn_ref-2) )*jacNRR;

/*

//apical
    Bp(all,range(0,2)) += ((lambda_ela*log(I3_gp_SR_ref)-mu_ela)/I3_gp_SR_ref*d_I3_gp_SR_ref+0.5*mu_ela*d_I1_gp_ref)*jacPRR;
    Bp(all,range(3,5)) += ((lambda_ela*log(I3_gp_SR_ref)-mu_ela)/I3_gp_SR_ref*dh_I3_gp_SR_ref+0.5*mu_ela*dh_I1_gp_ref)*jacPRR;
    //dd
    App(all,range(0,2),all,range(0,2))  +=((lambda_ela*log(I3_gp_SR_ref)-mu_ela)/I3_gp_SR_ref*dd_I3_gp_SR_ref+0.5*mu_ela*dd_I1_gp_ref)*jacPRR;
    App(all,range(0,2),all,range(0,2))  +=(lambda_ela*(1-log(I3_gp_SR_ref))+mu_ela)/(I3_gp_SR_ref*I3_gp_SR_ref)*outer(d_I3_gp_SR_ref,d_I3_gp_SR_ref)*jacPRR;// term 1
   //dh_d
    App(all,range(0,2),all,range(3,5))  +=((lambda_ela*log(I3_gp_SR_ref)-mu_ela)/I3_gp_SR_ref*dh_d_I3_gp_SR_ref.transpose({2,3,0,1})+0.5*mu_ela*dh_d_I1_gp_ref.transpose({2,3,0,1}))*jacPRR;
    App(all,range(0,2),all,range(3,5))  +=(lambda_ela*(1-log(I3_gp_SR_ref))+mu_ela)/(I3_gp_SR_ref*I3_gp_SR_ref)*outer(d_I3_gp_SR_ref,dh_I3_gp_SR_ref)*jacPRR;// term 1
    //d_dh
    App(all,range(3,5),all,range(0,2))  +=((lambda_ela*log(I3_gp_SR_ref)-mu_ela)/I3_gp_SR_ref*d_dh_I3_gp_SR_ref.transpose({2,3,0,1})+0.5*mu_ela*d_dh_I1_gp_ref.transpose({2,3,0,1}))*jacPRR;
    App(all,range(3,5),all,range(0,2))  +=(lambda_ela*(1-log(I3_gp_SR_ref))+mu_ela)/(I3_gp_SR_ref*I3_gp_SR_ref)*outer(dh_I3_gp_SR_ref,d_I3_gp_SR_ref)*jacPRR;// term 1
   //dh_dh
    App(all,range(3,5),all,range(3,5))  +=((lambda_ela*log(I3_gp_SR_ref)-mu_ela)/I3_gp_SR_ref*dh_dh_I3_gp_SR_ref+0.5*mu_ela*dh_dh_I1_gp_ref)*jacPRR;
    App(all,range(3,5),all,range(3,5))  +=(lambda_ela*(1-log(I3_gp_SR_ref))+mu_ela)/(I3_gp_SR_ref*I3_gp_SR_ref)*outer(dh_I3_gp_SR_ref,dh_I3_gp_SR_ref)*jacPRR;// term 1

//basal

    Bp(all,range(0,2)) +=((lambda_ela*log(I3_gn_SR_ref)-mu_ela)/I3_gn_SR_ref*d_I3_gn_SR_ref+0.5*mu_ela*d_I1_gn_ref)*jacNRR;
    Bp(all,range(3,5)) +=((lambda_ela*log(I3_gn_SR_ref)-mu_ela)/I3_gn_SR_ref*dh_I3_gn_SR_ref+0.5*mu_ela*dh_I1_gn_ref)*jacNRR;

    //dd
    App(all,range(0,2),all,range(0,2))  += ((lambda_ela*log(I3_gn_SR_ref)-mu_ela)/I3_gn_SR_ref*dd_I3_gn_SR_ref+0.5*mu_ela*dd_I1_gn_ref)*jacNRR;
    App(all,range(0,2),all,range(0,2))  +=(lambda_ela*(1-log(I3_gn_SR_ref))+mu_ela)/(I3_gn_SR_ref*I3_gn_SR_ref)*outer(d_I3_gn_SR_ref,d_I3_gn_SR_ref)*jacNRR;
    //dh_d
    App(all,range(0,2),all,range(3,5))  += ((lambda_ela*log(I3_gn_SR_ref)-mu_ela)/I3_gn_SR_ref*dh_d_I3_gn_SR_ref.transpose({2,3,0,1})+0.5*mu_ela*dh_d_I1_gn_ref.transpose({2,3,0,1}))*jacNRR;
    App(all,range(0,2),all,range(3,5))  +=(lambda_ela*(1-log(I3_gn_SR_ref))+mu_ela)/(I3_gn_SR_ref*I3_gn_SR_ref)*outer(d_I3_gn_SR_ref,dh_I3_gn_SR_ref)*jacNRR;
    //d_dh
    App(all,range(3,5),all,range(0,2))  += ((lambda_ela*log(I3_gn_SR_ref)-mu_ela)/I3_gn_SR_ref*d_dh_I3_gn_SR_ref.transpose({2,3,0,1})+0.5*mu_ela*d_dh_I1_gn_ref.transpose({2,3,0,1}))*jacNRR;
    App(all,range(3,5),all,range(0,2))  +=(lambda_ela*(1-log(I3_gn_SR_ref))+mu_ela)/(I3_gn_SR_ref*I3_gn_SR_ref)*outer(dh_I3_gn_SR_ref,d_I3_gn_SR_ref)*jacNRR;
    //dh_dh
    App(all,range(3,5),all,range(3,5))  += ((lambda_ela*log(I3_gn_SR_ref)-mu_ela)/I3_gn_SR_ref*dh_dh_I3_gn_SR_ref+0.5*mu_ela*dh_dh_I1_gn_ref)*jacNRR;
    App(all,range(3,5),all,range(3,5))  +=(lambda_ela*(1-log(I3_gn_SR_ref))+mu_ela)/(I3_gn_SR_ref*I3_gn_SR_ref)*outer(dh_I3_gn_SR_ref,dh_I3_gn_SR_ref)*jacNRR;
//////////
*/

    double EA                     = (0.5*lambda*log(I3_gp_SR)*log(I3_gp_SR)-mu*log(I3_gp_SR)+0.5*mu*(I1_gp-2) );
    double EB                    = (0.5*lambda*log(I3_gn_SR)*log(I3_gn_SR)-mu*log(I3_gn_SR)+0.5*mu*(I1_gn-2) );
    double EL                    = (0.5*lambda*log(I3L1_SR)*log(I3L1_SR)-mu*log(I3L1_SR)+0.5*mu*(I1L1-2) )+(0.5*lambda*log(I3L2_SR)*log(I3L2_SR)-mu*log(I3L2_SR)+0.5*mu*(I1L2-2) )+(0.5*lambda*log(I3L3_SR)*log(I3L3_SR)-mu*log(I3L3_SR)+0.5*mu*(I1L3-2) );
    //apical

    //VISCO-ELASTIC ENERGY
    double Evisco                     = (0.5*lambda*log(I3_gp_SR)*log(I3_gp_SR)-mu*log(I3_gp_SR)+0.5*mu*(I1_gp-2) )*jacPRR+ (0.5*lambda*log(I3_gn_SR)*log(I3_gn_SR)-mu*log(I3_gn_SR)+0.5*mu*(I1_gn-2) )*jacNRR;

    double Evisco1                     = (0.5*lambda*log(I3_gp_SR)*log(I3_gp_SR)-mu*log(I3_gp_SR)+0.5*mu*(I1_gp-2) )*jacPRR;
    double Evisco2                    = (0.5*lambda*log(I3_gn_SR)*log(I3_gn_SR)-mu*log(I3_gn_SR)+0.5*mu*(I1_gn-2) )*jacNRR;
    double Evisco3                    = (0.5*lambda*log(I3L1_SR)*log(I3L1_SR)-mu*log(I3L1_SR)+0.5*mu*(I1L1-2) )*f0/3*jacLRR1+(0.5*lambda*log(I3L2_SR)*log(I3L2_SR)-mu*log(I3L2_SR)+0.5*mu*(I1L2-2) )*f0/3*jacLRR2+(0.5*lambda*log(I3L3_SR)*log(I3L3_SR)-mu*log(I3L3_SR)+0.5*mu*(I1L3-2) )*f0/3*jacLRR3;
    //apical
    Bp(all,range(0,2)) += ((lambda*log(I3_gp_SR)-mu)/I3_gp_SR*d_I3_gp_SR+0.5*mu*d_I1_gp)*jacPRR;
    Bp(all,range(3,5)) += ((lambda*log(I3_gp_SR)-mu)/I3_gp_SR*dh_I3_gp_SR+0.5*mu*dh_I1_gp)*jacPRR;
    //dd_
    App(all,range(0,2),all,range(0,2))  +=((lambda*log(I3_gp_SR)-mu)/I3_gp_SR*dd_I3_gp_SR+0.5*mu*dd_I1_gp)*jacPRR;
    App(all,range(0,2),all,range(0,2))  +=(lambda*(1-log(I3_gp_SR))+mu)/(I3_gp_SR*I3_gp_SR)*outer(d_I3_gp_SR,d_I3_gp_SR)*jacPRR;// term 1
   //dh_d_
    App(all,range(0,2),all,range(3,5))  +=((lambda*log(I3_gp_SR)-mu)/I3_gp_SR*dh_d_I3_gp_SR.transpose({2,3,0,1})+0.5*mu*dh_d_I1_gp.transpose({2,3,0,1}))*jacPRR;
    App(all,range(0,2),all,range(3,5))  +=(lambda*(1-log(I3_gp_SR))+mu)/(I3_gp_SR*I3_gp_SR)*outer(d_I3_gp_SR,dh_I3_gp_SR)*jacPRR;// term 1
    //d_dh_
    App(all,range(3,5),all,range(0,2))  +=((lambda*log(I3_gp_SR)-mu)/I3_gp_SR*d_dh_I3_gp_SR.transpose({2,3,0,1})+0.5*mu*d_dh_I1_gp.transpose({2,3,0,1}))*jacPRR;
    App(all,range(3,5),all,range(0,2))  +=(lambda*(1-log(I3_gp_SR))+mu)/(I3_gp_SR*I3_gp_SR)*outer(dh_I3_gp_SR,d_I3_gp_SR)*jacPRR;// term 1
    //dh_dh_
    App(all,range(3,5),all,range(3,5))  +=((lambda*log(I3_gp_SR)-mu)/I3_gp_SR*dh_dh_I3_gp_SR+0.5*mu*dh_dh_I1_gp)*jacPRR;
    App(all,range(3,5),all,range(3,5))  +=(lambda*(1-log(I3_gp_SR))+mu)/(I3_gp_SR*I3_gp_SR)*outer(dh_I3_gp_SR,dh_I3_gp_SR)*jacPRR;// term 1

  //basal
     Bp(all,range(0,2)) +=((lambda*log(I3_gn_SR)-mu)/I3_gn_SR*d_I3_gn_SR+0.5*mu*d_I1_gn)*jacNRR;
     Bp(all,range(3,5)) +=((lambda*log(I3_gn_SR)-mu)/I3_gn_SR*dh_I3_gn_SR+0.5*mu*dh_I1_gn)*jacNRR;
     //dd_
     App(all,range(0,2),all,range(0,2))  += ((lambda*log(I3_gn_SR)-mu)/I3_gn_SR*dd_I3_gn_SR+0.5*mu*dd_I1_gn)*jacNRR;
     App(all,range(0,2),all,range(0,2))  +=(lambda*(1-log(I3_gn_SR))+mu)/(I3_gn_SR*I3_gn_SR)*outer(d_I3_gn_SR,d_I3_gn_SR)*jacNRR;

     //dh_d_
     App(all,range(0,2),all,range(3,5))  += ((lambda*log(I3_gn_SR)-mu)/I3_gn_SR*dh_d_I3_gn_SR.transpose({2,3,0,1})+0.5*mu*dh_d_I1_gn.transpose({2,3,0,1}))*jacNRR;
     App(all,range(0,2),all,range(3,5))  +=(lambda*(1-log(I3_gn_SR))+mu)/(I3_gn_SR*I3_gn_SR)*outer(d_I3_gn_SR,dh_I3_gn_SR)*jacNRR;
     //d_dh_
     App(all,range(3,5),all,range(0,2))  += ((lambda*log(I3_gn_SR)-mu)/I3_gn_SR*d_dh_I3_gn_SR.transpose({2,3,0,1})+0.5*mu*d_dh_I1_gn.transpose({2,3,0,1}))*jacNRR;
     App(all,range(3,5),all,range(0,2))  +=(lambda*(1-log(I3_gn_SR))+mu)/(I3_gn_SR*I3_gn_SR)*outer(dh_I3_gn_SR,d_I3_gn_SR)*jacNRR;
     //dh_dh_
     App(all,range(3,5),all,range(3,5))  += ((lambda*log(I3_gn_SR)-mu)/I3_gn_SR*dh_dh_I3_gn_SR+0.5*mu*dh_dh_I1_gn)*jacNRR;
     App(all,range(3,5),all,range(3,5))  +=(lambda*(1-log(I3_gn_SR))+mu)/(I3_gn_SR*I3_gn_SR)*outer(dh_I3_gn_SR,dh_I3_gn_SR)*jacNRR;

////////////////////////////////////////////////////


    //LATERAL
    //case 1
    Bp(all,range(0,2)) +=((lambda*log(I3L1_SR)-mu)/I3L1_SR*d_I3L1_SR+0.5*mu*d_I1L1)*f0/3*jacLRR1;
    Bp(all,range(3,5)) +=((lambda*log(I3L1_SR)-mu)/I3L1_SR*dh_I3L1_SR+0.5*mu*dh_I1L1)*f0/3*jacLRR1;

//dd_
    App(all,range(0,2),all,range(0,2))  +=((lambda*log(I3L1_SR)-mu)/I3L1_SR*dd_I3L1_SR+0.5*mu*dd_I1L1)*f0/3*jacLRR1;
    App(all,range(0,2),all,range(0,2))  +=(lambda*(1-log(I3L1_SR))+mu)/(I3L1_SR*I3L1_SR)*outer(d_I3L1_SR,d_I3L1_SR)*f0/3*jacLRR1;

    App(all,range(0,2),all,range(3,5))  +=((lambda*log(I3L1_SR)-mu)/I3L1_SR*dh_d_I3L1_SR.transpose({2,3,0,1})+0.5*mu*dh_d_I1L1.transpose({2,3,0,1}))*f0/3*jacLRR1;
    App(all,range(0,2),all,range(3,5))  +=(lambda*(1-log(I3L1_SR))+mu)/(I3L1_SR*I3L1_SR)*outer(d_I3L1_SR,dh_I3L1_SR)*f0/3*jacLRR1;
    //dh_d_
    App(all,range(3,5),all,range(0,2))  +=((lambda*log(I3L1_SR)-mu)/I3L1_SR*d_dh_I3L1_SR.transpose({2,3,0,1})+0.5*mu*d_dh_I1L1.transpose({2,3,0,1}))*f0/3*jacLRR1;
    App(all,range(3,5),all,range(0,2))  +=(lambda*(1-log(I3L1_SR))+mu)/(I3L1_SR*I3L1_SR)*outer(dh_I3L1_SR,d_I3L1_SR)*f0/3*jacLRR1;


//dh_dh_
    App(all,range(3,5),all,range(3,5))  +=((lambda*log(I3L1_SR)-mu)/I3L1_SR*dh_dh_I3L1_SR+0.5*mu*dh_dh_I1L1)*f0/3*jacLRR1;
    App(all,range(3,5),all,range(3,5))  +=(lambda*(1-log(I3L1_SR))+mu)/(I3L1_SR*I3L1_SR)*outer(dh_I3L1_SR,dh_I3L1_SR)*f0/3*jacLRR1;


//case2
    Bp(all,range(0,2)) +=((lambda*log(I3L2_SR)-mu)/I3L2_SR*d_I3L2_SR+0.5*mu*d_I1L2)*f0/3*jacLRR2;
    Bp(all,range(3,5)) +=((lambda*log(I3L2_SR)-mu)/I3L2_SR*dh_I3L2_SR+0.5*mu*dh_I1L2)*f0/3*jacLRR2;

    //dd_
    App(all,range(0,2),all,range(0,2))  +=((lambda*log(I3L2_SR)-mu)/I3L2_SR*dd_I3L2_SR+0.5*mu*dd_I1L2)*f0/3*jacLRR2;
    App(all,range(0,2),all,range(0,2))  +=(lambda*(1-log(I3L2_SR))+mu)/(I3L2_SR*I3L2_SR)*outer(d_I3L2_SR,d_I3L2_SR)*f0/3*jacLRR2;

    App(all,range(0,2),all,range(3,5))  +=((lambda*log(I3L2_SR)-mu)/I3L2_SR*dh_d_I3L2_SR.transpose({2,3,0,1})+0.5*mu*dh_d_I1L2.transpose({2,3,0,1}))*f0/3*jacLRR2;
    App(all,range(0,2),all,range(3,5))  +=(lambda*(1-log(I3L2_SR))+mu)/(I3L2_SR*I3L2_SR)*outer(d_I3L2_SR,dh_I3L2_SR)*f0/3*jacLRR2;

    //dh_d_
    App(all,range(3,5),all,range(0,2))  +=((lambda*log(I3L2_SR)-mu)/I3L2_SR*d_dh_I3L2_SR.transpose({2,3,0,1})+0.5*mu*d_dh_I1L2.transpose({2,3,0,1}))*f0/3*jacLRR2;
    App(all,range(3,5),all,range(0,2))  +=(lambda*(1-log(I3L2_SR))+mu)/(I3L2_SR*I3L2_SR)*outer(dh_I3L2_SR,d_I3L2_SR)*f0/3*jacLRR2;
    //dh_dh_
    App(all,range(3,5),all,range(3,5))  +=((lambda*log(I3L2_SR)-mu)/I3L2_SR*dh_dh_I3L2_SR+0.5*mu*dh_dh_I1L2)*f0/3*jacLRR2;
    App(all,range(3,5),all,range(3,5))  +=(lambda*(1-log(I3L2_SR))+mu)/(I3L2_SR*I3L2_SR)*outer(dh_I3L2_SR,dh_I3L2_SR)*f0/3*jacLRR2;

// case3
    Bp(all,range(0,2)) +=((lambda*log(I3L3_SR)-mu)/I3L3_SR*d_I3L3_SR+0.5*mu*d_I1L3)*f0/3*jacLRR3;
    Bp(all,range(3,5)) +=((lambda*log(I3L3_SR)-mu)/I3L3_SR*dh_I3L3_SR+0.5*mu*dh_I1L3)*f0/3*jacLRR3;

    //dd_
    App(all,range(0,2),all,range(0,2))  +=((lambda*log(I3L3_SR)-mu)/I3L3_SR*dd_I3L3_SR+0.5*mu*dd_I1L3)*f0/3*jacLRR3;
    App(all,range(0,2),all,range(0,2))  +=(lambda*(1-log(I3L3_SR))+mu)/(I3L3_SR*I3L3_SR)*outer(d_I3L3_SR,d_I3L3_SR)*f0/3*jacLRR3;

    App(all,range(0,2),all,range(3,5))  +=((lambda*log(I3L3_SR)-mu)/I3L3_SR*dh_d_I3L3_SR.transpose({2,3,0,1})+0.5*mu*dh_d_I1L3.transpose({2,3,0,1}))*f0/3*jacLRR3;
    App(all,range(0,2),all,range(3,5))  +=(lambda*(1-log(I3L3_SR))+mu)/(I3L3_SR*I3L3_SR)*outer(d_I3L3_SR,dh_I3L3_SR)*f0/3*jacLRR3;
    //dh_d_
    App(all,range(3,5),all,range(0,2))  +=((lambda*log(I3L3_SR)-mu)/I3L3_SR*d_dh_I3L3_SR.transpose({2,3,0,1})+0.5*mu*d_dh_I1L3.transpose({2,3,0,1}))*f0/3*jacLRR3;
    App(all,range(3,5),all,range(0,2))  +=(lambda*(1-log(I3L3_SR))+mu)/(I3L3_SR*I3L3_SR)*outer(dh_I3L3_SR,d_I3L3_SR)*f0/3*jacLRR3;

    //dh_dh_
    App(all,range(3,5),all,range(3,5))  +=((lambda*log(I3L3_SR)-mu)/I3L3_SR*dh_dh_I3L3_SR+0.5*mu*dh_dh_I1L3)*f0/3*jacLRR3;
    App(all,range(3,5),all,range(3,5))  +=(lambda*(1-log(I3L3_SR))+mu)/(I3L3_SR*I3L3_SR)*outer(dh_I3L3_SR,dh_I3L3_SR)*f0/3*jacLRR3;
    //[6.2] Friction


    double Efric                        = 0.5*fric*jac_n*product(x-x_n,x-x_n,{{0,0}});
    Lagrangian+=Efric;


    Bp(all,range(0,2))                 += fric*jac_n*outer(bf,x-x_n);
    App(all,range(0,2),all,range(0,2)) += fric*jac_n*outer(outer(bf,Id),bf).transpose({0,1,3,2});


  /*if (timestep > fric_start+control_fric)
    {
        Bp(all,range(0,2))                 += kspr_in*jacR*outer(bf,x-xR);
        App(all,range(0,2),all,range(0,2)) += kspr_in*jacR*outer(outer(bf,Id),bf).transpose({0,1,3,2});
   //cout<< ": "<<kspr_in<< ": "<<x_target<<    ": x" << x <<endl;

    }
*/







    //[6.3] Body force
    /*tensor<double,1> Id1(3);
    tensor<double,1> Id2(3);*/
    tensor<double,1> Id3(3);

    /* Id1(0)=1;
     Id1(1)=0;
     Id1(2)=0;

     Id2(0)=0;
     Id2(1)=1;
     Id2(2)=0;*/


    Id3(0)=0;
    Id3(1)=0;
    Id3(2)=1;


/*    tensor<double,2> d_x(eNN,numDOFs);
    tensor<double,2> d_y(eNN,numDOFs);*/
    tensor<double,2> d_z(eNN,numDOFs1);


    /*   d_x=outer(bf,Id1);
       d_y=outer(bf,Id2);*/
    d_z=outer(bf,Id3);


   // double  E_force = jacR * force * (x(2) - xR(2));
   // Bp(all, 2) += force * jacR * bf;

// bending


      //  Bp(all, range(0, 2)) += kappa*jacR * d_meancurva*meancurva;
      //  App(all, range(0, 2), all, range(0, 2)) += kappa*jacR *(dd_meancurva*meancurva+outer(d_meancurva,d_meancurva));

       

      

 if (timestep > fric_start+control_fric)
    {
    Bp(all,0)                 += kspr_in*jacR*bf*(x(0)-x_target(0));
    App(all,0,all,0) += kspr_in*jacR*outer(bf,bf);

    Bp(all,1)                 += kspr_in*jacR*bf*(x(1)-x_target(1));
    App(all,1,all,1) += kspr_in*jacR*outer(bf,bf);

    Bp(all,2)                 += kspr_in*jacR*bf*(x(2)-x_target(2));
    App(all,2,all,2) += kspr_in*jacR*outer(bf,bf);

    }





  if(spring>0.5)
{
    Lagrangian += jacR *0.5*kspr*product(x-xR,x-xR,{{0,0}});


    Bp(all,0)                 += kspr*jacR*bf*(x(0)-xR(0));
    App(all,0,all,0) += kspr*jacR*outer(bf,bf);

    Bp(all,1)                 += kspr*jacR*bf*(x(1)-xR(1));
    App(all,1,all,1) += kspr*jacR*outer(bf,bf);

    Bp(all,2)                 += kspr*jacR*bf*(x(2)-xR(2));
    App(all,2,all,2) += kspr*jacR*outer(bf,bf);


{


    Bp(all,range(3,5))                 += 0.25*kspr*jacR*outer(bf,H-HR);
    App(all,range(3,5),all,range(3,5)) += 0.25*kspr*jacR*outer(outer(bf,Id),bf).transpose({0,1,3,2});


     Bp(all,range(0,2))                 += 0.5*kspr*jacR*outer(bf,H-HR)*(-1);
     Bp(all,range(3,5))                 += 0.5*kspr*jacR*outer(bf,x-xR)*(-1);
     App(all,range(0,2),all,range(3,5)) += 0.5*kspr*(-1)*jacR*outer(outer(bf,Id),bf).transpose({0,1,3,2});
     App(all,range(3,5),all,range(0,2)) += 0.5*kspr*(-1)*jacR*outer(outer(bf,Id),bf).transpose({0,1,3,2});

}

   /* Bp(all,range(3,5))                 += kspr_out*jacR*outer(bf,H-HR);
    App(all,range(3,5),all,range(3,5)) += kspr_out*jacR*outer(outer(bf,Id),bf).transpose({0,1,3,2});
*/

 //   Bp(all,range(3,5))                 += kspr*jacR*outer(bf,H-HR);
   // App(all,range(3,5),all,range(3,5)) += kspr*jacR*outer(outer(bf,Id),bf).transpose({0,1,3,2});

}




    // [5.4]Adhesion Potential

    //Kconf=000000.0;
    Lagrangian += jac * ConfinementPotential(x, Kconf);

    Bp(all,range(0,2))                 += jac*d_z*dConfinementPotential(x, Kconf);
    Bp(all,range(0,2))        += d_jac*ConfinementPotential(x, Kconf);

    App(all,range(0,2),all,range(0,2)) +=outer(d_jac,d_z)*dConfinementPotential(x, Kconf);
    App(all,range(0,2),all,range(0,2)) +=jac*outer(d_z,d_z)*ddConfinementPotential(x, Kconf);

    App(all,range(0,2),all,range(0,2)) +=dd_jac*ConfinementPotential(x, Kconf);
    App(all,range(0,2),all,range(0,2)) +=outer(d_z,d_jac)*dConfinementPotential(x, Kconf);




//surface tension



        gamma_plus=gamma_plus*rho_api;
        gamma_minus=gamma_minus*rho_bas;

        Lagrangian += (jac_gp)*jacR * gamma_plus+(jac_gn)*jacR * gamma_minus+(gamma_l*f0/3)*((jacCL1)*rho_lat1+(jacCL2)*rho_lat2+(jacCL3)*rho_lat3)*jacR;

        double power1=gamma_plus*(jac_gp-jac_gp0)/(deltat*deltat);
        double power2=gamma_minus*(jac_gn-jac_gn0)/(deltat*deltat);
        double power3=(gamma_l)*((jacCL1-jacCL1_n)*rho_lat1+(jacCL2-jacCL2_n)*rho_lat2+(jacCL3-jacCL3_n)*rho_lat3);
        power3=power3/(deltat*deltat)*f0/3*jacR;



       double E_tens1=(jac_gp * gamma_plus*rho_api)/deltat;
       double E_tens2=(jac_gn * gamma_minus*rho_bas)/deltat;
       double E_tens3=(gamma_l*f0/3)*(jacCL1*rho_lat1+jacCL2*rho_lat2+jacCL3*rho_lat3)*jacR/deltat;


        Bp(all, range(0, 2)) += d_jac_gp * gamma_plus+d_jac_gn * gamma_minus;
        Bp(all, range(3, 5)) += dh_jac_gp * gamma_plus+dh_jac_gn * gamma_minus;


        App(all, range(0, 2), all, range(0, 2)) +=dd_jac_gp * gamma_plus+dd_jac_gn  * gamma_minus;
        App(all, range(0, 2), all, range(3, 5)) +=dh_d_jac_gp.transpose({2,3,0,1}) * gamma_plus+dh_d_jac_gn.transpose({2,3,0,1})  * gamma_minus;
        App(all, range(3, 5), all, range(0, 2)) +=d_dh_jac_gp.transpose({2,3,0,1}) * gamma_plus+d_dh_jac_gn.transpose({2,3,0,1})  * gamma_minus;

        App(all, range(3, 5), all, range(3, 5)) +=dh_dh_jac_gp * gamma_plus+dh_dh_jac_gn  * gamma_minus; //best


        Bp(all, range(0, 2)) += (gamma_l*f0/3)*(d_jacCL1*rho_lat1+d_jacCL2*rho_lat2+d_jacCL3*rho_lat3)*jacR;
        Bp(all, range(3, 5)) += (gamma_l*f0/3)*(dh_jacCL1*rho_lat1+dh_jacCL2*rho_lat2+dh_jacCL3*rho_lat3)*jacR;

        App(all, range(0, 2), all, range(0, 2)) +=(gamma_l*f0/3)*(dd_jacCL1*rho_lat1+dd_jacCL2*rho_lat2+dd_jacCL3*rho_lat3)*jacR;
        App(all, range(0, 2), all, range(3, 5)) +=(gamma_l*f0/3)*(dh_d_jacCL1.transpose({2,3,0,1})*rho_lat1+dh_d_jacCL2.transpose({2,3,0,1})*rho_lat2+dh_d_jacCL3.transpose({2,3,0,1})*rho_lat3)*jacR;

        App(all, range(3, 5), all, range(0, 2)) +=(gamma_l*f0/3)*(d_dh_jacCL1.transpose({2,3,0,1})*rho_lat1+d_dh_jacCL2.transpose({2,3,0,1})*rho_lat2+d_dh_jacCL3.transpose({2,3,0,1})*rho_lat3)*jacR;
        App(all, range(3, 5), all, range(3, 5)) +=(gamma_l*f0/3)*(dh_dh_jacCL1*rho_lat1+dh_dh_jacCL2*rho_lat2+dh_dh_jacCL3*rho_lat3)*jacR;






   // cout<<"print:jC "<< jacC_n<< " H: "<< H <<" lam: "<<lam   <<endl;
   // double Evol=0.0;


     Bp(all, range(0, 2)) += lam *(jacC*product(d_normal,H,{{2,0}})+d_jacC*product(H,normal,{{0,0}}))*jacR;
     Bp(all, range(3, 5)) += lam *(jacC*outer(bf, normal))*jacR;
     Bp(all, 6) += bf*(jacC*product(H,normal,{{0,0}})-height_in)*jacR; //checked

    App(all, range(0, 2), all, range(0, 2))+=lam *(jacC*product(dd_normal,H,{{4,0}})+dd_jacC*product(H,normal,{{0,0}})+outer(d_jacC,product(d_normal,H,{{2,0}}))+outer(product(d_normal,H,{{2,0}}),d_jacC))*jacR;
    App(all, range(0, 2), all, range(3, 5))+=lam*(jacC*outer(d_normal,bf).transpose({0,1,3,2})+outer(d_jacC,outer(bf,normal)))*jacR;
    App(all, range(0, 2), all, 6)+=outer((jacC*product(d_normal,H,{{2,0}})+d_jacC*product(H,normal,{{0,0}})),bf)*jacR;



    App(all, range(3, 5), all, range(0, 2))+=lam*(jacC*outer(bf,d_normal).transpose({0,3,1,2})+outer(outer(bf,normal),d_jacC))*jacR;
    App(all, range(3, 5), all, 6)+=jacC*outer(outer(bf, normal),bf)*jacR;




      App(all, 6, all, range(0, 2))+=jacC*outer(bf,(product(H,d_normal,{{0,2}})+d_jacC*product(H,normal,{{0,0}})))*jacR;
      App(all, 6, all, range(3, 5))+=jacC*outer(bf,outer(bf, normal))*jacR;



    //[6.4] Volume constraint


      double Evol=pressure/3.0 *(jac * xnormal - factor*jac_n*xnormal_n);
    Lagrangian+= Evol;


    Bp(all, range(0, 2)) += pressure / 3.0 * (d_jac * xnormal + jac * (d_normal * x) + jac * outer(bf, normal));


     Bg(0) += (jac * xnormal / 3) - factor * jacR * 1;
     App(all, range(0, 2), all, range(0, 2)) += pressure / 3.0 * (dd_jac * xnormal + outer(d_jac, d_normal * x) + outer(d_jac, outer(bf, normal)) +
                              outer(d_normal * x, d_jac) + outer(outer(bf, normal), d_jac) + jac * dd_normal * x +
                              jac * outer(d_normal, bf).transpose({0, 1, 3, 2}) +
                              jac * outer(bf, d_normal.transpose({2, 0, 1})));
     Apg(all, range(0, 2), 0) += 1.0 / 3.0 * (d_jac * xnormal + jac * d_normal * x + jac * outer(bf, normal));



    Agp = Apg.transpose({2,0,1});
    //[7] Add contributions to global integrals
    //[7] Add contributions to global integrals
    fillStr.addToGlobalIntegral("Dissipation",Efric/(deltat*deltat));

    fillStr.addToGlobalIntegral("E_tension1",x_target(0)*jacR);
    fillStr.addToGlobalIntegral("E_tension2",x_target(1)*jacR);
    fillStr.addToGlobalIntegral("E_tension3",x_target(2)*jacR);



    fillStr.addToGlobalIntegral("Energy",Evisco/deltat+Evisco3/deltat);

    fillStr.addToGlobalIntegral("E1",Evisco1/deltat);
    fillStr.addToGlobalIntegral("E2",Evisco2/deltat);
    fillStr.addToGlobalIntegral("E3",Evisco3/deltat);

    fillStr.addToGlobalIntegral("Lat_area", f0/3*(jacCL1+jacCL2+jacCL3)*jacR);
    fillStr.addToGlobalIntegral("Epress",Evol/deltat);

    // cout<< "king: " <<jacC<<endl;

    fillStr.addToGlobalIntegral("Trd",jac*jacR);
    fillStr.addToGlobalIntegral("area_n", jacR);
    fillStr.addToGlobalIntegral("area", jac);
    fillStr.addToGlobalIntegral("volume", (1.0/3.0) * jac*xnormal);
    fillStr.addToGlobalIntegral("Ela_small",E_small/deltat);

    fillStr.addToGlobalIntegral("area_gp",jac_gp);

    fillStr.addToGlobalIntegral("area_gn",jac_gn);


     fillStr.addToGlobalIntegral("P_tension1",power1);
     fillStr.addToGlobalIntegral("P_tension2",power2);
     fillStr.addToGlobalIntegral("P_tension3",power3);
    //This is to pass info to the other problem
    tensor<double,2> delI1_gp_delGP(pDim,pDim);//del I1 /del G ^bR
    tensor<double,2> delI1_gn_delGN(pDim,pDim);//


    tensor<double,2> delI1L1_delGL1(pDim,pDim);//
    tensor<double,2> delI1L2_delGL2(pDim,pDim);//
    tensor<double,2> delI1L3_delGL3(pDim,pDim);//

    delI1_gp_delGP=   (metric_gp)    ; //
    delI1_gn_delGN=   (metric_gn)    ; //_Br

    delI1L1_delGL1=(CL1);
    delI1L2_delGL2=(CL2);
    delI1L3_delGL3=(CL3);

    tensor<double,2> delI3_gp_delGP(pDim,pDim);//del I1 /del G ^bR
    tensor<double,2> delI3_gn_delGN(pDim,pDim);//

    tensor<double,2> delI3L1_delGL1(pDim,pDim);//
    tensor<double,2> delI3L2_delGL2(pDim,pDim);//
    tensor<double,2> delI3L3_delGL3(pDim,pDim);//


    delI3_gp_delGP=jac_gp*jac_gp*jacGP*iGP; //_Br
    delI3_gn_delGN=jac_gn*jac_gn*jacGN*iGN;

    delI3L1_delGL1=(jacCL1*jacCL1*jacGL1*iGL1);
    delI3L2_delGL2=(jacCL2*jacCL2*jacGL2*iGL2);
    delI3L3_delGL3=(jacCL3*jacCL3*jacGL3*iGL3);


    tensor<double,2> delWP_delGP(pDim,pDim);//del j+ /del G ^bR
    tensor<double,2> delWN_delGN(pDim,pDim);//

    tensor<double,2> delWL_delGL1(pDim,pDim);//
    tensor<double,2> delWL_delGL2(pDim,pDim);//
    tensor<double,2> delWL_delGL3(pDim,pDim);//


    delWP_delGP=0.5*(lambda*log(I3_gp_SR)-mu)/I3_gp*delI3_gp_delGP+0.5*mu*delI1_gp_delGP;
    delWN_delGN=0.5*(lambda*log(I3_gn_SR)-mu)/I3_gn*delI3_gn_delGN+0.5*mu*delI1_gn_delGN;

    delWL_delGL1=0.5*(lambda*log(I3L1_SR)-mu)/I3L1*delI3L1_delGL1+0.5*mu*delI1L1_delGL1;
    delWL_delGL2=0.5*(lambda*log(I3L2_SR)-mu)/I3L2*delI3L2_delGL2+0.5*mu*delI1L2_delGL2;
    delWL_delGL3=0.5*(lambda*log(I3L3_SR)-mu)/I3L3*delI3L3_delGL3+0.5*mu*delI1L3_delGL3;





    double* auxiliary_a = &fillStr.paramStr->a_aux[(subFill.loc_elemID*subFill.cubaInfo.iPts+subFill.kPt)*88];

    auxiliary_a[0]  = x(0);
    auxiliary_a[1]  = x(1);
    auxiliary_a[2]  = x(2);

    auxiliary_a[3]  =  metric_gp(0,0);
    auxiliary_a[4]  =  metric_gp(0,1);
    auxiliary_a[5]  =  metric_gp(1,0);
    auxiliary_a[6]  =  metric_gp(1,1);

    auxiliary_a[7]  =  metric_gn(0,0);
    auxiliary_a[8]  =  metric_gn(0,1);
    auxiliary_a[9]  =  metric_gn(1,0);
    auxiliary_a[10]  = metric_gn(1,1);


    auxiliary_a[11] = delWP_delGP(0,0);
    auxiliary_a[12] = delWP_delGP(0,1);
    auxiliary_a[13] = delWP_delGP(1,0);
    auxiliary_a[14] = delWP_delGP(1,1);

    auxiliary_a[15] = delWN_delGN(0,0);
    auxiliary_a[16] = delWN_delGN(0,1);
    auxiliary_a[17] = delWN_delGN(1,0);
    auxiliary_a[18] = delWN_delGN(1,1);



    auxiliary_a[19] = jacR;
    auxiliary_a[20] = jac;

    auxiliary_a[21] = Evisco1/deltat;
    auxiliary_a[22] =Evisco2/deltat;
    if (f0>0.001)
    {
        auxiliary_a[23] = 1 * delWL_delGL1(0, 0);
        auxiliary_a[24] = 1 * delWL_delGL1(1, 1);
        auxiliary_a[51] = 1 * delWL_delGL1(0, 1);

        auxiliary_a[31] = 1 * delWL_delGL2(0, 0);
        auxiliary_a[32] = 1 * delWL_delGL2(1, 1);
        auxiliary_a[52] = 1 * delWL_delGL2(0, 1);

        auxiliary_a[33] = 1 * delWL_delGL3(0, 0);
        auxiliary_a[34] = 1 * delWL_delGL3(1, 1);
        auxiliary_a[53] = 1 * delWL_delGL3(0, 1);


        auxiliary_a[25] = 1 * 0.000001 / deltat;
        auxiliary_a[26] = 1 * 0.000001 / deltat;
    }
    else
    {
        auxiliary_a[23] = 0.0 * delWL_delGL1(0, 0);
        auxiliary_a[24] = 0.0 * delWL_delGL1(1, 1);
        auxiliary_a[31] = 0.0 * delWL_delGL2(0, 0);
        auxiliary_a[32] = 0.0 * delWL_delGL2(1, 1);
        auxiliary_a[33] = 0.0 * delWL_delGL3(0, 0);
        auxiliary_a[34] = 0.0 * delWL_delGL3(1, 1);


        auxiliary_a[25] = 0.0 * 0.00001 / deltat;
        auxiliary_a[26] = 0.0 * 0.00001 / deltat;
    }

    auxiliary_a[27] =normal(0);
    auxiliary_a[28] =normal(1);
    auxiliary_a[29] = normal(2);
    auxiliary_a[30] = height0;

    auxiliary_a[35]  =  imetric_GPR(0,0);
    auxiliary_a[36]  =  imetric_GPR(0,1);
    auxiliary_a[37]  =  imetric_GPR(1,0);
    auxiliary_a[38]  = imetric_GPR(1,1);

    auxiliary_a[39]  =  imetric_GNR(0,0);
    auxiliary_a[40]  =  imetric_GNR(0,1);
    auxiliary_a[41]  =  imetric_GNR(1,0);
    auxiliary_a[42]  = imetric_GNR(1,1);



    auxiliary_a[43]  =  GP(0,0);
    auxiliary_a[44]  =  GP(1,0);
    auxiliary_a[45]  =  GP(0,1);
    auxiliary_a[46]  = GP(1,1);


    auxiliary_a[54]  =  CL1(0,0);
    auxiliary_a[55]  =  CL1(1,0);
    auxiliary_a[56]  =  CL1(0,1);
    auxiliary_a[57]  = CL1(1,1);

    auxiliary_a[58]  =  CL2(0,0);
    auxiliary_a[59]  =   CL2(1,0);
    auxiliary_a[60]  =   CL2(0,1);
    auxiliary_a[61]  = CL2(1,1);


    auxiliary_a[62]  =   CL3(0,0);
    auxiliary_a[63]  =  CL3(1,0);
    auxiliary_a[64]  =  CL3(0,1);
    auxiliary_a[65]  = CL3(1,1);


 double tr_rodt_plus = (jac_gp- jac_gp0) / jac_gp0;
 double tr_rodt_minus = (jac_gn- jac_gn0) / jac_gn0;
 double tr_rodt_lat1 = (jacCL1- jacCL1_n) / jacCL1_n;
 double tr_rodt_lat2 = (jacCL2- jacCL2_n) / jacCL2_n;
 double tr_rodt_lat3 = (jacCL3- jacCL3_n) / jacCL3_n;

    auxiliary_a[66]  =   tr_rodt_plus;
    auxiliary_a[67]  =  tr_rodt_minus;
    auxiliary_a[68]  =  tr_rodt_lat1;
    auxiliary_a[69]  = tr_rodt_lat2;
    auxiliary_a[70]  = tr_rodt_lat2;

    auxiliary_a[71]  =   jac_gp0;
    auxiliary_a[72]  =  jac_gn0;
    auxiliary_a[73]  =  jacCL1_n;
    auxiliary_a[74]  = jacCL2_n;
    auxiliary_a[75]  = jacCL3_n;

   auxiliary_a[76]  =   jacPRR;
    auxiliary_a[77]  =  jacNRR;
    auxiliary_a[78]  =  jacLRR1*f0/3.0;
    auxiliary_a[79]  = jacLRR2*f0/3.0;
    auxiliary_a[80]  = jacLRR3*f0/3.0;

auxiliary_a[81]  =EA;
auxiliary_a[82]  =EB;
auxiliary_a[83]  =EL;

auxiliary_a[84]  =power1;
auxiliary_a[85]  =power2;
auxiliary_a[86]  =power3;

auxiliary_a[87]  =jacC;

}


 void LS_ReacDif(hiperlife::FillStructure& fillStr)
{
    using namespace hiperlife;
    using namespace hiperlife::Tensor;

    // --------------------------------------------
    // [0] Initialize
    // --------------------------------------------
    auto& subFill = (fillStr)["viscoDHand"];

    //Parameters
    int numDOFs = subFill.numDOFs;
    int nDim    = subFill.nDim;
    int pDim    = subFill.pDim;
    int eNN     = subFill.eNN;

    //Coordinates and degrees of freedom
    wrapper<double,2> nborCoords(subFill.nborCoords.data(),eNN,nDim);
    wrapper<double,2> nborDOFs0(subFill.nborDOFs0.data(),eNN,numDOFs);
    wrapper<double,2> nborDOFs(subFill.nborDOFs.data(),eNN,numDOFs);

    //Basis functions and derivatives
    wrapper<double,1>    bf(subFill.nborBFsDers(0),eNN);
    wrapper<double,2>   Dbf(subFill.nborBFsDers(1),eNN,pDim);
    wrapper<double,3>  DDbf(subFill.nborBFsDers(2),eNN,pDim,pDim);

    //output
    wrapper<double,2>  Bk(fillStr.Bk(0).data(),eNN,numDOFs);
    wrapper<double,4>  Ak(fillStr.Ak(0,0).data(),eNN,numDOFs,eNN,numDOFs);


    // --------------------------------------------
    // [1] Model parameters
    // --------------------------------------------
    //[1.4] Parameters
 double deltat = fillStr.getRealParameter(MembParams::deltat);
 double fric2    = fillStr.getRealParameter(MembParams::fric2);
 int time_target    = fillStr.getIntParameter(MembParams::time_target);

    double young   = fillStr.getRealParameter(MembParams::young)* deltat;;
    double nu   = fillStr.getRealParameter(MembParams::poisson);

    int timestep =fillStr.getIntParameter(MembParams::timestep);//timestep
    
    double forward_new= fillStr.getRealParameter(MembParams::forward_new);





    double mu=0.5*young/(1+nu) ;
    double lambda=mu;
    // --------------------------------------------
    // [2] Auxiliary variables
    // --------------------------------------------
    //FIXME: this has to be loaded from some structure
    double* auxiliary_a = &fillStr.paramStr->a_aux[(subFill.loc_elemID*subFill.cubaInfo.iPts+subFill.kPt)*88];
    wrapper<double,1>      x(&auxiliary_a[0],3);
    wrapper<double,2> metric_gp(&auxiliary_a[3],2,2);
    wrapper<double,2> metric_gn(&auxiliary_a[7],2,2);
    wrapper<double,2> delWP_delGP(&auxiliary_a[11],2,2);
    wrapper<double,2> delWN_delGN(&auxiliary_a[15],2,2);

    wrapper<double,2> imetric_GPR(&auxiliary_a[35],2,2);
    wrapper<double,2>imetric_GNR(&auxiliary_a[39],2,2);
    tensor<double,2> metric_GPR=imetric_GPR.inv();
    tensor<double,2> metric_GNR=imetric_GNR.inv();

   wrapper<double,2>GP(&auxiliary_a[43],2,2);
   wrapper<double,2>GN(&auxiliary_a[47],2,2);

    wrapper<double,2>CL1(&auxiliary_a[54],2,2);
    wrapper<double,2>CL2(&auxiliary_a[58],2,2);
    wrapper<double,2>CL3(&auxiliary_a[62],2,2);

    //Jacobian
    double jacR= auxiliary_a[19];
    double jac= auxiliary_a[20];
    double E1_dt= auxiliary_a[21];
    double E2_dt= auxiliary_a[22];

    double jac_gp       = sqrt(metric_gp.det()); // J plus // sqrt(det gP)
    double jac_gn       = sqrt(metric_gn.det()); // J plus // sqrt(det gP)
    double jacCL1       = sqrt(CL1.det()); // J plus // sqrt(det gP)
    double jacCL2       = sqrt(CL2.det()); // J plus // sqrt(det gP)
    double jacCL3       = sqrt(CL3.det()); // J plus // sqrt(det gP)

    double delWL_delC1_dt= auxiliary_a[25];
    double delWL_delC2_dt= auxiliary_a[26];



    double delWL_delG1= auxiliary_a[23];
    double delWL_delG2= auxiliary_a[24];
    double delWL_delG12= auxiliary_a[51];

    double delWL_delG3= auxiliary_a[31];
    double delWL_delG4= auxiliary_a[32];
    double delWL_delG34= auxiliary_a[52];

    double delWL_delG5= auxiliary_a[33];
    double delWL_delG6= auxiliary_a[34];
    double delWL_delG56= auxiliary_a[53];


    //SUPG
    double G11P0{1.0}, G12P0{0.0},G22P0{1.0}, G11N0{1.0}, G12N0{0},G22N0{1.0} , G11P{1.0}, G12P{0.0},G22P{1.0}, G11N{1.0}, G12N{0},G22N{1.0}, G1L0{1.0},G2L0{1.0},G1L{1.0},G2L{1.0};
    double G3L0{1.0},G4L0{1.0},G5L0{1.0},G6L0{1.0},G3L{1.0},G4L{1.0},G5L{1.0},G6L{1.0},G12L{0.0},G34L{0.0},G56L{0.0},G12L0{0.0},G34L0{0.0},G56L0{0.0};


     G11P0  = nborDOFs0(all,0) * bf;
    G12P0  = nborDOFs0(all,1) * bf;
    G22P0  = nborDOFs0(all,2) * bf;
    G11N0  = nborDOFs0(all,3) * bf;
    G12N0  = nborDOFs0(all,4) * bf;
    G22N0  = nborDOFs0(all,5) * bf;

    G1L0  = nborDOFs0(all,6) * bf;
    G2L0  = nborDOFs0(all,7) * bf;

    G3L0  = nborDOFs0(all,8) * bf;
    G4L0  = nborDOFs0(all,9) * bf;
    G5L0  = nborDOFs0(all,10) * bf;
    G6L0  = nborDOFs0(all,11) * bf;


    G12L0  = nborDOFs0(all,12) * bf;
    G34L0  = nborDOFs0(all,13) * bf;
    G56L0  = nborDOFs0(all,14) * bf;
//cout<<"step: "<<timestep<<"check: "<<G11P0<<endl;

    // FIXME: We need to remove this

    tensor<double,2> delWP_delGP_CC(2,2);
    delWP_delGP_CC= product(imetric_GPR,product(delWP_delGP, imetric_GPR, {{1, 0}}), {{0, 0}});

    tensor<double,2> delWN_delGN_CC(2,2);
    delWN_delGN_CC= product(imetric_GNR,product(delWN_delGN, imetric_GNR, {{1, 0}}), {{0, 0}});




    // Compute reactions
    double dt_eta=-1/fric2;

    double enP11=dt_eta*delWP_delGP_CC(0,0);
    double enP12=dt_eta*delWP_delGP_CC(0,1);
    double enP22=dt_eta*delWP_delGP_CC(1,1);

    double enN11=dt_eta*delWN_delGN_CC(0,0);
    double enN12=dt_eta*delWN_delGN_CC(0,1);
    double enN22=dt_eta*delWN_delGN_CC(1,1);

    double enL1=dt_eta*delWL_delG1;
    double enL2=dt_eta*delWL_delG2;

    double enL3=dt_eta*delWL_delG3;
    double enL4=dt_eta*delWL_delG4;
    double enL5=dt_eta*delWL_delG5;
    double enL6=dt_eta*delWL_delG6;

    double enL12=dt_eta*delWL_delG12;
    double enL34=dt_eta*delWL_delG34;
    double enL56=dt_eta*delWL_delG56;

    //cout<<"check: "<<  " :"<<dt_eta*delWP_delGP0<<endl;


    // --------------------------------------------
    // [3] Fill rhs and matrix
    // --------------------------------------------
/*
    //Rhs
    Bk(all,0) = jacR * bf * (G11P0 + enP11);
    Bk(all,1) = jacR * bf * (G12P0 + enP12);
    Bk(all,2) = jacR * bf * (G22P0 + enP22);
    Bk(all,3) = jacR * bf * (G11N0 + enN11);
    Bk(all,4) = jacR * bf * (G12N0 + enN12);
    Bk(all,5) = jacR * bf * (G22N0 + enN22);

    Bk(all,6) = jacR * bf * (G1L0 + enL1);
    Bk(all,7) = jacR * bf * (G2L0 + enL2);

    Bk(all,8) = jacR * bf * (G3L0 + enL3);
    Bk(all,9) = jacR * bf * (G4L0 + enL4);
    Bk(all,10) = jacR * bf * (G5L0 + enL5);
    Bk(all,11) = jacR * bf * (G6L0 + enL6);

    Bk(all,12) = jacR * bf * (G12L0 + enL12);
    Bk(all,13) = jacR * bf * (G34L0 + enL34);
    Bk(all,14) = jacR * bf * (G56L0 + enL56);
    //Matrix
    Ak(all,0,all,0) = jacR * outer(bf, bf) ;//
    Ak(all,1,all,1) =jacR * outer(bf, bf);
    Ak(all,2,all,2) = jacR * outer(bf, bf);

    Ak(all,3,all,3) = jacR * outer(bf, bf);
    Ak(all,4,all,4) = jacR * outer(bf, bf);
    Ak(all,5,all,5) = jacR * outer(bf, bf);
    Ak(all,6,all,6) = jacR * outer(bf, bf);
    Ak(all,7,all,7) = jacR * outer(bf, bf);

    Ak(all,8,all,8) = jacR * outer(bf, bf);
    Ak(all,9,all,9) = jacR * outer(bf, bf);
    Ak(all,10,all,10) = jacR * outer(bf, bf);
    Ak(all,11,all,11) = jacR * outer(bf, bf);


    Ak(all,12,all,12) = jacR * outer(bf, bf);
    Ak(all,13,all,13) = jacR * outer(bf, bf);
    Ak(all,14,all,14) = jacR * outer(bf, bf);
*/

    // --------------------------------------------
    // [4] Extra
    // --------------------------------------------
    //cout<<"print: "<<G1L<< ": " <<G2L<< endl;
    // This is to pass info to the other problem


    double* auxiliary_b = &fillStr.paramStr->b_aux[(subFill.loc_elemID*subFill.cubaInfo.iPts+subFill.kPt)*18];



    if(forward_new>0.1)
    {
        if(timestep<time_target)
        {

            G11P =imetric_GPR(0,0);
            G12P = imetric_GPR(0,1);
            G22P = imetric_GPR(1,1);
            G11N = imetric_GNR(0,0);
            G12N= imetric_GNR(1,0);
            G22N = imetric_GNR(1,1);
            G1L = 1.0;
            G2L = 1.0;

            G3L = 1.0;
            G4L = 1.0;
            G5L = 1.0;
            G6L = 1.0;

            G12L = 0.0;
            G34L = 0.0;
            G56L = 0.0;
        }
        else
        {


            G11P = auxiliary_b[0]  + enP11;
            G12P = auxiliary_b[1]  + enP12;
            G22P = auxiliary_b[2]  + enP22;
            G11N = auxiliary_b[3]  + enN11;
            G12N = auxiliary_b[4]  + enN12;
            G22N= auxiliary_b[5]  + enN22;
            G1L = auxiliary_b[6]+enL1;
            G2L =auxiliary_b[7]+enL2;

            G3L = auxiliary_b[8]+enL3;
            G4L = auxiliary_b[9]+enL4;
            G5L = auxiliary_b[10]+enL5;
            G6L =auxiliary_b[11]+enL6;
            G12L = auxiliary_b[12]+enL12;
            G34L = auxiliary_b[13]+enL34;
            G56L =auxiliary_b[14]+enL56;

        }



    }

     tensor<double,2> GP_old(2,2);
     tensor<double,2> GN_old(2,2);
     tensor<double,2> GL1_old(2,2);
     tensor<double,2> GL2_old(2,2);
     tensor<double,2> GL3_old(2,2);




 if (forward_new<0.1)
 {
     tensor<double,2> GP_n(2,2);
     GP_n(0,0)=auxiliary_b[0];
     GP_n(0,1)=auxiliary_b[1];
     GP_n(1,0)=auxiliary_b[1];
     GP_n(1,1)=auxiliary_b[2];
     double jacGP_n       = sqrt(GP_n.det());
     tensor<double,2> iGP_n=GP_n.inv();

     tensor<double,2> GN_n(2,2);
     GN_n(0,0)=auxiliary_b[3];
     GN_n(0,1)=auxiliary_b[4];
     GN_n(1,0)=auxiliary_b[4];
     GN_n(1,1)=auxiliary_b[5];
     double jacGN_n       = sqrt(GN_n.det());
     tensor<double,2> iGN_n=GN_n.inv();

     tensor<double,2> GL1_n(2,2);
     GL1_n(0,0)=auxiliary_b[6];
     GL1_n(0,1)=auxiliary_b[12];
     GL1_n(1,0)=auxiliary_b[12];
     GL1_n(1,1)=auxiliary_b[7];
     double jacGL1_n       = sqrt(GL1_n.det());
     tensor<double,2> iGL1_n=GL1_n.inv();

     tensor<double,2> GL2_n(2,2);
     GL2_n(0,0)=auxiliary_b[8];
     GL2_n(0,1)=auxiliary_b[13];
     GL2_n(1,0)=auxiliary_b[13];
     GL2_n(1,1)=auxiliary_b[9];
     double jacGL2_n       = sqrt(GL2_n.det());
     tensor<double,2> iGL2_n=GL2_n.inv();

     tensor<double,2> GL3_n(2,2);
     GL3_n(0,0)=auxiliary_b[10];
     GL3_n(0,1)=auxiliary_b[14];
     GL3_n(1,0)=auxiliary_b[14];
     GL3_n(1,1)=auxiliary_b[11];
     double jacGL3_n       = sqrt(GL3_n.det());
     tensor<double,2> iGL3_n=GL3_n.inv();




//old

//old

     GP_old(0,0)=auxiliary_b[0];
     GP_old(0,1)=auxiliary_b[1];
     GP_old(1,0)=auxiliary_b[1];
     GP_old(1,1)=auxiliary_b[2];


     GN_old(0,0)=auxiliary_b[3];
     GN_old(0,1)=auxiliary_b[4];
     GN_old(1,0)=auxiliary_b[4];
     GN_old(1,1)=auxiliary_b[5];


     GL1_old(0,0)=auxiliary_b[6];
     GL1_old(0,1)=auxiliary_b[12];
     GL1_old(1,0)=auxiliary_b[12];
     GL1_old(1,1)=auxiliary_b[7];



     GL2_old(0,0)=auxiliary_b[8];
     GL2_old(0,1)=auxiliary_b[13];
     GL2_old(1,0)=auxiliary_b[13];
     GL2_old(1,1)=auxiliary_b[9];

     GL3_old(0,0)=auxiliary_b[10];
     GL3_old(0,1)=auxiliary_b[14];
     GL3_old(1,0)=auxiliary_b[14];
     GL3_old(1,1)=auxiliary_b[11];



     //nth step

     if(timestep<time_target)
     {

         G11P =imetric_GPR(0,0);
         G12P = imetric_GPR(0,1);
         G22P = imetric_GPR(1,1);
         G11N = imetric_GNR(0,0);
         G12N= imetric_GNR(1,0);
         G22N = imetric_GNR(1,1);
         G1L = 1.0;
         G2L = 1.0;

         G3L = 1.0;
         G4L = 1.0;
         G5L = 1.0;
         G6L = 1.0;

         G12L = 0.0;
         G34L = 0.0;
         G56L = 0.0;
     }
     else
     {

        //newton loop
         int nr_step=0;
         double GP_11_NR{0.0},GP_12_NR{0.0},GP_22_NR{0.0},GN_11_NR{0.0},GN_12_NR{0.0},GN_22_NR{0.0},GL1_11_NR{0.0},GL1_12_NR{0.0},GL1_22_NR{0.0},GL2_11_NR{0.0},GL2_12_NR{0.0},GL2_22_NR{0.0},GL3_11_NR{0.0},GL3_12_NR{0.0},GL3_22_NR{0.0};
         GP_11_NR= GP_n(0,0);
         GP_12_NR= GP_n(0,1);
         GP_22_NR= GP_n(1,1);

         GN_11_NR= GN_n(0,0);
         GN_12_NR= GN_n(0,1);
         GN_22_NR= GN_n(1,1);

         GL1_11_NR= GL1_n(0,0);
         GL1_12_NR= GL1_n(0,1);
         GL1_22_NR= GL1_n(1,1);

         GL2_11_NR= GL2_n(0,0);
         GL2_12_NR= GL2_n(0,1);
         GL2_22_NR= GL2_n(1,1);

         GL3_11_NR= GL3_n(0,0);
         GL3_12_NR= GL3_n(0,1);
         GL3_22_NR= GL3_n(1,1);

         double error=0.1;
         double error2=0.1;

         while ((error>0.0000001|| error2>0.0000001) && nr_step<30)
         {
              jacGP_n       = sqrt(GP_n.det());
              iGP_n=GP_n.inv();

             jacGN_n       = sqrt(GN_n.det());
             iGN_n=GN_n.inv();

              jacGL1_n       = sqrt(GL1_n.det());
              iGL1_n=GL1_n.inv();

             jacGL2_n       = sqrt(GL2_n.det());
             iGL2_n=GL2_n.inv();

             jacGL3_n       = sqrt(GL3_n.det());
             iGL3_n=GL3_n.inv();

             //This is to pass info to the other problem
             tensor<double,2> delI1_gp_delGP_n(pDim,pDim);//del I1 /del G ^bR
             tensor<double,2> delI1_gn_delGN_n(pDim,pDim);//


             tensor<double,2> delI1L1_delGL1_n(pDim,pDim);//
             tensor<double,2> delI1L2_delGL2_n(pDim,pDim);//
             tensor<double,2> delI1L3_delGL3_n(pDim,pDim);//

             delI1_gp_delGP_n=   (metric_gp)    ; //
             delI1_gn_delGN_n=   (metric_gn)    ; //_Br

             delI1L1_delGL1_n=(CL1);
             delI1L2_delGL2_n=(CL2);
             delI1L3_delGL3_n=(CL3);

             tensor<double,2> delI3_gp_delGP_n(pDim,pDim);//del I1 /del G ^bR
             tensor<double,2> delI3_gn_delGN_n(pDim,pDim);//

             tensor<double,2> delI3L1_delGL1_n(pDim,pDim);//
             tensor<double,2> delI3L2_delGL2_n(pDim,pDim);//
             tensor<double,2> delI3L3_delGL3_n(pDim,pDim);//


             delI3_gp_delGP_n=jac_gp*jac_gp*jacGP_n*jacGP_n*iGP_n; //_ab
             delI3_gn_delGN_n=jac_gn*jac_gn*jacGN_n*jacGN_n*iGN_n;

             delI3L1_delGL1_n=(jacCL1*jacCL1*jacGL1_n*jacGL1_n*iGL1_n);
             delI3L2_delGL2_n=(jacCL2*jacCL2*jacGL2_n*jacGL2_n*iGL2_n);
             delI3L3_delGL3_n=(jacCL3*jacCL3*jacGL3_n*jacGL3_n*iGL3_n);

             tensor<double,4> d_delI3_gp_delGP_n(pDim,pDim,pDim,pDim);//ab_cd
             d_delI3_gp_delGP_n=jac_gp*jac_gp*jacGP_n*jacGP_n*outer(iGP_n,iGP_n)-jac_gp*jac_gp*jacGP_n*jacGP_n*outer(iGP_n,iGP_n).transpose({0,2,3,1});//pqcd
             tensor<double,4> d_delI3_gn_delGN_n(pDim,pDim,pDim,pDim);//del I1 /del G ^bR
             d_delI3_gn_delGN_n=jac_gn*jac_gn*jacGN_n*jacGN_n*outer(iGN_n,iGN_n)-jac_gn*jac_gn*jacGN_n*jacGN_n*outer(iGN_n,iGN_n).transpose({0,2,3,1});
             tensor<double,4> d_delI3L1_delGL1_n(pDim,pDim,pDim,pDim);//del I1 /del G ^bR
             d_delI3L1_delGL1_n=jacCL1*jacCL1*jacGL1_n*jacGL1_n*outer(iGL1_n,iGL1_n)-jacCL1*jacCL1*jacGL1_n*jacGL1_n*outer(iGL1_n,iGL1_n).transpose({0,2,3,1});
             tensor<double,4> d_delI3L2_delGL2_n(pDim,pDim,pDim,pDim);//del I1 /del G ^bR
             d_delI3L2_delGL2_n=jacCL2*jacCL2*jacGL2_n*jacGL2_n*outer(iGL2_n,iGL2_n)-jacCL2*jacCL2*jacGL2_n*jacGL2_n*outer(iGL2_n,iGL2_n).transpose({0,2,3,1});
             tensor<double,4> d_delI3L3_delGL3_n(pDim,pDim,pDim,pDim);//del I1 /del G ^bR
             d_delI3L3_delGL3_n=jacCL3*jacCL3*jacGL3_n*jacGL3_n*outer(iGL3_n,iGL3_n)-jacCL3*jacCL3*jacGL3_n*jacGL3_n*outer(iGL3_n,iGL3_n).transpose({0,2,3,1});



             tensor<double,2> delWP_delGP_n(pDim,pDim);//del j+ /del G ^bR
             tensor<double,2> delWN_delGN_n(pDim,pDim);//

             tensor<double,2> delWL_delGL1_n(pDim,pDim);//
             tensor<double,2> delWL_delGL2_n(pDim,pDim);//
             tensor<double,2> delWL_delGL3_n(pDim,pDim);//

             tensor<double,4> d_delWP_delGP_n(pDim,pDim,pDim,pDim);//Jabcd* Fcd
             tensor<double,4> d_delWN_delGN_n(pDim,pDim,pDim,pDim);//

             tensor<double,4> d_delWL_delGL1_n(pDim,pDim,pDim,pDim);//
             tensor<double,4> d_delWL_delGL2_n(pDim,pDim,pDim,pDim);//
             tensor<double,4> d_delWL_delGL3_n(pDim,pDim,pDim,pDim);//


             double I3_gp_n=jac_gp*jac_gp* jacGP_n*jacGP_n;
             double I3_gn_n=jac_gn*jac_gn*jacGN_n*jacGN_n;


             double I3L1_n=jacCL1*jacCL1*jacGL1_n*jacGL1_n;
             double I3L2_n=jacCL2*jacCL2*jacGL2_n*jacGL2_n;
             double I3L3_n=jacCL3*jacCL3*jacGL3_n*jacGL3_n;



             double I3_gp_SR_n=(jacGP_n)*jac_gp;
             double I3_gn_SR_n=(jacGN_n)*jac_gn;

             double I3L1_SR_n=(jacGL1_n)*jacCL1;
             double I3L2_SR_n=(jacGL2_n)*jacCL2;
             double I3L3_SR_n=(jacGL3_n)*jacCL3;


             delWP_delGP_n=0.5*(lambda*log(I3_gp_SR_n)-mu)/I3_gp_n*delI3_gp_delGP_n+0.5*mu*delI1_gp_delGP_n;
             delWN_delGN_n=0.5*(lambda*log(I3_gn_SR_n)-mu)/I3_gn_n*delI3_gn_delGN_n+0.5*mu*delI1_gn_delGN_n;

             delWL_delGL1_n=0.5*(lambda*log(I3L1_SR_n)-mu)/I3L1_n*delI3L1_delGL1_n+0.5*mu*delI1L1_delGL1_n;
             delWL_delGL2_n=0.5*(lambda*log(I3L2_SR_n)-mu)/I3L2_n*delI3L2_delGL2_n+0.5*mu*delI1L2_delGL2_n;
             delWL_delGL3_n=0.5*(lambda*log(I3L3_SR_n)-mu)/I3L3_n*delI3L3_delGL3_n+0.5*mu*delI1L3_delGL3_n;//pq
          //pqcd
             d_delWP_delGP_n=(-2*lambda*log(I3_gp_SR_n)+2*mu+lambda)/(4*I3_gp_n*I3_gp_n)*outer(delI3_gp_delGP_n,delI3_gp_delGP_n)+0.5*(lambda*log(I3_gp_SR_n)-mu)/I3_gp_n*d_delI3_gp_delGP_n;
             d_delWN_delGN_n=(-2*lambda*log(I3_gn_SR_n)+2*mu+lambda)/(4*I3_gn_n*I3_gp_n)*outer(delI3_gn_delGN_n,delI3_gn_delGN_n)+0.5*(lambda*log(I3_gn_SR_n)-mu)/I3_gn_n*d_delI3_gn_delGN_n;

             d_delWL_delGL1_n=(-2*lambda*log(I3L1_SR_n)+2*mu+lambda)/(4*I3L1_n*I3L1_n)*outer(delI3L1_delGL1_n,delI3L1_delGL1_n)+0.5*(lambda*log(I3L1_SR_n)-mu)/I3L1_n*d_delI3L1_delGL1_n;
             d_delWL_delGL2_n=(-2*lambda*log(I3L2_SR_n)+2*mu+lambda)/(4*I3L2_n*I3L2_n)*outer(delI3L2_delGL2_n,delI3L2_delGL2_n)+0.5*(lambda*log(I3L2_SR_n)-mu)/I3L2_n*d_delI3L2_delGL2_n;
             d_delWL_delGL3_n=(-2*lambda*log(I3L3_SR_n)+2*mu+lambda)/(4*I3L3_n*I3L3_n)*outer(delI3L3_delGL3_n,delI3L3_delGL3_n)+0.5*(lambda*log(I3L3_SR_n)-mu)/I3L3_n*d_delI3L3_delGL3_n;


           //  tensor<double,2> delWP_delGP_CC(2,2);
           //  delWP_delGP_CC= product(imetric_GPR,product(delWP_delGP, imetric_GPR, {{1, 0}}), {{0, 0}});

           //  tensor<double,2> delWN_delGN_CC(2,2);
           //  delWN_delGN_CC= product(imetric_GNR,product(delWN_delGN, imetric_GNR, {{1, 0}}), {{0, 0}});

             tensor<double,2> delWP_delGP_n_CC(pDim,pDim);//del j+ /del G ^bR
             tensor<double,2> delWN_delGN_n_CC(pDim,pDim);//


             tensor<double,4> d_delWP_delGP_n_CC(pDim,pDim,pDim,pDim);//Jabcd* Fcd
             tensor<double,4> d_delWN_delGN_n_CC(pDim,pDim,pDim,pDim);//



             delWP_delGP_n_CC=product(imetric_GPR,product(delWP_delGP_n, imetric_GPR, {{1, 0}}), {{1, 0}});
             delWN_delGN_n_CC=product(imetric_GPR,product(delWN_delGN_n, imetric_GPR, {{1, 0}}), {{1, 0}});

             d_delWP_delGP_n_CC=product(imetric_GPR,product(d_delWP_delGP_n, imetric_GPR, {{1, 0}}), {{1, 0}}).transpose({0,3,1,2});//pqcd, ap*pqcd*qb, acdb
             d_delWN_delGN_n_CC=product(imetric_GNR,product(d_delWN_delGN_n, imetric_GNR, {{1, 0}}), {{1, 0}}).transpose({0,3,1,2});//ab_cd


             //FORCE
             double F_GP_11_NR{0.0},F_GP_12_NR{0.0},F_GP_22_NR{0.0},F_GN_11_NR{0.0},F_GN_12_NR{0.0},F_GN_22_NR{0.0},F_GL1_11_NR{0.0},F_GL1_12_NR{0.0},F_GL1_22_NR{0.0},F_GL2_11_NR{0.0},F_GL2_12_NR{0.0},F_GL2_22_NR{0.0},F_GL3_11_NR{0.0},F_GL3_12_NR{0.0},F_GL3_22_NR{0.0};

             F_GP_11_NR=GP_11_NR-auxiliary_b[0]-dt_eta*delWP_delGP_n_CC(0,0);
             F_GP_12_NR=GP_12_NR-auxiliary_b[1]-dt_eta*delWP_delGP_n_CC(0,1);
             F_GP_22_NR=GP_22_NR-auxiliary_b[2]-dt_eta*delWP_delGP_n_CC(1,1);

             F_GN_11_NR=GN_11_NR-auxiliary_b[3]-dt_eta*delWN_delGN_n_CC(0,0);
             F_GN_12_NR=GN_12_NR-auxiliary_b[4]-dt_eta*delWN_delGN_n_CC(0,1);
             F_GN_22_NR=GN_22_NR-auxiliary_b[5]-dt_eta*delWN_delGN_n_CC(1,1);

             F_GL1_11_NR=GL1_11_NR-auxiliary_b[6]-dt_eta*delWL_delGL1_n(0,0);
             F_GL1_12_NR=GL1_12_NR-auxiliary_b[12]-dt_eta*delWL_delGL1_n(0,1);
             F_GL1_22_NR=GL1_22_NR-auxiliary_b[7]-dt_eta*delWL_delGL1_n(1,1);

             F_GL2_11_NR=GL2_11_NR-auxiliary_b[8]-dt_eta*delWL_delGL2_n(0,0);
             F_GL2_12_NR=GL2_12_NR-auxiliary_b[13]-dt_eta*delWL_delGL2_n(0,1);
             F_GL2_22_NR=GL2_22_NR-auxiliary_b[9]-dt_eta*delWL_delGL2_n(1,1);

             F_GL3_11_NR=GL3_11_NR-auxiliary_b[10]-dt_eta*delWL_delGL3_n(0,0);
             F_GL3_12_NR=GL3_12_NR-auxiliary_b[14]-dt_eta*delWL_delGL3_n(0,1);
             F_GL3_22_NR=GL3_22_NR-auxiliary_b[11]-dt_eta*delWL_delGL3_n(1,1);



             std::vector<double> values = {
                 abs(F_GP_11_NR), abs(F_GP_12_NR), abs(F_GP_22_NR), abs(F_GN_11_NR), abs(F_GN_12_NR), abs(F_GN_22_NR),
                 abs(F_GL1_11_NR), abs(F_GL1_12_NR), abs(F_GL1_22_NR), abs(F_GL2_11_NR), abs(F_GL2_12_NR), abs(F_GL2_22_NR),
                 abs(F_GL3_11_NR), abs(F_GL3_12_NR), abs(F_GL3_22_NR)
             };

             error=*std::max_element(values.begin(), values.end());

             tensor<double,1> GP_voigt(3);//del j+ /del G ^bR
             GP_voigt(0)=GP_11_NR;
             GP_voigt(1)=GP_12_NR;
             GP_voigt(2)=GP_22_NR;
             tensor<double,1> GP_F_voigt(3);//del j+ /del G ^bR

             GP_F_voigt(0)=F_GP_11_NR;
             GP_F_voigt(1)=F_GP_12_NR;
             GP_F_voigt(2)=F_GP_22_NR;

             tensor<double,2> GP_jac_voigt(3,3);//del j+ /del G ^bR
             double  J_GP_11_NR_11=1.0-dt_eta*d_delWP_delGP_n_CC(0,0,0,0); //d_delWP_delGP_n_CC
             double  J_GP_11_NR_12=0.0-dt_eta*d_delWP_delGP_n_CC(0,0,0,1);
             double  J_GP_11_NR_22=0.0-dt_eta*d_delWP_delGP_n_CC(0,0,1,1);

             double  J_GP_12_NR_11=0.0-dt_eta*d_delWP_delGP_n_CC(0,1,0,0);
             double  J_GP_12_NR_12=1.0-dt_eta*d_delWP_delGP_n_CC(0,1,0,1);
             double  J_GP_12_NR_22=0.0-dt_eta*d_delWP_delGP_n_CC(0,1,1,1);


             double J_GP_22_NR_11=0.0-dt_eta*d_delWP_delGP_n_CC(1,1,0,0);
             double J_GP_22_NR_12=0.0-dt_eta*d_delWP_delGP_n_CC(1,1,0,1);
             double J_GP_22_NR_22=1.0-dt_eta*d_delWP_delGP_n_CC(1,1,1,1);

             GP_jac_voigt(0,0)=J_GP_11_NR_11;
             GP_jac_voigt(0,1)=J_GP_11_NR_12;
             GP_jac_voigt(0,2)=J_GP_11_NR_22;

             GP_jac_voigt(1,0)=J_GP_12_NR_11;
             GP_jac_voigt(1,1)=J_GP_12_NR_12;
             GP_jac_voigt(1,2)=J_GP_12_NR_22;

             GP_jac_voigt(2,0)=J_GP_22_NR_11;
             GP_jac_voigt(2,1)=J_GP_22_NR_12;
             GP_jac_voigt(2,2)=J_GP_22_NR_22;

             GP_voigt=GP_voigt-GP_jac_voigt.inv()*GP_F_voigt;

             GP_11_NR=GP_voigt(0);
             GP_12_NR=GP_voigt(1);
             GP_22_NR= GP_voigt(2);

             //GN

             tensor<double,1> GN_voigt(3);//del j+ /del G ^bR
             GN_voigt(0)=GN_11_NR;
             GN_voigt(1)=GN_12_NR;
             GN_voigt(2)=GN_22_NR;
             tensor<double,1> GN_F_voigt(3);//del j+ /del G ^bR

             GN_F_voigt(0)=F_GN_11_NR;
             GN_F_voigt(1)=F_GN_12_NR;
             GN_F_voigt(2)=F_GN_22_NR;

             tensor<double,2> GN_jac_voigt(3,3);//del j+ /del G ^bR
             double  J_GN_11_NR_11=1.0-dt_eta*d_delWN_delGN_n_CC(0,0,0,0);
             double  J_GN_11_NR_12=0.0-dt_eta*d_delWN_delGN_n_CC(0,0,0,1);
             double  J_GN_11_NR_22=0.0-dt_eta*d_delWN_delGN_n_CC(0,0,1,1);

             double  J_GN_12_NR_11=0.0-dt_eta*d_delWN_delGN_n_CC(0,1,0,0);
             double  J_GN_12_NR_12=1.0-dt_eta*d_delWN_delGN_n_CC(0,1,0,1);
             double  J_GN_12_NR_22=0.0-dt_eta*d_delWN_delGN_n_CC(0,1,1,1);


             double J_GN_22_NR_11=0.0-dt_eta*d_delWN_delGN_n_CC(1,1,0,0);
             double J_GN_22_NR_12=0.0-dt_eta*d_delWN_delGN_n_CC(1,1,0,1);
             double J_GN_22_NR_22=1.0-dt_eta*d_delWN_delGN_n_CC(1,1,1,1);

             GN_jac_voigt(0,0)=J_GN_11_NR_11;
             GN_jac_voigt(0,1)=J_GN_11_NR_12;
             GN_jac_voigt(0,2)=J_GN_11_NR_22;

             GN_jac_voigt(1,0)=J_GN_12_NR_11;
             GN_jac_voigt(1,1)=J_GN_12_NR_12;
             GN_jac_voigt(1,2)=J_GN_12_NR_22;

             GN_jac_voigt(2,0)=J_GN_22_NR_11;
             GN_jac_voigt(2,1)=J_GN_22_NR_12;
             GN_jac_voigt(2,2)=J_GN_22_NR_22;

             GN_voigt=GN_voigt-GN_jac_voigt.inv()*GN_F_voigt;

             GN_11_NR=GN_voigt(0);
             GN_12_NR=GN_voigt(1);
             GN_22_NR= GN_voigt(2);

             //GL1

             tensor<double,1> GL1_voigt(3);//del j+ /del G ^bR
             GL1_voigt(0)=GL1_11_NR;
             GL1_voigt(1)=GL1_12_NR;
             GL1_voigt(2)=GL1_22_NR;
             tensor<double,1> GL1_F_voigt(3);//del j+ /del G ^bR

             GL1_F_voigt(0)=F_GL1_11_NR;
             GL1_F_voigt(1)=F_GL1_12_NR;
             GL1_F_voigt(2)=F_GL1_22_NR;

             tensor<double,2> GL1_jac_voigt(3,3);//del j+ /del G ^bR
             double  J_GL1_11_NR_11=1.0-dt_eta*d_delWL_delGL1_n(0,0,0,0);
             double  J_GL1_11_NR_12=0.0-dt_eta*d_delWL_delGL1_n(0,0,0,1);
             double  J_GL1_11_NR_22=0.0-dt_eta*d_delWL_delGL1_n(0,0,1,1);

             double  J_GL1_12_NR_11=0.0-dt_eta*d_delWL_delGL1_n(0,1,0,0);
             double  J_GL1_12_NR_12=1.0-dt_eta*d_delWL_delGL1_n(0,1,0,1);
             double  J_GL1_12_NR_22=0.0-dt_eta*d_delWL_delGL1_n(0,1,1,1);


             double J_GL1_22_NR_11=0.0-dt_eta*d_delWL_delGL1_n(1,1,0,0);
             double J_GL1_22_NR_12=0.0-dt_eta*d_delWL_delGL1_n(1,1,0,1);
             double J_GL1_22_NR_22=1.0-dt_eta*d_delWL_delGL1_n(1,1,1,1);

             GL1_jac_voigt(0,0)=J_GL1_11_NR_11;
             GL1_jac_voigt(0,1)=J_GL1_11_NR_12;
             GL1_jac_voigt(0,2)=J_GL1_11_NR_22;

             GL1_jac_voigt(1,0)=J_GL1_12_NR_11;
             GL1_jac_voigt(1,1)=J_GL1_12_NR_12;
             GL1_jac_voigt(1,2)=J_GL1_12_NR_22;

             GL1_jac_voigt(2,0)=J_GL1_22_NR_11;
             GL1_jac_voigt(2,1)=J_GL1_22_NR_12;
             GL1_jac_voigt(2,2)=J_GL1_22_NR_22;

             GL1_voigt=GL1_voigt-GL1_jac_voigt.inv()*GL1_F_voigt;

             GL1_11_NR=GL1_voigt(0);
             GL1_12_NR=GL1_voigt(1);
             GL1_22_NR= GL1_voigt(2);


             //GL2

             tensor<double,1> GL2_voigt(3);//del j+ /del G ^bR
             GL2_voigt(0)=GL2_11_NR;
             GL2_voigt(1)=GL2_12_NR;
             GL2_voigt(2)=GL2_22_NR;
             tensor<double,1> GL2_F_voigt(3);//del j+ /del G ^bR

             GL2_F_voigt(0)=F_GL2_11_NR;
             GL2_F_voigt(1)=F_GL2_12_NR;
             GL2_F_voigt(2)=F_GL2_22_NR;

             tensor<double,2> GL2_jac_voigt(3,3);//del j+ /del G ^bR
             double  J_GL2_11_NR_11=1.0-dt_eta*d_delWL_delGL2_n(0,0,0,0);
             double  J_GL2_11_NR_12=0.0-dt_eta*d_delWL_delGL2_n(0,0,0,1);
             double  J_GL2_11_NR_22=0.0-dt_eta*d_delWL_delGL2_n(0,0,1,1);

             double  J_GL2_12_NR_11=0.0-dt_eta*d_delWL_delGL2_n(0,1,0,0);
             double  J_GL2_12_NR_12=1.0-dt_eta*d_delWL_delGL2_n(0,1,0,1);
             double  J_GL2_12_NR_22=0.0-dt_eta*d_delWL_delGL2_n(0,1,1,1);


             double J_GL2_22_NR_11=0.0-dt_eta*d_delWL_delGL2_n(1,1,0,0);
             double J_GL2_22_NR_12=0.0-dt_eta*d_delWL_delGL2_n(1,1,0,1);
             double J_GL2_22_NR_22=1.0-dt_eta*d_delWL_delGL2_n(1,1,1,1);

             GL2_jac_voigt(0,0)=J_GL2_11_NR_11;
             GL2_jac_voigt(0,1)=J_GL2_11_NR_12;
             GL2_jac_voigt(0,2)=J_GL2_11_NR_22;

             GL2_jac_voigt(1,0)=J_GL2_12_NR_11;
             GL2_jac_voigt(1,1)=J_GL2_12_NR_12;
             GL2_jac_voigt(1,2)=J_GL2_12_NR_22;

             GL2_jac_voigt(2,0)=J_GL2_22_NR_11;
             GL2_jac_voigt(2,1)=J_GL2_22_NR_12;
             GL2_jac_voigt(2,2)=J_GL2_22_NR_22;

             GL2_voigt=GL2_voigt-GL2_jac_voigt.inv()*GL2_F_voigt;

             GL2_11_NR=GL2_voigt(0);
             GL2_12_NR=GL2_voigt(1);
             GL2_22_NR= GL2_voigt(2);

             //GL3

             tensor<double,1> GL3_voigt(3);//del j+ /del G ^bR
             GL3_voigt(0)=GL3_11_NR;
             GL3_voigt(1)=GL3_12_NR;
             GL3_voigt(2)=GL3_22_NR;
             tensor<double,1> GL3_F_voigt(3);//del j+ /del G ^bR

             GL3_F_voigt(0)=F_GL3_11_NR;
             GL3_F_voigt(1)=F_GL3_12_NR;
             GL3_F_voigt(2)=F_GL3_22_NR;

             tensor<double,2> GL3_jac_voigt(3,3);//del j+ /del G ^bR
             double  J_GL3_11_NR_11=1.0-dt_eta*d_delWL_delGL3_n(0,0,0,0);
             double  J_GL3_11_NR_12=0.0-dt_eta*d_delWL_delGL3_n(0,0,0,1);
             double  J_GL3_11_NR_22=0.0-dt_eta*d_delWL_delGL3_n(0,0,1,1);

             double  J_GL3_12_NR_11=0.0-dt_eta*d_delWL_delGL3_n(0,1,0,0);
             double  J_GL3_12_NR_12=1.0-dt_eta*d_delWL_delGL3_n(0,1,0,1);
             double  J_GL3_12_NR_22=0.0-dt_eta*d_delWL_delGL3_n(0,1,1,1);


             double J_GL3_22_NR_11=0.0-dt_eta*d_delWL_delGL3_n(1,1,0,0);
             double J_GL3_22_NR_12=0.0-dt_eta*d_delWL_delGL3_n(1,1,0,1);
             double J_GL3_22_NR_22=1.0-dt_eta*d_delWL_delGL3_n(1,1,1,1);

             GL3_jac_voigt(0,0)=J_GL3_11_NR_11;
             GL3_jac_voigt(0,1)=J_GL3_11_NR_12;
             GL3_jac_voigt(0,2)=J_GL3_11_NR_22;

             GL3_jac_voigt(1,0)=J_GL3_12_NR_11;
             GL3_jac_voigt(1,1)=J_GL3_12_NR_12;
             GL3_jac_voigt(1,2)=J_GL3_12_NR_22;

             GL3_jac_voigt(2,0)=J_GL3_22_NR_11;
             GL3_jac_voigt(2,1)=J_GL3_22_NR_12;
             GL3_jac_voigt(2,2)=J_GL3_22_NR_22;

             GL3_voigt=GL3_voigt-GL3_jac_voigt.inv()*GL3_F_voigt;

             GL3_11_NR=GL3_voigt(0);
             GL3_12_NR=GL3_voigt(1);
             GL3_22_NR=GL3_voigt(2);

             double X_GP_11_NR{0.0},X_GP_12_NR{0.0},X_GP_22_NR{0.0},X_GN_11_NR{0.0},X_GN_12_NR{0.0},X_GN_22_NR{0.0},X_GL1_11_NR{0.0},X_GL1_12_NR{0.0},X_GL1_22_NR{0.0},X_GL2_11_NR{0.0},X_GL2_12_NR{0.0},X_GL2_22_NR{0.0},X_GL3_11_NR{0.0},X_GL3_12_NR{0.0},X_GL3_22_NR{0.0};

             tensor<double,1>delX_GP(3);//del j+ /del G ^bR
             delX_GP=-GP_jac_voigt.inv()*GP_F_voigt;
             X_GP_11_NR=delX_GP(0);
             X_GP_12_NR=delX_GP(1);
             X_GP_22_NR=delX_GP(2);

             tensor<double,1>delX_GN(3);//del j+ /del G ^bR
             delX_GN=-GN_jac_voigt.inv()*GN_F_voigt;
             X_GN_11_NR=delX_GN(0);
             X_GN_12_NR=delX_GN(1);
             X_GN_22_NR=delX_GN(2);

             tensor<double,1>delX_GL1(3);//del j+ /del G ^bR
             delX_GL1=-GL1_jac_voigt.inv()*GL1_F_voigt;
             X_GL1_11_NR=delX_GL1(0);
             X_GL1_12_NR=delX_GL1(1);
             X_GL1_22_NR=delX_GL1(2);

             tensor<double,1>delX_GL2(3);//del j+ /del G ^bR
             delX_GL2=-GL2_jac_voigt.inv()*GL2_F_voigt;
             X_GL2_11_NR=delX_GL2(0);
             X_GL2_12_NR=delX_GL2(1);
             X_GL2_22_NR=delX_GL2(2);

             tensor<double,1>delX_GL3(3);//del j+ /del G ^bR
             delX_GL3=-GL3_jac_voigt.inv()*GL3_F_voigt;
             X_GL3_11_NR=delX_GL3(0);
             X_GL3_12_NR=delX_GL3(1);
             X_GL3_22_NR=delX_GL3(2);

         double tol3=0.0000000001;
        std::vector<double> values2 = {
                 abs(X_GP_11_NR/(GP_11_NR+tol3)), abs(X_GP_12_NR/(GP_12_NR+tol3)), abs(X_GP_22_NR/(GP_22_NR+tol3)), abs(X_GN_11_NR/(GN_11_NR+tol3)), abs(X_GN_12_NR/(GN_12_NR+tol3)), abs(X_GN_22_NR/(GN_22_NR+tol3)),
                 abs(X_GL1_11_NR/(GL1_11_NR+tol3)), abs(X_GL1_12_NR/(GL1_12_NR+tol3)), abs(X_GL1_22_NR/(GL1_22_NR+tol3)), abs(X_GL2_11_NR/(GL2_11_NR+tol3)), abs(X_GL2_12_NR/(GL2_12_NR+tol3)), abs(X_GL2_22_NR/(GL2_22_NR+tol3)),
                 abs(X_GL3_11_NR/(GL3_11_NR+tol3)), abs(X_GL3_12_NR/(GL3_12_NR+tol3)), abs(X_GL3_22_NR/(GL3_22_NR+tol3))
             };

             error2=*std::max_element(values2.begin(), values2.end());





             //UPDATE
             GP_n(0,0)=GP_11_NR;
             GP_n(0,1)=GP_12_NR;
             GP_n(1,0)=GP_12_NR;
             GP_n(1,1)=GP_22_NR;

             GN_n(0,0)=GN_11_NR;
             GN_n(0,1)=GN_12_NR;
             GN_n(1,0)=GN_12_NR;
             GN_n(1,1)=GN_22_NR;

             GL1_n(0,0)=GL1_11_NR;
             GL1_n(0,1)=GL1_12_NR;
             GL1_n(1,0)=GL1_12_NR;
             GL1_n(1,1)=GL1_22_NR;

             GL2_n(0,0)=GL2_11_NR;
             GL2_n(0,1)=GL2_12_NR;
             GL2_n(1,0)=GL2_12_NR;
             GL2_n(1,1)=GL2_22_NR;

             GL3_n(0,0)=GL3_11_NR;
             GL3_n(0,1)=GL3_12_NR;
             GL3_n(1,0)=GL3_12_NR;
             GL3_n(1,1)=GL3_22_NR;





           //  if(timestep>450)
           //  cout<<"nr_step"<< nr_step <<   "residual G11: "<<   error<<endl;
            /* cout<<"nr_step"<< nr_step <<   "residual G12: "<<   F_GP_12_NR<<endl;
             cout<<"nr_step"<< nr_step <<   "residual G22: "<<   F_GP_22_NR <<endl;

             cout<<"nr_step"<< nr_step <<   "residual Gn11: "<<   F_GN_11_NR<<endl;
             cout<<"nr_step"<< nr_step <<   "residual Gn12: "<<   F_GN_12_NR<<endl;
             cout<<"nr_step"<< nr_step <<   "residual Gn22: "<<   F_GN_22_NR <<endl;

             cout<<"nr_step"<< nr_step <<   "residual GL1_11: "<<   F_GL1_11_NR<<endl;
             cout<<"nr_step"<< nr_step <<   "residual GL1_12: "<<   F_GL1_12_NR<<endl;
             cout<<"nr_step"<< nr_step <<   "residual GL1_22: "<<   F_GL1_22_NR <<endl;

             cout<<"nr_step"<< nr_step <<   "residual GL2_11: "<<   F_GL2_11_NR<<endl;
             cout<<"nr_step"<< nr_step <<   "residual GL2_12: "<<   F_GL2_12_NR<<endl;
             cout<<"nr_step"<< nr_step <<   "residual GL2_22: "<<   F_GL2_22_NR <<endl;*/

             if(nr_step>25)
               {
               cout<<"error warning: "<<endl;
               }

             nr_step=nr_step+1;
         }

         G11P =GP_11_NR;
         G12P = GP_12_NR;
         G22P = GP_22_NR;
         G11N = GN_11_NR;
         G12N= GN_12_NR;
         G22N =GN_22_NR;

         G1L = GL1_11_NR;
         G2L = GL1_22_NR;

         G3L = GL2_11_NR;
         G4L = GL2_22_NR;
         G5L = GL3_11_NR;
         G6L = GL3_22_NR;

         G12L =GL1_12_NR;
         G34L = GL2_12_NR;
         G56L = GL3_12_NR;
     }


 }




    auxiliary_b[0] = G11P;
    auxiliary_b[1] = G12P;
    auxiliary_b[2] = G22P;
    auxiliary_b[3] = G11N;
    auxiliary_b[4] = G12N;
    auxiliary_b[5] = G22N;

    auxiliary_b[6] = G1L;
    auxiliary_b[7] = G2L;
    auxiliary_b[8] = G3L;
    auxiliary_b[9] = G4L;
    auxiliary_b[10] = G5L;
    auxiliary_b[11] = G6L;

    auxiliary_b[12] = G12L;
    auxiliary_b[13] = G34L;
    auxiliary_b[14] = G56L;


//new
     tensor<double,2> GP_new(2,2);
     GP_new(0,0)=auxiliary_b[0];
     GP_new(0,1)=auxiliary_b[1];
     GP_new(1,0)=auxiliary_b[1];
     GP_new(1,1)=auxiliary_b[2];


     tensor<double,2> GN_new(2,2);
     GN_new(0,0)=auxiliary_b[3];
     GN_new(0,1)=auxiliary_b[4];
     GN_new(1,0)=auxiliary_b[4];
     GN_new(1,1)=auxiliary_b[5];


     tensor<double,2> GL1_new(2,2);
     GL1_new(0,0)=auxiliary_b[6];
     GL1_new(0,1)=auxiliary_b[12];
     GL1_new(1,0)=auxiliary_b[12];
     GL1_new(1,1)=auxiliary_b[7];

     tensor<double,2> GL2_new(2,2);
     GL2_new(0,0)=auxiliary_b[8];
     GL2_new(0,1)=auxiliary_b[13];
     GL2_new(1,0)=auxiliary_b[13];
     GL2_new(1,1)=auxiliary_b[9];

     tensor<double,2> GL3_new(2,2);
     GL3_new(0,0)=auxiliary_b[10];
     GL3_new(0,1)=auxiliary_b[14];
     GL3_new(1,0)=auxiliary_b[14];
     GL3_new(1,1)=auxiliary_b[11];

    double jacPRR=auxiliary_a[76];
    double  jacNRR=auxiliary_a[77];
    double  jacLRR1=auxiliary_a[78];
    double jacLRR2=auxiliary_a[79];
    double jacLRR3=auxiliary_a[80];

//new    product(GP_new-GP_old,metric_GPR,{{0,0}})
         double diss11=  product(product(GP_new-GP_old,metric_GPR,{{1,0}}),product(metric_GPR,GP_new-GP_old,{{1,0}}),{{0,0},{1,1}})/(deltat*deltat)*jacPRR;
         double diss22=  product(product(GN_new-GN_old,metric_GNR,{{1,0}}),product(metric_GNR,GN_new-GN_old,{{1,0}}),{{0,0},{1,1}})/(deltat*deltat)*jacNRR;


         double diss31=  product(GL1_new-GL1_old,GL1_new-GL1_old,{{0,0},{1,1}})/(deltat*deltat)*jacLRR1;
         double diss32=  product(GL2_new-GL2_old,GL2_new-GL2_old,{{0,0},{1,1}})/(deltat*deltat)*jacLRR2;
         double diss33=  product(GL3_new-GL3_old,GL3_new-GL3_old,{{0,0},{1,1}})/(deltat*deltat)*jacLRR3;

    // Global integrals
  //  fillStr.addToGlobalIntegral("totBonds",jacR);

 auxiliary_b[15] =diss11/jacPRR;
 auxiliary_b[16] =diss22/jacNRR;
 auxiliary_b[17] =diss31/jacLRR1+diss32/jacLRR2+diss33/jacLRR3;




    fillStr.addToGlobalIntegral("Diss1",diss11);
    fillStr.addToGlobalIntegral("Diss2",diss22);
    fillStr.addToGlobalIntegral("Diss3",diss31+diss32+diss33);


}



void LS_ED(hiperlife::FillStructure& fillStr)
{
    using namespace hiperlife;
    using namespace hiperlife::Tensor;

    // --------------------------------------------
    // [0] Initialize
    // --------------------------------------------
    auto& subFill = (fillStr)["EDHand"];

    //Parameters
    int numDOFs = subFill.numDOFs;
    int nDim    = subFill.nDim;
    int pDim    = subFill.pDim;
    int eNN     = subFill.eNN;
    double deltat = fillStr.paramStr->dparam[2];

    //Coordinates and degrees of freedom
    wrapper<double,2> nborCoords(subFill.nborCoords.data(),eNN,nDim);


    //Basis functions and derivatives
    wrapper<double,1>    bf(subFill.nborBFsDers(0),eNN);

    //output
    wrapper<double,2>  Bk1(fillStr.Bk(0).data(),eNN,numDOFs);
    wrapper<double,4>  Ak1(fillStr.Ak(0,0).data(),eNN,numDOFs,eNN,numDOFs);


    // --------------------------------------------
    // [1] Model parameters


    // --------------------------------------------
    // [2] Auxiliary variables
    // --------------------------------------------
    //FIXME: this has to be loaded from some structure
    double* auxiliary_a = &fillStr.paramStr->a_aux[(subFill.loc_elemID*subFill.cubaInfo.iPts+subFill.kPt)*88];



    //Jacobian
    double jacR= auxiliary_a[19];
    double jac= auxiliary_a[20];


    double nxx= auxiliary_a[27];
    double nyy= auxiliary_a[28];
    double nzz= auxiliary_a[29];
    double height1= auxiliary_a[30];


    double EA= auxiliary_a[81];
    double EB= auxiliary_a[82];
    double EL= auxiliary_a[83];
    double PA= auxiliary_a[84];
    double PB= auxiliary_a[85];
    double PL= auxiliary_a[86];
    double jacC=auxiliary_a[87];

    double *auxiliary_b = &fillStr.paramStr->b_aux[(subFill.loc_elemID * subFill.cubaInfo.iPts + subFill.kPt) * 18];

    double DA= auxiliary_b[15];
    double DB= auxiliary_b[16];
    double DL= auxiliary_b[17];


    // --------------------------------------------
    // [3] Fill rhs and matrix
    // --------------------------------------------

//Rhs


    Bk1(all,0) = jacR * bf * (nxx);
    Bk1(all,1) = jacR * bf * (nyy);
    Bk1(all,2) = jacR * bf * (nzz);
    Bk1(all,3) = jacR * bf * (height1);

  Bk1(all,4) = jacR * bf * (EA);
    Bk1(all,5) = jacR * bf * (EB);
    Bk1(all,6) = jacR * bf * (EL);
    Bk1(all,7) = jacR * bf * (DA);
  Bk1(all,8) = jacR * bf * (DB);
    Bk1(all,9) = jacR * bf * (DL);
    Bk1(all,10) = jacR * bf * (PA);
    Bk1(all,11) = jacR * bf * (PB);
    Bk1(all,12) = jacR * bf * (PL);
    Bk1(all,13) = jacR * bf * (jacC);

    //Matrix
    Ak1(all,0,all,0) = jacR * outer(bf, bf) ;//
    Ak1(all,1,all,1) =jacR * outer(bf, bf);
    Ak1(all,2,all,2) = jacR * outer(bf, bf) ;//
    Ak1(all,3,all,3) =jacR * outer(bf, bf);
    Ak1(all,4,all,4) = jacR * outer(bf, bf) ;//
    Ak1(all,5,all,5) =jacR * outer(bf, bf);
    Ak1(all,6,all,6) = jacR * outer(bf, bf) ;//
    Ak1(all,7,all,7) =jacR * outer(bf, bf);
    Ak1(all,8,all,8) = jacR * outer(bf, bf) ;//
    Ak1(all,9,all,9) =jacR * outer(bf, bf);
    Ak1(all,10,all,10) = jacR * outer(bf, bf) ;//
    Ak1(all,11,all,11) =jacR * outer(bf, bf);
    Ak1(all,12,all,12) =jacR * outer(bf, bf);

    Ak1(all,13,all,13) =jacR * outer(bf, bf);





}

void LS_Rho(hiperlife::FillStructure& fillStr)
{
    using namespace hiperlife;
    using namespace hiperlife::Tensor;

    // --------------------------------------------
    // [0] Initialize
    // --------------------------------------------
    auto& subFill = (fillStr)["RhoDHand"];

    //Parameters
    int numDOFs = subFill.numDOFs;
    int nDim    = subFill.nDim;
    int pDim    = subFill.pDim;
    int eNN     = subFill.eNN;
    double deltat = fillStr.getRealParameter(MembParams::deltat);

    double k_p=fillStr.getRealParameter(MembParams::k_p)*deltat;
    double k_d=fillStr.getRealParameter(MembParams::k_d)*deltat;

    //Coordinates and degrees of freedom
    wrapper<double,2> nborCoords(subFill.nborCoords.data(),eNN,nDim);
    wrapper<double,2> nborDOFs0(subFill.nborDOFs0.data(),eNN,numDOFs);
    wrapper<double,2> nborDOFs(subFill.nborDOFs.data(),eNN,numDOFs);

    //Basis functions and derivatives
    wrapper<double,1>    bf(subFill.nborBFsDers(0),eNN);
   // tensor<double,2>   Dbf(subFill.nborBFsDers(1),eNN,pDim);
   // tensor<double,3>  DDbf(subFill.nborBFsDers(2),eNN,pDim,pDim);

    //output
    wrapper<double,2>  Bk2(fillStr.Bk(0).data(),eNN,numDOFs);
    wrapper<double,4>  Ak2(fillStr.Ak(0,0).data(),eNN,numDOFs,eNN,numDOFs);


    // --------------------------------------------
    // [1] Model parameters


    // --------------------------------------------
    // [2] Auxiliary variables
    // --------------------------------------------
    //FIXME: this has to be loaded from some structure
    double* auxiliary_a = &fillStr.paramStr->a_aux[(subFill.loc_elemID*subFill.cubaInfo.iPts+subFill.kPt)*88];



    //Jacobian
    double jacR= auxiliary_a[19];
    double jac= auxiliary_a[20];

    double f0 = fillStr.getRealParameter(MembParams::f0);

    double  tr_rodt_api= auxiliary_a[66];
    double  tr_rodt_bas= auxiliary_a[67];
    double  tr_rodt_lat1= auxiliary_a[68];
    double   tr_rodt_lat2= auxiliary_a[69];
    double  tr_rodt_lat3= auxiliary_a[70];

    double   jac_gp0= auxiliary_a[71];
    double  jac_gn0= auxiliary_a[72];
    double  jacCL1_n= auxiliary_a[73]*f0/3.0;
    double jacCL2_n= auxiliary_a[74]*f0/3.0;
    double jacCL3_n= auxiliary_a[75]*f0/3.0;

  //Concentration values
    double rho_api  = nborDOFs(all,0) * bf;
    double rho_bas  = nborDOFs(all,1) * bf;
    double rho_lat1  = nborDOFs(all,2) * bf;
    double rho_lat2 = nborDOFs(all,3) * bf;
    double rho_lat3 = nborDOFs(all,4) * bf;

    double rho0_api  = nborDOFs0(all,0) * bf;
    double rho0_bas  = nborDOFs0(all,1) * bf;
    double rho0_lat1  = nborDOFs0(all,2) * bf;
    double rho0_lat2 = nborDOFs0(all,3) * bf;
    double rho0_lat3 = nborDOFs0(all,4) * bf;








    // --------------------------------------------
    // [3] Fill rhs and matrix
    // --------------------------------------------

//Rhs


    Bk2(all,0) = jac_gp0 * bf * (rho0_api+k_p);
    Bk2(all,1) = jac_gn0 * bf * (rho0_bas+k_p);
    Bk2(all,2) = jacCL1_n * jacR*bf * (rho0_lat1+k_p);
    Bk2(all,3) = jacCL2_n * jacR*bf * (rho0_lat2+k_p);
    Bk2(all,4) = jacCL3_n * jacR*bf * (rho0_lat3+k_p);

    //Matrix
    Ak2(all,0,all,0) = jac_gp0 * (outer(bf, bf) * (1.0 + tr_rodt_api+k_d));
    Ak2(all,1,all,1) = jac_gn0 * (outer(bf, bf) * (1.0 + tr_rodt_bas+k_d));
    Ak2(all,2,all,2) = jacCL1_n *jacR* (outer(bf, bf) * (1.0 + tr_rodt_lat1+k_d));
    Ak2(all,3,all,3) = jacCL2_n * jacR*(outer(bf, bf) * (1.0 + tr_rodt_lat2+k_d));
    Ak2(all,4,all,4) = jacCL3_n *jacR* (outer(bf, bf) * (1.0 + tr_rodt_lat3+k_d));

double* auxiliary_c = &fillStr.paramStr->c_aux[(subFill.loc_elemID*subFill.cubaInfo.iPts+subFill.kPt)*9];

auxiliary_c[0]=rho_api;
auxiliary_c[1]=rho_bas;
auxiliary_c[2]=rho_lat1;
auxiliary_c[3]=rho_lat2;
auxiliary_c[4]=rho_lat3;


}


void LS_node_end(hiperlife::FillStructure& fillStr)
{
    using namespace hiperlife;
    using namespace hiperlife::Tensor;

    // --------------------------------------------
    // [0] Initialize
    // --------------------------------------------
    auto& subFill = (fillStr)["xyzDHand"];

    //Parameters
    int numDOFs = subFill.numDOFs;
    int nDim    = subFill.nDim;
    int pDim    = subFill.pDim;
    int eNN     = subFill.eNN;
    double deltat = fillStr.paramStr->dparam[2];



    //Coordinates and degrees of freedom
    wrapper<double,2> nborCoords(subFill.nborCoords.data(),eNN,nDim);
    wrapper<double,2> nborDOFs0(subFill.nborDOFs0.data(),eNN,numDOFs);
    wrapper<double,2> nborDOFs(subFill.nborDOFs.data(),eNN,numDOFs);

    //Basis functions and derivatives
    wrapper<double,1>    bf(subFill.nborBFsDers(0),eNN);
   // tensor<double,2>   Dbf(subFill.nborBFsDers(1),eNN,pDim);
   // tensor<double,3>  DDbf(subFill.nborBFsDers(2),eNN,pDim,pDim);

    //output

    // --------------------------------------------
    // [1] Model parameters


    // --------------------------------------------
    // [2] Auxiliary variables
    // --------------------------------------------
    //FIXME: this has to be loaded from some structure
    double* auxiliary_a = &fillStr.paramStr->a_aux[(subFill.loc_elemID*subFill.cubaInfo.iPts+subFill.kPt)*88];



    //Jacobian
    double jacR= auxiliary_a[19];
    double jac= auxiliary_a[20];









    // --------------------------------------------
    // [3] Fill rhs and matrix
    // --------------------------------------------

//Rhs

    wrapper<double,1>      x(&auxiliary_a[0],3);


double* auxiliary_c = &fillStr.paramStr->c_aux[(subFill.loc_elemID*subFill.cubaInfo.iPts+subFill.kPt)*9];

auxiliary_c[5]=x(0);
auxiliary_c[6]=x(1);
auxiliary_c[7]=x(2);
auxiliary_c[8]=jac;


}
