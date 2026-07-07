

//#include "Amesos.h"


/// hiperlife headers
#include "hl_FillStructure.h"
#include "hl_Geometry.h"
#include "hl_SurfLagrParam.h"
#include "hl_Tensor.h"
//#include "hl_LinearSolver_Direct_Amesos2.h"
#include <hl_LinearSolver_Direct_MUMPS.h>
#include <fstream>

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


    tensor<double, 2,false> nborDOFs0(subFill.nborDOFs0.data(), eNN, numDOFs);
    tensor<double, 2,false> nborDOFs(subFill.nborDOFs.data(), eNN, numDOFs);
    tensor<double, 2,false> nborCoords(subFill.nborCoords.data(), eNN, nDim);
    tensor<double,2,false> nborAuxF(subFill.nborAuxF.data(),eNN,numAuxF);

    tensor<double, 1,false> bf(subFill.nborBFs(), eNN);
    tensor<double, 2,false> Dbf(subFill.nborBFsGrads(), eNN, pDim);
    tensor<double, 3,false> DDbf(subFill.nborBFsHess(), eNN, pDim, pDim);

    //[1.3] Global constraints
    auto &g_subFill = (fillStr)["gloDHand"];
    int g_numDOFs = g_subFill.numDOFs;

    tensor<double, 1,false> gDOFs(g_subFill.nborDOFs.data(), g_numDOFs);
    double pressure = gDOFs(0);
    tensor<double, 1> F = gDOFs(range(1, 3));
    tensor<double, 1> L = gDOFs(range(4, 6));

    //[1.4] Parameters
    double deltat = fillStr.getRealParameter(MembParams::deltat);
    double fric = fillStr.getRealParameter(MembParams::fric);
    double factor = fillStr.getRealParameter(MembParams::factor);

    double thick = fillStr.getRealParameter(MembParams::thick);
    double Kconf = fillStr.getRealParameter(MembParams::Kconf) * deltat;
    int timestep = fillStr.getIntParameter(MembParams::timestep);//timestep

 double R    = fillStr.getRealParameter(MembParams::R);
 double width   = fillStr.getRealParameter(MembParams::width);




     double force=fillStr.getRealParameter(MembParams::force)*deltat;
    double fric_fact = fillStr.getRealParameter(MembParams::fric_fact);

    double fric2 = fillStr.getRealParameter(MembParams::fric2);
    double height_in = fillStr.getRealParameter(MembParams::height_in);
    int control_fric =  fillStr.getIntParameter(MembParams::control_fric);


    double Lagrangian{};

    double kspr= fillStr.getRealParameter(MembParams::kspr) * deltat;
    double sp_gap= fillStr.getRealParameter(MembParams::sp_gap);

    double spring= fillStr.getRealParameter(MembParams::spring);
 int gap = fillStr.getIntParameter(MembParams::gap);

 int vnstep = fillStr.getIntParameter(MembParams::vnstep);


    int fric_start = gap+vnstep;

    double gamma_minus=fillStr.getRealParameter(MembParams::gamma_minus) * deltat;
    double gamma_plus=fillStr.getRealParameter(MembParams::gamma_plus) * deltat;
    double fact_elastic=fillStr.getRealParameter(MembParams::fact_elastic);
 double gamma=fillStr.getRealParameter(MembParams::gamma) * deltat;



 double lambda   = fillStr.getRealParameter(MembParams::lambda)*deltat;
 double mu   = fillStr.getRealParameter(MembParams::mu)*deltat;






    //OUTPUTS
    tensor<double, 2,false> Bp(fillStr.Bk(0).data(), eNN, numDOFs);
    tensor<double, 1,false> Bg(fillStr.Bk(1).data(), g_numDOFs);

    tensor<double, 4,false> App(fillStr.Ak(0, 0).data(), eNN, numDOFs, eNN, numDOFs);
    tensor<double, 3,false> Apg(fillStr.Ak(0, 1).data(), eNN, numDOFs, g_numDOFs);
    tensor<double, 3,false> Agp(fillStr.Ak(1, 0).data(), g_numDOFs, eNN, numDOFs);

     double eps=0.0000001;
    //------------------------------------------------------------------
    // [2] Compute variables
    // [2.1] Compute geometry in the reference
    tensor<double,1> xRv(nDim);
    tensor<double,1> xRu(nDim);
    tensor<double,1> xRuu(nDim);
    tensor<double,1> xRuv(nDim);
    tensor<double,1> xRvv(nDim);
    tensor<double,1> normalR(nDim);
    tensor<double,2> metricR(pDim,pDim);
    tensor<double,2> curvatureR(pDim,pDim);
    tensor<double,3> cristSymR(pDim,pDim,pDim);
    SurfLagrParam::ChristoffelSymbols(cristSymR, curvatureR, metricR, normalR, xRu, xRv, xRuu, xRvv, xRuv, eNN, nborCoords, Dbf, DDbf);

      if(timestep<1)
        {
         curvatureR(0,0)=curvatureR(0,0)+eps;
         curvatureR(1,1)=curvatureR(1,1)+eps;

        }


    //[2.1] Previous time-step
    tensor<double,1> xu_n(nDim);
    tensor<double,1> xv_n(nDim);
    tensor<double,1> normal_n(nDim);
    tensor<double,2> metric_n(pDim,pDim);
    SurfLagrParam::Normal(normal_n, metric_n, xu_n, xv_n, eNN, nborDOFs0, Dbf);

    tensor<double,1> xuu_n(nDim);
    tensor<double,1> xuv_n(nDim);
    tensor<double,1> xvv_n(nDim);
    tensor<double,2> curva_n(pDim,pDim);
    tensor<double,3> christsym_n(pDim,pDim,pDim);
    SurfLagrParam::ChristoffelSymbols(christsym_n, curva_n, metric_n, normal_n, xu_n, xv_n, xuu_n, xvv_n, xuv_n, eNN, nborDOFs0, Dbf, DDbf);


     if(timestep<1)
        {
         curva_n(0,0)=curva_n(0,0)+eps;
         curva_n(1,1)=curva_n(1,1)+eps;

        }

    tensor<double,1> x_n = bf * nborDOFs0;
    tensor<double,2> imetric_n = metric_n.inv();

    tensor<double,1> xR = bf * nborCoords;
    tensor<double,2> imetricR = metricR.inv();
    double jacR = sqrt(metricR.det());
    double jac_n = sqrt(metric_n.det());
    double xnormal_n = x_n * normal_n;
    tensor<double,2>    icurvatureR= curvatureR.inv();




    //[2.2] Current time-step
    tensor<double,1> xu(nDim);
    tensor<double,1> xv(nDim);
    tensor<double,1> xuu(nDim);
    tensor<double,1> xuv(nDim);
    tensor<double,1> xvv(nDim);
    tensor<double,1> normal(nDim);
    tensor<double,2> metric(pDim,pDim);
    tensor<double,2> curva(pDim,pDim);
    tensor<double,3> christsym(pDim,pDim,pDim);
    SurfLagrParam::ChristoffelSymbols(christsym, curva, metric, normal, xu, xv, xuu, xvv, xuv, eNN, nborDOFs, Dbf, DDbf);


     if(timestep<1)
        {
         curva(0,0)=curva(0,0)+eps;
         curva(1,1)=curva(1,1)+eps;

        }


    tensor<double,1> x =  bf * nborDOFs;


    tensor<double,2> imetric = metric.inv();
    double jac       = sqrt(metric.det());
    double meancurva = product(imetric,curva,{{0,0},{1,1}});


   double kspr_out=kspr/1.0;
    double rad=sqrt(xR(0)*xR(0)+xR(1)*xR(1));


    if(rad<R-sp_gap)
    {
        kspr=0.0;
    }



 if (timestep > fric_start+control_fric)
 {
  fric = fric * fric_fact + (fric * fric_fact - fric) * tanh(width* (thick - x_n(2)));
 }




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

    //------------------------------------------------------------------
    //------------------------------------------------------------------
    // [3.1] First derivatives

    tensor<double,4> d_metric(eNN,numDOFs,pDim,pDim);
    SurfLagrParam::d_Metric(d_metric,Dbf,xu,xv);

    tensor<double,4> d_metric_CC  = product(product(d_metric,imetric,{{2,0}}),imetric,{{2,0}});
    tensor<double,4> d_metric_CnCn = product(product(d_metric,imetric_n,{{2,0}}),imetric_n,{{2,0}});

    tensor<double,4> d_imetric = -1.0 * d_metric_CC;

    tensor<double,2> d_jac(eNN,numDOFs);
    SurfLagrParam::d_Jac(d_jac,Dbf,xu,xv,normal);


    tensor<double,3> d_normal(eNN,numDOFs,nDim);
    SurfLagrParam::d_Normal(d_normal,Dbf,xu,xv,normal,jac,d_jac);

    tensor<double,4> d_curva(eNN,numDOFs,pDim,pDim);
    SurfLagrParam::d_Curva(d_curva,DDbf,xuu,xuv,xvv,normal,d_normal);

    tensor<double,2> d_meancurva = product(imetric,d_curva,{{0,2},{1,3}}) + product(curva,d_imetric,{{0,2},{1,3}});



    //------------------------------------------------------------------
    // [4] second derivatives

    tensor<double,6> dd_metric(eNN,numDOFs,eNN,numDOFs,pDim,pDim);
    SurfLagrParam::dd_Metric(dd_metric,Dbf);

    tensor<double,6> dd_metric_CC  = product(product(dd_metric,imetric,{{4,0}}),imetric,{{4,0}});

    tensor<double,6> dd_imetric = -1.0 * dd_metric_CC - product(product(d_metric,d_imetric,{{2,2}}).transpose({0,1,3,4,5,2}),imetric,{{5,0}}) - product(product(d_metric,imetric,{{2,0}}),d_imetric,{{2,2}}).transpose({0,1,3,4,2,5});

    tensor<double,4> dd_jac(eNN,numDOFs,eNN,numDOFs);
    SurfLagrParam::dd_Jac(dd_jac,Dbf,xu,xv,normal,d_normal);

    tensor<double,5> dd_normal(eNN,numDOFs,eNN,numDOFs,nDim);
    SurfLagrParam::dd_Normal (dd_normal, Dbf, jac, d_jac, dd_jac, normal, d_normal);

    tensor<double,6> dd_curva(eNN,numDOFs,eNN,numDOFs,pDim,pDim);
    SurfLagrParam::dd_Curva(dd_curva, DDbf, xuu, xuv, xvv, d_normal, dd_normal);


    tensor<double,4> dd_meancurva = product(imetric,dd_curva,{{0,4},{1,5}}) + product(d_imetric,d_curva,{{2,2},{3,3}})+ product(d_curva,d_imetric,{{2,2},{3,3}}) + product(curva,dd_imetric,{{0,4},{1,5}});



   // tensor<double, 4> dd_meancurva = product(imetric, dd_curva, {{0, 4},{1, 5}}) + product(d_imetric, d_curva, {{2, 2},{3, 3}}) +product(d_curva, d_imetric, {{2, 2},{3, 3}}) + product(curva, dd_imetric, {{0, 4},{1, 5}});



  
   // cout<<"print: "<<NR-normalR   <<endl;
    // STRETCH CALCULATION  AND EVOLVE HEIGHT
    double jacC = jac/jacR; //
    double jacC_n = jac_n/jacR; //
    double thickness=height_in/jacC_n;

    //double mu=0.5*young/(1+nu) ;
   // double lambda=mu ;







 tensor<double, 2> Id2(pDim, pDim); //C_ab
    Id2(0,0)=1.0;
    Id2(0,1)=0.0;
    Id2(1,0)=0.0;
    Id2(1,1)=1.0;


    tensor<double, 2> metric_GPR(pDim, pDim); //gp ^A_b

    metric_GPR = metricR; //Cp_ab

    tensor<double,2> imetric_GPR(pDim,pDim); //C+^ab

    imetric_GPR=metric_GPR.inv();

    double jac_GPR     = sqrt(metric_GPR.det()); // J plus // sqrt(det gP)


    //REAL EVOLVING

    double *auxiliary_b = &fillStr.paramStr->b_aux[(subFill.loc_elemID * subFill.cubaInfo.iPts + subFill.kPt) * 11];

    double G11PG = auxiliary_b[0];
    double G12PG = auxiliary_b[1];
    double G22PG = auxiliary_b[2];

    double G11NG = auxiliary_b[3];
    double G12NG = auxiliary_b[4];
    double G22NG = auxiliary_b[5];

    double G11PG0 = auxiliary_b[6];
    double G12PG0 = auxiliary_b[7];
    double G22PG0 = auxiliary_b[8];

    // cout<<"print: "<< "height:11p: "<<GL11<<  " 12p: " <<GL22 << " :22p: " << GL33<< endl;

    tensor<double, 2> GP(pDim, pDim);
    tensor<double, 2> GN(pDim, pDim);
    tensor<double, 2> GP0(pDim, pDim);


    GP(0, 0) = G11PG; // G ^AB // imetricR(0,0)
    GP(0, 1) = G12PG;
    GP(1, 0) = G12PG;
    GP(1, 1) = G22PG;

    GN(0, 0) = G11NG;
    GN(0, 1) = G12NG;
    GN(1, 0) = G12NG;
    GN(1, 1) = G22NG;


    GP0(0, 0) = G11PG0; // G ^AB // imetricR(0,0)
    GP0(0, 1) = G12PG0;
    GP0(1, 0) = G12PG0;
    GP0(1, 1) = G22PG0;

    //GN=GN*0.0;



      /*GP=imetricR;
     GP0=imetricR;
       GN=curvatureR;*/

    tensor<double, 2> iGP = GP.inv(); //G_AB
    tensor<double, 2> iGN = GN.inv();


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


   

    tensor<double, 1> G1(nDim);
    G1 = xRu(all);
    tensor<double, 1> G2(nDim);
    G2 = xRv(all);



    tensor<double,2> d_jacC=d_jac/jacR;
    tensor<double,4> dd_jacC=dd_jac/jacR;
//till here updated


    // jacobian of viscoelastic metric

    double jacGP       = (GP.det()); // J of GPP (G11P*G22P-G12P*G12P)
//first and third invariants of visco elastic tensor
    double I1=product(GP,metric ,{{0,0},{1,1}}) ;
    double I3=jac*jac* jacGP;
    double I3_SR=sqrt(jacGP)*jac;


//d_I1, d_I3
    tensor<double,2> d_I1(eNN,numDOFs);
    tensor<double,2> d_I3(eNN,numDOFs);

    tensor<double,2> d_I3_SR(eNN,numDOFs);
    d_I1=product(GP,d_metric ,{{0,2},{1,3}}) ;
    d_I3=2*jac*d_jac* jacGP;
    d_I3_SR=sqrt(jacGP)*d_jac;


//dd_I1,dd_I3
    tensor<double,4> dd_I1(eNN,numDOFs,eNN,numDOFs);

    tensor<double,4> dd_I3(eNN,numDOFs,eNN,numDOFs);
    tensor<double,4> dd_I3_SR(eNN,numDOFs,eNN,numDOFs);


    dd_I1=product(GP,dd_metric ,{{0,4},{1,5}}) ;
    dd_I3=(2*jac*dd_jac+2*outer(d_jac,d_jac))* jacGP;
    dd_I3_SR=sqrt(jacGP)*dd_jac;



    double jacPRR=jacR;

    //  viscoelastic metric

    double jacGP_ref       = 1/(jacR*jacR); // J of GPP (G11P*G22P-G12P*G12P)


//first and third invariants
    double I1_ref=product(imetricR,metric ,{{0,0},{1,1}}) ;

    double I3_ref=jac*jac* jacGP_ref;

    double I3_SR_ref=sqrt(jacGP_ref)*jac;



    tensor<double,2> d_I1_ref(eNN,numDOFs);
    tensor<double,2> d_I3_ref(eNN,numDOFs);
    tensor<double,2> d_I3_SR_ref(eNN,numDOFs);


    d_I1_ref=product(imetricR,d_metric ,{{0,2},{1,3}}) ;
    d_I3_ref=2*jac*d_jac* jacGP_ref;
    d_I3_SR_ref=sqrt(jacGP_ref)*d_jac;


    //dd_I1_ref
    tensor<double,4> dd_I1_ref(eNN,numDOFs,eNN,numDOFs);

    tensor<double,4> dd_I3_ref(eNN,numDOFs,eNN,numDOFs);

    tensor<double,4> dd_I3_SR_ref(eNN,numDOFs,eNN,numDOFs);

    dd_I1_ref=product(imetricR,dd_metric ,{{0,4},{1,5}}) ;
    dd_I3_ref=(2*jac*dd_jac+2*outer(d_jac,d_jac))* jacGP_ref;
    dd_I3_SR_ref=sqrt(jacGP_ref)*dd_jac;



    double EA                     = (0.5*lambda*log(I3_SR)*log(I3_SR)-mu*log(I3_SR)+0.5*mu*(I1-2) )*2.0;


    //VISCO-ELASTIC ENERGY
    double Evisco                     = (0.5*lambda*log(I3_SR)*log(I3_SR)-mu*log(I3_SR)+0.5*mu*(I1-2) )*jacPRR*2.0;

    //apical
    Bp(all,range(0,2)) += ((lambda*log(I3_SR)-mu)/I3_SR*d_I3_SR+0.5*mu*d_I1)*jacPRR*2.0;
    //dd_
    App(all,range(0,2),all,range(0,2))  +=((lambda*log(I3_SR)-mu)/I3_SR*dd_I3_SR+0.5*mu*dd_I1)*jacPRR*2.0;
    App(all,range(0,2),all,range(0,2))  +=(lambda*(1-log(I3_SR))+mu)/(I3_SR*I3_SR)*outer(d_I3_SR,d_I3_SR)*jacPRR*2.0;// term 1




 double mu_ela=1.5797*deltat*fact_elastic ;
 double lambda_ela=1.2147*fact_elastic*deltat;

 //VISCO-ELAST

 //IC ENERGY
 double E_small                     = (0.5*lambda_ela*log(I3_SR_ref)*log(I3_SR_ref)-mu_ela*log(I3_SR_ref)+0.5*mu_ela*(I1_ref-2) )*jacPRR;

 //apical
 Bp(all,range(0,2)) += ((lambda_ela*log(I3_SR_ref)-mu_ela)/I3_SR_ref*d_I3_SR_ref+0.5*mu_ela*d_I1_ref)*jacPRR;
 //dd_
 App(all,range(0,2),all,range(0,2))  +=((lambda_ela*log(I3_SR_ref)-mu_ela)/I3_SR_ref*dd_I3_SR_ref+0.5*mu_ela*dd_I1_ref)*jacPRR;
 App(all,range(0,2),all,range(0,2))  +=(lambda_ela*(1-log(I3_SR_ref))+mu_ela)/(I3_SR_ref*I3_SR_ref)*outer(d_I3_SR_ref,d_I3_SR_ref)*jacPRR;// term 1






    // SHELL ANALYSIS, BENDING

    //[5.0] Tangent modulus 4th order
    tensor<double,2> a0(pDim,pDim);//INVERSE METRIC TENSOR

    a0(0,0)=imetricR(0,0);
    a0(0,1)=imetricR(0,1);
    a0(1,0)=imetricR(1,0);
    a0(1,1)=imetricR(1,1);
    tensor<double,4> BIG_C(pDim,pDim,pDim,pDim);


    //Cabcd=E/(1-nu²)*[nu g^ab g^cd+0.5*(1-nu) (g^ad g^bc+g^ac g^bd]
    BIG_C(0,0,0,0)=(lambda*a0(0,0)*a0(0,0)+mu*(a0(0,0)*a0(0,0)+a0(0,0)*a0(0,0)));
    BIG_C(0,0,0,1)=(lambda*a0(0,0)*a0(0,1)+mu*(a0(0,1)*a0(0,0)+a0(0,0)*a0(0,1)));
    BIG_C(0,0,1,0)=(lambda*a0(0,0)*a0(1,0)+mu*(a0(0,0)*a0(0,1)+a0(0,1)*a0(0,0)));
    BIG_C(0,0,1,1)=(lambda*a0(0,0)*a0(1,1)+mu*(a0(0,1)*a0(0,1)+a0(0,1)*a0(0,1)));
    BIG_C(0,1,0,0)=(lambda*a0(0,1)*a0(0,0)+mu*(a0(0,0)*a0(1,0)+a0(0,0)*a0(1,0)));
    BIG_C(0,1,0,1)=(lambda*a0(0,1)*a0(0,1)+mu*(a0(0,1)*a0(1,0)+a0(0,0)*a0(1,1)));
    BIG_C(0,1,1,0)=(lambda*a0(0,1)*a0(1,0)+mu*(a0(0,0)*a0(1,1)+a0(0,1)*a0(1,0)));
    BIG_C(0,1,1,1)=(lambda*a0(0,1)*a0(1,1)+mu*(a0(0,1)*a0(1,1)+a0(0,1)*a0(1,1)));

    BIG_C(1,0,0,0)=(lambda*a0(1,0)*a0(0,0)+mu*(a0(1,0)*a0(0,0)+a0(1,0)*a0(0,0)));
    BIG_C(1,0,0,1)=(lambda*a0(1,0)*a0(0,1)+mu*(a0(1,1)*a0(0,0)+a0(1,0)*a0(0,1)));
    BIG_C(1,0,1,0)=(lambda*a0(1,0)*a0(1,0)+mu*(a0(1,0)*a0(0,1)+a0(1,1)*a0(0,0)));
    BIG_C(1,0,1,1)=(lambda*a0(1,0)*a0(1,1)+mu*(a0(1,1)*a0(0,1)+a0(1,1)*a0(0,1)));
    BIG_C(1,1,0,0)=(lambda*a0(1,1)*a0(0,0)+mu*(a0(1,0)*a0(1,0)+a0(1,0)*a0(1,0)));
    BIG_C(1,1,0,1)=(lambda*a0(1,1)*a0(0,1)+mu*(a0(1,1)*a0(1,0)+a0(1,0)*a0(1,1)));
    BIG_C(1,1,1,0)=(lambda*a0(1,1)*a0(1,0)+mu*(a0(1,0)*a0(1,1)+a0(1,1)*a0(1,0)));
    BIG_C(1,1,1,1)=(lambda*a0(1,1)*a0(1,1)+mu*(a0(1,1)*a0(1,1)+a0(1,1)*a0(1,1)));


    // membrane and bending strains
   // tensor<double,2> mem_str=0.5*(metric-metricR);//epsilon_ab
    tensor<double,2> bend_str=curva-GN;//rho



    // derivative of membrane and bending stress
  //  tensor<double,4> d_mem_str=0.5*d_metric;// del_epsilon_pq/del_(x_I ^a)
    tensor<double,4> d_bend_str=d_curva;
    // second derivatives
    //tensor<double,6> dd_mem_str=0.5*dd_metric;
    tensor<double,6> dd_bend_str=dd_curva;



    // fill_matrix

    // [6.1] Energies and dissipations
    double const_bend=thickness*thickness/4.0;


    tensor<double,2> en_m(2,2);//
    tensor<double,2> en_b(2,2);//

    en_b=product(BIG_C,bend_str,{{0,0},{1,1}}); // C_abcd*rho_ab

    // STRAIN ENERGY DENSITY
    double E_bend{};

    E_bend+=jacR*const_bend*(product(en_b,bend_str,{{0,0},{1,1}}));// C_abcd*bend_ab*bend_cd

  double EB=E_bend;
    // derivatives for force MATRIX

    tensor<double,4> d_bend1(eNN,numDOFs,2,2);//(eNN,numDOFs,2,2)
    tensor<double,4> d_bend2(eNN,numDOFs,2,2);//(eNN,numDOFs,2,2)

    d_bend1=product(d_bend_str,BIG_C,{{2,0},{3,1}});//drho_Napq*C^pqrs
    d_bend2=product(d_bend_str,BIG_C,{{2,2},{3,3}});//drho_Nars*C^pqrs


    Bp(all,range(0,2))  +=jacR*const_bend*product(d_bend1,bend_str,{{2,0},{3,1}}) ;//(drho_Napq*C^pqrs)*rho_rs
    Bp(all,range(0,2))  +=jacR*const_bend*product(d_bend2,bend_str,{{2,0},{3,1}}) ;//(drho_Nars*C^pqrs)*rho_pq


    // double derivative for HESSIAN

    // tensor<double,6> dd_curva(eNN,numDOFs,eNN,numDOFs,pDim,pDim);
    tensor<double,6> dd_bend1(eNN,numDOFs,eNN,numDOFs,2,2);//(eNN,numDOFs,2,2)
    tensor<double,6> dd_bend2(eNN,numDOFs,eNN,numDOFs,2,2);//(eNN,numDOFs,2,2)

    dd_bend1=product(dd_bend_str,BIG_C,{{4,0},{5,1}});//dd_rho_JbIapq*C^pqrs
    dd_bend2=product(dd_bend_str,BIG_C,{{4,2},{5,3}});//dd_rho_JbIars*C^pqrs

//from the derivative of memebrane and bending strains

    App(all,range(0,2),all,range(0,2)) +=jacR*const_bend*product(dd_bend1,bend_str,{{4,0},{5,1}}) ;
    App(all,range(0,2),all,range(0,2)) +=jacR*const_bend*product(dd_bend2,bend_str,{{4,0},{5,1}}) ;

    // cross diodic terms

    //from the derivative of memebrane and bending strains
    tensor<double,8> d_bend_d_bend(eNN,numDOFs,2,2,eNN,numDOFs,2,2);//(drho diodioc drho)_JbrsIapq
    d_bend_d_bend=outer(d_bend_str,d_bend_str);

    tensor<double,4> d_bend_d_bend_Cpqrs(eNN,numDOFs,eNN,numDOFs);//(drho diodioc drho)_JbrsIapq*C^pqrs
    tensor<double,4> d_bend_d_bend_Crspq(eNN,numDOFs,eNN,numDOFs);//(drho diodioc drho)_JbpqIars*C^pqrs

    d_bend_d_bend_Cpqrs= product( d_bend_d_bend,BIG_C,{{2,2},{3,3},{6,0},{7,1}});//(drho diodioc drho)_JbrsIapq*C^pqrs
    d_bend_d_bend_Crspq= product( d_bend_d_bend,BIG_C,{{2,0},{3,1},{6,2},{7,3}});//(drho diodioc drho)_JbpqIars*C^pqrs


    App(all,range(0,2),all,range(0,2)) +=jacR*const_bend*d_bend_d_bend_Cpqrs;
    App(all,range(0,2),all,range(0,2)) +=jacR*const_bend*d_bend_d_bend_Crspq ;





    double Efric                        = 0.5*fric*jac_n*product(x-x_n,x-x_n,{{0,0}});
    Lagrangian+=Efric;


    Bp(all,range(0,2))                 += fric*jac_n*outer(bf,x-x_n);
    App(all,range(0,2),all,range(0,2)) += fric*jac_n*outer(outer(bf,Id),bf).transpose({0,1,3,2});



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
    tensor<double,2> d_z(eNN,numDOFs);

    double  E_force = jacR * force * (x(2) - xR(2));
     //Bp(all, 2) += -force * jacR * bf;

    /*   d_x=outer(bf,Id1);
       d_y=outer(bf,Id2);*/
    d_z=outer(bf,Id3);

  if(spring>0.5)
{
    Lagrangian += jacR *0.5*kspr*product(x-xR,x-xR,{{0,0}});


    Bp(all,0)                 += kspr*jacR*bf*(x(0)-xR(0));
    App(all,0,all,0) += kspr*jacR*outer(bf,bf);

    Bp(all,1)                 += kspr*jacR*bf*(x(1)-xR(1));
    App(all,1,all,1) += kspr*jacR*outer(bf,bf);

    Bp(all,2)                 += kspr*jacR*bf*(x(2)-xR(2));
    App(all,2,all,2) += kspr*jacR*outer(bf,bf);


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



        Lagrangian +=0.0;

       // double power1=-0.5*gamma*product(iGP,(GP-GP0),{{0,0},{1,1}})/(deltat*deltat)*jac;

    //double power1=-0.5*gamma*product(iGP,(GP-GP0),{{0,0},{1,1}})*jac/(deltat*deltat);
      double power1=gamma*(jac-jac_n)*1.0/(deltat*deltat);

       double E_tens1=0.0;


       // Bp(all, range(0, 2)) += -0.5*gamma*product(iGP,(GP-GP0),{{0,0},{1,1}})*d_jac;
       // App(all, range(0, 2), all, range(0, 2)) +=-0.5*gamma*product(iGP,(GP-GP0),{{0,0},{1,1}})*dd_jac;


         Bp(all, range(0, 2)) += gamma*d_jac;
         App(all, range(0, 2), all, range(0, 2)) +=gamma*dd_jac;




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
    fillStr.addToGlobalIntegral("Dissipation",Efric/(deltat*deltat));

    fillStr.addToGlobalIntegral("E_tension1",E_tens1);
   fillStr.addToGlobalIntegral("E_elastic",E_small/deltat);



    fillStr.addToGlobalIntegral("Energy",Evisco/deltat);

    fillStr.addToGlobalIntegral("E1",E_bend/deltat);

    fillStr.addToGlobalIntegral("Epress",Evol/deltat);

    // cout<< "king: " <<jacC<<endl;

    fillStr.addToGlobalIntegral("Trd",jac*jacR);
    fillStr.addToGlobalIntegral("area_n", jacR);
    fillStr.addToGlobalIntegral("area", jac);
    fillStr.addToGlobalIntegral("volume", (1.0/3.0) * jac*xnormal);

     fillStr.addToGlobalIntegral("P_tension1",power1);

    //This is to pass info to the other problem
    tensor<double,2> delI1_delGP(pDim,pDim);//del I1 /del G ^bR
    tensor<double,2> delI3_delGP(pDim,pDim);//del I1 /del G ^bR


    delI1_delGP=metric;

    delI3_delGP=jac*jac*jacGP*iGP; //_Br


    tensor<double,2> delWP_delGP(pDim,pDim);//del j+ /del G ^bR
    tensor<double,2> delWN_delGN(pDim,pDim);//



    delWP_delGP=0.5*(lambda*log(I3_SR)-mu)/I3*delI3_delGP+0.5*mu*delI1_delGP;





    delWN_delGN=-1*product(BIG_C,bend_str,{{0,0},{1,1}})*const_bend*0.0;

 tensor<double,2> delI1_delC(pDim,pDim);//del j+ /del G ^bR
 tensor<double,2> delI3_delC(pDim,pDim);//
 tensor<double,2> delWP_delC(pDim,pDim);//del j+ /del G ^bR


 delI1_delC=GP ;
 delI3_delC=jac*jac*jacGP*imetric; //_Br

 delWP_delC=0.5*(lambda*log(I3_SR)-mu)/I3*delI3_delC+0.5*mu*delI1_delC;

 tensor<double,2> sigma(pDim,pDim);//del j+ /del G ^bR
 sigma=2*delWP_delC/jac*jacR;

  double trace_sigma=product(sigma,metric,{{0,0},{1,1}});
 tensor<double,2> dev_sigma1(pDim,pDim);//del j+ /del G ^bR
 dev_sigma1=sigma-0.5*trace_sigma*imetric;
 double dev_sigma=product(dev_sigma1*metric,metric*dev_sigma1,{{0,0},{1,1}});

 trace_sigma=trace_sigma/deltat;
 dev_sigma=sqrt(dev_sigma)/deltat;


    double* auxiliary_a = &fillStr.paramStr->a_aux[(subFill.loc_elemID*subFill.cubaInfo.iPts+subFill.kPt)*58];

    auxiliary_a[0]  = x(0);
    auxiliary_a[1]  = x(1);
    auxiliary_a[2]  = x(2);

    auxiliary_a[3]  =  metric(0,0);
    auxiliary_a[4]  =  metric(0,1);
    auxiliary_a[5]  =  metric(1,0);
    auxiliary_a[6]  =  metric(1,1);

   auxiliary_a[7]  =  metric(0,0);
    auxiliary_a[8]  =  metric(0,1);
    auxiliary_a[9]  =  metric(1,0);
    auxiliary_a[10]  =  metric(1,1);



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

    auxiliary_a[21] = Evisco/deltat;
    auxiliary_a[22] =0.0/deltat;

    auxiliary_a[23] =0.0;
    auxiliary_a[24] =0.0;
    auxiliary_a[25] = 0.0;
    auxiliary_a[26] = 0.0;

    auxiliary_a[27] =normal(0);
    auxiliary_a[28] =normal(1);
    auxiliary_a[29] = normal(2);


    auxiliary_a[30] = thickness;

    auxiliary_a[31] =0.0;
    auxiliary_a[32] =0.0;
    auxiliary_a[33] = 0.0;
    auxiliary_a[34] = 0.0;

    auxiliary_a[35]  =  imetricR(0,0);
    auxiliary_a[36]  =  imetricR(0,1);
    auxiliary_a[37]  =  imetricR(1,0);
    auxiliary_a[38]  = imetricR(1,1);

    auxiliary_a[39]  =  icurvatureR(0,0);
    auxiliary_a[40]  =  icurvatureR(0,1);
    auxiliary_a[41]  =  icurvatureR(1,0);
    auxiliary_a[42]  = icurvatureR(1,1);



    auxiliary_a[43]  =  GP(0,0);
    auxiliary_a[44]  =  GP(1,0);
    auxiliary_a[45]  =  GP(0,1);
    auxiliary_a[46]  = GP(1,1);



   auxiliary_a[47]  =   jacPRR;


auxiliary_a[48]  =EA/deltat;
auxiliary_a[49]  =EB/deltat;

auxiliary_a[50]  =power1;

    auxiliary_a[51]  =  curva(0,0);
    auxiliary_a[52]  =  curva(1,0);
    auxiliary_a[53]  =  curva(0,1);
    auxiliary_a[54]  = curva(1,1);

    auxiliary_a[55]  =const_bend;

 auxiliary_a[56]  =trace_sigma;
 auxiliary_a[57]  =dev_sigma;


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
    tensor<double,2,false> nborCoords(subFill.nborCoords.data(),eNN,nDim);
    tensor<double,2,false> nborDOFs0(subFill.nborDOFs0.data(),eNN,numDOFs);
    tensor<double,2,false> nborDOFs(subFill.nborDOFs.data(),eNN,numDOFs);

    //Basis functions and derivatives
    tensor<double,1,false>    bf(subFill.nborBFsDers(0),eNN);
    tensor<double,2,false>   Dbf(subFill.nborBFsDers(1),eNN,pDim);
    tensor<double,3,false>  DDbf(subFill.nborBFsDers(2),eNN,pDim,pDim);

    //output
    tensor<double,2,false>  Bk(fillStr.Bk(0).data(),eNN,numDOFs);
    tensor<double,4,false>  Ak(fillStr.Ak(0,0).data(),eNN,numDOFs,eNN,numDOFs);


    // --------------------------------------------
    // [1] Model parameters
    // --------------------------------------------
    //[1.4] Parameters
    double deltat = fillStr.getRealParameter(MembParams::deltat);
    double fric2    = fillStr.getRealParameter(MembParams::fric2);
    int time_target    = fillStr.getIntParameter(MembParams::time_target);
    double fric3    =  fillStr.getRealParameter(MembParams::fric3);
    double fric    = fillStr.getRealParameter(MembParams::fric);

 double factor    =  fillStr.getRealParameter(MembParams::factor);

 double height_in    =  fillStr.getRealParameter(MembParams::height_in);

 double fric_fact    = fillStr.getRealParameter(MembParams::fric_fact);

    int timestep = fillStr.getIntParameter(MembParams::timestep);//timestep



    double forward_new= fillStr.getRealParameter(MembParams::forward_new);

     double lambda   = fillStr.getRealParameter(MembParams::lambda)*deltat;
    double mu   =fillStr.getRealParameter(MembParams::mu) *deltat;

 double gamma=fillStr.getRealParameter(MembParams::gamma) * deltat*0.0;



 double gamma_minus=fillStr.getRealParameter(MembParams::gamma_minus) * deltat;
 double gamma_plus=fillStr.getRealParameter(MembParams::gamma_plus) * deltat;



    //double mu=0.5*young/(1+nu) ;
    //double lambda=mu ;
    // --------------------------------------------
    // [2] Auxiliary variables
    // --------------------------------------------
    //FIXME: this has to be loaded from some structure
    double* auxiliary_a = &fillStr.paramStr->a_aux[(subFill.loc_elemID*subFill.cubaInfo.iPts+subFill.kPt)*58];
    tensor<double,1,false>      x(&auxiliary_a[0],3);
    tensor<double,2,false> metric(&auxiliary_a[3],2,2);
    //tensor<double,2> metric_gn(&auxiliary_a[7],2,2);
    tensor<double,2,false> delWP_delGP(&auxiliary_a[11],2,2);
   // tensor<double,2> delWN_delGN(&auxiliary_a[15],2,2);
    tensor<double,2> delWB_delK(pDim,pDim);//

    tensor<double,2,false> imetricR(&auxiliary_a[35],2,2);
    tensor<double,2,false>icurvatureR(&auxiliary_a[39],2,2);
    tensor<double,2,false>curva(&auxiliary_a[51],2,2);
    tensor<double,2>curvatureR=icurvatureR.inv();

    tensor<double,2> metricR=imetricR.inv();
   // tensor<double,2> metric_GNR=imetric_GNR.inv();


    //Jacobian
    double jacR= auxiliary_a[19];
    double jac= auxiliary_a[20];
      double const_bend= auxiliary_a[55];

//[5.0] Tangent modulus 4th order
    tensor<double,2> a0(pDim,pDim);//INVERSE METRIC TENSOR

    a0(0,0)=imetricR(0,0);
    a0(0,1)=imetricR(0,1);
    a0(1,0)=imetricR(1,0);
    a0(1,1)=imetricR(1,1);
    tensor<double,4> BIG_C(pDim,pDim,pDim,pDim);


    //Cabcd=E/(1-nu²)*[nu g^ab g^cd+0.5*(1-nu) (g^ad g^bc+g^ac g^bd]
    BIG_C(0,0,0,0)=(lambda*a0(0,0)*a0(0,0)+mu*(a0(0,0)*a0(0,0)+a0(0,0)*a0(0,0)));
    BIG_C(0,0,0,1)=(lambda*a0(0,0)*a0(0,1)+mu*(a0(0,1)*a0(0,0)+a0(0,0)*a0(0,1)));
    BIG_C(0,0,1,0)=(lambda*a0(0,0)*a0(1,0)+mu*(a0(0,0)*a0(0,1)+a0(0,1)*a0(0,0)));
    BIG_C(0,0,1,1)=(lambda*a0(0,0)*a0(1,1)+mu*(a0(0,1)*a0(0,1)+a0(0,1)*a0(0,1)));
    BIG_C(0,1,0,0)=(lambda*a0(0,1)*a0(0,0)+mu*(a0(0,0)*a0(1,0)+a0(0,0)*a0(1,0)));
    BIG_C(0,1,0,1)=(lambda*a0(0,1)*a0(0,1)+mu*(a0(0,1)*a0(1,0)+a0(0,0)*a0(1,1)));
    BIG_C(0,1,1,0)=(lambda*a0(0,1)*a0(1,0)+mu*(a0(0,0)*a0(1,1)+a0(0,1)*a0(1,0)));
    BIG_C(0,1,1,1)=(lambda*a0(0,1)*a0(1,1)+mu*(a0(0,1)*a0(1,1)+a0(0,1)*a0(1,1)));

    BIG_C(1,0,0,0)=(lambda*a0(1,0)*a0(0,0)+mu*(a0(1,0)*a0(0,0)+a0(1,0)*a0(0,0)));
    BIG_C(1,0,0,1)=(lambda*a0(1,0)*a0(0,1)+mu*(a0(1,1)*a0(0,0)+a0(1,0)*a0(0,1)));
    BIG_C(1,0,1,0)=(lambda*a0(1,0)*a0(1,0)+mu*(a0(1,0)*a0(0,1)+a0(1,1)*a0(0,0)));
    BIG_C(1,0,1,1)=(lambda*a0(1,0)*a0(1,1)+mu*(a0(1,1)*a0(0,1)+a0(1,1)*a0(0,1)));
    BIG_C(1,1,0,0)=(lambda*a0(1,1)*a0(0,0)+mu*(a0(1,0)*a0(1,0)+a0(1,0)*a0(1,0)));
    BIG_C(1,1,0,1)=(lambda*a0(1,1)*a0(0,1)+mu*(a0(1,1)*a0(1,0)+a0(1,0)*a0(1,1)));
    BIG_C(1,1,1,0)=(lambda*a0(1,1)*a0(1,0)+mu*(a0(1,0)*a0(1,1)+a0(1,1)*a0(1,0)));
    BIG_C(1,1,1,1)=(lambda*a0(1,1)*a0(1,1)+mu*(a0(1,1)*a0(1,1)+a0(1,1)*a0(1,1)));


    // membrane and bending strains




    //SUPG
    double G11P0{1.0}, G12P0{0.0},G22P0{1.0}, G11N0{1.0}, G12N0{0},G22N0{1.0} , G11P{1.0}, G12P{0.0},G22P{1.0}, G11N{1.0}, G12N{0},G22N{1.0};

/*
     G11P0  = nborDOFs0(all,0) * bf;
    G12P0  = nborDOFs0(all,1) * bf;
    G22P0  = nborDOFs0(all,2) * bf;
    G11N0  = nborDOFs0(all,3) * bf;
    G12N0  = nborDOFs0(all,4) * bf;
    G22N0  = nborDOFs0(all,5) * bf;
*/

//cout<<"step: "<<timestep<<"check: "<<G11P0<<endl;

    // FIXME: We need to remove this




    //cout<<"check: "<<  " :"<<dt_eta*delWP_delGP0<<endl;


    // --------------------------------------------
    // [3] Fill rhs and matrix
    // --------------------------------------------

    // --------------------------------------------
    // [4] Extra
    // --------------------------------------------
    //cout<<"print: "<<G1L<< ": " <<G2L<< endl;
    // This is to pass info to the other problem


    double* auxiliary_b = &fillStr.paramStr->b_aux[(subFill.loc_elemID*subFill.cubaInfo.iPts+subFill.kPt)*11];
     tensor<double,2> GP_old(2,2);
     tensor<double,2> GN_old(2,2);

 //cout<<"curva R<<"<< curvatureR<<endl;
// cout<<"icurva R<<"<< icurvatureR<<endl;


      if(timestep<time_target)
        {

             GP_old(0,0) =1.0*imetricR(0,0);
             GP_old(1,0) = 1.0*imetricR(0,1);
             GP_old(0,1) = 1.0*imetricR(0,1);
             GP_old(1,1) =1.0* imetricR(1,1);


            GN_old(0,0) = curvatureR(0,0);
            GN_old(0,1)=  curvatureR(1,0);
            GN_old(1,0)=  curvatureR(1,0);
            GN_old(1,1) = curvatureR(1,1);

        }
        else
        {

              GP_old(0,0)=auxiliary_b[0];
             GP_old(1,0)=auxiliary_b[1];
             GP_old(0,1)=auxiliary_b[1];
                 GP_old(1,1)=auxiliary_b[2];

               GN_old(0,0)=auxiliary_b[3];
              GN_old(1,0)=auxiliary_b[4];
               GN_old(0,1)=auxiliary_b[4];
              GN_old(1,1)=auxiliary_b[5];


        }

   double dt_eta=-1/fric2;
    double dt_eta_k=-1/fric3;


   if(forward_new>0.1)
   {

     tensor<double,2> delWP_delGP_n=delWP_delGP;
     tensor<double,2> delWP_delGP_CC(2,2);
     delWP_delGP_CC= product(metricR.inv(),product(delWP_delGP_n, metricR.inv(), {{1, 0}}), {{0, 0}});


    tensor<double,2> bend_str=curva-GN_old;//rho

    tensor<double,2>  en_b=product(BIG_C,bend_str,{{0,0},{1,1}}); // C_abcd*rho_ab
    tensor<double,2> delWN_delGN=const_bend*(-2.0)*en_b;

     tensor<double,2> delWN_delGN_CC(2,2);
     delWN_delGN_CC= product(metricR,product(delWN_delGN,metricR, {{1, 0}}), {{0, 0}});


    // Compute reactions


    double enP11=dt_eta*delWP_delGP_CC(0,0);
    double enP12=dt_eta*delWP_delGP_CC(0,1);
    double enP22=dt_eta*delWP_delGP_CC(1,1);

    double enN11=dt_eta_k*delWN_delGN_CC(0,0);
    double enN12=dt_eta_k*delWN_delGN_CC(0,1);
    double enN22=dt_eta_k*delWN_delGN_CC(1,1);



        if(timestep<time_target)
        {

            G11P =1.0*imetricR(0,0);
            G12P = 1.0*imetricR(0,1);
            G22P = 1.0*imetricR(1,1);
            G11N = curvatureR(0,0);
            G12N=  curvatureR(1,0);
            G22N = curvatureR(1,1);

        }
        else
        {


            G11P = auxiliary_b[0]  + enP11;
            G12P = auxiliary_b[1]  + enP12;
            G22P = auxiliary_b[2]  + enP22;
            G11N = auxiliary_b[3]  + enN11;
            G12N = auxiliary_b[4]  + enN12;
            G22N= auxiliary_b[5]  + enN22;


        }



    }





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



     //nth step

     if(timestep<time_target)
     {

            G11P =1.0*imetricR(0,0);
            G12P = 1.0*imetricR(0,1);
            G22P = 1.0*imetricR(1,1);
            G11N = curvatureR(0,0);
            G12N=  curvatureR(1,0);
            G22N = curvatureR(1,1);

     }
     else
     {

        //newton loop
         int nr_step=0;
         double GP_11_NR{0.0},GP_12_NR{0.0},GP_22_NR{0.0},GN_11_NR{0.0},GN_12_NR{0.0},GN_22_NR{0.0};
         GP_11_NR= GP_n(0,0);
         GP_12_NR= GP_n(0,1);
         GP_22_NR= GP_n(1,1);

         GN_11_NR= GN_n(0,0);
         GN_12_NR= GN_n(0,1);
         GN_22_NR= GN_n(1,1);



         double error=0.1;
         double error2=0.1;

         while ((error>0.0000001|| error2>0.0000001) && nr_step<30)
         {
              jacGP_n       = sqrt(GP_n.det());
              iGP_n=GP_n.inv();

             jacGN_n       = sqrt(GN_n.det());
             iGN_n=GN_n.inv();


             //This is to pass info to the other problem
             tensor<double,2> delI1_delGP_n(pDim,pDim);//del I1 /del G ^bR
             tensor<double,2> delI1_gn_delGN_n(pDim,pDim);//


             delI1_delGP_n=   (metric)    ; //
            // delI1_gn_delGN_n=   (metric_gn)    ; //_Br


             tensor<double,2> delI3_delGP_n(pDim,pDim);//del I1 /del G ^bR
             tensor<double,2> delI3_gn_delGN_n(pDim,pDim);//



             delI3_delGP_n=jac*jac*jacGP_n*jacGP_n*iGP_n; //_ab
           //  delI3_gn_delGN_n=jac_gn*jac_gn*jacGN_n*jacGN_n*iGN_n;


             tensor<double,4> d_delI3_delGP_n(pDim,pDim,pDim,pDim);//ab_cd
             d_delI3_delGP_n=jac*jac*jacGP_n*jacGP_n*outer(iGP_n,iGP_n)-jac*jac*jacGP_n*jacGP_n*outer(iGP_n,iGP_n).transpose({0,2,3,1});//pqcd
            // tensor<double,4> d_delI3_gn_delGN_n(pDim,pDim,pDim,pDim);//del I1 /del G ^bR
            // d_delI3_gn_delGN_n=jac_gn*jac_gn*jacGN_n*jacGN_n*outer(iGN_n,iGN_n)-jac_gn*jac_gn*jacGN_n*jacGN_n*outer(iGN_n,iGN_n).transpose({0,2,3,1});


             tensor<double,2> delWP_delGP_n(pDim,pDim);//del j+ /del G ^bR
             tensor<double,2> delWN_delGN_n(pDim,pDim);//

             tensor<double,4> d_delWP_delGP_n(pDim,pDim,pDim,pDim);//Jabcd* Fcd
             tensor<double,4> d_delWN_delGN_n(pDim,pDim,pDim,pDim);//



             double I3_n=jac*jac* jacGP_n*jacGP_n;
             double I3_SR_n=(jacGP_n)*jac;
             delWP_delGP_n=0.5*(lambda*log(I3_SR_n)-mu)/I3_n*delI3_delGP_n+0.5*mu*delI1_delGP_n;

                tensor<double,2> bend_str=curva-GN_n;//rho
               tensor<double,2>  en_b=product(BIG_C,bend_str,{{0,0},{1,1}}); // C_abcd*rho_ab
              delWN_delGN_n=const_bend*(-2.0)*en_b;


             d_delWP_delGP_n=(-2*lambda*log(I3_SR_n)+2*mu+lambda)/(4*I3_n*I3_n)*outer(delI3_delGP_n,delI3_delGP_n)+0.5*(lambda*log(I3_SR_n)-mu)/I3_n*d_delI3_delGP_n;
             d_delWN_delGN_n=const_bend*(2.0)*BIG_C;

             tensor<double,2> delWP_delGP_n_CC(pDim,pDim);//del j+ /del G ^bR
             tensor<double,2> delWN_delGN_n_CC(pDim,pDim);//


             tensor<double,4> d_delWP_delGP_n_CC(pDim,pDim,pDim,pDim);//Jabcd* Fcd
             tensor<double,4> d_delWN_delGN_n_CC(pDim,pDim,pDim,pDim);//



             delWP_delGP_n_CC=product(metricR.inv(),product(delWP_delGP_n, metricR.inv(), {{1, 0}}), {{1, 0}});
             delWN_delGN_n_CC=product(metricR,product(delWN_delGN_n, metricR, {{1, 0}}), {{1, 0}});

             d_delWP_delGP_n_CC=product(metricR.inv(),product(d_delWP_delGP_n, metricR.inv(), {{1, 0}}), {{1, 0}}).transpose({0,3,1,2});//pqcd, ap*pqcd*qb, acdb
             d_delWN_delGN_n_CC=product(metricR,product(d_delWN_delGN_n, metricR, {{1, 0}}), {{1, 0}}).transpose({0,3,1,2});//ab_cd


             //FORCE
             double F_11_NR{0.0},F_12_NR{0.0},F_22_NR{0.0},F_GN_11_NR{0.0},F_GN_12_NR{0.0},F_GN_22_NR{0.0};

             F_11_NR=GP_11_NR-auxiliary_b[0]-dt_eta*delWP_delGP_n_CC(0,0);
             F_12_NR=GP_12_NR-auxiliary_b[1]-dt_eta*delWP_delGP_n_CC(0,1);
             F_22_NR=GP_22_NR-auxiliary_b[2]-dt_eta*delWP_delGP_n_CC(1,1);

             F_GN_11_NR=GN_11_NR-auxiliary_b[3]-dt_eta_k*delWN_delGN_n_CC(0,0);
             F_GN_12_NR=GN_12_NR-auxiliary_b[4]-dt_eta_k*delWN_delGN_n_CC(0,1);
             F_GN_22_NR=GN_22_NR-auxiliary_b[5]-dt_eta_k*delWN_delGN_n_CC(1,1);





             std::vector<double> values = {
                 abs(F_11_NR), abs(F_12_NR), abs(F_22_NR), abs(F_GN_11_NR), abs(F_GN_12_NR), abs(F_GN_22_NR)

             };

             error=*std::max_element(values.begin(), values.end());

             tensor<double,1> GP_voigt(3);//del j+ /del G ^bR
             GP_voigt(0)=GP_11_NR;
             GP_voigt(1)=GP_12_NR;
             GP_voigt(2)=GP_22_NR;
             tensor<double,1> GP_F_voigt(3);//del j+ /del G ^bR

             GP_F_voigt(0)=F_11_NR;
             GP_F_voigt(1)=F_12_NR;
             GP_F_voigt(2)=F_22_NR;

             tensor<double,2> GP_jac_voigt(3,3);//del j+ /del G ^bR
             double  J_11_NR_11=1.0-dt_eta*d_delWP_delGP_n_CC(0,0,0,0); //d_delWP_delGP_n_CC
             double  J_11_NR_12=0.0-dt_eta*d_delWP_delGP_n_CC(0,0,0,1);
             double  J_11_NR_22=0.0-dt_eta*d_delWP_delGP_n_CC(0,0,1,1);

             double  J_12_NR_11=0.0-dt_eta*d_delWP_delGP_n_CC(0,1,0,0);
             double  J_12_NR_12=1.0-dt_eta*d_delWP_delGP_n_CC(0,1,0,1);
             double  J_12_NR_22=0.0-dt_eta*d_delWP_delGP_n_CC(0,1,1,1);


             double J_22_NR_11=0.0-dt_eta*d_delWP_delGP_n_CC(1,1,0,0);
             double J_22_NR_12=0.0-dt_eta*d_delWP_delGP_n_CC(1,1,0,1);
             double J_22_NR_22=1.0-dt_eta*d_delWP_delGP_n_CC(1,1,1,1);

             GP_jac_voigt(0,0)=J_11_NR_11;
             GP_jac_voigt(0,1)=J_11_NR_12;
             GP_jac_voigt(0,2)=J_11_NR_22;

             GP_jac_voigt(1,0)=J_12_NR_11;
             GP_jac_voigt(1,1)=J_12_NR_12;
             GP_jac_voigt(1,2)=J_12_NR_22;

             GP_jac_voigt(2,0)=J_22_NR_11;
             GP_jac_voigt(2,1)=J_22_NR_12;
             GP_jac_voigt(2,2)=J_22_NR_22;

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
             double  J_GN_11_NR_11=1.0-dt_eta_k*d_delWN_delGN_n_CC(0,0,0,0);
             double  J_GN_11_NR_12=0.0-dt_eta_k*d_delWN_delGN_n_CC(0,0,0,1);
             double  J_GN_11_NR_22=0.0-dt_eta_k*d_delWN_delGN_n_CC(0,0,1,1);

             double  J_GN_12_NR_11=0.0-dt_eta_k*d_delWN_delGN_n_CC(0,1,0,0);
             double  J_GN_12_NR_12=1.0-dt_eta_k*d_delWN_delGN_n_CC(0,1,0,1);
             double  J_GN_12_NR_22=0.0-dt_eta_k*d_delWN_delGN_n_CC(0,1,1,1);


             double J_GN_22_NR_11=0.0-dt_eta_k*d_delWN_delGN_n_CC(1,1,0,0);
             double J_GN_22_NR_12=0.0-dt_eta_k*d_delWN_delGN_n_CC(1,1,0,1);
             double J_GN_22_NR_22=1.0-dt_eta_k*d_delWN_delGN_n_CC(1,1,1,1);

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




             double X_11_NR{0.0},X_12_NR{0.0},X_22_NR{0.0},X_GN_11_NR{0.0},X_GN_12_NR{0.0},X_GN_22_NR{0.0},X_GL1_11_NR{0.0},X_GL1_12_NR{0.0},X_GL1_22_NR{0.0},X_GL2_11_NR{0.0},X_GL2_12_NR{0.0},X_GL2_22_NR{0.0},X_GL3_11_NR{0.0},X_GL3_12_NR{0.0},X_GL3_22_NR{0.0};

             tensor<double,1>delX(3);//del j+ /del G ^bR
             delX=-GP_jac_voigt.inv()*GP_F_voigt;
             X_11_NR=delX(0);
             X_12_NR=delX(1);
             X_22_NR=delX(2);

             tensor<double,1>delX_GN(3);//del j+ /del G ^bR
             delX_GN=-GN_jac_voigt.inv()*GN_F_voigt;
             X_GN_11_NR=delX_GN(0);
             X_GN_12_NR=delX_GN(1);
             X_GN_22_NR=delX_GN(2);



         double tol3=0.0000000001;
        std::vector<double> values2 = {
                 abs(X_11_NR/(GP_11_NR+tol3)), abs(X_12_NR/(GP_12_NR+tol3)), abs(X_22_NR/(GP_22_NR+tol3)), abs(X_GN_11_NR/(GN_11_NR+tol3)), abs(X_GN_12_NR/(GN_12_NR+tol3)), abs(X_GN_22_NR/(GN_22_NR+tol3))
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

     }


 }




    auxiliary_b[0] = G11P;
    auxiliary_b[1] = G12P;
    auxiliary_b[2] = G22P;
    auxiliary_b[3] = G11N;
    auxiliary_b[4] = G12N;
    auxiliary_b[5] = G22N;

    auxiliary_b[6] =    GP_old(0,0);
     auxiliary_b[7] =    GP_old(0,1);
    auxiliary_b[8] =    GP_old(1,1);





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


    double jacPRR=auxiliary_a[47];


//new    product(GP_new-GP_old,metricR,{{0,0}})
         double diss11=fric2*  product(product(GP_new-GP_old,metricR,{{1,0}}),product(metricR,GP_new-GP_old,{{1,0}}),{{0,0},{1,1}})/(deltat*deltat)*jacPRR;

         double diss22= fric3* product(product(GN_new-GN_old,imetricR,{{1,0}}),product(imetricR,GN_new-GN_old,{{1,0}}),{{0,0},{1,1}})/(deltat*deltat)*jacPRR;



    // Global integrals
  //  fillStr->addToGlobalIntegral("totBonds",jacR);

 auxiliary_b[9] =diss11/jacPRR;
 auxiliary_b[10] =diss22/jacPRR;




    fillStr.addToGlobalIntegral("Diss1",diss11);
    fillStr.addToGlobalIntegral("Diss2",diss22);



}



void LS_ED(hiperlife::FillStructure & fillStr)
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

    //Coordinates and degrees of freedom
    tensor<double,2,false> nborCoords(subFill.nborCoords.data(),eNN,nDim);
   // tensor<double,2> nborDOFs0(subFill.nborDOFs0.data(),eNN,numDOFs);
   // tensor<double,2> nborDOFs(subFill.nborDOFs.data(),eNN,numDOFs);

    //Basis functions and derivatives
    tensor<double,1,false>    bf(subFill.nborBFsDers(0),eNN);
   // tensor<double,2>   Dbf(subFill.nborBFsDers(1),eNN,pDim);
   // tensor<double,3>  DDbf(subFill.nborBFsDers(2),eNN,pDim,pDim);

    //output
    tensor<double,2,false>  Bk1(fillStr.Bk(0).data(),eNN,numDOFs);
    tensor<double,4,false>  Ak1(fillStr.Ak(0,0).data(),eNN,numDOFs,eNN,numDOFs);


    // --------------------------------------------
    // [1] Model parameters


    // --------------------------------------------
    // [2] Auxiliary variables
    // --------------------------------------------
    //FIXME: this has to be loaded from some structure
    double* auxiliary_a = &fillStr.paramStr->a_aux[(subFill.loc_elemID*subFill.cubaInfo.iPts+subFill.kPt)*58];

     double trace_sigma= auxiliary_a[56]  ;
      double dev_sigma=auxiliary_a[57];

    //Jacobian
    double jacR= auxiliary_a[19];
    double jac= auxiliary_a[20];


    double nxx= auxiliary_a[27];
    double nyy= auxiliary_a[28];
    double nzz= auxiliary_a[29];
    double height1= auxiliary_a[30];


    double EA= auxiliary_a[48];
    double EB= auxiliary_a[49];


    double PA= auxiliary_a[50];


    double *auxiliary_b = &fillStr.paramStr->b_aux[(subFill.loc_elemID * subFill.cubaInfo.iPts + subFill.kPt) * 11];

    double DA= auxiliary_b[9];
    double DB= auxiliary_b[10];


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


    Bk1(all,6) = jacR * bf * (DA);
    Bk1(all,7) = jacR * bf * (DB);

    Bk1(all,8) = jacR * bf * (PA);

    Bk1(all,9) = jacR * bf * (trace_sigma);

    Bk1(all,10) = jacR * bf * (dev_sigma);

    //Matrix
    Ak1(all,0,all,0) = jacR * outer(bf, bf) ;//
    Ak1(all,1,all,1) =jacR * outer(bf, bf);
    Ak1(all,2,all,2) = jacR * outer(bf, bf) ;//
    Ak1(all,3,all,3) =jacR * outer(bf, bf);
    Ak1(all,4,all,4) = jacR * outer(bf, bf) ;//
    Ak1(all,5,all,5) =jacR * outer(bf, bf);
    Ak1(all,6,all,6) = jacR * outer(bf, bf) ;//

    Ak1(all,7,all,7) = jacR * outer(bf, bf) ;//
    Ak1(all,8,all,8) = jacR * outer(bf, bf) ;//

 Ak1(all,9,all,9) = jacR * outer(bf, bf) ;//
 Ak1(all,10,all,10) = jacR * outer(bf, bf) ;//


}


