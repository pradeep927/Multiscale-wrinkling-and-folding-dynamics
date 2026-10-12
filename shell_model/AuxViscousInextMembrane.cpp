/*
    Element fillings of shell_model (see AuxViscousInextMembrane.h):
      LS          mechanics of the shell + enclosed-volume constraint
      LS_ReacDif  relaxation of the material metric G and curvature kbar
                  (local, per Gauss point; writes b_aux)
      LS_ED       L2 projection of Gauss-point fields to the nodes (output)
    Same arithmetic as adam_shell (unused code removed).
*/
#include "hl_FillStructure.h"
#include "hl_Geometry.h"
#include "hl_SurfLagrParam.h"
#include "hl_Tensor.h"
#include <hl_LinearSolver_Direct_MUMPS.h>
#include <fstream>
#include "hl_Parser.h"
#include "hl_Core.h"
#include "AuxViscousInextMembrane.h"

void LS(hiperlife::FillStructure& fillStr)
{
    using namespace hiperlife;
    using namespace hiperlife::Tensor;
    auto &subFill = (fillStr)["posDHand"];
    int nDim = subFill.nDim;
    int pDim = subFill.pDim;
    int eNN = subFill.eNN;
    int numDOFs = subFill.numDOFs;
    wrapper<double, 2> nborDOFs0(subFill.nborDOFs0.data(), eNN, numDOFs);
    wrapper<double, 2> nborDOFs(subFill.nborDOFs.data(), eNN, numDOFs);
    wrapper<double, 2> nborCoords(subFill.nborCoords.data(), eNN, nDim);
    wrapper<double, 1> bf(subFill.nborBFs(), eNN);
    wrapper<double, 2> Dbf(subFill.nborBFsGrads(), eNN, pDim);
    wrapper<double, 3> DDbf(subFill.nborBFsHess(), eNN, pDim, pDim);
    auto &g_subFill = (fillStr)["gloDHand"];
    int g_numDOFs = g_subFill.numDOFs;
    wrapper<double, 1> gDOFs(g_subFill.nborDOFs.data(), g_numDOFs);
    double pressure = gDOFs(0);
    double deltat = fillStr.getRealParameter(MembParams::deltat);
    double fric = fillStr.getRealParameter(MembParams::fric);
    double factor = fillStr.getRealParameter(MembParams::factor);
    double thick = fillStr.getRealParameter(MembParams::thick);
    double Kconf =fillStr.getRealParameter(MembParams::Kconf) * deltat;
    double gamma =fillStr.getRealParameter(MembParams::gamma) * deltat;
    int timestep =fillStr.getIntParameter(MembParams::timestep);
    double R = fillStr.getRealParameter(MembParams::R);
    double width =fillStr.getRealParameter(MembParams::width);
    double fric_fact =fillStr.getRealParameter(MembParams::fric_fact);
    double height_in =fillStr.getRealParameter(MembParams::height_in);
    int control_fric =fillStr.getIntParameter(MembParams::control_fric);
    double kspr =fillStr.getRealParameter(MembParams::kspr) * deltat;
    double sp_gap =fillStr.getRealParameter(MembParams::sp_gap);
    int spring = fillStr.getIntParameter(MembParams::spring);
    int fric_start = fillStr.getIntParameter(MembParams::fric_start);
    double alpha =fillStr.getRealParameter(MembParams::alpha);
    double fact_elastic =fillStr.getRealParameter(MembParams::fact_elastic);
    double lambda =fillStr.getRealParameter(MembParams::lambda)* deltat;
    double mu =fillStr.getRealParameter(MembParams::mu)* deltat;
    wrapper<double, 2> Bp(fillStr.Bk(0).data(), eNN, numDOFs);
    wrapper<double, 1> Bg(fillStr.Bk(1).data(), g_numDOFs);
    wrapper<double, 4> App(fillStr.Ak(0, 0).data(), eNN, numDOFs, eNN, numDOFs);
    wrapper<double, 3> Apg(fillStr.Ak(0, 1).data(), eNN, numDOFs, g_numDOFs);
    wrapper<double, 3> Agp(fillStr.Ak(1, 0).data(), g_numDOFs, eNN, numDOFs);
    double eps=0.000001;
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
    if(timestep<2)
    {
        curvatureR(0,0)=curvatureR(0,0)+eps;
        curvatureR(1,1)=curvatureR(1,1)+eps;
    }
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
    if(timestep<2)
    {
        curva_n(0,0)=curva_n(0,0)+eps;
        curva_n(1,1)=curva_n(1,1)+eps;
    }
    tensor<double,1> x_n = bf * nborDOFs0;
    tensor<double,1> xR = bf * nborCoords;
    tensor<double,2> imetricR = metricR.inv();
    double jacR = sqrt(metricR.det());
    double jac_n = sqrt(metric_n.det());
    double xnormal_n = x_n * normal_n;
    tensor<double,2>    icurvatureR= curvatureR.inv();
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
    if(timestep<2)
    {
        curva(0,0)=curva(0,0)+eps;
        curva(1,1)=curva(1,1)+eps;
    }
    tensor<double,1> x =  bf * nborDOFs;
    tensor<double,2> imetric = metric.inv();
    double jac       = sqrt(metric.det());
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
    tensor<double,4> d_metric(eNN,numDOFs,pDim,pDim);
    SurfLagrParam::d_Metric(d_metric,Dbf,xu,xv);
    tensor<double,2> d_jac(eNN,numDOFs);
    SurfLagrParam::d_Jac(d_jac,Dbf,xu,xv,normal);
    tensor<double,3> d_normal(eNN,numDOFs,nDim);
    SurfLagrParam::d_Normal(d_normal,Dbf,xu,xv,normal,jac,d_jac);
    tensor<double,4> d_curva(eNN,numDOFs,pDim,pDim);
    SurfLagrParam::d_Curva(d_curva,DDbf,xuu,xuv,xvv,normal,d_normal);
    tensor<double,6> dd_metric(eNN,numDOFs,eNN,numDOFs,pDim,pDim);
    SurfLagrParam::dd_Metric(dd_metric,Dbf);
    tensor<double,4> dd_jac(eNN,numDOFs,eNN,numDOFs);
    SurfLagrParam::dd_Jac(dd_jac,Dbf,xu,xv,normal,d_normal);
    tensor<double,5> dd_normal(eNN,numDOFs,eNN,numDOFs,nDim);
    SurfLagrParam::dd_Normal (dd_normal, Dbf, jac, d_jac, dd_jac, normal, d_normal);
    tensor<double,6> dd_curva(eNN,numDOFs,eNN,numDOFs,pDim,pDim);
    SurfLagrParam::dd_Curva(dd_curva, DDbf, xuu, xuv, xvv, d_normal, dd_normal);
    double jacC_n = jac_n/jacR;
    double thickness=height_in/jacC_n;
    double *auxiliary_b = &fillStr.paramStr->b_aux[(subFill.loc_elemID * subFill.cubaInfo.iPts + subFill.kPt)*kNB];
    double G11PG = auxiliary_b[0];
    double G12PG = auxiliary_b[1];
    double G22PG = auxiliary_b[2];
    double G11NG = auxiliary_b[3];
    double G12NG = auxiliary_b[4];
    double G22NG = auxiliary_b[5];
    tensor<double, 2> GP(pDim, pDim);
    tensor<double, 2> GN(pDim, pDim);
    GP(0, 0) = G11PG;
    GP(0, 1) = G12PG;
    GP(1, 0) = G12PG;
    GP(1, 1) = G22PG;
    GN(0, 0) = G11NG;
    GN(0, 1) = G12NG;
    GN(1, 0) = G12NG;
    GN(1, 1) = G22NG;
    tensor<double, 2> iGP = GP.inv();
    double jacGP       = (GP.det());
    double I1=product(GP,metric ,{{0,0},{1,1}}) ;
    double I3=jac*jac* jacGP;
    double I3_SR=sqrt(jacGP)*jac;
    tensor<double,2> d_I1(eNN,numDOFs);
    tensor<double,2> d_I3_SR(eNN,numDOFs);
    d_I1=product(GP,d_metric ,{{0,2},{1,3}}) ;
    d_I3_SR=sqrt(jacGP)*d_jac;
    tensor<double,4> dd_I1(eNN,numDOFs,eNN,numDOFs);
    tensor<double,4> dd_I3_SR(eNN,numDOFs,eNN,numDOFs);
    dd_I1=product(GP,dd_metric ,{{0,4},{1,5}}) ;
    dd_I3_SR=sqrt(jacGP)*dd_jac;
    double jacPRR=jacR;
    double jacGP_ref       = 1/(jacR*jacR);
    double I1_ref=product(imetricR,metric ,{{0,0},{1,1}}) ;
    double I3_SR_ref=sqrt(jacGP_ref)*jac;
    tensor<double,2> d_I1_ref(eNN,numDOFs);
    tensor<double,2> d_I3_SR_ref(eNN,numDOFs);
    d_I1_ref=product(imetricR,d_metric ,{{0,2},{1,3}}) ;
    d_I3_SR_ref=sqrt(jacGP_ref)*d_jac;
    tensor<double,4> dd_I1_ref(eNN,numDOFs,eNN,numDOFs);
    tensor<double,4> dd_I3_SR_ref(eNN,numDOFs,eNN,numDOFs);
    dd_I1_ref=product(imetricR,dd_metric ,{{0,4},{1,5}}) ;
    dd_I3_SR_ref=sqrt(jacGP_ref)*dd_jac;
    double EA                     = (0.5*lambda*log(I3_SR)*log(I3_SR)-mu*log(I3_SR)+0.5*mu*(I1-2) )*2.0;
    double Evisco                     = (0.5*lambda*log(I3_SR)*log(I3_SR)-mu*log(I3_SR)+0.5*mu*(I1-2) )*jacPRR*2.0;
    Bp(all,range(0,2)) += ((lambda*log(I3_SR)-mu)/I3_SR*d_I3_SR+0.5*mu*d_I1)*jacPRR*2.0;
    App(all,range(0,2),all,range(0,2))  +=((lambda*log(I3_SR)-mu)/I3_SR*dd_I3_SR+0.5*mu*dd_I1)*jacPRR*2.0;
    App(all,range(0,2),all,range(0,2))  +=(lambda*(1-log(I3_SR))+mu)/(I3_SR*I3_SR)*outer(d_I3_SR,d_I3_SR)*jacPRR*2.0;
    double mu_ela=fillStr.getRealParameter(MembParams::mu_pre)*deltat*fact_elastic ;
    double lambda_ela=fillStr.getRealParameter(MembParams::lambda_pre)*fact_elastic*deltat;
    double E_small                     = (0.5*lambda_ela*log(I3_SR_ref)*log(I3_SR_ref)-mu_ela*log(I3_SR_ref)+0.5*mu_ela*(I1_ref-2) )*jacPRR;
    Bp(all,range(0,2)) += ((lambda_ela*log(I3_SR_ref)-mu_ela)/I3_SR_ref*d_I3_SR_ref+0.5*mu_ela*d_I1_ref)*jacPRR;
    App(all,range(0,2),all,range(0,2))  +=((lambda_ela*log(I3_SR_ref)-mu_ela)/I3_SR_ref*dd_I3_SR_ref+0.5*mu_ela*dd_I1_ref)*jacPRR;
    App(all,range(0,2),all,range(0,2))  +=(lambda_ela*(1-log(I3_SR_ref))+mu_ela)/(I3_SR_ref*I3_SR_ref)*outer(d_I3_SR_ref,d_I3_SR_ref)*jacPRR;
    tensor<double,2> a0(pDim,pDim);
    a0(0,0)=imetricR(0,0);
    a0(0,1)=imetricR(0,1);
    a0(1,0)=imetricR(1,0);
    a0(1,1)=imetricR(1,1);
    tensor<double,4> BIG_C(pDim,pDim,pDim,pDim);
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
    tensor<double,2> bend_str=curva-GN;
    tensor<double,4> d_bend_str=d_curva;
    tensor<double,6> dd_bend_str=dd_curva;
    double const_bend=thickness*thickness/4.0*alpha;
    tensor<double,2> en_b(2,2);
    en_b=product(BIG_C,bend_str,{{0,0},{1,1}});
    double E_bend{};
    E_bend+=jacR*const_bend*(product(en_b,bend_str,{{0,0},{1,1}}));
    double EB=E_bend;
    tensor<double,4> d_bend1(eNN,numDOFs,2,2);
    tensor<double,4> d_bend2(eNN,numDOFs,2,2);
    d_bend1=product(d_bend_str,BIG_C,{{2,0},{3,1}});
    d_bend2=product(d_bend_str,BIG_C,{{2,2},{3,3}});
    Bp(all,range(0,2))  +=jacR*const_bend*product(d_bend1,bend_str,{{2,0},{3,1}}) ;
    Bp(all,range(0,2))  +=jacR*const_bend*product(d_bend2,bend_str,{{2,0},{3,1}}) ;
    tensor<double,6> dd_bend1(eNN,numDOFs,eNN,numDOFs,2,2);
    tensor<double,6> dd_bend2(eNN,numDOFs,eNN,numDOFs,2,2);
    dd_bend1=product(dd_bend_str,BIG_C,{{4,0},{5,1}});
    dd_bend2=product(dd_bend_str,BIG_C,{{4,2},{5,3}});
    App(all,range(0,2),all,range(0,2)) +=jacR*const_bend*product(dd_bend1,bend_str,{{4,0},{5,1}}) ;
    App(all,range(0,2),all,range(0,2)) +=jacR*const_bend*product(dd_bend2,bend_str,{{4,0},{5,1}}) ;
    tensor<double,8> d_bend_d_bend(eNN,numDOFs,2,2,eNN,numDOFs,2,2);
    d_bend_d_bend=outer(d_bend_str,d_bend_str);
    tensor<double,4> d_bend_d_bend_Cpqrs(eNN,numDOFs,eNN,numDOFs);
    tensor<double,4> d_bend_d_bend_Crspq(eNN,numDOFs,eNN,numDOFs);
    d_bend_d_bend_Cpqrs= product( d_bend_d_bend,BIG_C,{{2,2},{3,3},{6,0},{7,1}});
    d_bend_d_bend_Crspq= product( d_bend_d_bend,BIG_C,{{2,0},{3,1},{6,2},{7,3}});
    App(all,range(0,2),all,range(0,2)) +=jacR*const_bend*d_bend_d_bend_Cpqrs;
    App(all,range(0,2),all,range(0,2)) +=jacR*const_bend*d_bend_d_bend_Crspq ;
    double Efric                        = 0.5*fric*jac_n*product(x-x_n,x-x_n,{{0,0}});
    Bp(all,range(0,2))                 += fric*jac_n*outer(bf,x-x_n);
    App(all,range(0,2),all,range(0,2)) += fric*jac_n*outer(outer(bf,Id),bf).transpose({0,1,3,2});
    tensor<double,1> Id3(3);
    Id3(0)=0;
    Id3(1)=0;
    Id3(2)=1;
    tensor<double,2> d_z(eNN,numDOFs);
    d_z=outer(bf,Id3);
    if(spring>0.5)
    {
        Bp(all,0)                 += kspr*jacR*bf*(x(0)-xR(0));
        App(all,0,all,0) += kspr*jacR*outer(bf,bf);
        Bp(all,1)                 += kspr*jacR*bf*(x(1)-xR(1));
        App(all,1,all,1) += kspr*jacR*outer(bf,bf);
        Bp(all,2)                 += kspr*jacR*bf*(x(2)-xR(2));
        App(all,2,all,2) += kspr*jacR*outer(bf,bf);
    }
    Bp(all,range(0,2))                 += jac*d_z*dConfinementPotential(x, Kconf);
    Bp(all,range(0,2))        += d_jac*ConfinementPotential(x, Kconf);
    App(all,range(0,2),all,range(0,2)) +=outer(d_jac,d_z)*dConfinementPotential(x, Kconf);
    App(all,range(0,2),all,range(0,2)) +=jac*outer(d_z,d_z)*ddConfinementPotential(x, Kconf);
    App(all,range(0,2),all,range(0,2)) +=dd_jac*ConfinementPotential(x, Kconf);
    App(all,range(0,2),all,range(0,2)) +=outer(d_z,d_jac)*dConfinementPotential(x, Kconf);
    double power1=gamma*(jac-jac_n)*1.0/(deltat*deltat);
    double E_tens1=0.0;
    Bp(all, range(0, 2)) += gamma*d_jac;
    App(all, range(0, 2), all, range(0, 2)) +=gamma*dd_jac;
    double Evol=pressure/3.0 *(jac * xnormal - factor*jac_n*xnormal_n);
    Bp(all, range(0, 2)) += pressure / 3.0 * (d_jac * xnormal + jac * (d_normal * x) + jac * outer(bf, normal));
    Bg(0) += (jac * xnormal / 3) - factor * jacR * 1;
    App(all, range(0, 2), all, range(0, 2)) += pressure / 3.0 * (dd_jac * xnormal + outer(d_jac, d_normal * x) + outer(d_jac, outer(bf, normal)) +
        outer(d_normal * x, d_jac) + outer(outer(bf, normal), d_jac) + jac * dd_normal * x +
        jac * outer(d_normal, bf).transpose({0, 1, 3, 2}) +
        jac * outer(bf, d_normal.transpose({2, 0, 1})));
    Apg(all, range(0, 2), 0) += 1.0 / 3.0 * (d_jac * xnormal + jac * d_normal * x + jac * outer(bf, normal));
    Agp = Apg.transpose({2,0,1});
    fillStr.addToGlobalIntegral("Dissipation",Efric/(deltat*deltat));
    fillStr.addToGlobalIntegral("E_tension1",E_tens1);
    fillStr.addToGlobalIntegral("E_elastic",E_small/deltat);
    fillStr.addToGlobalIntegral("Energy",Evisco/deltat);
    fillStr.addToGlobalIntegral("E1",E_bend/deltat);
    fillStr.addToGlobalIntegral("Epress",Evol/deltat);
    fillStr.addToGlobalIntegral("Trd",jac*jacR);
    fillStr.addToGlobalIntegral("area_n", jacR);
    fillStr.addToGlobalIntegral("area", jac);
    fillStr.addToGlobalIntegral("volume", (1.0/3.0) * jac*xnormal);
    fillStr.addToGlobalIntegral("P_tension1",power1);
    tensor<double,2> delI1_delGP(pDim,pDim);
    tensor<double,2> delI3_delGP(pDim,pDim);
    delI1_delGP=metric;
    delI3_delGP=jac*jac*jacGP*iGP;
    tensor<double,2> delWP_delGP(pDim,pDim);
    tensor<double,2> delWN_delGN(pDim,pDim);
    delWP_delGP=0.5*(lambda*log(I3_SR)-mu)/I3*delI3_delGP+0.5*mu*delI1_delGP;
    delWN_delGN=-1*product(BIG_C,bend_str,{{0,0},{1,1}})*const_bend*0.0;
    tensor<double,2> delI1_delC(pDim,pDim);
    tensor<double,2> delI3_delC(pDim,pDim);
    tensor<double,2> delWP_delC(pDim,pDim);
    delI1_delC=GP ;
    delI3_delC=jac*jac*jacGP*imetric;
    delWP_delC=0.5*(lambda*log(I3_SR)-mu)/I3*delI3_delC+0.5*mu*delI1_delC;
    tensor<double,2> sigma(pDim,pDim);
    sigma=2*delWP_delC/jac*jacR;
    double trace_sigma=product(sigma,metric,{{0,0},{1,1}});
    tensor<double,2> dev_sigma1(pDim,pDim);
    dev_sigma1=sigma-0.5*trace_sigma*imetric;
    double dev_sigma=product(dev_sigma1*metric,metric*dev_sigma1,{{0,0},{1,1}});
    trace_sigma=trace_sigma/deltat;
    dev_sigma=sqrt(dev_sigma)/deltat;
    double* auxiliary_a = &fillStr.paramStr->a_aux[(subFill.loc_elemID*subFill.cubaInfo.iPts+subFill.kPt)*kNA];
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
    auto& subFill = (fillStr)["viscoDHand"];
    int numDOFs = subFill.numDOFs;
    int pDim    = subFill.pDim;
    int eNN     = subFill.eNN;
    double deltat =fillStr.getRealParameter(MembParams::deltat);
    double fric2 =fillStr.getRealParameter(MembParams::fric2);
    int time_target =fillStr.getIntParameter(MembParams::time_target);
    double fric3 =fillStr.getRealParameter(MembParams::fric3);
    int timestep =fillStr.getIntParameter(MembParams::timestep);
    double alpha =1.0;
    double forward_new =fillStr.getRealParameter(MembParams::forward_new);
    double lambda =fillStr.getRealParameter(MembParams::lambda)* deltat;
    double mu =fillStr.getRealParameter(MembParams::mu)* deltat;
    double* auxiliary_a = &fillStr.paramStr->a_aux[(subFill.loc_elemID*subFill.cubaInfo.iPts+subFill.kPt)*kNA];
    wrapper<double,2> metric(&auxiliary_a[3],2,2);
    wrapper<double,2> delWP_delGP(&auxiliary_a[11],2,2);
    wrapper<double,2> imetricR(&auxiliary_a[35],2,2);
    wrapper<double,2>icurvatureR(&auxiliary_a[39],2,2);
    wrapper<double,2>curva(&auxiliary_a[51],2,2);
    tensor<double,2>curvatureR=icurvatureR.inv();
    tensor<double,2> metricR=imetricR.inv();
    double jac= auxiliary_a[20];
    double const_bend= auxiliary_a[55];
    tensor<double,2> a0(pDim,pDim);
    a0(0,0)=imetricR(0,0);
    a0(0,1)=imetricR(0,1);
    a0(1,0)=imetricR(1,0);
    a0(1,1)=imetricR(1,1);
    tensor<double,4> BIG_C(pDim,pDim,pDim,pDim);
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
    double G11P{1.0};
    double G12P{0.0};
    double G22P{1.0};
    double G11N{1.0};
    double G12N{0};
    double G22N{1.0};
    double* auxiliary_b = &fillStr.paramStr->b_aux[(subFill.loc_elemID*subFill.cubaInfo.iPts+subFill.kPt)*kNB];
    tensor<double,2> GP_old(2,2);
    tensor<double,2> GN_old(2,2);
    if(timestep<time_target)
    {
        GP_old(0,0) =1.0/alpha*imetricR(0,0);
        GP_old(1,0) = 1.0/alpha*imetricR(0,1);
        GP_old(0,1) = 1.0/alpha*imetricR(0,1);
        GP_old(1,1) =1.0/alpha* imetricR(1,1);
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
        tensor<double,2> bend_str=curva-GN_old;
        tensor<double,2>  en_b=product(BIG_C,bend_str,{{0,0},{1,1}});
        tensor<double,2> delWN_delGN=const_bend*(-2.0)*en_b;
        tensor<double,2> delWN_delGN_CC(2,2);
        delWN_delGN_CC= product(metricR,product(delWN_delGN,metricR, {{1, 0}}), {{0, 0}});
        double enP11=dt_eta*delWP_delGP_CC(0,0);
        double enP12=dt_eta*delWP_delGP_CC(0,1);
        double enP22=dt_eta*delWP_delGP_CC(1,1);
        double enN11=dt_eta_k*delWN_delGN_CC(0,0);
        double enN12=dt_eta_k*delWN_delGN_CC(0,1);
        double enN22=dt_eta_k*delWN_delGN_CC(1,1);
        if(timestep<time_target)
        {
            G11P =1.0/alpha*imetricR(0,0);
            G12P = 1.0/alpha*imetricR(0,1);
            G22P = 1.0/alpha*imetricR(1,1);
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
        if(timestep<time_target)
        {
            G11P =1.0/alpha*imetricR(0,0);
            G12P = 1.0/alpha*imetricR(0,1);
            G22P = 1.0/alpha*imetricR(1,1);
            G11N = curvatureR(0,0);
            G12N=  curvatureR(1,0);
            G22N = curvatureR(1,1);
        }
        else
        {
            int nr_step=0;
            double GP_11_NR{0.0};
            double GP_12_NR{0.0};
            double GP_22_NR{0.0};
            double GN_11_NR{0.0};
            double GN_12_NR{0.0};
            double GN_22_NR{0.0};
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
                tensor<double,2> delI1_delGP_n(pDim,pDim);
                delI1_delGP_n=   (metric)    ;
                tensor<double,2> delI3_delGP_n(pDim,pDim);
                delI3_delGP_n=jac*jac*jacGP_n*jacGP_n*iGP_n;
                tensor<double,4> d_delI3_delGP_n(pDim,pDim,pDim,pDim);
                d_delI3_delGP_n=jac*jac*jacGP_n*jacGP_n*outer(iGP_n,iGP_n)-jac*jac*jacGP_n*jacGP_n*outer(iGP_n,iGP_n).transpose({0,2,3,1});
                tensor<double,2> delWP_delGP_n(pDim,pDim);
                tensor<double,2> delWN_delGN_n(pDim,pDim);
                tensor<double,4> d_delWP_delGP_n(pDim,pDim,pDim,pDim);
                tensor<double,4> d_delWN_delGN_n(pDim,pDim,pDim,pDim);
                double I3_n=jac*jac* jacGP_n*jacGP_n;
                double I3_SR_n=(jacGP_n)*jac;
                delWP_delGP_n=0.5*(lambda*log(I3_SR_n)-mu)/I3_n*delI3_delGP_n+0.5*mu*delI1_delGP_n;
                tensor<double,2> bend_str=curva-GN_n;
                tensor<double,2>  en_b=product(BIG_C,bend_str,{{0,0},{1,1}});
                delWN_delGN_n=const_bend*(-2.0)*en_b;
                d_delWP_delGP_n=(-2*lambda*log(I3_SR_n)+2*mu+lambda)/(4*I3_n*I3_n)*outer(delI3_delGP_n,delI3_delGP_n)+0.5*(lambda*log(I3_SR_n)-mu)/I3_n*d_delI3_delGP_n;
                d_delWN_delGN_n=const_bend*(2.0)*BIG_C;
                tensor<double,2> delWP_delGP_n_CC(pDim,pDim);
                tensor<double,2> delWN_delGN_n_CC(pDim,pDim);
                tensor<double,4> d_delWP_delGP_n_CC(pDim,pDim,pDim,pDim);
                tensor<double,4> d_delWN_delGN_n_CC(pDim,pDim,pDim,pDim);
                delWP_delGP_n_CC=product(metricR.inv(),product(delWP_delGP_n, metricR.inv(), {{1, 0}}), {{1, 0}});
                delWN_delGN_n_CC=product(metricR,product(delWN_delGN_n, metricR, {{1, 0}}), {{1, 0}});
                d_delWP_delGP_n_CC=product(metricR.inv(),product(d_delWP_delGP_n, metricR.inv(), {{1, 0}}), {{1, 0}}).transpose({0,3,1,2});
                d_delWN_delGN_n_CC=product(metricR,product(d_delWN_delGN_n, metricR, {{1, 0}}), {{1, 0}}).transpose({0,3,1,2});
                double F_11_NR{0.0};
                double F_12_NR{0.0};
                double F_22_NR{0.0};
                double F_GN_11_NR{0.0};
                double F_GN_12_NR{0.0};
                double F_GN_22_NR{0.0};
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
                tensor<double,1> GP_voigt(3);
                GP_voigt(0)=GP_11_NR;
                GP_voigt(1)=GP_12_NR;
                GP_voigt(2)=GP_22_NR;
                tensor<double,1> GP_F_voigt(3);
                GP_F_voigt(0)=F_11_NR;
                GP_F_voigt(1)=F_12_NR;
                GP_F_voigt(2)=F_22_NR;
                tensor<double,2> GP_jac_voigt(3,3);
                double  J_11_NR_11=1.0-dt_eta*d_delWP_delGP_n_CC(0,0,0,0);
                double  J_11_NR_12=0.0-dt_eta*(d_delWP_delGP_n_CC(0,0,0,1)+d_delWP_delGP_n_CC(0,0,1,0));
                double  J_11_NR_22=0.0-dt_eta*d_delWP_delGP_n_CC(0,0,1,1);
                double  J_12_NR_11=0.0-dt_eta*d_delWP_delGP_n_CC(0,1,0,0);
                double  J_12_NR_12=1.0-dt_eta*(d_delWP_delGP_n_CC(0,1,0,1)+d_delWP_delGP_n_CC(0,1,1,0));
                double  J_12_NR_22=0.0-dt_eta*d_delWP_delGP_n_CC(0,1,1,1);
                double J_22_NR_11=0.0-dt_eta*d_delWP_delGP_n_CC(1,1,0,0);
                double J_22_NR_12=0.0-dt_eta*(d_delWP_delGP_n_CC(1,1,0,1)+d_delWP_delGP_n_CC(1,1,1,0));
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
                tensor<double,1> GN_voigt(3);
                GN_voigt(0)=GN_11_NR;
                GN_voigt(1)=GN_12_NR;
                GN_voigt(2)=GN_22_NR;
                tensor<double,1> GN_F_voigt(3);
                GN_F_voigt(0)=F_GN_11_NR;
                GN_F_voigt(1)=F_GN_12_NR;
                GN_F_voigt(2)=F_GN_22_NR;
                tensor<double,2> GN_jac_voigt(3,3);
                double  J_GN_11_NR_11=1.0-dt_eta_k*d_delWN_delGN_n_CC(0,0,0,0);
                double  J_GN_11_NR_12=0.0-dt_eta_k*(d_delWN_delGN_n_CC(0,0,0,1)+d_delWN_delGN_n_CC(0,0,1,0));
                double  J_GN_11_NR_22=0.0-dt_eta_k*d_delWN_delGN_n_CC(0,0,1,1);
                double  J_GN_12_NR_11=0.0-dt_eta_k*d_delWN_delGN_n_CC(0,1,0,0);
                double  J_GN_12_NR_12=1.0-dt_eta_k*(d_delWN_delGN_n_CC(0,1,0,1)+d_delWN_delGN_n_CC(0,1,1,0));
                double  J_GN_12_NR_22=0.0-dt_eta_k*d_delWN_delGN_n_CC(0,1,1,1);
                double J_GN_22_NR_11=0.0-dt_eta_k*d_delWN_delGN_n_CC(1,1,0,0);
                double J_GN_22_NR_12=0.0-dt_eta_k*(d_delWN_delGN_n_CC(1,1,0,1)+d_delWN_delGN_n_CC(1,1,1,0));
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
                double X_11_NR{0.0};
                double X_12_NR{0.0};
                double X_22_NR{0.0};
                double X_GN_11_NR{0.0};
                double X_GN_12_NR{0.0};
                double X_GN_22_NR{0.0};
                tensor<double,1>delX(3);
                delX=-GP_jac_voigt.inv()*GP_F_voigt;
                X_11_NR=delX(0);
                X_12_NR=delX(1);
                X_22_NR=delX(2);
                tensor<double,1>delX_GN(3);
                delX_GN=-GN_jac_voigt.inv()*GN_F_voigt;
                X_GN_11_NR=delX_GN(0);
                X_GN_12_NR=delX_GN(1);
                X_GN_22_NR=delX_GN(2);
                // update relative to the size of each tensor (shear components and the curvatures
                // can be ~0 or negative): metric GP by its diagonal, curvature GN by max(|kbar|,|k|)
                auto rel3 = [](double x1, double x2, double x3, double sc)
                { sc = std::max(sc, 1.0e-12); return std::max(fabs(x1), std::max(fabs(x2), fabs(x3))) / sc; };
                const double scP = std::max(fabs(GP_11_NR), fabs(GP_22_NR));
                const double scN = std::max(std::max(fabs(GN_11_NR), fabs(GN_22_NR)), std::max(fabs(curva(0,0)), fabs(curva(1,1))));
                std::vector<double> values2 = { rel3(X_11_NR, X_12_NR, X_22_NR, scP), rel3(X_GN_11_NR, X_GN_12_NR, X_GN_22_NR, scN) };
                error2=*std::max_element(values2.begin(), values2.end());
                GP_n(0,0)=GP_11_NR;
                GP_n(0,1)=GP_12_NR;
                GP_n(1,0)=GP_12_NR;
                GP_n(1,1)=GP_22_NR;
                GN_n(0,0)=GN_11_NR;
                GN_n(0,1)=GN_12_NR;
                GN_n(1,0)=GN_12_NR;
                GN_n(1,1)=GN_22_NR;
                if (nr_step == 29 && (error > 0.0000001 || error2 > 0.0000001))
                cout << "WARNING LS_ReacDif: local Newton for G, kbar not converged in 30 iterations (residual " << error << ", relative update " << error2 << ")" << endl;
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
    double diss11=fric2*  product(product(GP_new-GP_old,metricR,{{1,0}}),product(metricR,GP_new-GP_old,{{1,0}}),{{0,0},{1,1}})/(deltat*deltat)*jacPRR;
    double diss22= fric3* product(product(GN_new-GN_old,imetricR,{{1,0}}),product(imetricR,GN_new-GN_old,{{1,0}}),{{0,0},{1,1}})/(deltat*deltat)*jacPRR;
    auxiliary_b[9] =diss11/jacPRR;
    auxiliary_b[10] =diss22/jacPRR;
    fillStr.addToGlobalIntegral("Diss1",diss11);
    fillStr.addToGlobalIntegral("Diss2",diss22);
}

void LS_ED(hiperlife::FillStructure& fillStr)
{
    using namespace hiperlife;
    using namespace hiperlife::Tensor;
    auto& subFill = (fillStr)["EDHand"];
    int numDOFs = subFill.numDOFs;
    int eNN     = subFill.eNN;
    wrapper<double,1>    bf(subFill.nborBFsDers(0),eNN);
    wrapper<double,2>  Bk1(fillStr.Bk(0).data(),eNN,numDOFs);
    wrapper<double,4>  Ak1(fillStr.Ak(0,0).data(),eNN,numDOFs,eNN,numDOFs);
    double* auxiliary_a = &fillStr.paramStr->a_aux[(subFill.loc_elemID*subFill.cubaInfo.iPts+subFill.kPt)*kNA];
    double trace_sigma= auxiliary_a[56]  ;
    double dev_sigma=auxiliary_a[57];
    double jacR= auxiliary_a[19];
    double nxx= auxiliary_a[27];
    double nyy= auxiliary_a[28];
    double nzz= auxiliary_a[29];
    double height1= auxiliary_a[30];
    double EA= auxiliary_a[48];
    double EB= auxiliary_a[49];
    double PA= auxiliary_a[50];
    double *auxiliary_b = &fillStr.paramStr->b_aux[(subFill.loc_elemID * subFill.cubaInfo.iPts + subFill.kPt)*kNB];
    double DA= auxiliary_b[9];
    double DB= auxiliary_b[10];
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
    Ak1(all,0,all,0) = jacR * outer(bf, bf) ;
    Ak1(all,1,all,1) =jacR * outer(bf, bf);
    Ak1(all,2,all,2) = jacR * outer(bf, bf) ;
    Ak1(all,3,all,3) =jacR * outer(bf, bf);
    Ak1(all,4,all,4) = jacR * outer(bf, bf) ;
    Ak1(all,5,all,5) =jacR * outer(bf, bf);
    Ak1(all,6,all,6) = jacR * outer(bf, bf) ;
    Ak1(all,7,all,7) = jacR * outer(bf, bf) ;
    Ak1(all,8,all,8) = jacR * outer(bf, bf) ;
    Ak1(all,9,all,9) = jacR * outer(bf, bf) ;
    Ak1(all,10,all,10) = jacR * outer(bf, bf) ;
}


// ============================================================================
// restart files of the Gauss-point state b_aux (G, kbar and output densities)
//   <prefix>.<rank>.gp, binary: int64 nElem, gPts, kNB; then per element
//   centroid (3 doubles, mean of its node coordinates) + gPts*kNB doubles.
// On reading, every saved element is matched to the local element with the
// nearest centroid (within tol): the text restart files of the mesh round the
// coordinates, and the local element order need not be the same.
// ============================================================================
#include <array>
#include <map>
#include <cstdio>
#include <cmath>
namespace
{
    std::array<double,3> ElemCentroid(hiperlife::SmartPtr<hiperlife::DistributedMesh>& mesh, int e)
    {
        using namespace hiperlife;
        std::vector<double> c = mesh->elemNborNodeCoords(e, IndexType::Local);
        const int n = int(c.size()) / 3;
        std::array<double,3> s = {0.0, 0.0, 0.0};
        for (int i = 0; i < n; i++) for (int d = 0; d < 3; d++) s[d] += c[3*i+d];
        for (int d = 0; d < 3; d++) s[d] /= n;
        return s;
    }
}

void GaussPointStateWrite(hiperlife::SmartPtr<hiperlife::DistributedMesh>& mesh,
                          const std::vector<double>& b_aux, int gPts, const std::string& prefix)
{
    using namespace hiperlife;
    const int nE = mesh->loc_nElem();
    std::string fname = prefix + "." + std::to_string(MyRank()) + ".gp";
    FILE* f = std::fopen(fname.c_str(), "wb");
    if (!f) { std::cout << "WARNING: cannot write " << fname << std::endl; return; }
    long long hdr[3] = {nE, gPts, kNB};
    std::fwrite(hdr, sizeof(long long), 3, f);
    const int blk = gPts * kNB;
    for (int e = 0; e < nE; e++)
    {
        auto c = ElemCentroid(mesh, e);
        std::fwrite(c.data(), sizeof(double), 3, f);
        std::fwrite(&b_aux[size_t(e)*blk], sizeof(double), blk, f);
    }
    std::fclose(f);
}

bool GaussPointStateRead(hiperlife::SmartPtr<hiperlife::DistributedMesh>& mesh,
                         std::vector<double>& b_aux, int gPts, const std::string& prefix, double tol)
{
    using namespace hiperlife;
    const int nE = mesh->loc_nElem();
    const int blk = gPts * kNB;
    int ok = 1;
    std::string fname = prefix + "." + std::to_string(MyRank()) + ".gp";
    FILE* f = std::fopen(fname.c_str(), "rb");
    long long hdr[3] = {0,0,0};
    if (!f || std::fread(hdr, sizeof(long long), 3, f) != 3 || hdr[0] != nE || hdr[1] != gPts || hdr[2] != kNB)
        ok = 0;
    if (ok)
    {
        const double cell = tol;     // hash grid, neighbouring cells searched
        auto key = [&](double x, double y, double z) {
            return std::array<long long,3>{(long long)std::floor(x/cell), (long long)std::floor(y/cell), (long long)std::floor(z/cell)}; };
        std::vector<std::array<double,3>> cen(nE);
        std::map<std::array<long long,3>, std::vector<int>> grid;
        for (int e = 0; e < nE; e++) { cen[e] = ElemCentroid(mesh, e); grid[key(cen[e][0], cen[e][1], cen[e][2])].push_back(e); }
        std::vector<double> rec(3 + blk);
        std::vector<int> seen(nE, 0);
        if (b_aux.size() < size_t(nE)*blk) b_aux.assign(size_t(nE)*blk, 0.0);
        for (int r = 0; r < nE && ok; r++)
        {
            if (std::fread(rec.data(), sizeof(double), rec.size(), f) != rec.size()) { ok = 0; break; }
            auto k = key(rec[0], rec[1], rec[2]);
            int best = -1; double bestD2 = tol*tol;
            for (long long dx = -1; dx <= 1; dx++) for (long long dy = -1; dy <= 1; dy++) for (long long dz = -1; dz <= 1; dz++)
            {
                auto it = grid.find({k[0]+dx, k[1]+dy, k[2]+dz});
                if (it == grid.end()) continue;
                for (int e : it->second)
                {
                    double d2 = 0.0;
                    for (int d = 0; d < 3; d++) d2 += (cen[e][d]-rec[d])*(cen[e][d]-rec[d]);
                    if (d2 < bestD2) { bestD2 = d2; best = e; }
                }
            }
            if (best < 0 || seen[best]) { ok = 0; break; }
            std::copy(rec.begin() + 3, rec.end(), b_aux.begin() + size_t(best)*blk);
            seen[best] = 1;
        }
        for (int e = 0; e < nE; e++) if (!seen[e]) ok = 0;
    }
    if (f) std::fclose(f);
    int okAll = 0;
    MPI_Allreduce(&ok, &okAll, 1, MPI_INT, MPI_MIN, MpiComm());
    return okAll == 1;
}
