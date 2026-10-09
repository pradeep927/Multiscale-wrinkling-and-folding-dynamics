// Element fillings of the phenomenological viscoelastic shell (see header).

#include "hl_FillStructure.h"
#include "hl_Geometry.h"
#include "hl_SurfLagrParam.h"
#include "hl_Tensor.h"
#include <hl_LinearSolver_Direct_MUMPS.h>
#include <algorithm>

#include "AuxViscousInextMembrane.h"

// Isotropic 4th-order tensor C^{abcd} = lambda a^ab a^cd + mu (a^ac a^bd + a^ad a^bc)
// built on the reference inverse metric a0 = g_R^{-1} (bending, Eq. 47/56).
static void IsotropicTensor(hiperlife::Tensor::tensor<double,4>& C,
                            hiperlife::Tensor::tensor<double,2>& a0,
                            double lambda, double mu)
{
    C(0,0,0,0)=(lambda*a0(0,0)*a0(0,0)+mu*(a0(0,0)*a0(0,0)+a0(0,0)*a0(0,0)));
    C(0,0,0,1)=(lambda*a0(0,0)*a0(0,1)+mu*(a0(0,1)*a0(0,0)+a0(0,0)*a0(0,1)));
    C(0,0,1,0)=(lambda*a0(0,0)*a0(1,0)+mu*(a0(0,0)*a0(0,1)+a0(0,1)*a0(0,0)));
    C(0,0,1,1)=(lambda*a0(0,0)*a0(1,1)+mu*(a0(0,1)*a0(0,1)+a0(0,1)*a0(0,1)));
    C(0,1,0,0)=(lambda*a0(0,1)*a0(0,0)+mu*(a0(0,0)*a0(1,0)+a0(0,0)*a0(1,0)));
    C(0,1,0,1)=(lambda*a0(0,1)*a0(0,1)+mu*(a0(0,1)*a0(1,0)+a0(0,0)*a0(1,1)));
    C(0,1,1,0)=(lambda*a0(0,1)*a0(1,0)+mu*(a0(0,0)*a0(1,1)+a0(0,1)*a0(1,0)));
    C(0,1,1,1)=(lambda*a0(0,1)*a0(1,1)+mu*(a0(0,1)*a0(1,1)+a0(0,1)*a0(1,1)));

    C(1,0,0,0)=(lambda*a0(1,0)*a0(0,0)+mu*(a0(1,0)*a0(0,0)+a0(1,0)*a0(0,0)));
    C(1,0,0,1)=(lambda*a0(1,0)*a0(0,1)+mu*(a0(1,1)*a0(0,0)+a0(1,0)*a0(0,1)));
    C(1,0,1,0)=(lambda*a0(1,0)*a0(1,0)+mu*(a0(1,0)*a0(0,1)+a0(1,1)*a0(0,0)));
    C(1,0,1,1)=(lambda*a0(1,0)*a0(1,1)+mu*(a0(1,1)*a0(0,1)+a0(1,1)*a0(0,1)));
    C(1,1,0,0)=(lambda*a0(1,1)*a0(0,0)+mu*(a0(1,0)*a0(1,0)+a0(1,0)*a0(1,0)));
    C(1,1,0,1)=(lambda*a0(1,1)*a0(0,1)+mu*(a0(1,1)*a0(1,0)+a0(1,0)*a0(1,1)));
    C(1,1,1,0)=(lambda*a0(1,1)*a0(1,0)+mu*(a0(1,0)*a0(1,1)+a0(1,1)*a0(1,0)));
    C(1,1,1,1)=(lambda*a0(1,1)*a0(1,1)+mu*(a0(1,1)*a0(1,1)+a0(1,1)*a0(1,1)));
}

// =============================================================================
// LS: mechanics of the shell + enclosed-volume constraint (Newton problem)
// =============================================================================
void LS(hiperlife::FillStructure& fillStr)
{
    using namespace hiperlife;
    using namespace hiperlife::Tensor;

    // ---- [1] inputs ---------------------------------------------------------
    auto &subFill = (fillStr)["posDHand"];
    int nDim = subFill.nDim;
    int pDim = subFill.pDim;
    int eNN = subFill.eNN;
    int numDOFs = subFill.numDOFs;

    tensor<double, 2,false> nborDOFs0(subFill.nborDOFs0.data(), eNN, numDOFs);
    tensor<double, 2,false> nborDOFs(subFill.nborDOFs.data(), eNN, numDOFs);
    tensor<double, 2,false> nborCoords(subFill.nborCoords.data(), eNN, nDim);

    tensor<double, 1,false> bf(subFill.nborBFs(), eNN);
    tensor<double, 2,false> Dbf(subFill.nborBFsGrads(), eNN, pDim);
    tensor<double, 3,false> DDbf(subFill.nborBFsHess(), eNN, pDim, pDim);

    // global unknown: volume Lagrange multiplier
    auto &g_subFill = (fillStr)["gloDHand"];
    int g_numDOFs = g_subFill.numDOFs;
    tensor<double, 1,false> gDOFs(g_subFill.nborDOFs.data(), g_numDOFs);
    double pressure = gDOFs(0);

    // parameters (energy-like coefficients are multiplied by deltat)
    double deltat    = fillStr.getRealParameter(MembParams::deltat);
    double fric      = fillStr.getRealParameter(MembParams::fric);
    double factor    = fillStr.getRealParameter(MembParams::factor);
    double thick     = fillStr.getRealParameter(MembParams::thick);
    double Kconf     = fillStr.getRealParameter(MembParams::Kconf) * deltat;
    int    timestep  = fillStr.getIntParameter(MembParams::timestep);
    double width     = fillStr.getRealParameter(MembParams::width);
    double fric_fact = fillStr.getRealParameter(MembParams::fric_fact);
    double height_in = fillStr.getRealParameter(MembParams::height_in);
    int control_fric = fillStr.getIntParameter(MembParams::control_fric);
    int gap          = fillStr.getIntParameter(MembParams::gap);
    int vnstep       = fillStr.getIntParameter(MembParams::vnstep);
    int fric_start   = gap+vnstep;

    double fact_elastic = fillStr.getRealParameter(MembParams::fact_elastic);
    double gamma  = fillStr.getRealParameter(MembParams::gamma) * deltat;
    double lambda = fillStr.getRealParameter(MembParams::lambda)*deltat;
    double mu     = fillStr.getRealParameter(MembParams::mu)*deltat;

    // outputs
    tensor<double, 2,false> Bp(fillStr.Bk(0).data(), eNN, numDOFs);
    tensor<double, 1,false> Bg(fillStr.Bk(1).data(), g_numDOFs);
    tensor<double, 4,false> App(fillStr.Ak(0, 0).data(), eNN, numDOFs, eNN, numDOFs);
    tensor<double, 3,false> Apg(fillStr.Ak(0, 1).data(), eNN, numDOFs, g_numDOFs);
    tensor<double, 3,false> Agp(fillStr.Ak(1, 0).data(), g_numDOFs, eNN, numDOFs);

    // tiny curvature offset at the first step so that the (flat) reference
    // curvature tensor is invertible
    double eps=0.0000001;

    // ---- [2] geometry ---------------------------------------------------------
    // reference (flat disc)
    tensor<double,1> xRv(nDim), xRu(nDim), xRuu(nDim), xRuv(nDim), xRvv(nDim), normalR(nDim);
    tensor<double,2> metricR(pDim,pDim), curvatureR(pDim,pDim);
    tensor<double,3> cristSymR(pDim,pDim,pDim);
    SurfLagrParam::ChristoffelSymbols(cristSymR, curvatureR, metricR, normalR, xRu, xRv, xRuu, xRvv, xRuv, eNN, nborCoords, Dbf, DDbf);
    if(timestep<1)
    {
        curvatureR(0,0)=curvatureR(0,0)+eps;
        curvatureR(1,1)=curvatureR(1,1)+eps;
    }

    // previous time step (normal and metric are lagged in drag and thickness)
    tensor<double,1> xu_n(nDim), xv_n(nDim), xuu_n(nDim), xuv_n(nDim), xvv_n(nDim), normal_n(nDim);
    tensor<double,2> metric_n(pDim,pDim), curva_n(pDim,pDim);
    tensor<double,3> christsym_n(pDim,pDim,pDim);
    SurfLagrParam::ChristoffelSymbols(christsym_n, curva_n, metric_n, normal_n, xu_n, xv_n, xuu_n, xvv_n, xuv_n, eNN, nborDOFs0, Dbf, DDbf);

    tensor<double,1> x_n = bf * nborDOFs0;
    tensor<double,1> xR = bf * nborCoords;
    tensor<double,2> imetricR = metricR.inv();
    double jacR = sqrt(metricR.det());
    double jac_n = sqrt(metric_n.det());
    double xnormal_n = x_n * normal_n;
    tensor<double,2> icurvatureR = curvatureR.inv();

    // current (unknown) configuration
    tensor<double,1> xu(nDim), xv(nDim), xuu(nDim), xuv(nDim), xvv(nDim), normal(nDim);
    tensor<double,2> metric(pDim,pDim), curva(pDim,pDim);
    tensor<double,3> christsym(pDim,pDim,pDim);
    SurfLagrParam::ChristoffelSymbols(christsym, curva, metric, normal, xu, xv, xuu, xvv, xuv, eNN, nborDOFs, Dbf, DDbf);
    if(timestep<1)
    {
        curva(0,0)=curva(0,0)+eps;
        curva(1,1)=curva(1,1)+eps;
    }

    tensor<double,1> x = bf * nborDOFs;
    tensor<double,2> imetric = metric.inv();
    double jac = sqrt(metric.det());
    double xnormal = x * normal;

    // optional: drag increases near the substrate (z < thick) late in deflation
    if (timestep > fric_start+control_fric)
    {
        fric = fric * fric_fact + (fric * fric_fact - fric) * tanh(width* (thick - x_n(2)));
    }

    tensor<double, 2> Id(3, 3);
    Id(0,0)=1; Id(0,1)=0; Id(0,2)=0;
    Id(1,0)=0; Id(1,1)=1; Id(1,2)=0;
    Id(2,0)=0; Id(2,1)=0; Id(2,2)=1;

    // first and second derivatives w.r.t. nodal positions
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

    // thickness used for bending: volume conservation of the tissue, h = h0 / J_n
    double jacC_n = jac_n/jacR;
    double thickness=height_in/jacC_n;

    // ---- [3] internal variables: material metric G^{IJ}, material curvature kbar
    double *auxiliary_b = &fillStr.paramStr->b_aux[(subFill.loc_elemID * subFill.cubaInfo.iPts + subFill.kPt) * N_AUXB];
    tensor<double, 2> GP(pDim, pDim);
    tensor<double, 2> GN(pDim, pDim);
    GP(0,0) = auxiliary_b[B_GP+0]; GP(0,1) = auxiliary_b[B_GP+1];
    GP(1,0) = auxiliary_b[B_GP+1]; GP(1,1) = auxiliary_b[B_GP+2];
    GN(0,0) = auxiliary_b[B_GN+0]; GN(0,1) = auxiliary_b[B_GN+1];
    GN(1,0) = auxiliary_b[B_GN+1]; GN(1,1) = auxiliary_b[B_GN+2];
    tensor<double, 2> iGP = GP.inv();

    // ---- [4] relaxing membrane energy w_mem(C,G), Eq. (55) ----------------------
    //   I1 = G^{IJ} g_IJ,  I3_SR = sqrt(det G) * J
    double jacGP = (GP.det());
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

    double EA     = (0.5*lambda*log(I3_SR)*log(I3_SR)-mu*log(I3_SR)+0.5*mu*(I1-2) )*2.0;
    double Evisco = (0.5*lambda*log(I3_SR)*log(I3_SR)-mu*log(I3_SR)+0.5*mu*(I1-2) )*jacPRR*2.0;
    Bp(all,range(0,2)) += ((lambda*log(I3_SR)-mu)/I3_SR*d_I3_SR+0.5*mu*d_I1)*jacPRR*2.0;
    App(all,range(0,2),all,range(0,2)) += ((lambda*log(I3_SR)-mu)/I3_SR*dd_I3_SR+0.5*mu*dd_I1)*jacPRR*2.0;
    App(all,range(0,2),all,range(0,2)) += (lambda*(1-log(I3_SR))+mu)/(I3_SR*I3_SR)*outer(d_I3_SR,d_I3_SR)*jacPRR*2.0;

    // ---- [5] non-relaxing pre-tension elasticity w_pre(C), Eq. (54) -------------
    //   same neo-Hookean form relative to the flat reference (G = g_R^{-1})
    double jacGP_ref = 1/(jacR*jacR);
    double I1_ref    = product(imetricR,metric ,{{0,0},{1,1}}) ;
    double I3_SR_ref = sqrt(jacGP_ref)*jac;

    tensor<double,2> d_I1_ref(eNN,numDOFs);
    tensor<double,2> d_I3_SR_ref(eNN,numDOFs);
    d_I1_ref=product(imetricR,d_metric ,{{0,2},{1,3}}) ;
    d_I3_SR_ref=sqrt(jacGP_ref)*d_jac;

    tensor<double,4> dd_I1_ref(eNN,numDOFs,eNN,numDOFs);
    tensor<double,4> dd_I3_SR_ref(eNN,numDOFs,eNN,numDOFs);
    dd_I1_ref=product(imetricR,dd_metric ,{{0,4},{1,5}}) ;
    dd_I3_SR_ref=sqrt(jacGP_ref)*dd_jac;

    double mu_ela     = fillStr.getRealParameter(MembParams::mu_pre)*deltat*fact_elastic;
    double lambda_ela = fillStr.getRealParameter(MembParams::lambda_pre)*fact_elastic*deltat;

    double E_small = (0.5*lambda_ela*log(I3_SR_ref)*log(I3_SR_ref)-mu_ela*log(I3_SR_ref)+0.5*mu_ela*(I1_ref-2) )*jacPRR;
    Bp(all,range(0,2)) += ((lambda_ela*log(I3_SR_ref)-mu_ela)/I3_SR_ref*d_I3_SR_ref+0.5*mu_ela*d_I1_ref)*jacPRR;
    App(all,range(0,2),all,range(0,2)) += ((lambda_ela*log(I3_SR_ref)-mu_ela)/I3_SR_ref*dd_I3_SR_ref+0.5*mu_ela*dd_I1_ref)*jacPRR;
    App(all,range(0,2),all,range(0,2)) += (lambda_ela*(1-log(I3_SR_ref))+mu_ela)/(I3_SR_ref*I3_SR_ref)*outer(d_I3_SR_ref,d_I3_SR_ref)*jacPRR;

    // ---- [6] relaxing bending energy w_bend(k,kbar), Eq. (56) -------------------
    //   E = jacR h^2/4 C^{abcd} (k-kbar)_ab (k-kbar)_cd, C built with (lambda, mu)
    tensor<double,2> a0(pDim,pDim);
    a0(0,0)=imetricR(0,0); a0(0,1)=imetricR(0,1);
    a0(1,0)=imetricR(1,0); a0(1,1)=imetricR(1,1);
    tensor<double,4> BIG_C(pDim,pDim,pDim,pDim);
    IsotropicTensor(BIG_C, a0, lambda, mu);

    tensor<double,2> bend_str=curva-GN;
    tensor<double,4> d_bend_str=d_curva;
    tensor<double,6> dd_bend_str=dd_curva;

    double const_bend=thickness*thickness/4.0;

    tensor<double,2> en_b(2,2);
    en_b=product(BIG_C,bend_str,{{0,0},{1,1}});
    double E_bend{};
    E_bend+=jacR*const_bend*(product(en_b,bend_str,{{0,0},{1,1}}));
    double EB=E_bend;

    tensor<double,4> d_bend1(eNN,numDOFs,2,2);
    tensor<double,4> d_bend2(eNN,numDOFs,2,2);
    d_bend1=product(d_bend_str,BIG_C,{{2,0},{3,1}});
    d_bend2=product(d_bend_str,BIG_C,{{2,2},{3,3}});
    Bp(all,range(0,2)) += jacR*const_bend*product(d_bend1,bend_str,{{2,0},{3,1}}) ;
    Bp(all,range(0,2)) += jacR*const_bend*product(d_bend2,bend_str,{{2,0},{3,1}}) ;

    tensor<double,6> dd_bend1(eNN,numDOFs,eNN,numDOFs,2,2);
    tensor<double,6> dd_bend2(eNN,numDOFs,eNN,numDOFs,2,2);
    dd_bend1=product(dd_bend_str,BIG_C,{{4,0},{5,1}});
    dd_bend2=product(dd_bend_str,BIG_C,{{4,2},{5,3}});
    App(all,range(0,2),all,range(0,2)) += jacR*const_bend*product(dd_bend1,bend_str,{{4,0},{5,1}}) ;
    App(all,range(0,2),all,range(0,2)) += jacR*const_bend*product(dd_bend2,bend_str,{{4,0},{5,1}}) ;

    tensor<double,8> d_bend_d_bend(eNN,numDOFs,2,2,eNN,numDOFs,2,2);
    d_bend_d_bend=outer(d_bend_str,d_bend_str);
    tensor<double,4> d_bend_d_bend_Cpqrs(eNN,numDOFs,eNN,numDOFs);
    tensor<double,4> d_bend_d_bend_Crspq(eNN,numDOFs,eNN,numDOFs);
    d_bend_d_bend_Cpqrs= product( d_bend_d_bend,BIG_C,{{2,2},{3,3},{6,0},{7,1}});
    d_bend_d_bend_Crspq= product( d_bend_d_bend,BIG_C,{{2,0},{3,1},{6,2},{7,3}});
    App(all,range(0,2),all,range(0,2)) += jacR*const_bend*d_bend_d_bend_Cpqrs;
    App(all,range(0,2),all,range(0,2)) += jacR*const_bend*d_bend_d_bend_Crspq ;

    // ---- [7] drag: fric/2 |x - x_n|^2 J_n ---------------------------------------
    double Efric = 0.5*fric*jac_n*product(x-x_n,x-x_n,{{0,0}});
    Bp(all,range(0,2))                 += fric*jac_n*outer(bf,x-x_n);
    App(all,range(0,2),all,range(0,2)) += fric*jac_n*outer(outer(bf,Id),bf).transpose({0,1,3,2});

    // ---- [8] confinement above z = 0 ----------------------------------------------
    tensor<double,1> Id3(3);
    Id3(0)=0; Id3(1)=0; Id3(2)=1;
    tensor<double,2> d_z(eNN,numDOFs);
    d_z=outer(bf,Id3);

    Bp(all,range(0,2))                 += jac*d_z*dConfinementPotential(x, Kconf);
    Bp(all,range(0,2))                 += d_jac*ConfinementPotential(x, Kconf);
    App(all,range(0,2),all,range(0,2)) += outer(d_jac,d_z)*dConfinementPotential(x, Kconf);
    App(all,range(0,2),all,range(0,2)) += jac*outer(d_z,d_z)*ddConfinementPotential(x, Kconf);
    App(all,range(0,2),all,range(0,2)) += dd_jac*ConfinementPotential(x, Kconf);
    App(all,range(0,2),all,range(0,2)) += outer(d_z,d_jac)*dConfinementPotential(x, Kconf);

    // ---- [9] active pre-tension sigma_a J, Eq. (53) ---------------------------------
    double power1=gamma*(jac-jac_n)*1.0/(deltat*deltat);   // power of the tension (output)
    double E_tens1=0.0;                                      // kept as a column of globalIntegrals.dat
    Bp(all, range(0, 2)) += gamma*d_jac;
    App(all, range(0, 2), all, range(0, 2)) += gamma*dd_jac;

    // ---- [10] enclosed volume: (1/3) int x.n J = factor * reference area -------------
    //   (Lagrange multiplier "pressure"; the volume is measured relative to z = 0)
    double Evol=pressure/3.0 *(jac * xnormal - factor*jac_n*xnormal_n);
    Bp(all, range(0, 2)) += pressure / 3.0 * (d_jac * xnormal + jac * (d_normal * x) + jac * outer(bf, normal));
    Bg(0) += (jac * xnormal / 3) - factor * jacR * 1;
    App(all, range(0, 2), all, range(0, 2)) += pressure / 3.0 * (dd_jac * xnormal + outer(d_jac, d_normal * x) + outer(d_jac, outer(bf, normal)) +
                              outer(d_normal * x, d_jac) + outer(outer(bf, normal), d_jac) + jac * dd_normal * x +
                              jac * outer(d_normal, bf).transpose({0, 1, 3, 2}) +
                              jac * outer(bf, d_normal.transpose({2, 0, 1})));
    Apg(all, range(0, 2), 0) += 1.0 / 3.0 * (d_jac * xnormal + jac * d_normal * x + jac * outer(bf, normal));
    Agp = Apg.transpose({2,0,1});

    // ---- [11] global integrals (rates: energies / deltat) --------------------------
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

    // ---- [12] data for LS_ReacDif (relaxation) and LS_ED (output) ------------------
    // dW_mem/dG^{IJ}
    tensor<double,2> delI1_delGP(pDim,pDim);
    tensor<double,2> delI3_delGP(pDim,pDim);
    delI1_delGP=metric;
    delI3_delGP=jac*jac*jacGP*iGP;
    tensor<double,2> delWP_delGP(pDim,pDim);
    delWP_delGP=0.5*(lambda*log(I3_SR)-mu)/I3*delI3_delGP+0.5*mu*delI1_delGP;

    // membrane stress (output): sigma = 2 dW/dC / J
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

    double* aa = &fillStr.paramStr->a_aux[(subFill.loc_elemID*subFill.cubaInfo.iPts+subFill.kPt)*N_AUXA];
    aa[A_METRIC+0] = metric(0,0);   aa[A_METRIC+1] = metric(0,1);
    aa[A_METRIC+2] = metric(1,0);   aa[A_METRIC+3] = metric(1,1);
    aa[A_DWPDGP+0] = delWP_delGP(0,0); aa[A_DWPDGP+1] = delWP_delGP(0,1);
    aa[A_DWPDGP+2] = delWP_delGP(1,0); aa[A_DWPDGP+3] = delWP_delGP(1,1);
    aa[A_IMETRICR+0] = imetricR(0,0);  aa[A_IMETRICR+1] = imetricR(0,1);
    aa[A_IMETRICR+2] = imetricR(1,0);  aa[A_IMETRICR+3] = imetricR(1,1);
    aa[A_ICURVR+0] = icurvatureR(0,0); aa[A_ICURVR+1] = icurvatureR(0,1);
    aa[A_ICURVR+2] = icurvatureR(1,0); aa[A_ICURVR+3] = icurvatureR(1,1);
    aa[A_CURVA+0] = curva(0,0);   aa[A_CURVA+1] = curva(1,0);   // (symmetric)
    aa[A_CURVA+2] = curva(0,1);   aa[A_CURVA+3] = curva(1,1);
    aa[A_JACR]    = jacPRR;
    aa[A_JAC]     = jac;
    aa[A_CONSTB]  = const_bend;
    aa[A_NORMAL+0] = normal(0); aa[A_NORMAL+1] = normal(1); aa[A_NORMAL+2] = normal(2);
    aa[A_THICK]   = thickness;
    aa[A_EA]      = EA/deltat;
    aa[A_EB]      = EB/deltat;
    aa[A_PA]      = power1;
    aa[A_TRSIG]   = trace_sigma;
    aa[A_DEVSIG]  = dev_sigma;
}


// =============================================================================
// LS_ReacDif: relaxation of the material metric G and curvature kbar at every
// Gauss point (no global unknowns; run through FillGlobalIntegrals once per
// converged step). Eqs. (58)-(59):
//     dG^{IJ}/dt   = -(1/fric2) [dW_mem/dG]       (raised with g_R^{-1})
//     dkbar_IJ/dt  = -(1/fric3) [dW_bend/dkbar]   (lowered with g_R)
// forward_new = 1: explicit Euler; 0: implicit Euler with a local Newton.
// While timestep < time_target, G = g_R^{-1} and kbar = k_R (pure elasticity).
// =============================================================================
void LS_ReacDif(hiperlife::FillStructure& fillStr)
{
    using namespace hiperlife;
    using namespace hiperlife::Tensor;

    auto& subFill = (fillStr)["viscoDHand"];
    int pDim = subFill.pDim;

    double deltat      = fillStr.getRealParameter(MembParams::deltat);
    double fric2       = fillStr.getRealParameter(MembParams::fric2);
    int    time_target = fillStr.getIntParameter(MembParams::time_target);
    double fric3       = fillStr.getRealParameter(MembParams::fric3);
    int    timestep    = fillStr.getIntParameter(MembParams::timestep);
    double forward_new = fillStr.getRealParameter(MembParams::forward_new);
    double lambda      = fillStr.getRealParameter(MembParams::lambda)*deltat;
    double mu          = fillStr.getRealParameter(MembParams::mu) *deltat;

    // data of the converged mechanics (written by LS)
    double* aa = &fillStr.paramStr->a_aux[(subFill.loc_elemID*subFill.cubaInfo.iPts+subFill.kPt)*N_AUXA];
    tensor<double,2,false> metric(&aa[A_METRIC],2,2);
    tensor<double,2,false> delWP_delGP(&aa[A_DWPDGP],2,2);
    tensor<double,2,false> imetricR(&aa[A_IMETRICR],2,2);
    tensor<double,2,false> icurvatureR(&aa[A_ICURVR],2,2);
    tensor<double,2,false> curva(&aa[A_CURVA],2,2);
    tensor<double,2> curvatureR=icurvatureR.inv();
    tensor<double,2> metricR=imetricR.inv();
    double jac        = aa[A_JAC];
    double const_bend = aa[A_CONSTB];
    double jacPRR     = aa[A_JACR];

    tensor<double,2> a0(pDim,pDim);
    a0(0,0)=imetricR(0,0); a0(0,1)=imetricR(0,1);
    a0(1,0)=imetricR(1,0); a0(1,1)=imetricR(1,1);
    tensor<double,4> BIG_C(pDim,pDim,pDim,pDim);
    IsotropicTensor(BIG_C, a0, lambda, mu);

    double G11P{1.0}, G12P{0.0}, G22P{1.0}, G11N{1.0}, G12N{0}, G22N{1.0};

    double* auxiliary_b = &fillStr.paramStr->b_aux[(subFill.loc_elemID*subFill.cubaInfo.iPts+subFill.kPt)*N_AUXB];

    // state at the beginning of the step (for the dissipation)
    tensor<double,2> GP_old(2,2);
    tensor<double,2> GN_old(2,2);
    if(timestep<time_target)
    {
        GP_old(0,0) = 1.0*imetricR(0,0);
        GP_old(1,0) = 1.0*imetricR(0,1);
        GP_old(0,1) = 1.0*imetricR(0,1);
        GP_old(1,1) = 1.0*imetricR(1,1);
        GN_old(0,0) = curvatureR(0,0);
        GN_old(0,1) = curvatureR(1,0);
        GN_old(1,0) = curvatureR(1,0);
        GN_old(1,1) = curvatureR(1,1);
    }
    else
    {
        GP_old(0,0) = auxiliary_b[B_GP+0];
        GP_old(1,0) = auxiliary_b[B_GP+1];
        GP_old(0,1) = auxiliary_b[B_GP+1];
        GP_old(1,1) = auxiliary_b[B_GP+2];
        GN_old(0,0) = auxiliary_b[B_GN+0];
        GN_old(1,0) = auxiliary_b[B_GN+1];
        GN_old(0,1) = auxiliary_b[B_GN+1];
        GN_old(1,1) = auxiliary_b[B_GN+2];
    }

    double dt_eta   = -1/fric2;
    double dt_eta_k = -1/fric3;

    if (timestep<time_target)
    {
        G11P = 1.0*imetricR(0,0);
        G12P = 1.0*imetricR(0,1);
        G22P = 1.0*imetricR(1,1);
        G11N = curvatureR(0,0);
        G12N = curvatureR(1,0);
        G22N = curvatureR(1,1);
    }
    else if (forward_new>0.1)
    {
        // ---- explicit Euler ----
        tensor<double,2> delWP_delGP_n=delWP_delGP;
        tensor<double,2> delWP_delGP_CC(2,2);
        delWP_delGP_CC= product(metricR.inv(),product(delWP_delGP_n, metricR.inv(), {{1, 0}}), {{0, 0}});

        tensor<double,2> bend_str=curva-GN_old;
        tensor<double,2> en_b=product(BIG_C,bend_str,{{0,0},{1,1}});
        tensor<double,2> delWN_delGN=const_bend*(-2.0)*en_b;
        tensor<double,2> delWN_delGN_CC(2,2);
        delWN_delGN_CC= product(metricR,product(delWN_delGN,metricR, {{1, 0}}), {{0, 0}});

        G11P = auxiliary_b[B_GP+0] + dt_eta*delWP_delGP_CC(0,0);
        G12P = auxiliary_b[B_GP+1] + dt_eta*delWP_delGP_CC(0,1);
        G22P = auxiliary_b[B_GP+2] + dt_eta*delWP_delGP_CC(1,1);
        G11N = auxiliary_b[B_GN+0] + dt_eta_k*delWN_delGN_CC(0,0);
        G12N = auxiliary_b[B_GN+1] + dt_eta_k*delWN_delGN_CC(0,1);
        G22N = auxiliary_b[B_GN+2] + dt_eta_k*delWN_delGN_CC(1,1);
    }
    else
    {
        // ---- implicit Euler: local Newton for (G^11,G^12,G^22) and (kbar_11,kbar_12,kbar_22)
        tensor<double,2> GP_n(2,2);
        GP_n(0,0)=auxiliary_b[B_GP+0]; GP_n(0,1)=auxiliary_b[B_GP+1];
        GP_n(1,0)=auxiliary_b[B_GP+1]; GP_n(1,1)=auxiliary_b[B_GP+2];
        tensor<double,2> GN_n(2,2);
        GN_n(0,0)=auxiliary_b[B_GN+0]; GN_n(0,1)=auxiliary_b[B_GN+1];
        GN_n(1,0)=auxiliary_b[B_GN+1]; GN_n(1,1)=auxiliary_b[B_GN+2];

        int nr_step=0;
        double GP_11_NR = GP_n(0,0), GP_12_NR = GP_n(0,1), GP_22_NR = GP_n(1,1);
        double GN_11_NR = GN_n(0,0), GN_12_NR = GN_n(0,1), GN_22_NR = GN_n(1,1);

        double error=0.1;
        double error2=0.1;
        while ((error>0.0000001|| error2>0.0000001) && nr_step<30)
        {
            double jacGP_n = sqrt(GP_n.det());
            tensor<double,2> iGP_n = GP_n.inv();

            // dW_mem/dG (at the new G) and its derivative
            tensor<double,2> delI1_delGP_n(pDim,pDim);
            delI1_delGP_n = (metric);
            tensor<double,2> delI3_delGP_n(pDim,pDim);
            delI3_delGP_n=jac*jac*jacGP_n*jacGP_n*iGP_n;
            tensor<double,4> d_delI3_delGP_n(pDim,pDim,pDim,pDim);
            d_delI3_delGP_n=jac*jac*jacGP_n*jacGP_n*outer(iGP_n,iGP_n)-jac*jac*jacGP_n*jacGP_n*outer(iGP_n,iGP_n).transpose({0,2,3,1});

            double I3_n=jac*jac* jacGP_n*jacGP_n;
            double I3_SR_n=(jacGP_n)*jac;
            tensor<double,2> delWP_delGP_n(pDim,pDim);
            delWP_delGP_n=0.5*(lambda*log(I3_SR_n)-mu)/I3_n*delI3_delGP_n+0.5*mu*delI1_delGP_n;
            tensor<double,4> d_delWP_delGP_n(pDim,pDim,pDim,pDim);
            d_delWP_delGP_n=(-2*lambda*log(I3_SR_n)+2*mu+lambda)/(4*I3_n*I3_n)*outer(delI3_delGP_n,delI3_delGP_n)+0.5*(lambda*log(I3_SR_n)-mu)/I3_n*d_delI3_delGP_n;

            // dW_bend/dkbar (at the new kbar) and its derivative
            tensor<double,2> bend_str=curva-GN_n;
            tensor<double,2> en_b=product(BIG_C,bend_str,{{0,0},{1,1}});
            tensor<double,2> delWN_delGN_n(pDim,pDim);
            delWN_delGN_n=const_bend*(-2.0)*en_b;
            tensor<double,4> d_delWN_delGN_n(pDim,pDim,pDim,pDim);
            d_delWN_delGN_n=const_bend*(2.0)*BIG_C;

            // index raising (G) / lowering (kbar) with the reference metric
            tensor<double,2> delWP_delGP_n_CC(pDim,pDim);
            tensor<double,2> delWN_delGN_n_CC(pDim,pDim);
            tensor<double,4> d_delWP_delGP_n_CC(pDim,pDim,pDim,pDim);
            tensor<double,4> d_delWN_delGN_n_CC(pDim,pDim,pDim,pDim);
            delWP_delGP_n_CC=product(metricR.inv(),product(delWP_delGP_n, metricR.inv(), {{1, 0}}), {{1, 0}});
            delWN_delGN_n_CC=product(metricR,product(delWN_delGN_n, metricR, {{1, 0}}), {{1, 0}});
            d_delWP_delGP_n_CC=product(metricR.inv(),product(d_delWP_delGP_n, metricR.inv(), {{1, 0}}), {{1, 0}}).transpose({0,3,1,2});
            d_delWN_delGN_n_CC=product(metricR,product(d_delWN_delGN_n, metricR, {{1, 0}}), {{1, 0}}).transpose({0,3,1,2});

            // residuals  F = G - G_n - dt/eta * dW/dG
            double F_11_NR=GP_11_NR-auxiliary_b[B_GP+0]-dt_eta*delWP_delGP_n_CC(0,0);
            double F_12_NR=GP_12_NR-auxiliary_b[B_GP+1]-dt_eta*delWP_delGP_n_CC(0,1);
            double F_22_NR=GP_22_NR-auxiliary_b[B_GP+2]-dt_eta*delWP_delGP_n_CC(1,1);
            double F_GN_11_NR=GN_11_NR-auxiliary_b[B_GN+0]-dt_eta_k*delWN_delGN_n_CC(0,0);
            double F_GN_12_NR=GN_12_NR-auxiliary_b[B_GN+1]-dt_eta_k*delWN_delGN_n_CC(0,1);
            double F_GN_22_NR=GN_22_NR-auxiliary_b[B_GN+2]-dt_eta_k*delWN_delGN_n_CC(1,1);

            std::vector<double> values = {abs(F_11_NR), abs(F_12_NR), abs(F_22_NR), abs(F_GN_11_NR), abs(F_GN_12_NR), abs(F_GN_22_NR)};
            error=*std::max_element(values.begin(), values.end());

            // Newton update for G (Voigt 11,12,22)
            tensor<double,1> GP_voigt(3);
            GP_voigt(0)=GP_11_NR; GP_voigt(1)=GP_12_NR; GP_voigt(2)=GP_22_NR;
            tensor<double,1> GP_F_voigt(3);
            GP_F_voigt(0)=F_11_NR; GP_F_voigt(1)=F_12_NR; GP_F_voigt(2)=F_22_NR;
            tensor<double,2> GP_jac_voigt(3,3);
            GP_jac_voigt(0,0)=1.0-dt_eta*d_delWP_delGP_n_CC(0,0,0,0);
            GP_jac_voigt(0,1)=0.0-dt_eta*d_delWP_delGP_n_CC(0,0,0,1);
            GP_jac_voigt(0,2)=0.0-dt_eta*d_delWP_delGP_n_CC(0,0,1,1);
            GP_jac_voigt(1,0)=0.0-dt_eta*d_delWP_delGP_n_CC(0,1,0,0);
            GP_jac_voigt(1,1)=1.0-dt_eta*d_delWP_delGP_n_CC(0,1,0,1);
            GP_jac_voigt(1,2)=0.0-dt_eta*d_delWP_delGP_n_CC(0,1,1,1);
            GP_jac_voigt(2,0)=0.0-dt_eta*d_delWP_delGP_n_CC(1,1,0,0);
            GP_jac_voigt(2,1)=0.0-dt_eta*d_delWP_delGP_n_CC(1,1,0,1);
            GP_jac_voigt(2,2)=1.0-dt_eta*d_delWP_delGP_n_CC(1,1,1,1);
            GP_voigt=GP_voigt-GP_jac_voigt.inv()*GP_F_voigt;
            GP_11_NR=GP_voigt(0); GP_12_NR=GP_voigt(1); GP_22_NR=GP_voigt(2);

            // Newton update for kbar
            tensor<double,1> GN_voigt(3);
            GN_voigt(0)=GN_11_NR; GN_voigt(1)=GN_12_NR; GN_voigt(2)=GN_22_NR;
            tensor<double,1> GN_F_voigt(3);
            GN_F_voigt(0)=F_GN_11_NR; GN_F_voigt(1)=F_GN_12_NR; GN_F_voigt(2)=F_GN_22_NR;
            tensor<double,2> GN_jac_voigt(3,3);
            GN_jac_voigt(0,0)=1.0-dt_eta_k*d_delWN_delGN_n_CC(0,0,0,0);
            GN_jac_voigt(0,1)=0.0-dt_eta_k*d_delWN_delGN_n_CC(0,0,0,1);
            GN_jac_voigt(0,2)=0.0-dt_eta_k*d_delWN_delGN_n_CC(0,0,1,1);
            GN_jac_voigt(1,0)=0.0-dt_eta_k*d_delWN_delGN_n_CC(0,1,0,0);
            GN_jac_voigt(1,1)=1.0-dt_eta_k*d_delWN_delGN_n_CC(0,1,0,1);
            GN_jac_voigt(1,2)=0.0-dt_eta_k*d_delWN_delGN_n_CC(0,1,1,1);
            GN_jac_voigt(2,0)=0.0-dt_eta_k*d_delWN_delGN_n_CC(1,1,0,0);
            GN_jac_voigt(2,1)=0.0-dt_eta_k*d_delWN_delGN_n_CC(1,1,0,1);
            GN_jac_voigt(2,2)=1.0-dt_eta_k*d_delWN_delGN_n_CC(1,1,1,1);
            GN_voigt=GN_voigt-GN_jac_voigt.inv()*GN_F_voigt;
            GN_11_NR=GN_voigt(0); GN_12_NR=GN_voigt(1); GN_22_NR=GN_voigt(2);

            // relative size of the Newton increments
            tensor<double,1> delX(3);
            delX=-GP_jac_voigt.inv()*GP_F_voigt;
            tensor<double,1> delX_GN(3);
            delX_GN=-GN_jac_voigt.inv()*GN_F_voigt;
            double tol3=0.0000000001;
            std::vector<double> values2 = {
                abs(delX(0)/(GP_11_NR+tol3)), abs(delX(1)/(GP_12_NR+tol3)), abs(delX(2)/(GP_22_NR+tol3)),
                abs(delX_GN(0)/(GN_11_NR+tol3)), abs(delX_GN(1)/(GN_12_NR+tol3)), abs(delX_GN(2)/(GN_22_NR+tol3))};
            error2=*std::max_element(values2.begin(), values2.end());

            GP_n(0,0)=GP_11_NR; GP_n(0,1)=GP_12_NR; GP_n(1,0)=GP_12_NR; GP_n(1,1)=GP_22_NR;
            GN_n(0,0)=GN_11_NR; GN_n(0,1)=GN_12_NR; GN_n(1,0)=GN_12_NR; GN_n(1,1)=GN_22_NR;

            if(nr_step>25)
                cout<<"error warning: "<<endl;   // local Newton not converging
            nr_step=nr_step+1;
        }
        G11P = GP_11_NR; G12P = GP_12_NR; G22P = GP_22_NR;
        G11N = GN_11_NR; G12N = GN_12_NR; G22N = GN_22_NR;
    }

    auxiliary_b[B_GP+0] = G11P;
    auxiliary_b[B_GP+1] = G12P;
    auxiliary_b[B_GP+2] = G22P;
    auxiliary_b[B_GN+0] = G11N;
    auxiliary_b[B_GN+1] = G12N;
    auxiliary_b[B_GN+2] = G22N;

    // viscous dissipation of the two relaxations (output)
    tensor<double,2> GP_new(2,2);
    GP_new(0,0)=auxiliary_b[B_GP+0]; GP_new(0,1)=auxiliary_b[B_GP+1];
    GP_new(1,0)=auxiliary_b[B_GP+1]; GP_new(1,1)=auxiliary_b[B_GP+2];
    tensor<double,2> GN_new(2,2);
    GN_new(0,0)=auxiliary_b[B_GN+0]; GN_new(0,1)=auxiliary_b[B_GN+1];
    GN_new(1,0)=auxiliary_b[B_GN+1]; GN_new(1,1)=auxiliary_b[B_GN+2];

    double diss11=fric2* product(product(GP_new-GP_old,metricR,{{1,0}}),product(metricR,GP_new-GP_old,{{1,0}}),{{0,0},{1,1}})/(deltat*deltat)*jacPRR;
    double diss22=fric3* product(product(GN_new-GN_old,imetricR,{{1,0}}),product(imetricR,GN_new-GN_old,{{1,0}}),{{0,0},{1,1}})/(deltat*deltat)*jacPRR;

    auxiliary_b[B_DA] = diss11/jacPRR;
    auxiliary_b[B_DB] = diss22/jacPRR;

    fillStr.addToGlobalIntegral("Diss1",diss11);
    fillStr.addToGlobalIntegral("Diss2",diss22);
}

// =============================================================================
// LS_ED: L2 projection of Gauss-point quantities onto the nodes (output only)
// fields: nx ny nz height EA EB DA DB PA trace_sigma dev_sigma
// =============================================================================
void LS_ED(hiperlife::FillStructure & fillStr)
{
    using namespace hiperlife;
    using namespace hiperlife::Tensor;

    auto& subFill = (fillStr)["EDHand"];
    int numDOFs = subFill.numDOFs;
    int eNN     = subFill.eNN;

    tensor<double,1,false> bf(subFill.nborBFsDers(0),eNN);
    tensor<double,2,false> Bk1(fillStr.Bk(0).data(),eNN,numDOFs);
    tensor<double,4,false> Ak1(fillStr.Ak(0,0).data(),eNN,numDOFs,eNN,numDOFs);

    double* aa = &fillStr.paramStr->a_aux[(subFill.loc_elemID*subFill.cubaInfo.iPts+subFill.kPt)*N_AUXA];
    double* ab = &fillStr.paramStr->b_aux[(subFill.loc_elemID*subFill.cubaInfo.iPts+subFill.kPt)*N_AUXB];
    double jacR = aa[A_JACR];

    const double vals[11] = {aa[A_NORMAL+0], aa[A_NORMAL+1], aa[A_NORMAL+2], aa[A_THICK],
                             aa[A_EA], aa[A_EB], ab[B_DA], ab[B_DB], aa[A_PA],
                             aa[A_TRSIG], aa[A_DEVSIG]};
    for (int k = 0; k < 11; k++)
    {
        Bk1(all,k) = jacR * bf * (vals[k]);
        Ak1(all,k,all,k) = jacR * outer(bf, bf);
    }
}
