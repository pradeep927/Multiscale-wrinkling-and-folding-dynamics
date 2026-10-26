// AssemblyKernels.h -- three small helpers for the metric algebra of the shell element (LS).
//
// They replace generic tensor contractions (temporary tensors and index bookkeeping) by direct
// loops over contiguous storage, with the same arithmetic order as the original expressions.
//
// Storage assumptions
//   * hiperlife tensors are row-major (last index fastest) and contiguous;
//   * in every metric-type tensor the LAST TWO indices are the 2x2 surface indices (a,b);
//     the leading indices (e,k) or (e,k,f,l) are node / dof indices.
//
// Helpers (I, J = leading index groups)
//   contract(A, D)       C[I]   = sum_ab A_ab D[I,a,b]                       (metric contraction)
//   dInverse(Gi, dG)     R[I]   = - Gi dG[I] Gi                               (inverse-metric derivative)
//   jacHessian(...)      H[I,J] = 0.5 j (g^-1 : d_IJ g) + 0.5 j (d_I g^-1 : d_J g)
//                                 + 0.5 d_I j (g^-1 : d_J g)                  (Jacobian-Hessian)
//     with j = sqrt(det g); this is the second derivative of the area Jacobian.
//
// Tested against the original tensor expressions in test_AssemblyKernels.cpp.
#ifndef ASSEMBLY_KERNELS_H
#define ASSEMBLY_KERNELS_H
#include "hl_Tensor.h"

namespace AssemblyKernels
{
    using hiperlife::Tensor::tensor;

    // C[e,k] = sum_ab A_ab D[e,k,a,b]  == product(A, D, {{0,2},{1,3}})
    inline tensor<double,2> contract(const tensor<double,2>& A, const tensor<double,4>& D)
    {
        auto s = D.shape(); tensor<double,2> C(s[0], s[1]);
        const double* a = A.data(); const double* d = D.data(); double* c = C.data();
        for (int i = 0, n = s[0]*s[1]; i < n; i++, d += 4) c[i] = a[0]*d[0] + a[1]*d[1] + a[2]*d[2] + a[3]*d[3];
        return C;
    }
    // C[e,k,f,l] = sum_ab A_ab D[e,k,f,l,a,b]  == product(A, D, {{0,4},{1,5}})
    inline tensor<double,4> contract(const tensor<double,2>& A, const tensor<double,6>& D)
    {
        auto s = D.shape(); tensor<double,4> C(s[0], s[1], s[2], s[3]);
        const double* a = A.data(); const double* d = D.data(); double* c = C.data();
        for (int i = 0, n = s[0]*s[1]*s[2]*s[3]; i < n; i++, d += 4) c[i] = a[0]*d[0] + a[1]*d[1] + a[2]*d[2] + a[3]*d[3];
        return C;
    }
    // R[e,k] = -Gi dG[e,k] Gi   (Gi = G^{-1})  == -product(product(Gi, dG, {{1,2}}).transpose({1,2,0,3}), Gi, {{3,0}})
    //   evaluated in the same two stages as the tensor expression: T = Gi dG, then R = -(T Gi)
    inline tensor<double,4> dInverse(const tensor<double,2>& Gi, const tensor<double,4>& dG)
    {
        auto s = dG.shape(); tensor<double,4> R(s[0], s[1], 2, 2);
        const double* g = Gi.data(); const double* d = dG.data(); double* r = R.data();
        for (int i = 0, n = s[0]*s[1]; i < n; i++, d += 4, r += 4)
        {
            double t[4];
            for (int a = 0; a < 2; a++)
                for (int e = 0; e < 2; e++) t[2*a + e] = g[2*a]*d[e] + g[2*a + 1]*d[2 + e];
            for (int a = 0; a < 2; a++)
                for (int b = 0; b < 2; b++) r[2*a + b] = -(t[2*a]*g[b] + t[2*a + 1]*g[2 + b]);
        }
        return R;
    }
    // H[e,k,f,l] = 0.5 j (Gi : ddG[e,k,f,l]) + 0.5 j (dGiA[e,k] : dGB[f,l]) + 0.5 dJA[e,k] (Gi : dGB[f,l])
    //   == 0.5*j*product(Gi,ddG,{{0,4},{1,5}}) + 0.5*j*product(dGiA,dGB,{{2,2},{3,3}})
    //      + 0.5*outer(dJA, product(Gi,dGB,{{0,2},{1,3}}))
    //   A and B are the two differentiation variables (x or h); dGiA = dInverse(Gi, dG_A).
    inline tensor<double,4> jacHessian(double j, const tensor<double,2>& Gi, const tensor<double,6>& ddG,
                                       const tensor<double,4>& dGiA, const tensor<double,4>& dGB,
                                       const tensor<double,2>& dJA)
    {
        auto sa = dGiA.shape(), sb = dGB.shape();
        const int na = sa[0]*sa[1], nb = sb[0]*sb[1];
        tensor<double,4> H(sa[0], sa[1], sb[0], sb[1]);
        tensor<double,2> cB = contract(Gi, dGB);
        const double* g = Gi.data(); const double* dd = ddG.data(); const double* x = dGiA.data();
        const double* y = dGB.data(); const double* dj = dJA.data(); const double* cb = cB.data(); double* h = H.data();
        const double hj = 0.5*j;
        for (int i = 0; i < na; i++)
            for (int k = 0; k < nb; k++, dd += 4)
            {
                const double c = g[0]*dd[0] + g[1]*dd[1] + g[2]*dd[2] + g[3]*dd[3];
                const double p = x[4*i]*y[4*k] + x[4*i+1]*y[4*k+1] + x[4*i+2]*y[4*k+2] + x[4*i+3]*y[4*k+3];
                h[i*nb + k] = hj*c + hj*p + 0.5*(dj[i]*cb[k]);
            }
        return H;
    }
}
#endif
