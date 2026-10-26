// Independent test of AssemblyKernels.h against the original hiperlife tensor expressions (random data).
// Compile separately with the application's HiPerLife compiler flags and -I.. from this directory.
// Link LAPACK/BLAS as required by the HiPerLife tensor headers. Exit code 0 = agreement to round-off.
#include "AssemblyKernels.h"
#include <random>
#include <cstdio>
#include <cmath>
using namespace hiperlife::Tensor;
int main()
{
    std::mt19937 gen(1); std::uniform_real_distribution<double> U(-1, 1);
    auto fill = [&](double* p, int n) { for (int i = 0; i < n; i++) p[i] = U(gen); };
    auto err  = [](const double* a, const double* b, int n) { double m = 0; for (int i = 0; i < n; i++) m = std::max(m, std::fabs(a[i]-b[i])); return m; };
    const int e = 3, k = 3;                                   // linear triangle, 3 displacement dofs
    tensor<double,2> G(2,2); fill(G.data(), 4); G(0,1) = G(1,0); G(0,0) += 3; G(1,1) += 3;   // SPD metric
    tensor<double,2> Gi = G.inv();
    double j = std::sqrt(G.det());
    tensor<double,4> dGA(e,k,2,2), dGB(e,k,2,2); fill(dGA.data(), e*k*4); fill(dGB.data(), e*k*4);
    tensor<double,6> ddG(e,k,e,k,2,2); fill(ddG.data(), e*k*e*k*4);
    tensor<double,2> dJA(e,k); fill(dJA.data(), e*k);

    // 1) metric contraction
    tensor<double,2> c1 = product(Gi, dGA, {{0,2},{1,3}}), c1n = AssemblyKernels::contract(Gi, dGA);
    tensor<double,4> c2 = product(Gi, ddG, {{0,4},{1,5}}), c2n = AssemblyKernels::contract(Gi, ddG);
    // 2) inverse-metric derivative
    tensor<double,4> s1 = product(product(Gi, dGA, {{1,2}}).transpose({1,2,0,3}), Gi, {{3,0}}); s1 = -1.0*s1;
    tensor<double,4> s1n = AssemblyKernels::dInverse(Gi, dGA);
    // 3) Jacobian-Hessian (original expression of LS)
    tensor<double,2> cB = product(Gi, dGB, {{0,2},{1,3}});
    tensor<double,4> h1 = 0.5*j*product(Gi, ddG, {{0,4},{1,5}}) + 0.5*j*product(s1, dGB, {{2,2},{3,3}}) + 0.5*outer(dJA, cB);
    tensor<double,4> h1n = AssemblyKernels::jacHessian(j, Gi, ddG, s1n, dGB, dJA);

    double e1 = err(c1.data(), c1n.data(), e*k), e2 = err(c2.data(), c2n.data(), e*k*e*k);
    double e3 = err(s1.data(), s1n.data(), e*k*4), e4 = err(h1.data(), h1n.data(), e*k*e*k);
    std::printf("contract(2,4) %.2e  contract(2,6) %.2e  dInverse %.2e  jacHessian %.2e\n", e1, e2, e3, e4);
    bool ok = e1 < 1e-14 && e2 < 1e-14 && e3 < 1e-14 && e4 < 1e-14;   // round-off level
    std::printf("%s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
