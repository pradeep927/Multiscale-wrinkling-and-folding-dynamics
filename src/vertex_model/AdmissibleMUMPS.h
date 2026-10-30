#pragma once
#include "InvalidGeometry.h"
#include "hl_LinearSolver_Direct_MUMPS.h"

// Keep model-specific admissibility out of the generic linear solver.
class AdmissibleMUMPS : public hiperlife::MUMPSDirectLinearSolver {
public:
    bool solve() override {
        vertex_admissibility::invalid = false;
        _timer->start("filling");
        _hiperProbl->FillLinearSystem();
        _timer->end("filling");
        if (!vertex_admissibility::validEverywhere()) {
            _converged = false;
            _linearProblem->GetLHS()->PutScalar(0.0);
            if (myRank() == 0) std::cerr << "Newton rejected: inadmissible trial state" << std::endl;
            return false;
        }
        // Already assembled above; preserve the retained single-assembly optimization.
        const bool result = MUMPSDirectLinearSolver::solve();
        return result;
    }
};
