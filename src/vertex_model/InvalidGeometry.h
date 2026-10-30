#pragma once
#include <cmath>
#include <mpi.h>
// No collectives inside element callbacks: some ranks have different workloads.
#include <atomic>
namespace vertex_admissibility {
inline std::atomic<bool> invalid{false};
inline bool positive(double value) { return std::isfinite(value) && value > 0.0; }
inline bool validEverywhere() {
    int local = invalid.load() ? 1 : 0, global = 0;
    MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    return global == 0;
}
}
