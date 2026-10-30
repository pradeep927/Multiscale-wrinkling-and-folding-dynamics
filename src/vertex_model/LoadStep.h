#pragma once
#include <algorithm>

namespace vertex_loading {
struct Step {
    double dt;
    double end;
    double rate;
};

// Times are absolute. Each attempted interval stays within one loading phase.
// The onset cap is applied once, so a rejected attempt keeps its reduced dt.
inline Step prepare(double time, double proposedDt, double inflationStart,
                    double inflationEnd, double deflationStart,
                    double deflationEnd, double inflationRate,
                    double deflationRate, bool& resetDeflationDt)
{
    if (resetDeflationDt && time >= deflationStart) {
        proposedDt = std::min(proposedDt, 0.05);
        resetDeflationDt = false;
    }
    double end = time + proposedDt;
    for (double event : {inflationStart, inflationEnd, deflationStart, deflationEnd})
        if (time < event && event < end)
            end = event;
    // Select by the left endpoint: an interval ending at an event belongs
    // to the preceding phase; an interval starting there belongs to the next.
    const double rate = time >= inflationStart && time < inflationEnd
        ? inflationRate
        : (time >= deflationStart && time < deflationEnd ? deflationRate : 0.0);
    // Preserve the proposed dt exactly unless a boundary actually shortened it.
    const double dt = end == time + proposedDt ? proposedDt : end - time;
    return {dt, end, rate};
}
} // namespace vertex_loading
