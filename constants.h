#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <limits>
#include <math.h>

namespace Constants {
    const double PI =  2.*acos(0);
    const double INF = std::numeric_limits<double>::max();
    const double EPS = (pow(std::numeric_limits<double>::epsilon(), 0.5));
    const double TST_THRSH = 0.05; // Percent error for test comparison
}

#endif
