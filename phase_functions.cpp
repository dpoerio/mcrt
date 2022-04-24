#include "phase_functions.h"
#include "constants.h"
#include <cmath>
#include <stdexcept>
#include <string>


double Phase_Functions::hg(double g, double mu) {
    if (std::abs(g) > 1. or std::abs(mu) > 1.) {
        throw std::domain_error(std::string("abs(g) and abs(mu) must be below \
            1 with the Henyey-Greenstein function"));
    }
    return (0.5)*(1-g*g)/pow((1+g*g-2*g*mu), 1.5);
}

double Phase_Functions::ihg(double g, double rn) {
    if (std::abs(g) > 1. or std::abs(rn) > 1.) {
        throw std::domain_error(std::string("std::abs(g) and std::abs(rn) must be below \
            1 with the Henyey-Greenstein function"));
    }
    if (std::abs(g) < Constants::EPS){
        return 2.*rn-1.;
    } else {
        return (1./(2.*g))*(1+g*g-pow((1-g*g)/(1+g*(2.*rn-1)),2));
    }
}

double Phase_Functions::rlgh(double mu) {
    if (std::abs(mu) > 1.) {
        throw std::domain_error(std::string("abs(mu) must be below \
            1 with the Rayleigh function"));
    }
    return 0.75*(1.0 + pow(mu, 2.0));
}

double Phase_Functions::irlgh(double rn) {
    if (std::abs(rn) > 1.) {
        throw std::domain_error(std::string("abs(rn) must be below \
            1 with the Rayleigh function"));
    }
    double q = 4.0*rn-2.0;
    double u = pow(-q + pow(1 + pow(q, 2.0), 0.5), 1.0/3.0);
    return (u - (1/u));
}
