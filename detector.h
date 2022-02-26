#ifndef DETECTOR_H
#define DETECTOR_H

#include <math.h>
#include <array>

struct Detector {
    const double tht_angl;  // Angle the detector is oriented relative to normal
    const double ph_angl;   // Phi angle of the detector
    const double vrt_dpth;
    const double aprtr;
    const double y_angl;
    const double nrmlzr;
    double d_vec[3];
    const double rot_mat[3][3];
    double intnsty_dtctd = 0.;
    double diffuse = 0.;

    void detect_photon(double wght,double trns_frc, double tau_to_dtctr, 
        double p_scat);

    void detect_photon(double wght, std::array<double, 3> phtn_vec);

    void detect_diffuse(double wght);

    Detector(const double tht_angl, const double ph_angl, 
        const double vrt_dpth = 0., const double aprtr = 5.0);

};

#endif
