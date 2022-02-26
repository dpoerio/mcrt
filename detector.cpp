#include "detector.h"
#include "constants.h"
#include <array>

Detector::Detector(const double tht_angl, const double ph_angl, 
    const double vrt_dpth, const double aprtr) : 
    tht_angl(tht_angl), 
    ph_angl(ph_angl), 
    vrt_dpth(vrt_dpth), 
    aprtr(aprtr), 
    y_angl(Constants::PI / 2 - tht_angl), 
    nrmlzr(pow(cos(tht_angl), 2.0) + pow(cos(Constants::PI / 2 - tht_angl), 
        2.0)),
    d_vec{cos(tht_angl) / nrmlzr, cos(y_angl) / nrmlzr, 0},
    rot_mat{{1, 0, 0},
        {0, cos(ph_angl), -sin(ph_angl)},
        {0, sin(ph_angl), cos(ph_angl)}} 
    {
        if (ph_angl > Constants::EPS){
            double d_vec_rot[3];
            d_vec_rot[0] = rot_mat[0][0] * d_vec[0] + rot_mat[0][1] * 
                d_vec[1] + rot_mat[0][2] * d_vec[2];
            d_vec_rot[1] = rot_mat[1][0] * d_vec[0] + rot_mat[1][1] * 
                d_vec[1] + rot_mat[1][2] * d_vec[2];
            d_vec_rot[2] = rot_mat[2][0] * d_vec[0] + rot_mat[2][1] * 
                d_vec[1] + rot_mat[2][2] * d_vec[2];
            d_vec[0] = d_vec_rot[0];
            d_vec[1] = d_vec_rot[1];
            d_vec[2] = d_vec_rot[2];
        }
    }

void Detector::detect_photon(double wght, double trns_frc, double tau_to_dtctr,
    double p_scat) {
    // Detect the photon using the local/directional estimate technique
    // TODO: Error check for invalid inputs and write tests
#pragma omp atomic
    intnsty_dtctd += (wght * trns_frc * p_scat / 
        (2 * Constants::PI) * exp(-tau_to_dtctr) / abs(d_vec[0]));
}

void Detector::detect_photon(double wght, std::array<double, 3> phtn_vec) {
    double ang = acos(-phtn_vec[0]*d_vec[0] + phtn_vec[1]*d_vec[1] + 
        phtn_vec[2]*d_vec[2])*180./Constants::PI;
    if (ang <= aprtr) {
#pragma omp atomic
        intnsty_dtctd += wght / (2.0 * Constants::PI * (1.0 - 
        cos(aprtr*(Constants::PI/180.0)))) / abs(phtn_vec[0]);
    }
}

void Detector::detect_diffuse(double wght){
#pragma omp atomic
    diffuse += wght;
}

