#ifndef PHOTON_H
#define PHOTON_H

#include <random>
#include "slab.h"
#include "detector.h"
#include "constants.h"
#include <array>

struct Photon {
    Slab& slab;
    Detector& detector;
    Detector& bottom_detector;
    std::mt19937& gen;
    double inc_angl;
    double ph_angl;
    double nrmlzr;
    double dstnc_thrgh_slb = Constants::INF;
    std::array<double, 3> phtn_vec;
    double rot_mat[3][3];    
    double phtn_vec_rot[3];
    bool rchd_scat_evnt = false;  // Variable for while loop to determine when 
        //the event is reached
    double vrt_dpth = 0;  // Vertical depth of photon in the slab
    double intnsty = 1;   // Intensity of the photon at the present time
    const double intnsty_thrshld =
        1e-4;            // Intensity at which Russian Roulette kicks in
    const int rlt_cnstnt = 6;  // The Roulette constant
    bool in_slb = true;  // Whether or not the photon is in the slab
    int cr_lyr = 0;      // The current layer in the slab
    double tau_tos = 0;  // The optical distance to the top of the slab
    double dt_sct_ang;   // The needed scattering angle to hit the detector
    double scat_evnt_dstnc;
    double scat_dstnc_to_trvl;  // This is the distance left to travel -
                                // initially it is the
    // same as the attenuation event distance but can decrease in intermediate
    // steps if the photon has to travel through an interface
    double opt_dst_tos;
    double opt_dst_bos;

    double phys_dstnc_thrgh_lyr;
    double scat_dstnc_thrgh_lyr;


    enum class Evnt{
        Initiation, /* For the purpose of getting the downward direct flux. 
            TODO: Need to explicitly set photons as collimated or diffuse since
            only the collimated beam contributes to downward direct flux */
        Refraction,
        Reflection,
        Absorption,
        Scattering,
        Annihilation
    };

    enum class Trvl_Dir{
        Up = 0,
        Down = 1,
        Horizontal = 2
    };

    Trvl_Dir trvl_dir = Photon::Trvl_Dir::Down; // Always start traveling down

    Evnt last_evnt = Photon::Evnt::Initiation;


    std::array<double, 3> nec_scat_dir;

    void updt_trvl_dir();
    
    void lcl_intnsty_estmt();

    double clc_dir_of_inc(std::array<double, 3> dir_vec, 
        std::array<double, 3> pnt_vec);

    double clc_dir_of_rfrc(double dir_of_inc, double n);

    std::array<double, 3> rfrctn(std::array<double, 3> dir_vec, 
        std::array<double, 3> pnt_vec, double n, double cos_th1, 
        double cos_th2);

    std::array<double, 3> rflctn(std::array<double, 3> dir_vec, 
        std::array<double, 3> pnt_vec, double cos_th1);

    double clct_rfrc_angs_tos(std::vector<std::array<double, 3>> 
        &rfrc_cos_angs);

    void clct_opt_dst_tos(std::vector<std::array<double, 3>> &rfrc_cos_angs);

    void clct_opt_dst_toinb();

    double clct_rfrc_angs_bos(std::vector<std::array<double, 3>> 
        &rfrc_cos_angs);

    void clct_opt_dst_bos(std::vector<std::array<double, 3>> &rfrc_cos_angs);

    void smpl_scat_evnt_dstnc();

    void trvl_to_scat_evnt();

    void absrb(double init_dpth);

    void scttr();

    void updt_dir_cos(double sct_tht, double sct_phi);

    void get_cur_lyr();

    void updt_scat_dst_to_trvl(double init_dpth);
    
    void clct_dstnc_thrgh_lyr();

    void clct_dstnc_thrgh_slb();

    void escp();
    
    void trvl_to_annhltn();

    void dtct();

    void rltt();

    void intrfc_intrctn(double n1, double n2);

    double frsnl(double cos_th1, double cos_th2, double n1, double n2);

    Photon(Slab& slb, std::vector<Detector> &dtctrs, std::mt19937& prng, 
        double inc_angl=Constants::PI/4., double ph_angl=0);

};

#endif
