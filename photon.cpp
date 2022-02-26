#include "photon.h"
#include "slab.h"
#include "detector.h"
#include <random>
#include "constants.h"
#include <vector>
#include <array>
#include <cassert>
#include "math.h"
#include <algorithm>

static bool lcl_estmt = false;

Photon::Photon(Slab& slb, std::vector<Detector> &dtctrs, std::mt19937& prng, 
    double inc_angl, double ph_angl) : 
    slab(slb), 
    detector(dtctrs[0]), 
    bottom_detector(dtctrs[1]), 
    gen(prng),
    inc_angl(inc_angl), 
    ph_angl(ph_angl), 
    nrmlzr(pow(cos(inc_angl), 2.0) + pow(cos(Constants::PI / 2 - inc_angl), 
        2.0)),
    phtn_vec{cos(inc_angl) / nrmlzr, cos(Constants::PI / 2 - inc_angl) / 
        nrmlzr, 0},
    rot_mat{{1, 0, 0},
        {0, cos(ph_angl), -sin(ph_angl)},
        {0, sin(ph_angl), cos(ph_angl)}}
    {
        if (ph_angl > Constants::EPS){
            double phtn_vec_rot[3];
            phtn_vec_rot[0] = rot_mat[0][0] * phtn_vec[0] + rot_mat[0][1] * 
                phtn_vec[1] + rot_mat[0][2] * phtn_vec[2];
            phtn_vec_rot[1] = rot_mat[1][0] * phtn_vec[0] + rot_mat[1][1] * 
                phtn_vec[1] + rot_mat[1][2] * phtn_vec[2];
            phtn_vec_rot[2] = rot_mat[2][0] * phtn_vec[0] + rot_mat[2][1] * 
                phtn_vec[1] + rot_mat[2][2] * phtn_vec[2];
            phtn_vec[0] = phtn_vec_rot[0];
            phtn_vec[1] = phtn_vec_rot[1];
            phtn_vec[2] = phtn_vec_rot[2];
        }
    }

double Photon::clct_rfrc_angs_tos(std::vector<std::array<double, 3>> 
    &rfrc_cos_angs){
    // Trace from the detector to the photon position
    // TODO: Account for arbitrary detector position
    std::array<double, 3> pnt_vec = {-1., 0., 0.};
    rfrc_cos_angs.push_back({detector.d_vec[0], detector.d_vec[1], 
        detector.d_vec[2]});
    double trnsff = 1.;
    for (int i = 0; i < cr_lyr; i++) {
        double n1 = slab.layers[i].n;
        double n2 = slab.layers[i+1].n;
        double n = n1 / n2;
        if(n1 > n2){
            double theta_crit = asin(1/n);
            if (acos(abs(rfrc_cos_angs[i][0])) > theta_crit){
                // The photon cannot hit the detector
                trnsff = 0.;
                return trnsff;
            }
        }
        double cos_th1 = clc_dir_of_inc(rfrc_cos_angs[i], pnt_vec);
        double cos_th2 = clc_dir_of_rfrc(cos_th1, n);
        trnsff *= (1. - frsnl(cos_th1, cos_th2, n1, n2))*(pow(n,2));
        std::array<double, 3> v_rfrc_i = rfrctn(rfrc_cos_angs[i], pnt_vec, n, 
            cos_th1, cos_th2);
        rfrc_cos_angs.push_back(v_rfrc_i);
    }
    return trnsff;
}

void Photon::clct_opt_dst_toinb(){
    // Determine if a refracting interface is present between photon and
    // detector
    std::array<double, 3> dir_vec = {detector.d_vec[0], detector.d_vec[1],
        detector.d_vec[2]};
    std::array<double, 3> pnt_vec = {-1., 0., 0.};
    for (int i = cr_lyr; i > 0; i--) {
        double n1 = slab.layers[i].n;
        double n2 = slab.layers[i-1].n;
        if (std::abs(n1 - n2) > Constants::EPS) {
            return;
        }
    }
    // No refracting interface between photon and detector. Check if there 
    // is an interface away from the detector to reflect from. Tally the optical
    // distance along the way in case we are reflecting back for the local
    // estimate
    opt_dst_tos = 0.0;
    double cr_dpth = vrt_dpth;
    double dst_to_tol;
    int intrfc_lyr;
    bool intrfc_exists = false;
    double rflctn_coef;
    for (int i = cr_lyr; i <= slab.n_layers-1; i++) {
        double n1 = slab.layers[i].n;
        if (i + 1 >= slab.n_layers) {
            break;
        }
        double n2 = slab.layers[i+1].n;
        dst_to_tol = abs(cr_dpth-slab.layer_interfaces[i][1]) / 
            detector.d_vec[0];
        cr_dpth = slab.layer_interfaces[i][1];
        opt_dst_tos = opt_dst_tos +
                      (slab.layers[i].s + slab.layers[i].k) * dst_to_tol;
        if (std::abs(n1 - n2) > Constants::EPS) {
            double n = n1/n2;
            // There is a refracting interface away from the detector. Determine
            // the direction under which we would totally internally reflect
            // from this interface to the detector
            double theta_crit = asin(1/n);
            if (acos(abs(dir_vec[0])) > theta_crit){
                rflctn_coef = 1.;
            } else {
                double cos_th1 = clc_dir_of_inc(dir_vec, pnt_vec);
                double cos_th2 = clc_dir_of_rfrc(cos_th1, n);
                rflctn_coef = frsnl(cos_th1, cos_th2, n1, n2);
            }
            intrfc_lyr = i;
            intrfc_exists = true;
            break;
        }
    }
    if (!intrfc_exists) {
        return;
    }
    // Now calculate the optical distance back to the detector
    for (int i = intrfc_lyr; i >= 0; i--) {
        // We already know the RI doesn't change
        dst_to_tol = abs(cr_dpth-slab.layer_interfaces[i][0]) / 
            detector.d_vec[0];
        cr_dpth = slab.layer_interfaces[i][0];
        opt_dst_tos = opt_dst_tos +
                      (slab.layers[i].s + slab.layers[i].k) * dst_to_tol;
    }
    double scat_cos = phtn_vec[0]*dir_vec[0] + 
        phtn_vec[1]*dir_vec[1] +
        phtn_vec[2]*dir_vec[2];
    double p_scat = slab.layers[cr_lyr].p[0](scat_cos);
    detector.detect_photon(intnsty, rflctn_coef, opt_dst_tos, p_scat);
    return;
}

void Photon::clct_opt_dst_tos(std::vector<std::array<double, 3>> 
    &rfrc_cos_angs) {
    /* Calculate the optical distance from the current position to the top of
    the slab.This is solely to make the local radiance estimate, thus we have to 
    consider the directional angle of the detector.*/
    // TODO: Account for arbitrary detector position
    opt_dst_tos = 0.0;
    double cr_dpth = vrt_dpth;
    double dst_to_tol;
    nec_scat_dir = rfrc_cos_angs[cr_lyr];
    nec_scat_dir[0] *= -1.;
    for (int i = cr_lyr; i >= 0; i--) {
        // Is this the top layer?
        if (i == 0) {
            // We are calculating from the top layer
            dst_to_tol = abs(cr_dpth) / rfrc_cos_angs[i][0];
            cr_dpth = 0.0;
        } else {
            // Interface[0] is smaller than interface[1] ([1] is deeper)
            dst_to_tol = abs(cr_dpth - slab.layer_interfaces[i][0]) /
                         rfrc_cos_angs[i][0];
            cr_dpth = slab.layer_interfaces[i][0];
        }
        opt_dst_tos = opt_dst_tos +
                      (slab.layers[i].s + slab.layers[i].k) * dst_to_tol;
    }
}

double Photon::clct_rfrc_angs_bos(std::vector<std::array<double, 3>> 
    &rfrc_cos_angs){
    // Trace from the bottom detector to the photon position
    // TODO: Account for arbitrary detector position
    std::array<double, 3> pnt_vec = {1., 0., 0.};
    // TODO: Need a bottom detector vector
    rfrc_cos_angs.push_back({bottom_detector.d_vec[0], bottom_detector.d_vec[1], 
        bottom_detector.d_vec[2]});
    double trnsff = 1.;
    int n_layers = slab.n_layers;
    int cntr = 0;
    for (int i = n_layers-1; i > cr_lyr; i--) {
        double n1 = slab.layers[i].n;
        double n2 = slab.layers[i-1].n;
        double n = n1 / n2;
        if(n1 > n2){
            double theta_crit = asin(1/n);
            if (acos(abs(rfrc_cos_angs[cntr][0])) > theta_crit){
                // The photon cannot hit the detector
                trnsff = 0.;
                return trnsff;
            }
        }
        double cos_th1 = clc_dir_of_inc(rfrc_cos_angs[cntr], pnt_vec);
        double cos_th2 = clc_dir_of_rfrc(cos_th1, n);
        trnsff *= (1. - frsnl(cos_th1, cos_th2, n1, n2))*(pow(n,2));
        std::array<double, 3> v_rfrc_i = rfrctn(rfrc_cos_angs[cntr], pnt_vec, n, 
            cos_th1, cos_th2);
        rfrc_cos_angs.push_back(v_rfrc_i);
        cntr++;
    }
    return trnsff;
}

void Photon::clct_opt_dst_bos(std::vector<std::array<double, 3>> 
    &rfrc_cos_angs) {
    /* Calculate the optical distance from the current position to the bottom of
    the slab.This is solely to make the local radiance estimate, thus we have to 
    consider the directional angle of the detector.*/
    // TODO: Account for arbitrary detector position
    opt_dst_bos = 0.0;
    double cr_dpth = vrt_dpth;
    double dst_to_tol;
    nec_scat_dir = rfrc_cos_angs[0]; // Measured from current layer (ind 0) to bot
    nec_scat_dir[0] *= -1.;
    int n_layers = slab.n_layers - 1;
    // Map between loop index and actual layer
    std::vector<int> ind_map;
    for (int i = cr_lyr; i <= n_layers; i++) {
        ind_map.push_back(i);
    }

    for (int i = 0; i <= (n_layers - cr_lyr); i++) {
        // Is this the bottom layer?
        if (i == (n_layers-cr_lyr)) {
            // We are calculating from the bottom layer
            dst_to_tol = abs(cr_dpth - 
                slab.layer_interfaces[ind_map[i]][1]) / 
                abs(rfrc_cos_angs[i][0]);
            cr_dpth = slab.layer_interfaces[ind_map[i]][1];
        } else {
            // Interface[0] is smaller than interface[1] ([1] is deeper)
            dst_to_tol = abs(cr_dpth - 
                slab.layer_interfaces[ind_map[i]][1]) /
                         abs(rfrc_cos_angs[i][0]);
            cr_dpth = slab.layer_interfaces[ind_map[i]][1];
        }
        opt_dst_bos = opt_dst_bos +
                      (slab.layers[ind_map[i]].s + 
                       slab.layers[ind_map[i]].k) * 
                      dst_to_tol;
    }
}

void Photon::smpl_scat_evnt_dstnc() {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    double rn = dist(gen);
    scat_evnt_dstnc = -log(1 - (1. - exp(-dstnc_thrgh_slb))*rn);
    // Need to set this distance as the "distance left to travel" - this is
    // the strategy for handling hitting interfaces in the simulation
    scat_dstnc_to_trvl = scat_evnt_dstnc;
    return;
}

void Photon::trvl_to_scat_evnt() {
    double init_dpth = vrt_dpth;
    get_cur_lyr();
    //clct_dstnc_thrgh_slb();
    smpl_scat_evnt_dstnc();
    while (!rchd_scat_evnt) {
        clct_dstnc_thrgh_lyr();
        if (scat_dstnc_to_trvl > scat_dstnc_thrgh_lyr) {
            // In this case, we travel to the interface
            if (trvl_dir == Photon::Trvl_Dir::Up) {
                // Travel to top of current layer
                vrt_dpth = slab.layer_interfaces[cr_lyr][0];
                // Can also static_cast<int>(trvl_dir) to convert scoped enum
            } else if (trvl_dir == Photon::Trvl_Dir::Down) {
                // Travel to bottom of current layer
                vrt_dpth = slab.layer_interfaces[cr_lyr][1];
            } else {
                throw std::logic_error(std::string("Horizontal travel not"
                    " yet properly accounted for"));
            }
            updt_scat_dst_to_trvl(init_dpth);
            absrb(init_dpth);
            dtct();
            rltt();
            if (!in_slb) {
                return;
            }
            if (trvl_dir == Photon::Trvl_Dir::Up) {
                intrfc_intrctn(slab.layers[cr_lyr].n, 
                        slab.layers[cr_lyr - 1].n);
            } else if (trvl_dir == Photon::Trvl_Dir::Down) {
                intrfc_intrctn(slab.layers[cr_lyr].n, 
                        slab.layers[cr_lyr + 1].n);
            } else {
                throw std::logic_error(std::string("Horizontal travel not"
                    " yet properly accounted for"));
            }
            updt_trvl_dir();
            init_dpth = vrt_dpth;
        } else {
        // This is now the layer we reach the attenuation event in
            vrt_dpth +=
                phtn_vec[0] * scat_dstnc_to_trvl / slab.layers[cr_lyr].s;
            absrb(init_dpth);
            /* intnsty *= (1. - exp(-dstnc_thrgh_slb)); //  ref Kattawar, remove 
             //bias for thin layers */
            dtct();
            rltt();
            if (!in_slb) {
                return;
            }
            updt_scat_dst_to_trvl(init_dpth);
            rchd_scat_evnt = true;
            if (lcl_estmt) {
                lcl_intnsty_estmt();
            }
            scttr();
            updt_trvl_dir();
        }
    }
    rchd_scat_evnt = false;
    return;
}

void Photon::updt_trvl_dir(){
    if (phtn_vec[0] > 0) {
        trvl_dir = Photon::Trvl_Dir::Down;
    } else if (phtn_vec[0] < 0) {
        trvl_dir = Photon::Trvl_Dir::Up;
    } else {
        trvl_dir = Photon::Trvl_Dir::Horizontal;
        throw std::logic_error(std::string("Horizontal travel not"
            " yet properly accounted for"));
    }
}

void Photon::dtct(){
    if (vrt_dpth < Constants::EPS && in_slb) {
        assert(trvl_dir == Photon::Trvl_Dir::Up);
        if (!lcl_estmt) {
            detector.detect_photon(intnsty, phtn_vec);
        }
        slab.accum_flux(intnsty, cr_lyr, 1);
        escp();
    } else if (abs(vrt_dpth - slab.layer_interfaces[slab.n_layers-1][1]) < 
        Constants::EPS && in_slb) {
        assert(trvl_dir == Photon::Trvl_Dir::Down);
        if (!lcl_estmt) {
            bottom_detector.detect_photon(intnsty, phtn_vec);
        }
        slab.accum_flux(intnsty, cr_lyr, 0);
        escp();
    }
}

void Photon::rltt(){
    if (intnsty < intnsty_thrshld && in_slb) {
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        double rn = dist(gen);
        if (rn <= 1./rlt_cnstnt){
            escp();
        } else {
            intnsty *= rlt_cnstnt;
        }
    }
}

void Photon::absrb(double init_dpth) {
    intnsty *= exp(-slab.layers[cr_lyr].k * abs(init_dpth - vrt_dpth) /
                   abs(phtn_vec[0]));
    last_evnt = Photon::Evnt::Absorption;
    return;
}

void Photon::lcl_intnsty_estmt(){
    /* Note: There are TWO ways a photon can scatter and hit a detector. In all
       cases, the photon can scatter directly into the detector (accounting for
       refraction if there are interfaces between the photon and detector).
       However, IFF the photon is on the same side of a refracting interface as
       a detector, the photon can also scatter in the direction OF THE
       INTERFACE, AWAY FROM THE DETECTOR, such that it is totally internally
       reflected into the direction of the detector. In this case we must
       multiply by the reflection coefficient instead of the transmission 
       coefficient */
    std::vector<std::array<double, 3>> rfrc_cos_angs;
    std::vector<std::array<double, 3>> rfrc_cos_angs_bot;
    double trnsff = clct_rfrc_angs_tos(rfrc_cos_angs);
    if (trnsff > Constants::EPS) {
        clct_opt_dst_toinb();
        clct_opt_dst_tos(rfrc_cos_angs);
        double scat_cos = phtn_vec[0]*nec_scat_dir[0] + 
            phtn_vec[1]*nec_scat_dir[1] +
            phtn_vec[2]*nec_scat_dir[2];
        double p_scat = slab.layers[cr_lyr].p[0](scat_cos);
        detector.detect_photon(intnsty, trnsff, opt_dst_tos,
            p_scat);
    }
    double trnsff_bot = clct_rfrc_angs_bos(rfrc_cos_angs_bot);
    if (trnsff_bot > Constants::EPS) {
        clct_opt_dst_bos(rfrc_cos_angs_bot);
        double scat_cos = phtn_vec[0]*nec_scat_dir[0] + 
            phtn_vec[1]*nec_scat_dir[1] +
            phtn_vec[2]*nec_scat_dir[2];
        double p_scat = slab.layers[cr_lyr].p[0](scat_cos);
        bottom_detector.detect_photon(intnsty, trnsff_bot, opt_dst_bos,
            p_scat);
    }
    return;
}

void Photon::scttr() {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    double rn1 = dist(gen);
    double rn2 = dist(gen);
    double sct_cos_tht = slab.layers[cr_lyr].p[1](rn1);
    double sct_tht = acos(sct_cos_tht);
    double sct_phi = 2 * rn2 * Constants::PI;
    // TODO: Always assuming azimuthal symmetry at the moment
    updt_dir_cos(sct_tht, sct_phi);
    return;
}

void Photon::updt_dir_cos(double sct_tht, double sct_phi) {
    // Rodrigues rotation formulas
    // https://en.wikipedia.org/wiki/Monte_Carlo_method_for_photon_transport
    // and
    // https://www.scratchapixel.com/lessons/mathematics-physics-for-computer-graphics/monte-carlo-methods-in-practice/monte-carlo-simulation
    // Update the direction cosines based on the sampled scattering angles
    double sct_sin_tht = sin(sct_tht);
    double sct_cos_phi = cos(sct_phi);
    double sct_sin_phi = sin(sct_phi);
    double sct_cos_tht = cos(sct_tht);
    double new_cos_tht;
    double new_cos_y;
    double new_cos_phi;
    // Implement the special case formulas
    if (abs(phtn_vec[0] - 1) < Constants::EPS) {
        phtn_vec[0] = sct_cos_tht;
        phtn_vec[1] = sct_sin_tht * sct_cos_phi;
        phtn_vec[2] = sct_sin_tht * sct_sin_phi;
    } else if (abs(1 + phtn_vec[0]) < Constants::EPS) {
        phtn_vec[0] = -sct_cos_tht;
        phtn_vec[1] = -sct_sin_tht * sct_sin_phi;
        phtn_vec[2] = sct_sin_tht * sct_cos_phi;
    } else {
        // order is x, y, phi
        new_cos_tht = -pow(1.0 - pow(phtn_vec[0], 2.0), 0.5) * sct_sin_tht *
                          sct_cos_phi +
                      phtn_vec[0] * sct_cos_tht;
        new_cos_y = sct_sin_tht *
                        (phtn_vec[1] * phtn_vec[0] * sct_cos_phi +
                         phtn_vec[2] * sct_sin_phi) /
                        pow(1.0 - pow(phtn_vec[0], 2.0), 0.5) +
                    phtn_vec[1] * sct_cos_tht;
        new_cos_phi = sct_sin_tht *
                          (phtn_vec[2] * phtn_vec[0] * sct_cos_phi -
                           phtn_vec[1] * sct_sin_phi) /
                          pow(1.0 - pow(phtn_vec[0], 2.0), 0.5) +
                      phtn_vec[2] * sct_cos_tht;
        phtn_vec[0] = new_cos_tht;
        phtn_vec[1] = new_cos_y;
        phtn_vec[2] = new_cos_phi;
    }
    last_evnt = Photon::Evnt::Scattering;
    return;
}

void Photon::get_cur_lyr() {
    for (int i = 0; i < slab.n_layers; i++) {
        if (vrt_dpth >= slab.layer_interfaces[i][0] &&
            vrt_dpth < slab.layer_interfaces[i][1]) {
            cr_lyr = i;
            break;
        }
    }
    return;
}

void Photon::updt_scat_dst_to_trvl(double init_dpth) {
    double dst_trvld = slab.layers[cr_lyr].s / abs(phtn_vec[0]) *
                       abs(init_dpth - vrt_dpth);
    scat_dstnc_to_trvl -= dst_trvld;
    return;
}

void Photon::clct_dstnc_thrgh_slb() {
    dstnc_thrgh_slb = 0;
    if (trvl_dir == Photon::Trvl_Dir::Down) { 
        double cr_dpth = vrt_dpth;
        std::vector<std::array<double, 3>> rfrc_cos_angs_bot;
        double trnsff_bot = clct_rfrc_angs_bos(rfrc_cos_angs_bot);
        if (trnsff_bot > Constants::EPS) {
            std::reverse(rfrc_cos_angs_bot.begin(), rfrc_cos_angs_bot.end());
            for (int i = cr_lyr; i <= slab.n_layers-1; i++) {
                double dstnc_thrgh_lyr = abs(cr_dpth - slab.layer_interfaces[i][1]) /
                    abs(rfrc_cos_angs_bot[i][0]);
                dstnc_thrgh_slb += dstnc_thrgh_lyr*(slab.layers[i].k + 
                    slab.layers[i].s);
                cr_dpth = slab.layer_interfaces[i][1];
            }
        }
    } else if (trvl_dir == Photon::Trvl_Dir::Up) {
        double cr_dpth = 0;
        std::vector<std::array<double, 3>> rfrc_cos_angs;
        double trnsff = clct_rfrc_angs_bos(rfrc_cos_angs);
        if (trnsff > Constants::EPS) {
            for (int i = 0; i <= cr_lyr; i++) {
                double dstnc_thrgh_lyr = abs(cr_dpth - slab.layer_interfaces[i][0]) /
                    rfrc_cos_angs[i][0];
                dstnc_thrgh_slb += dstnc_thrgh_lyr*(slab.layers[i].k + 
                    slab.layers[i].s);
                cr_dpth = slab.layer_interfaces[i][0];
            }
        }
    }
}

void Photon::clct_dstnc_thrgh_lyr() {
    // Calculate the physical distance and scattering optical path through the
    // layer, based on the current direction of travel
    if (trvl_dir == Photon::Trvl_Dir::Down) {
        // Traveling down - need to get distance to cur_lyr + 1
        phys_dstnc_thrgh_lyr =
            (slab.layer_interfaces[cr_lyr][1] - vrt_dpth) /
            abs(phtn_vec[0]);
    } else if (trvl_dir == Photon::Trvl_Dir::Up) {
        // Traveling up - need to get distance to cur_lyr - 0
        phys_dstnc_thrgh_lyr =
            (vrt_dpth - slab.layer_interfaces[cr_lyr][0]) /
            abs(phtn_vec[0]);
    } else {
        // Traveling horizontally - we cannot reach the next VERTICAL layer
        // however, for a horizontally inhomogenous layer we will need to
        // break this up
        // In practice this should basically never happen.
        // TODO: Implement horizontal gridding
        phys_dstnc_thrgh_lyr = Constants::INF;
        throw(std::logic_error(std::string("Horizontal travel direction not "
            "yet accounted for!")));
    }
    scat_dstnc_thrgh_lyr = phys_dstnc_thrgh_lyr * slab.layers[cr_lyr].s;
    return;
}

double Photon::frsnl(double cos_th1, double cos_th2, double n1, double n2) {
    double rs = pow(abs((n1*cos_th1 - n2*cos_th2)/((n1*cos_th1 + n2*cos_th2))),
        2);
    double rp = pow(abs((n1*cos_th2 - n2*cos_th1)/((n1*cos_th2 + n2*cos_th1))),
        2);
    return (0.5 * (rs + rp));
}

std::array<double, 3> Photon::rflctn(std::array<double, 3> dir_vec, 
    std::array<double, 3> pnt_vec, double cos_th1){
    // Total internal reflection
    std::array<double ,3> v_rflct;
    double tmp[3] = {2*cos_th1*pnt_vec[0], 2*cos_th1*pnt_vec[1], 
        2*cos_th1*pnt_vec[2]};
    v_rflct[0] = dir_vec[0] + tmp[0];
    v_rflct[1] = dir_vec[1] + tmp[1];
    v_rflct[2] = dir_vec[2] + tmp[2];
    return v_rflct;
}

std::array<double, 3> Photon::rfrctn(std::array<double, 3> dir_vec, 
    std::array<double, 3> pnt_vec, double n, double cos_th1, double cos_th2){
    // Refraction
    double tmp[3] = {n*dir_vec[0], n*dir_vec[1],n*dir_vec[2]}; 
    double tmp2 = n*cos_th1-cos_th2;
    double tmp3[3] = {pnt_vec[0]*tmp2, pnt_vec[1]*tmp2, pnt_vec[2]*tmp2};
    std::array<double, 3> v_rfrct = {tmp3[0]+tmp[0], tmp3[1]+tmp[1], 
        tmp3[2]+tmp[2]};
    return v_rfrct;
}

double Photon::clc_dir_of_inc(std::array<double, 3> dir_vec,
    std::array<double, 3> pnt_vec){
    double cos_th1 = -(pnt_vec[0]*dir_vec[0] + pnt_vec[1]*dir_vec[1] + 
        pnt_vec[2]*dir_vec[2]);
    return cos_th1;
}

double Photon::clc_dir_of_rfrc(double dir_of_inc, double n){
    double cos_th2 = pow(1. - pow(n, 2)*(1-pow(dir_of_inc,2)),0.5);
    return cos_th2;
}

void Photon::intrfc_intrctn(double n1, double n2) {
    double n = n1/n2;
    int sgn;
    if (phtn_vec[0] > 0){
        sgn = 1;
    } else {
        sgn = -1;
    }
    std::array<double, 3> pnt_vec = {-1.*sgn, 0., 0.};
    double cos_th1 = clc_dir_of_inc(phtn_vec, pnt_vec);
    if(n1 > n2){
        double theta_crit = asin(1/n);
        if (acos(abs(phtn_vec[0])) > theta_crit){
            std::array<double, 3> v_rflct = rflctn(phtn_vec, pnt_vec, cos_th1);
            phtn_vec[0] = v_rflct[0];
            phtn_vec[1] = v_rflct[1];
            phtn_vec[2] = v_rflct[2];
            last_evnt = Photon::Evnt::Reflection;
            return;
        }
    }
    double cos_th2 = clc_dir_of_rfrc(cos_th1, n);
    std::array<double, 3> v_rfrct = rfrctn(phtn_vec, pnt_vec, n, cos_th1, 
        cos_th2);
    double reff = frsnl(phtn_vec[0], v_rfrct[0], n1, n2);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    double rn = dist(gen);
    if (rn < reff){
        // Total internal reflection
        std::array<double, 3> v_rflct = rflctn(phtn_vec, pnt_vec, cos_th1);
        phtn_vec[0] = v_rflct[0];
        phtn_vec[1] = v_rflct[1];
        phtn_vec[2] = v_rflct[2];
        last_evnt = Photon::Evnt::Reflection;
    } else {
        // Refraction
        phtn_vec[0] = v_rfrct[0];
        phtn_vec[1] = v_rfrct[1];
        phtn_vec[2] = v_rfrct[2];
        last_evnt = Photon::Evnt::Refraction;
        if (trvl_dir == Photon::Trvl_Dir::Up) {
            slab.accum_flux(intnsty, cr_lyr, 1);
            cr_lyr--;
        } else if (trvl_dir == Photon::Trvl_Dir::Down) {
            slab.accum_flux(intnsty, cr_lyr, 0);
            cr_lyr++;
        } else {
            throw std::logic_error(std::string("Horizontal travel during "
                "refraction!"));
        }
    }
    return;
}

void Photon::escp() {
    // TODO: Account for escaping out the bottom of the slab for non infinite
    // slabs
    in_slb = false;
    last_evnt = Photon::Evnt::Annihilation;
    return;
}

void Photon::trvl_to_annhltn() {
    if (in_slb) {
        while (in_slb) {
            trvl_to_scat_evnt();
        }
    } else {
        throw std::logic_error(std::string("Call to trvl_to_annhltn on already"
            " annihilated photon!"));
    }
    return;
}

