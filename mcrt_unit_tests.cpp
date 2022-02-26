#include "detector.h"
#include "layer.h"
#include "photon.h"
#include "slab.h"
#include "constants.h"
#include <random>
#include <string>
#include <functional>
#include <stdexcept>
#include <cmath>
#include "test_helpers.h"

int test_layer_phase_calc_lut(){
    // To be implemented, not tested currently
    Layer::phass calc_phase = Layer::Phase(std::function<double(double)>(
        [](double mu){return 0.5;}));
    return EXIT_SUCCESS;
}

int test_layer_explicit_phase(){
    std::function<double(double)> pdf = std::function<double(double)>(
        [](double mu){return 0.5;});
    std::function<double(double)> icdf = std::function<double(double)>(
        [](double rn){return 2.*rn-1.;});
    Layer::phass phas = Layer::Phase(pdf, icdf);
    return EXIT_SUCCESS;
}

int test_make_slab(){
    Slab Slab1 = Slab();
    CHECK(Slab1.n_layers == 0);
    CHECK(Slab1.layers.size() == 0);
    CHECK(Slab1.layer_interfaces.size() == 0);
    return EXIT_SUCCESS;
}

int test_make_slab_one_layer(){
    Slab Slab1 = Slab();
    std::function<double(double)> pdf = std::function<double(double)>(
        [](double mu){return 0.5;});
    std::function<double(double)> icdf = std::function<double(double)>(
        [](double rn){return 2.*rn-1.;});
    double k = 1.;
    double s = 1.;
    double t = 1.;
    double n = 1.;
    Slab1.add_layer(k, s, n, t, Layer::Phase(pdf, icdf));
    CHECK(Slab1.layer_interfaces[0][0] <= Constants::EPS);
    CHECK(abs(1. - Slab1.layer_interfaces[0][1]) <= Constants::EPS);
    CHECK(Slab1.n_layers == 1);
    CHECK(Slab1.layers.size() == 1);
    CHECK(Slab1.layer_interfaces.size() == 1);
    return EXIT_SUCCESS;
}

int test_make_slab_multi_layer(){
    Slab Slab1 = Slab();
    std::function<double(double)> pdf = std::function<double(double)>(
        [](double mu){return 0.5;});
    std::function<double(double)> icdf = std::function<double(double)>(
        [](double rn){return 2.*rn-1.;});
    double k = 1.;
    double s = 1.;
    double t = 1.;
    double n = 1.;
    Slab1.add_layer(k, s, n, t, Layer::Phase(pdf, icdf));
    Slab1.add_layer(k, s, n, t, Layer::Phase(pdf, icdf));
    Slab1.add_layer(k, s, n, t, Layer::Phase(pdf, icdf));
    CHECK(Slab1.n_layers == 3);
    CHECK(Slab1.layers.size() == 3);
    CHECK(Slab1.layer_interfaces.size() == 3);
    CHECK(abs(Slab1.layer_interfaces[0][1] - Slab1.layer_interfaces[1][0])
        <= Constants::EPS);
    CHECK(abs(Slab1.layer_interfaces[1][0] + Slab1.layer_interfaces[0][1] - 
        Slab1.layer_interfaces[2][0]) <= Constants::EPS);
    return EXIT_SUCCESS;
}

int test_make_detector_2d(){
    double det_ang = 45.;
    double det_ph_ang = 0.;
    // All angles within the classes themselves are in radians
    Detector Detector1 = Detector(det_ang*Constants::PI/180, 
        det_ph_ang*Constants::PI/180);
    CHECK(abs(1. - (pow(Detector1.d_vec[0], 2.) + pow(Detector1.d_vec[1], 2.)
        + pow(Detector1.d_vec[2], 2.))) < Constants::EPS);
    return EXIT_SUCCESS;
}

int test_make_detector_3d(){
    double det_ang = 45.;
    double det_ph_ang = 45.;
    Detector Detector1 = Detector(det_ang*Constants::PI/180, 
        det_ph_ang*Constants::PI/180);
    CHECK(abs(1. - (pow(Detector1.d_vec[0], 2.) + pow(Detector1.d_vec[1], 2.)
        + pow(Detector1.d_vec[2], 2.))) < Constants::EPS);
    return EXIT_SUCCESS;
}

int test_detect_photon_intensity(){
    // Simply checking that it runs
    // TODO: Need better test here
    double det_ang = 0.;
    double det_ph_ang = 0.;
    Detector Detector1 = Detector(det_ang*Constants::PI/180, 
        det_ph_ang*Constants::PI/180);
    double wght = 1.;
    double trns_frc = 1.;
    double tau_to_dtctr = 1.;
    double p_scat = 1.;
    Detector1.detect_photon(wght, trns_frc, tau_to_dtctr, p_scat);
    return EXIT_SUCCESS;
}

int test_detect_photon_flux(){
    double det_ang = 0.;
    double det_ph_ang = 0.;
    Detector Detector1 = Detector(det_ang*Constants::PI/180, 
        det_ph_ang*Constants::PI/180);
    double wght = 1.;
    Detector1.detect_diffuse(wght);
    CHECK(Detector1.diffuse == wght);
    return EXIT_SUCCESS;
}

int test_make_photon_2d(){
    double inc_ang = 45.;
    double inc_ph_ang = 0.;
    double det_ang = 45.;
    double det_ph_ang = 0.;
    Detector Detector1 = Detector(det_ang*Constants::PI/180, 
        det_ph_ang*Constants::PI/180);
    std::vector<Detector> Detectors{Detector1, Detector1}; // Dummy bot det
    Slab Slab1 = Slab();
    std::function<double(double)> pdf = std::function<double(double)>(
        [](double mu){return 0.5;});
    std::function<double(double)> icdf = std::function<double(double)>(
        [](double rn){return 2.*rn-1.;});
    double k = 1.;
    double s = 1.;
    double t = 1e8;
    double n = 1.;
    Slab1.add_layer(k, s, n, t, Layer::Phase(pdf, icdf));
    std::mt19937 gen(0);
    Photon Photon1 = Photon(Slab1, Detectors, gen, inc_ang*Constants::PI/4.,
        inc_ph_ang*Constants::PI/4.);
    CHECK(abs(1. - (pow(Photon1.phtn_vec[0], 2.) + pow(Photon1.phtn_vec[1], 2.)
        + pow(Photon1.phtn_vec[2], 2.))) < Constants::EPS);
    return EXIT_SUCCESS;
}

int test_make_photon_3d(){
    double inc_ang = 45.;
    double inc_ph_ang = 45.;
    double det_ang = 45.;
    double det_ph_ang = 0.;
    Detector Detector1 = Detector(det_ang*Constants::PI/180, 
        det_ph_ang*Constants::PI/180);
    std::vector<Detector> Detectors{Detector1, Detector1}; // Dummy bot det
    Slab Slab1 = Slab();
    std::function<double(double)> pdf = std::function<double(double)>(
        [](double mu){return 0.5;});
    std::function<double(double)> icdf = std::function<double(double)>(
        [](double rn){return 2.*rn-1.;});
    double k = 1.;
    double s = 1.;
    double t = 1e8;
    double n = 1.;
    Slab1.add_layer(k, s, n, t, Layer::Phase(pdf, icdf));
    std::mt19937 gen(0);
    Photon Photon1 = Photon(Slab1, Detectors, gen, inc_ang*Constants::PI/4.,
        inc_ph_ang*Constants::PI/4.);
    CHECK(abs(1. - (pow(Photon1.phtn_vec[0], 2.) + pow(Photon1.phtn_vec[1], 2.)
        + pow(Photon1.phtn_vec[2], 2.))) < Constants::EPS);
    return EXIT_SUCCESS;
}

int test_photon_escp(){
    double det_ang = 45.;
    double det_ph_ang = 0.;
    Detector Detector1 = Detector(det_ang*Constants::PI/180, 
        det_ph_ang*Constants::PI/180);
    std::vector<Detector> Detectors{Detector1, Detector1}; // Dummy bot det
    Slab Slab1 = Slab();
    std::function<double(double)> pdf = std::function<double(double)>(
        [](double mu){return 0.5;});
    std::function<double(double)> icdf = std::function<double(double)>(
        [](double rn){return 2.*rn-1.;});
    double k = 1.;
    double s = 1.;
    double t = 1e8;
    double n = 1.;
    Slab1.add_layer(k, s, n, t, Layer::Phase(pdf, icdf));
    std::mt19937 gen(0);
    Photon Photon1 = Photon(Slab1, Detectors, gen);
    Photon1.escp();
    CHECK(!Photon1.in_slb);
    bool cant_trvl = false;
    try{
        Photon1.trvl_to_annhltn();
    } catch(const std::logic_error& e) {
        cant_trvl = true;
    }
    CHECK(cant_trvl);
    return EXIT_SUCCESS;
}

int main(int argc, char *argv[]){
    if (argc == 1) {
        throw(std::logic_error("mcrt_tests must be called with command line "
            "argument specifying test number!"));
        return 0;
    }
    int test_num;
    try {
        test_num = std::stoi(argv[1]);
    } catch (...) {
        throw std::invalid_argument("command line argument to mcrt_tests must "
            "be convertible to integer!");
    }
    switch (test_num){
        case 1:
            return test_layer_explicit_phase();
        case 2:
            return test_make_slab();
        case 3:
            return test_make_slab_one_layer();
        case 4:
            return test_make_slab_multi_layer();
        case 5:
            return test_make_detector_2d();
        case 6:
            return test_make_detector_3d();
        case 7:
            return test_detect_photon_intensity();
        case 8:
            return test_detect_photon_flux();
        case 9:
            return test_make_photon_2d();
        case 10:
            return test_make_photon_3d();
        case 11:
            return test_photon_escp();
        default:
            return EXIT_FAILURE;
    }
}

