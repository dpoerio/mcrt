#include <math.h>
#include <chrono>
#include <exception>
#include <iostream>
#include <random>
#include <stdexcept>
#include <functional>
#include "detector.h"
#include "layer.h"
#include "slab.h"
#include "photon.h"
#include "constants.h"
#include "phase_functions.h"
#ifdef _OPENMP
#include <omp.h>
#endif

int main() {
    // Run a Monte Carlo radiative transfer simulation
    int nphotons = 5000000;
    double pht_ang = 20;
    double top_det_ang = 30.;
    double bot_det_ang = 150.;
    double top_det_ph_ang = 50.;
    double bot_det_ph_ang = 10.;
    double pht_ph_ang = 40.;
    std::vector<Detector> Detectors;
    Detectors.push_back(Detector(top_det_ang*Constants::PI/180, 
        top_det_ph_ang*Constants::PI/180));


    Slab Slab1 = Slab();
    double k = 0.2;
    double s = 0.1;
    double n = 1.0;
    double t = 2.0;
    double g = 0.7;
    std::function<double(double)> hgl = std::function<double(double)>(
        [g](double mu){return Phase_Functions::hg(g, mu);});
    std::function<double(double)> hgl2 = std::function<double(double)>(
        [g](double rn){return Phase_Functions::ihg(g, rn);});
    Slab1.add_layer(k, s, n, t, Layer::Phase(hgl, hgl2));


    k = 0.1;
    s = 0.8;
    n = 1.0;
    t = 1.0;
    g = -0.7;
    hgl = std::function<double(double)>([g](double mu){
        return Phase_Functions::hg(g, mu);});
    hgl2 = std::function<double(double)>([g](double rn){
        return Phase_Functions::ihg(g, rn);});
    Slab1.add_layer(k, s, n, t, Layer::Phase(hgl, hgl2));

    k = 0.3;
    s = 0.3;
    n = 1.0;
    t = 1.0;
    hgl = std::function<double(double)>([](double mu){
        return Phase_Functions::rlgh(mu);});
    hgl2 = std::function<double(double)>([](double rn){
        return Phase_Functions::irlgh(rn);});
    Slab1.add_layer(k, s, n, t, Layer::Phase(hgl, hgl2));



    k = 0.2;
    s = 0.5;
    n = 1.0;
    t = 1.0;
    hgl = std::function<double(double)>([](double mu){
        return Phase_Functions::rlgh(mu);});
    hgl2 = std::function<double(double)>([](double rn){
        return Phase_Functions::irlgh(rn);});
    Slab1.add_layer(k, s, n, t, Layer::Phase(hgl, hgl2));


    k = 0.5;
    s = 0.1;
    n = 1.0;
    t = 1.0;
    g = 0.0;
    hgl = std::function<double(double)>([g](double mu){
        return Phase_Functions::hg(g, mu);});
    hgl2 = std::function<double(double)>([g](double rn){
        return Phase_Functions::ihg(g, rn);});
    Slab1.add_layer(k, s, n, t, Layer::Phase(hgl, hgl2));


    const double slab_thickness = Slab1.layer_interfaces[Slab1.n_layers-1][1];
    Detectors.push_back(Detector(bot_det_ang*Constants::PI/180, 
        bot_det_ph_ang*Constants::PI/180, slab_thickness));
#ifdef _OPENMP
    int thrds = omp_get_max_threads();
#else
    int thrds = 1;
#endif
    std::vector<std::mt19937> prngs;
    for (int i = 0; i < thrds; i++){
        // Create PRNGs for each thread
        prngs.push_back(std::mt19937(i));
    }
#pragma omp parallel for shared (Slab1, Detectors)
    for (int i = 0; i < nphotons; i++) {
#ifdef _OPENMP
        int thrdi = omp_get_thread_num();
#else
        int thrdi = 0;
#endif
        Photon Photon_i(Slab1, Detectors, prngs[thrdi], 
            pht_ang*Constants::PI/180, pht_ph_ang*Constants::PI/180);
        Photon_i.trvl_to_annhltn();
    }
    std::cout << "Final top intensity is " << Detectors[0].intnsty_dtctd / nphotons 
        << std::endl;
    std::cout << "Final bot intensity is " << Detectors[1].intnsty_dtctd / nphotons
        << std::endl;
    for (int i = 0; i < Slab1.n_layers; i++) {
        std::cout << "Final layer " << i << " upward flux is " 
            << Slab1.layer_upwd_fluxes[i] / nphotons << std::endl;
        std::cout << "Final layer " << i << " downward flux is " 
            << Slab1.layer_dwwd_fluxes[i] / nphotons << std::endl;
    }
    return 0;
}

