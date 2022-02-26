#include "detector.h"
#include "layer.h"
#include "photon.h"
#include "slab.h"
#include "constants.h"
#include "phase_functions.h"
#include <random>
#include <string>
#include <functional>
#include <stdexcept>
#include <cmath>
#include "test_helpers.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/* NOTES:
DISORT references using disort4.0.99, 200 discrete ordinate streams, 299 phase
function moments. In some cases, intensities don't seem to match with highly
peaked phased functions.

mcml is used for some tests of flux/irradiance with different RIs.

Tests with RI differences at arbitrary detector location and different phase
function across layers are only tested against previous MCRT results, as I've
not used any other software with this feature.
*/

inline double percent_error(double target_val, double mcrt_val){
    double pcte = std::abs(target_val - mcrt_val)/target_val;
    return pcte;
}

static int nphotons = 5000000;

int test_two_layer_angld_pht_mtch_ri(){
    double pht_ang = 30.;
    double top_det_ang = 0.;
    double det_ph_ang = 0.;
    double pht_ph_ang = 0.;
    Detector Detector1 = Detector(top_det_ang*Constants::PI/180, 
        det_ph_ang*Constants::PI/180);
    std::vector<Detector>  Detectors{Detector1, Detector1}; // Just a dummy
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
    t = 1e10;
    g = -0.7;
    hgl = std::function<double(double)>([g](double mu){
        return Phase_Functions::hg(g, mu);});
    hgl2 = std::function<double(double)>([g](double rn){
        return Phase_Functions::ihg(g, rn);});
    Slab1.add_layer(k, s, n, t, Layer::Phase(hgl, hgl2));
#ifdef _OPENMP
    int thrds = omp_get_max_threads();
#else
    int thrds = 1;
#endif
    std::vector<std::mt19937> prngs;
    for (int i = 0; i < thrds; i++){
        // Create PRNGs for each thread
        std::mt19937 gen(i);
        prngs.push_back(gen);
    }
#pragma omp parallel for shared (Slab1, Detector1)
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
    double disort_intnsty = 7.57978e-02;
    double disort_upwd_flux = 1.79753e-01;
    double mcrt_intnsty = Detectors[0].intnsty_dtctd / nphotons;
    double mcrt_upwd_flux = Slab1.layer_upwd_fluxes[0] / nphotons;
    double pcte_intnsty = percent_error(disort_intnsty, mcrt_intnsty);
    double pcte_upwd_flux = percent_error(disort_upwd_flux, mcrt_upwd_flux);
    // TODO: Downward direct intensity is 1 in this case since the detector is
    // at the top of the slab
    CHECK(pcte_intnsty < Constants::TST_THRSH);
    CHECK(pcte_upwd_flux < Constants::TST_THRSH);
    return EXIT_SUCCESS;
}

int test_two_layer_angld_det_mtch_ri_flip(){
    // Same as previous test, pht_ang and det_ang switched
    double pht_ang = 0.;
    double top_det_ang = 30.;
    double det_ph_ang = 0.;
    double pht_ph_ang = 0.;
    Detector Detector1 = Detector(top_det_ang*Constants::PI/180, 
        det_ph_ang*Constants::PI/180);
    std::vector<Detector>  Detectors{Detector1, Detector1}; // Just a dummy
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
    t = 1e10;
    g = -0.7;
    hgl = std::function<double(double)>([g](double mu){
        return Phase_Functions::hg(g, mu);});
    hgl2 = std::function<double(double)>([g](double rn){
        return Phase_Functions::ihg(g, rn);});
    Slab1.add_layer(k, s, n, t, Layer::Phase(hgl, hgl2));
#ifdef _OPENMP
    int thrds = omp_get_max_threads();
#else
    int thrds = 1;
#endif
    std::vector<std::mt19937> prngs;
    for (int i = 0; i < thrds; i++){
        // Create PRNGs for each thread
        std::mt19937 gen(i);
        prngs.push_back(gen);
    }
#pragma omp parallel for shared (Slab1, Detector1)
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
    double disort_intnsty = 7.57978e-02;
    double disort_upwd_flux = 1.95076e-01; // The diffuse flux is not recip.
    double mcrt_intnsty = Detectors[0].intnsty_dtctd / nphotons;
    double mcrt_upwd_flux = Slab1.layer_upwd_fluxes[0] / nphotons;
    double pcte_intnsty = percent_error(disort_intnsty, mcrt_intnsty);
    double pcte_upwd_flux = percent_error(disort_upwd_flux, mcrt_upwd_flux);
    // TODO: Downward direct intensity is 1 in this case since the detector is
    // at the top of the slab
    CHECK(pcte_intnsty < Constants::TST_THRSH);
    CHECK(pcte_upwd_flux < Constants::TST_THRSH);
    return EXIT_SUCCESS;
}

int test_two_layer_angld_pht_diff_ri(){
    // Same as previous test, pht_ang and det_ang switched
    double pht_ang = 30.;
    double top_det_ang = 0.;
    double det_ph_ang = 0.;
    double pht_ph_ang = 0.;
    Detector Detector1 = Detector(top_det_ang*Constants::PI/180, 
        det_ph_ang*Constants::PI/180);
    std::vector<Detector>  Detectors{Detector1, Detector1}; // Just a dummy
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
    n = 1.5;
    t = 1e10;
    g = -0.7;
    hgl = std::function<double(double)>([g](double mu){
        return Phase_Functions::hg(g, mu);});
    hgl2 = std::function<double(double)>([g](double rn){
        return Phase_Functions::ihg(g, rn);});
    Slab1.add_layer(k, s, n, t, Layer::Phase(hgl, hgl2));
#ifdef _OPENMP
    int thrds = omp_get_max_threads();
#else
    int thrds = 1;
#endif
    std::vector<std::mt19937> prngs;
    for (int i = 0; i < thrds; i++){
        // Create PRNGs for each thread
        std::mt19937 gen(i);
        prngs.push_back(gen);
    }
#pragma omp parallel for shared (Slab1, Detector1)
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
    double mcml_upwd_flux = 0.13921757;
    double mcrt_upwd_flux = Slab1.layer_upwd_fluxes[0] / nphotons;
    double pcte_upwd_flux = percent_error(mcml_upwd_flux, mcrt_upwd_flux);
    // TODO: Downward direct intensity is 1 in this case since the detector is
    // at the top of the slab
    CHECK(pcte_upwd_flux < Constants::TST_THRSH);
    return EXIT_SUCCESS;
}

int test_two_layer_angld_det_diff_ri_flip(){
    // Same as previous test, pht_ang and det_ang switched
    double pht_ang = 0.;
    double top_det_ang = 30.;
    double det_ph_ang = 0.;
    double pht_ph_ang = 0.;
    Detector Detector1 = Detector(top_det_ang*Constants::PI/180,
        det_ph_ang*Constants::PI/180);
    std::vector<Detector>  Detectors{Detector1, Detector1}; // Just a dummy
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
    n = 1.5;
    t = 1e10;
    g = -0.7;
    hgl = std::function<double(double)>([g](double mu){
        return Phase_Functions::hg(g, mu);});
    hgl2 = std::function<double(double)>([g](double rn){
        return Phase_Functions::ihg(g, rn);});
    Slab1.add_layer(k, s, n, t, Layer::Phase(hgl, hgl2));
#ifdef _OPENMP
    int thrds = omp_get_max_threads();
#else
    int thrds = 1;
#endif
    std::vector<std::mt19937> prngs;
    for (int i = 0; i < thrds; i++){
        // Create PRNGs for each thread
        std::mt19937 gen(i);
        prngs.push_back(gen);
    }
#pragma omp parallel for shared (Slab1, Detector1)
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
    double mcml_upwd_flux = 0.15545996;
    double mcrt_upwd_flux = Slab1.layer_upwd_fluxes[0] / nphotons;
    double pcte_upwd_flux = percent_error(mcml_upwd_flux, mcrt_upwd_flux);
    // TODO: Downward direct intensity is 1 in this case since the detector is
    // at the top of the slab
    CHECK(pcte_upwd_flux < Constants::TST_THRSH);
    return EXIT_SUCCESS;
}

int test_two_layer_angld_pht_and_det_phi_rot0(){
    double pht_ang = 20.;
    double top_det_ang = 60.;
    double det_ph_ang = 0.;
    double pht_ph_ang = 0.;
    Detector Detector1 = Detector(top_det_ang*Constants::PI/180, 
        det_ph_ang*Constants::PI/180);
    std::vector<Detector>  Detectors{Detector1, Detector1}; // Just a dummy
    Slab Slab1 = Slab();
    double k = 0.2;
    double s = 0.1;
    double n = 1.0;
    double t = 2.0;
    double g = -0.9;
    std::function<double(double)> hgl = std::function<double(double)>(
        [g](double mu){return Phase_Functions::hg(g, mu);});
    std::function<double(double)> hgl2 = std::function<double(double)>(
        [g](double rn){return Phase_Functions::ihg(g, rn);});
    Slab1.add_layer(k, s, n, t, Layer::Phase(hgl, hgl2));
    k = 0.1;
    s = 0.8;
    n = 1.0;
    t = 1e10;
    g = 0.15;
    hgl = std::function<double(double)>([g](double mu){
        return Phase_Functions::hg(g, mu);});
    hgl2 = std::function<double(double)>([g](double rn){
        return Phase_Functions::ihg(g, rn);});
    Slab1.add_layer(k, s, n, t, Layer::Phase(hgl, hgl2));
#ifdef _OPENMP
    int thrds = omp_get_max_threads();
#else
    int thrds = 1;
#endif
    std::vector<std::mt19937> prngs;
    for (int i = 0; i < thrds; i++){
        // Create PRNGs for each thread
        std::mt19937 gen(i);
        prngs.push_back(gen);
    }
#pragma omp parallel for shared (Slab1, Detector1)
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
    double disort_intnsty = 2.53700e-02;
    double disort_upwd_flux = 1.97156e-01;
    double mcrt_intnsty = Detectors[0].intnsty_dtctd / nphotons;
    double mcrt_upwd_flux = Slab1.layer_upwd_fluxes[0] / nphotons;
    double pcte_intnsty = percent_error(disort_intnsty, mcrt_intnsty);
    double pcte_upwd_flux = percent_error(disort_upwd_flux, mcrt_upwd_flux);
    // TODO: Downward direct intensity is 1 in this case since the detector is
    // at the top of the slab
    CHECK(pcte_intnsty < Constants::TST_THRSH);
    CHECK(pcte_upwd_flux < Constants::TST_THRSH);
    return EXIT_SUCCESS;
}

int test_two_layer_angld_pht_and_det_phi_rot180(){
    double pht_ang = 20.;
    double top_det_ang = 60.;
    double det_ph_ang = 180.;
    double pht_ph_ang = 0.;
    Detector Detector1 = Detector(top_det_ang*Constants::PI/180, 
        det_ph_ang*Constants::PI/180);
    std::vector<Detector>  Detectors{Detector1, Detector1}; // Just a dummy
    Slab Slab1 = Slab();
    double k = 0.2;
    double s = 0.1;
    double n = 1.0;
    double t = 2.0;
    double g = -0.9;
    std::function<double(double)> hgl = std::function<double(double)>(
        [g](double mu){return Phase_Functions::hg(g, mu);});
    std::function<double(double)> hgl2 = std::function<double(double)>(
        [g](double rn){return Phase_Functions::ihg(g, rn);});
    Slab1.add_layer(k, s, n, t, Layer::Phase(hgl, hgl2));
    k = 0.1;
    s = 0.8;
    n = 1.0;
    t = 1e10;
    g = 0.15;
    hgl = std::function<double(double)>([g](double mu){
        return Phase_Functions::hg(g, mu);});
    hgl2 = std::function<double(double)>([g](double rn){
        return Phase_Functions::ihg(g, rn);});
    Slab1.add_layer(k, s, n, t, Layer::Phase(hgl, hgl2));
#ifdef _OPENMP
    int thrds = omp_get_max_threads();
#else
    int thrds = 1;
#endif
    std::vector<std::mt19937> prngs;
    for (int i = 0; i < thrds; i++){
        // Create PRNGs for each thread
        std::mt19937 gen(i);
        prngs.push_back(gen);
    }
#pragma omp parallel for shared (Slab1, Detector1)
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
    double disort_intnsty = 3.26817e-02;
    double disort_upwd_flux = 1.97156e-01;
    double mcrt_intnsty = Detectors[0].intnsty_dtctd / nphotons;
    double mcrt_upwd_flux = Slab1.layer_upwd_fluxes[0] / nphotons;
    double pcte_intnsty = percent_error(disort_intnsty, mcrt_intnsty);
    double pcte_upwd_flux = percent_error(disort_upwd_flux, mcrt_upwd_flux);
    // TODO: Downward direct intensity is 1 in this case since the detector is
    // at the top of the slab
    CHECK(pcte_intnsty < Constants::TST_THRSH);
    CHECK(pcte_upwd_flux < Constants::TST_THRSH);
    return EXIT_SUCCESS;
}

int test_two_layer_angld_det_and_pht_phi_rot180(){
    double pht_ang = 20.;
    double top_det_ang = 60.;
    double det_ph_ang = 0.;
    double pht_ph_ang = 180.;
    Detector Detector1 = Detector(top_det_ang*Constants::PI/180, 
        det_ph_ang*Constants::PI/180);
    std::vector<Detector>  Detectors{Detector1, Detector1}; // Just a dummy
    Slab Slab1 = Slab();
    double k = 0.2;
    double s = 0.1;
    double n = 1.0;
    double t = 2.0;
    double g = -0.9;
    std::function<double(double)> hgl = std::function<double(double)>(
        [g](double mu){return Phase_Functions::hg(g, mu);});
    std::function<double(double)> hgl2 = std::function<double(double)>(
        [g](double rn){return Phase_Functions::ihg(g, rn);});
    Slab1.add_layer(k, s, n, t, Layer::Phase(hgl, hgl2));
    k = 0.1;
    s = 0.8;
    n = 1.0;
    t = 1e10;
    g = 0.15;
    hgl = std::function<double(double)>([g](double mu){
        return Phase_Functions::hg(g, mu);});
    hgl2 = std::function<double(double)>([g](double rn){
        return Phase_Functions::ihg(g, rn);});
    Slab1.add_layer(k, s, n, t, Layer::Phase(hgl, hgl2));
#ifdef _OPENMP
    int thrds = omp_get_max_threads();
#else
    int thrds = 1;
#endif
    std::vector<std::mt19937> prngs;
    for (int i = 0; i < thrds; i++){
        // Create PRNGs for each thread
        std::mt19937 gen(i);
        prngs.push_back(gen);
    }
#pragma omp parallel for shared (Slab1, Detector1)
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
    double disort_intnsty = 3.26817e-02;
    double disort_upwd_flux = 1.97156e-01;
    double mcrt_intnsty = Detectors[0].intnsty_dtctd / nphotons;
    double mcrt_upwd_flux = Slab1.layer_upwd_fluxes[0] / nphotons;
    double pcte_intnsty = percent_error(disort_intnsty, mcrt_intnsty);
    double pcte_upwd_flux = percent_error(disort_upwd_flux, mcrt_upwd_flux);
    // TODO: Downward direct intensity is 1 in this case since the detector is
    // at the top of the slab
    CHECK(pcte_intnsty < Constants::TST_THRSH);
    CHECK(pcte_upwd_flux < Constants::TST_THRSH);
    return EXIT_SUCCESS;
}

int test_two_layer_thin_slab_intrfc_fluxes_trnsmn(){
    double pht_ang = 0;
    double top_det_ang = 0.;
    double det_ph_ang = 0.;
    double pht_ph_ang = 0.;
    std::vector<Detector> Detectors;
    Detectors.push_back(Detector(top_det_ang*Constants::PI/180,
        det_ph_ang*Constants::PI/180));
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
    t = 1;
    g = -0.7;
    hgl = std::function<double(double)>([g](double mu){
        return Phase_Functions::hg(g, mu);});
    hgl2 = std::function<double(double)>([g](double rn){
        return Phase_Functions::ihg(g, rn);});
    Slab1.add_layer(k, s, n, t, Layer::Phase(hgl, hgl2));
    const double slab_thickness = Slab1.layer_interfaces[Slab1.n_layers-1][1];
    double det_ang = 150.;
    Detectors.push_back(Detector(det_ang*Constants::PI/180, 
        det_ph_ang*Constants::PI/180, slab_thickness));
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
    double disort_top_intnsty = 2.00249e-01;
    double disort_bot_intnsty = 3.60812e-02;
    // double disort_intrfc_down_intsnty = 3.37606e-01;
    // double disort_intrfc_up_intsnty = 1.71054e-01;
    double disort_top_upwd_flux = 1.43781e-01;
    double disort_intrfc_upwd_flux = 2.46037e-01;
    // Downward direct + downward diffuse
    double disort_intrfc_dwwd_flux = 5.48812e-01 + 1.04391e-01;
    double disort_bot_dwwd_flux = 2.23130e-01 +  1.03818e-01;
    double mcrt_top_intnsty = Detectors[0].intnsty_dtctd / nphotons;
    double mcrt_bot_intnsty = Detectors[1].intnsty_dtctd / nphotons;
    double mcrt_top_upwd_flux = Slab1.layer_upwd_fluxes[0] / nphotons;
    double mcrt_intrfc_dwwd_flux = Slab1.layer_dwwd_fluxes[0] / nphotons;
    double mcrt_intrfc_upwd_flux = Slab1.layer_upwd_fluxes[1] / nphotons;
    double mcrt_bot_dwwd_flux = Slab1.layer_dwwd_fluxes[1] / nphotons;
    double pcte_top_intnsty = percent_error(disort_top_intnsty, 
        mcrt_top_intnsty);
    double pcte_bot_intnsty = percent_error(disort_bot_intnsty, 
        mcrt_bot_intnsty);
    double pcte_top_upwd_flux = percent_error(disort_top_upwd_flux, 
        mcrt_top_upwd_flux);
    double pcte_bot_dwwd_flux = percent_error(disort_bot_dwwd_flux, 
        mcrt_bot_dwwd_flux);
    double pcte_intrfc_dwwd_flux = percent_error(disort_intrfc_dwwd_flux, 
        mcrt_intrfc_dwwd_flux);
    double pcte_intrfc_upwd_flux = percent_error(disort_intrfc_upwd_flux, 
        mcrt_intrfc_upwd_flux);
    CHECK(pcte_bot_intnsty < Constants::TST_THRSH);
    CHECK(pcte_top_intnsty < Constants::TST_THRSH);
    CHECK(pcte_top_upwd_flux < Constants::TST_THRSH);
    CHECK(pcte_bot_dwwd_flux < Constants::TST_THRSH);
    CHECK(pcte_intrfc_upwd_flux < Constants::TST_THRSH);
    CHECK(pcte_intrfc_dwwd_flux < Constants::TST_THRSH);
    return EXIT_SUCCESS;
}

int test_two_layer_thin_slab_intrfc_edge_fluxes_tos_rad_diff_ri(){
    double pht_ang = 10;
    double top_det_ang = 0.;
    double det_ph_ang = 0.;
    double pht_ph_ang = 0.;
    std::vector<Detector> Detectors;
    Detectors.push_back(Detector(top_det_ang*Constants::PI/180, 
        det_ph_ang*Constants::PI/180));
    Slab Slab1 = Slab();
    double k = 0.2;
    double s = 0.1;
    double n = 2.7;
    double t = 1.0;
    double g = -0.7;
    std::function<double(double)> hgl = std::function<double(double)>(
        [g](double mu){return Phase_Functions::hg(g, mu);});
    std::function<double(double)> hgl2 = std::function<double(double)>(
        [g](double rn){return Phase_Functions::ihg(g, rn);});
    Slab1.add_layer(k, s, n, t, Layer::Phase(hgl, hgl2));
    k = 0.1;
    s = 0.8;
    n = 1.0;
    t = 2.0;
    g = 0.7;
    hgl = std::function<double(double)>([g](double mu){
        return Phase_Functions::hg(g, mu);});
    hgl2 = std::function<double(double)>([g](double rn){
        return Phase_Functions::ihg(g, rn);});
    Slab1.add_layer(k, s, n, t, Layer::Phase(hgl, hgl2));
    const double slab_thickness = Slab1.layer_interfaces[Slab1.n_layers-1][1];
    double det_ang = 180.;
    Detectors.push_back(Detector(det_ang*Constants::PI/180, 
        det_ph_ang*Constants::PI/180, slab_thickness));
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
    double mcml_top_upwd_flux = 0.23860179;
    double mcml_bot_dwwd_flux = 0.36221261;
    double mcrt_top_upwd_flux = Slab1.layer_upwd_fluxes[0] / nphotons;
    double mcrt_bot_dwwd_flux = Slab1.layer_dwwd_fluxes[1] / nphotons;
    double pcte_top_upwd_flux = percent_error(mcml_top_upwd_flux, 
        mcrt_top_upwd_flux);
    double pcte_bot_dwwd_flux = percent_error(mcml_bot_dwwd_flux, 
        mcrt_bot_dwwd_flux);
    CHECK(pcte_top_upwd_flux < Constants::TST_THRSH);
    CHECK(pcte_bot_dwwd_flux < Constants::TST_THRSH);
    return EXIT_SUCCESS;
}

int test_two_layer_diff_ri_reciprocity(){
    double pht_ang = 0;
    double det_ang = 70.;
    double det_ph_ang = 0.;
    double pht_ph_ang = 0.;
    std::vector<Detector> Detectors;
    Detectors.push_back(Detector(det_ang*Constants::PI/180, 
        det_ph_ang*Constants::PI/180));
    Slab Slab1 = Slab();
    double k = 0.2;
    double s = 0.1;
    double n = 2.7;
    double t = 1.0;
    double g = -0.7;
    std::function<double(double)> hgl = std::function<double(double)>(
        [g](double mu){return Phase_Functions::hg(g, mu);});
    std::function<double(double)> hgl2 = std::function<double(double)>(
        [g](double rn){return Phase_Functions::ihg(g, rn);});
    Slab1.add_layer(k, s, n, t, Layer::Phase(hgl, hgl2));
    k = 0.1;
    s = 0.8;
    n = 1.0;
    t = 2.0;
    g = 0.7;
    hgl = std::function<double(double)>([g](double mu){
        return Phase_Functions::hg(g, mu);});
    hgl2 = std::function<double(double)>([g](double rn){
        return Phase_Functions::ihg(g, rn);});
    Slab1.add_layer(k, s, n, t, Layer::Phase(hgl, hgl2));
    const double slab_thickness = Slab1.layer_interfaces[Slab1.n_layers-1][1];
    det_ang = 180.;
    Detectors.push_back(Detector(det_ang*Constants::PI/180, 
        det_ph_ang*Constants::PI/180, slab_thickness));
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
    double tos_intnsty_angled_det = Detectors[0].intnsty_dtctd/nphotons;
    pht_ang = 70;
    det_ang = 0.;
    Detectors.clear();
    Detectors.push_back(Detector(det_ang*Constants::PI/180, 
        det_ph_ang*Constants::PI/180));
    det_ang = 180.;
    Detectors.push_back(Detector(det_ang*Constants::PI/180, 
        det_ph_ang*Constants::PI/180, slab_thickness));
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
    double tos_intnsty_angled_pht = Detectors[0].intnsty_dtctd/nphotons;
    double pct_diff = percent_error(tos_intnsty_angled_det, 
        tos_intnsty_angled_pht);
    CHECK(pct_diff < Constants::TST_THRSH);
    return EXIT_SUCCESS;
}

int test_mltylyr_same_ri(){
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
    double disort_top_intnsty = 3.11075e-02;
    double disort_bot_intnsty = 1.52039e-02;
    double mcrt_top_intnsty = Detectors[0].intnsty_dtctd / nphotons;
    double mcrt_bot_intnsty = Detectors[1].intnsty_dtctd / nphotons;
    double top_pcte = percent_error(disort_top_intnsty, mcrt_top_intnsty);
    double bot_pcte = percent_error(disort_bot_intnsty, mcrt_bot_intnsty);
    CHECK(top_pcte < Constants::TST_THRSH);
    CHECK(bot_pcte < Constants::TST_THRSH);
    return EXIT_SUCCESS;
}


// TODO: Integration test - test for equality with equivalent phi rotations of
// the photons and the detector

int main(int argc, char *argv[]){
    if (argc == 1) {
        throw(std::logic_error("mcrt_tests must be called with command line "
            "argument specifying test number!"));
        return EXIT_SUCCESS;
    }
    int test_num;
    try {
        test_num = std::stoi(argv[1]);
    } catch (...) {
        throw std::invalid_argument("command line argument to "
            "mcrt_integration_tests must be convertible to integer!");
    }
    switch (test_num){
        case 1:
            return test_two_layer_angld_pht_mtch_ri();
        case 2:
            return test_two_layer_angld_det_mtch_ri_flip();
        case 3:
            return test_two_layer_angld_pht_diff_ri();
        case 4:
            return test_two_layer_angld_det_diff_ri_flip();
        case 5:
            return test_two_layer_angld_pht_and_det_phi_rot0();
        case 6:
            return test_two_layer_angld_pht_and_det_phi_rot180();
        case 7:
            return test_two_layer_angld_det_and_pht_phi_rot180();
        case 8:
            return test_two_layer_thin_slab_intrfc_fluxes_trnsmn();
        case 9:
            return test_two_layer_thin_slab_intrfc_edge_fluxes_tos_rad_diff_ri();
        case 10:
            return test_two_layer_diff_ri_reciprocity();
        case 11:
            return test_mltylyr_same_ri();
        default:
            return EXIT_FAILURE;
    }
}
