#ifndef SLAB_H
#define SLAB_H

#include "layer.h"
#include <vector>

struct Slab {
    // The slab is the full plane parallel system, composed of many layers
    std::vector<Layer> layers;
    std::vector<std::vector<double>> layer_interfaces;
    std::vector<double> layer_upwd_fluxes;
    std::vector<double> layer_dwwd_fluxes;
    int n_layers = 0;

    void add_layer(double k, double s, double n, double t, Layer::phass p);

    void accum_flux(double intnsty, int lyr, int dir);
};

#endif
