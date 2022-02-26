#include "slab.h"
#include "layer.h"
#include <string>
#include <stdexcept>

void Slab::add_layer(double k, double s, double n, double t, Layer::phass p){
    if ((k <= 0) | (s <= 0) | (t <= 0) | (n <= 0)){
        throw std::domain_error(std::string("Optical properties and "
            "thickness must be positive!"));
    }
    Layer NewLayer(k, s, n, t, p);
    layers.push_back(NewLayer);
    std::vector<double> layer_interface;
    if (n_layers == 0) {
        layer_interface.push_back(0);
        layer_interface.push_back(t);
        layer_interfaces.push_back(layer_interface);
        n_layers++;
    } else {
        layer_interface.push_back(layer_interfaces[n_layers - 1][1]);
        layer_interface.push_back(layer_interfaces[n_layers - 1][1] + t);
        layer_interfaces.push_back(layer_interface);
        n_layers++;
    }
    layer_dwwd_fluxes.push_back(0.);
    layer_upwd_fluxes.push_back(0.);
};

void Slab::accum_flux(double intnsty, int lyr, int dir){
    if (dir == 0) {
#pragma omp atomic
        layer_dwwd_fluxes[lyr] += intnsty;
    } else if (dir == 1) {
#pragma omp atomic
        layer_upwd_fluxes[lyr] += intnsty;
    }
};
