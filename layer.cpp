#include "layer.h"

Layer::Layer(double absr, double scat, double rfrc, double thck, phass phase) : 
        k(absr), 
        s(scat), 
        n(rfrc), 
        t(thck), 
        p(phase) {}
        // TODO: Add a validation of the phase functions by numerical
        // integration

