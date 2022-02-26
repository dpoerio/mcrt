#ifndef LAYER_H
#define LAYER_H

#include <memory>
#include <functional>
#include <vector>
#include <stdexcept>
#include <string>

struct Layer {

    typedef std::vector<std::function<double(double)>> phass;
    double k;
    double s;
    double n;
    double t;
    phass p;

    Layer(double absr, double scat, double rfrc, double thck, phass phase);
        /*This program uses the convention that the phase function is a f(cos
        theta) and is normalized such that the integral from -1 to 1 is 1.
        Inverse phase function CDFs are define cos theta = f(U[0,1]).*/

    static phass Phase(std::function<double(double)> pdf, 
        std::function<double(double)> icdf){
        phass pf;
        pf.push_back(pdf);
        pf.push_back(icdf);
        // TODO: Validate PDF and ICDF
        return pf;
    }


    static phass Phase(std::function<double(double)> pdf) {
        phass pf;
        pf.push_back(pdf);
        // TODO: Calculate the ICDF from the PDF
        throw std::logic_error(std::string("Automatic calculation of inverse \
            CDF from arbitrary phase function not yet implemented"));
        //return pf;
    }

};

#endif
