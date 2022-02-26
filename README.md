# MCRT

MCRT is a C++11 program for performing Monte Carlo radiative transfer
simulations of multi-layer, horizontally infinite plane parallel geometries. The
program is modeled in terms of well understood constructs in radiative transfer:
there are self contained classes for layers, slabs, detectors, and photons (the
workhorse). Simulations are run by writing a driver describing the system being
modeled (a slab consisting of an arbitrary number of layers, a detector, and as
many photons as desired being propagated through the system), and linking
against a library containing the key classes.

MCRT includes the following features:

* Complete handling of refractive index as an inherent optical property (IOP) of
  the layer.
* Intensity/radiance detectors at the top and bottom of the system, as well as
  flux/irradiance detectors at all interfaces, including the top and bottom of
the system.
* Arbitrary photon initiation direction and intensity detector direction.
* Parallelization via OpenMP.
* A test suite with integration tests comparing results to
  [DISORT](http://www.rtatmocn.com/disort/) and
[MCML](http://coilab.caltech.edu/mc.html), two public domain radiative transfer
solvers.

In the future, we would like to add the following features:

* Support for arbitrary position of intensity/radiance detectors. This can be
  achieved presently for flux/irradiance by creating an "interface" (even if the
two layers are exactly the same) where the detection is desired, since flux is
detected at layer interfaces.
* Support for arbitrary, user defined scattering phase functions, without
  needing to specify an inverse cumulative distribution function (ICDF) for
sampling. We intend to validate the phase function, and then compute an
efficient ICDF via numerical integration to generate a look up table (LUT), with
exact values estimated via linear interpolation from the LUT during sampling.
* Polarization effects
* Support for horizontal inhomogeneities (3D gridding of the system)
* Validation of results by comparison against other packages.  Currently, we
  validate against [DISORT](http://www.rtatmocn.com/disort/) for radiances and
fluxes in multilayer, non-refracting systems, and
[MCML](http://coilab.caltech.edu/mc.html) for fluxes in refracting systems. We
are confident in the results, but if you know of a suitable package for
comparison or can help us validate our implementation against experimental data,
feel free to reach out.

# Compilation

To compile MCRT, run the CMake based build script, `build.sh`. It will create
a build directory, `build`, and execute an out of tree build. Ultimately, it
will compile the main classes in a library archive, `libmcrt_classes.a`, and
finally an executable, `mcrt`, for the simulation defined in `mcrt.cpp`. You can
also look at `mcrt_integration_tests.cpp` for examples of other drivers.

# Testing

To run the tests, compile MCRT by running the build script, `cd` to the `build`
directory, and run `ctest`.

# Acknowledgments

The main methods implemented in MCRT are described in:

* [Radiative transfer in the cloudy atmosphere B.  Mayer EPJ Web of Conferences
  1 75-99 (2009) DOI:
10.1140/epjconf/e2009-00912-1](https://doi.org/10.1140/epjconf/e2009-00912-1)
* [H. Hatcher Tynes, George W. Kattawar, Eleonora P. Zege, Iosif L. Katsev,
  Alexander S. Prikhach, and Ludmila I. Chaikovskaya, "Monte Carlo and
multicomponent approximation methods for vector radiative transfer by use of
effective Mueller matrix calculations," Appl. Opt. 40, 400-412 (2001), DOI:
10.1364/AO.40.000400](https://doi.org/10.1364/AO.40.000400)

