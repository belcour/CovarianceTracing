# Covariance tracing

This repository contains the source code to use [Covariance tracing](https://hal.inria.fr/hal-00814164) in a rendering engine. This C++ code is header-only and has no dependency except of STL.

To use covariance tracing in your software, simply include the corresponding header file available in `include` and specify the two arguments of the template class (floatting point representation and Vector class):

      #include <Covariance/Covariance4D.hpp>

      Covariance::Covariance4D<float, Vector> cov;

See the different tutorials and source code documentation to get a better view on how to use this class in your code. To build the tutorials, remember to load the `tinyexr` git submodule.
