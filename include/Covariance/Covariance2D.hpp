#pragma once

namespace Covariance {

   /*! The 2D isotropic version of Covariance Tracing.
    *  This code assumes that light transport in 2 dimensions and only
    *  store a 2D covariance matrix for space and angles.
    */
   struct Covariance2D {

      float matrix[3];

   };
}
