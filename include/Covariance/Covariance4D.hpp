#pragma once

// STL includes
#include <array>
#include <cmath>
#include <limits>

// Local includes
#include "Matrix.hpp"

#define COV_MAX_FLOAT 1.0E+5
#define COV_MIN_FLOAT 1.0E-5

#define USE_WOODBURY_IDENTITY

namespace Covariance {

   /* The 4D version of Covariance Tracing.
    * This is the timeless version of Covariance Tracing. It allows to compute
    * covariance information of the local lightfield around the central ray
    * position (Belcour 2013).
    *
    * The covariance matrix represent the covariance of the local light field
    * in the (X, Y, U, V) local coordinate frame.
    *
    * The covariance matrix is a 4x4 symetric Matrix and we only deal with the
    * upper triangle for all the operations. Variance terms are at indices
    * 0, 2, 5, 9.
    *
    * The matrix is indexed the following way:
    *
    * C =  ( 0  1  3  6)
    *      ( *  2  4  7)
    *      ( *  *  5  8)
    *      ( *  *  *  9)
    *
    * Thus the full indexing of the matrix is the following:
    *
    * C =  ( 0  1  3  6)
    *      ( 1  2  4  7)
    *      ( 3  4  5  8)
    *      ( 6  7  8  9)
    */
   template<class Vector, typename Float>
   struct Covariance4D {

      std::array<Float, 10> matrix;
      Vector x, y, z;


      ////////////////////////
      //  Atomic operators  //
      ////////////////////////

      /* Travel operator
       *
       * 'd' distance of travel along the central ray
       */
      inline void Travel(Float d) {
         ShearAngleSpace(d, d);
      }

      /* Curvature operator
       *
       * 'kx' curvature along the X direction
       * 'ky' curvature along the Y direction
       */
      inline void Curvature(Float kx, Float ky) {
         ShearSpaceAngle(kx, ky);
      }

      /* Cosine operator
       *
       * 'wz' The incident direction's elevation in the local frame
       */
      inline void Cosine(Float wz) {
         const Float theta = acos(wz);
         const Float dist  = std::abs(0.5*M_PI-theta);
         const Float frequ = 2.0 / M_PI;
         const Float freqv = 1.0 / fmax(dist, 1.0E-10);
         matrix[5] += frequ*frequ;
         matrix[9] += freqv*freqv;
      }

      /* Reflection operator
       *
       * 'suu' the covariance of the BRDF along the X axis.
       * 'svv' the covariance of the BRDF along the Y axis.
       */
      inline void Reflection(Float suu, Float svv) {
         if(suu < COV_MAX_FLOAT && svv < COV_MAX_FLOAT) {
            ProductUV(fmax(suu, COV_MIN_FLOAT), fmax(svv, COV_MIN_FLOAT));
         }
      }


      /////////////////////////////
      //  Local Frame alignment  //
      /////////////////////////////

      /* Perform the projection of the incomming lightfield on the surface with
       * normal n. This function assumes that the surface normal is in the
       * opposite direction to the main vector of the lightfield.
       *
       * 'n' the surface normal.
       */
      inline void Projection(const Vector& n) {

         const auto cx = Vector::Dot(x, n);
         const auto cy = Vector::Dot(y, n);

         // Rotate the Frame to be aligned with plane.
         const Float alpha = (cx != 0.0) ? atan2(cx, cy) : 0.0;
         const Float c = cos(alpha), s = -sin(alpha);
         Rotate(c, s);

         // Scale the componnent that project by the cosine of the ray direction
         // and the normal.
         const Float cosine = Vector::Dot(z, n);
         ScaleY(std::abs(cosine));

         // Update direction vectors.
         x = c*x + s*y;
         z = (cosine < 0.0f) ? -n : n;
         y = (cosine < 0.0f) ?  Vector::Cross(x, z) : Vector::Cross(z, x);
      }

      /* Perform the projection of a lightfield defined on a surface to an
       * outgoing direction. This function assumes that the lightfield main
       * vector is the surface normal and that the outgoing vector is in the
       * same direction.
       *
       * 'd' the outgoing direction.
       */
      inline void InverseProjection(const Vector& d) {

         const auto cx = Vector::Dot(x, d);
         const auto cy = Vector::Dot(y, d);

         // Rotate the Frame to be aligned with plane.
         const Float alpha = (cx != 0.0) ? atan2(cx, cy) : 0.0;
         const Float c = cos(alpha), s = -sin(alpha);
         Rotate(c, s); // Rotate of -alpha

         // Scale the componnent that project by the inverse cosine of the ray
         // direction and the normal.
         const Float cosine = Vector::Dot(z, d);
         if(cosine < 0.0f) {
            ScaleV(-1.0f);
            ScaleU(-1.0f);
         }
         ScaleY(1.0/fmax(fabs(cosine), COV_MIN_FLOAT));

         // Update direction vectors.
         x = c*x + s*y;
         z = d;
         y = Vector::Cross(z, x);
      }


      /* Note: for the case of tracking the local frame of the light field
       * performing the symmetry makes not difference since the local frame
       * is ajusted with respect to symmetry.
       */
      inline void Symmetry() {
         matrix[3] = -matrix[3];
         matrix[4] = -matrix[4];
         matrix[6] = -matrix[6];
         matrix[7] = -matrix[7];
      }


      ////////////////////////
      //   Matrix scaling   //
      ////////////////////////

      inline void ScaleX(Float alpha) {
         matrix[0] *= alpha*alpha;
         matrix[1] *= alpha;
         matrix[3] *= alpha;
         matrix[6] *= alpha;
      }

      inline void ScaleY(Float alpha) {
         matrix[1] *= alpha;
         matrix[2] *= alpha*alpha;
         matrix[4] *= alpha;
         matrix[7] *= alpha;
      }

      inline void ScaleU(Float alpha) {
         matrix[3] *= alpha;
         matrix[4] *= alpha;
         matrix[5] *= alpha*alpha;
         matrix[8] *= alpha;
      }

      inline void ScaleV(Float alpha) {
         matrix[6] *= alpha;
         matrix[7] *= alpha;
         matrix[8] *= alpha;
         matrix[9] *= alpha*alpha;
      }


      ////////////////////////
      //   Matrix shearing  //
      ////////////////////////

      // Shear the Spatial (x,y) domain by the Angular (u,v).
      // \param cx amount of shear along the x direction.
      // \param cy amount of shear along the y direction.
      inline void ShearSpaceAngle(Float cx, Float cy) {
         matrix[0] += (matrix[5]*cx - 2*matrix[3])*cx;
         matrix[1] +=  matrix[8]*cx*cy - (matrix[4]*cy + matrix[6]*cx);
         matrix[2] += (matrix[9]*cy - 2*matrix[7])*cy;
         matrix[3] -=  matrix[5]*cx;
         matrix[4] -=  matrix[8]*cy;
         matrix[6] -=  matrix[8]*cx;
         matrix[7] -=  matrix[9]*cy;
      }

      // Shear the angular (U, V) domain by the spatial (X, Y) domain.
      // \param cu amount of shear along the U direction.
      // \param cy amount of shear along the V direction.
      inline void ShearAngleSpace(Float cu, Float cv) {
         matrix[5] += (matrix[0]*cu - 2*matrix[3])*cu;
         matrix[3] -=  matrix[0]*cu;
         matrix[8] +=  matrix[1]*cu*cv - (matrix[4]*cv + matrix[6]*cu);
         matrix[4] -=  matrix[1]*cv;
         matrix[6] -=  matrix[1]*cu;
         matrix[9] += (matrix[2]*cv - 2*matrix[7])*cv;
         matrix[7] -=  matrix[2]*cv;
      }


      ////////////////////////
      //   Matrix rotation  //
      ////////////////////////

      // alpha rotation angle
      inline void Rotate(Float alpha) {
         const Float c = cos(alpha);
         const Float s = sin(alpha);
         Rotate(c, s);
      }

      // 'c' the cosine of the rotation angle
      // 's' the sine of the rotation angle
      inline void Rotate(Float c, Float s) {
         const Float cs = c*s;
         const Float c2 = c*c;
         const Float s2 = s*s;

         const Float cov_xx = matrix[0];
         const Float cov_xy = matrix[1];
         const Float cov_yy = matrix[2];
         const Float cov_xu = matrix[3];
         const Float cov_yu = matrix[4];
         const Float cov_uu = matrix[5];
         const Float cov_xv = matrix[6];
         const Float cov_yv = matrix[7];
         const Float cov_uv = matrix[8];
         const Float cov_vv = matrix[9];

         // Rotation of the space
         matrix[0] = c2 * cov_xx + 2*cs * cov_xy + s2 * cov_yy;
         matrix[1] = (c2-s2) * cov_xy + cs * (cov_yy - cov_xx);
         matrix[2] = c2 * cov_yy - 2*cs * cov_xy + s2 * cov_xx;

         // Rotation of the angle
         matrix[5] = c2 * cov_uu + 2*cs * cov_uv + s2 * cov_vv;
         matrix[8] = (c2-s2) * cov_uv + cs * (cov_vv - cov_uu);
         matrix[9] = c2 * cov_vv - 2*cs * cov_uv + s2 * cov_uu;

         // Covariances
         matrix[3] = c2 * cov_xu + cs * (cov_xv + cov_yu) + s2 * cov_yv;
         matrix[4] = c2 * cov_yu + cs * (cov_yv - cov_xu) - s2 * cov_xv;
         matrix[6] = c2 * cov_xv + cs * (cov_yv - cov_xu) - s2 * cov_yu;
         matrix[7] = c2 * cov_yv - cs * (cov_xv + cov_yu) + s2 * cov_xu;
      }


      ////////////////////////
      //    Matrix inverse  //
      ////////////////////////

      /* Compute the inverse of the covariance matrix.
       *
       * We add an epsilon to the diagonal in order to ensure that the matrix
       * can be inverted.
       *
       * 'inverse' needs to be a 4x4 preallocated matrix (16 Floats)..
       */
      void InverseMatrix(Float* inverse) const {
         // Compute the inverse matrix. We need to add an epsilon to the
         // diagonal in order to ensure that the matrix can be inverted.
         inverse[ 0] = matrix[ 0] + COV_MIN_FLOAT;
         inverse[ 1] = matrix[ 1]; inverse[ 4] = matrix[ 1];
         inverse[ 5] = matrix[ 2] + COV_MIN_FLOAT;
         inverse[ 2] = matrix[ 3]; inverse[ 8] = matrix[ 3];
         inverse[ 6] = matrix[ 4]; inverse[ 9] = matrix[ 4];
         inverse[10] = matrix[ 5] + COV_MIN_FLOAT;
         inverse[ 3] = matrix[ 6]; inverse[12] = matrix[ 6];
         inverse[ 7] = matrix[ 7]; inverse[13] = matrix[ 7];
         inverse[11] = matrix[ 8]; inverse[14] = matrix[ 8];
         inverse[15] = matrix[ 9] + COV_MIN_FLOAT;

         if(!Inverse<Float>(inverse, 4)) { throw 1; }
      }



      ////////////////////////
      // Product of signals //
      ////////////////////////

      // Evaluate the covariance matrix of the product of the local lightfield
      // and a angularly varying only signal (like a BSDF).
      //
      // 'su' is the sigma u of the inverse angular signal's covariance matrix.
      // 'sv' the sigma v of the inverse angular signal's covariance matrix.
      inline void ProductUV(Float su, Float sv) {

         if(su == 0.0 && sv == 0.0) {
            return;
         }

#ifdef USE_WOODBURY_IDENTITY
         const Float cov_xu = matrix[3];
         const Float cov_yu = matrix[4];
         const Float cov_uu = matrix[5];
         const Float cov_xv = matrix[6];
         const Float cov_yv = matrix[7];
         const Float cov_uv = matrix[8];
         const Float cov_vv = matrix[9];

         // The following is an application of the Woodbury matrix identity for
         // a rank 2 C matrix. See:
         //         http://en.wikipedia.org/wiki/Woodbury_matrix_identity
         const Float sig_u = cov_uu+su, sig_v = cov_vv+sv;
         const Float og = (sig_u*sig_v - cov_uv*cov_uv);
         const Float  g = 1.0f / og;

         matrix[0] -=  g*(cov_xu*(sig_v*cov_xu-cov_uv*cov_xv) +
                          cov_xv*(sig_u*cov_xv-cov_uv*cov_xu));

         matrix[1] -=  g*(cov_yu*(sig_v*cov_xu-cov_uv*cov_xv) +
                          cov_yv*(sig_u*cov_xv-cov_uv*cov_xu));
         matrix[2] -=  g*(cov_yu*(sig_v*cov_yu-cov_uv*cov_yv) +
                          cov_yv*(sig_u*cov_yv-cov_uv*cov_yu));

         matrix[3] -=  g*(cov_uu*(sig_v*cov_xu-cov_uv*cov_xv) +
                          cov_uv*(sig_u*cov_xv-cov_uv*cov_xu));
         matrix[4] -=  g*(cov_uu*(sig_v*cov_yu-cov_uv*cov_yv) +
                          cov_uv*(sig_u*cov_yv-cov_uv*cov_yu));
         matrix[5] -=  g*(cov_uu*(sig_v*cov_uu-cov_uv*cov_uv) +
                          cov_uv*(sig_u*cov_uv-cov_uv*cov_uu));

         matrix[6] -=  g*(cov_uv*(sig_v*cov_xu-cov_uv*cov_xv) +
                          cov_vv*(sig_u*cov_xv-cov_uv*cov_xu));
         matrix[7] -=  g*(cov_uv*(sig_v*cov_yu-cov_uv*cov_yv) +
                          cov_vv*(sig_u*cov_yv-cov_uv*cov_yu));
         matrix[8] -=  g*(cov_uv*(sig_v*cov_uu-cov_uv*cov_uv) +
                          cov_vv*(sig_u*cov_uv-cov_uv*cov_uu));
         matrix[9] -=  g*(cov_uv*(sig_v*cov_uv-cov_uv*cov_vv) +
                          cov_vv*(sig_u*cov_vv-cov_uv*cov_uv));
#else
         // This method computes the inverse of the covariance matrix
         // explicitely, adds the inverse of the BRDF matrix to the
         // angular component, and then invert the matrix again to obtain
         // the covariance matrix. This is very slow method to check
         // the Woodbury formula.

         Float inverse[16];
         inverse[ 0] = matrix[ 0] + COV_MIN_FLOAT;
         inverse[ 1] = matrix[ 1]; inverse[ 4] = matrix[ 1];
         inverse[ 5] = matrix[ 2] + COV_MIN_FLOAT;
         inverse[ 2] = matrix[ 3]; inverse[ 8] = matrix[ 3];
         inverse[ 6] = matrix[ 4]; inverse[ 9] = matrix[ 4];
         inverse[10] = matrix[ 5] + COV_MIN_FLOAT;
         inverse[ 3] = matrix[ 6]; inverse[12] = matrix[ 6];
         inverse[ 7] = matrix[ 7]; inverse[13] = matrix[ 7];
         inverse[11] = matrix[ 8]; inverse[14] = matrix[ 8];
         inverse[15] = matrix[ 9] + COV_MIN_FLOAT;

         if(!Inverse<Float>(inverse, 4)) { throw 1; }
         inverse[10] += 1.0f/std::max<Float>(sv, COV_MIN_FLOAT);
         inverse[15] += 1.0f/std::max<Float>(su, COV_MIN_FLOAT);

         if(!Inverse<Float>(inverse, 4)) { throw 1; }
         matrix[ 0] = inverse[ 0];
         matrix[ 1] = inverse[ 1];
         matrix[ 2] = inverse[ 5];
         matrix[ 3] = inverse[ 2];
         matrix[ 4] = inverse[ 6];
         matrix[ 5] = inverse[10];
         matrix[ 6] = inverse[ 3];
         matrix[ 7] = inverse[ 7];
         matrix[ 8] = inverse[11];
         matrix[ 9] = inverse[15];
#endif
      }

      /////////////////////////////
      // Add covariance together //
      /////////////////////////////

      void Add(const Covariance4D& cov, Float L1=1.0f, Float L2=1.0f) {
         const Float L = L1+L2;
         if(L <= 0.0f) return;

         for(unsigned short i=0; i<10; ++i) {
            matrix[i] = (L1*matrix[i] + L2*cov.matrix[i]) / L;
         }
      }

      ////////////////////////////
      // Spatio-angular Filters //
      ////////////////////////////

      /*  Compute the spatio-angular extent of the space related to the
       *  covariance matrix's filter. The extent is provided as vectors
       *  'Dx', 'Dy', and 'Du', 'Dv'. Those vectors are the main axis of
       *  the filter's fooprint.

       *  Vector like 'Dx' and 'Dy' are expressed in the local tangent frame
       *  x, y using the first two components:
       *         Dx = [x, y, 0]

       *  To extract Dx, Dy, Du and Dv, we do an eigen-decomposition of the
       *  covariance matrix and use the normalized eigen-vectors as the axis
       *  of the extent and the eigen-values are the squared extent of the
       *  polygonal shape.
       *
       *  This spatio-angular polygonal shape can be used to specify an
       *  equivalent ray differential [Igehy 1999].
       */
      void Extent(Vector& Dx, Vector& Dy, Vector& Du, Vector& Dv) const {
         // Compute the inverse matrix. 
         Float inverse[16];
         this->InverseMatrix(inverse);

          // T = trace and D = det of the spatial submatrix
          Float T = inverse[0]+inverse[5];
          Float D = inverse[0]*inverse[5] - inverse[1]*inverse[1];

          // Solve the 2nd order polynomial roots of p(l) = l^2 - l T + D.
          // This gives us the eigen values.
          Float d  = 0.25*T*T - D;
          if(d < 0.0) { throw 1; } // No solution exists
          Float l1 = 0.5*T + sqrt(d);
          Float l2 = 0.5*T - sqrt(d);

          if(abs(inverse[1]) > COV_MIN_FLOAT) {
            Dx.x = l1 - inverse[5];
            Dx.y = inverse[1];
            Dx.z = 0.0;
            Dy.x = l2 - inverse[5];
            Dy.y = inverse[1];
            Dy.z = 0.0;

            Dx.Normalize();
            Dy.Normalize();

            Dx = sqrt(l1)/(2.0*M_PI) * Dx;
            Dy = sqrt(l2)/(2.0*M_PI) * Dy;
          } else {
            Dx.x = sqrt(inverse[0])/(2.0*M_PI);
            Dx.y = 0.0;
            Dx.z = 0.0;
            Dy.x = 0.0;
            Dy.y = sqrt(inverse[5])/(2.0*M_PI);
            Dy.z = 0.0;
          }

          // T = trace and D = det of the angluar submatrix
          T = inverse[10]+inverse[15];
          D = inverse[10]*inverse[15] - inverse[11]*inverse[11];

          // Solve the 2nd order polynomial roots of p(l) = l^2 - l T + D.
          // This gives us the eigen values.
          d  = 0.25*T*T - D;
          if(d < 0.0) { throw 1; } // No solution exists
          l1 = 0.5*T + sqrt(d);
          l2 = 0.5*T - sqrt(d);

          if(abs(inverse[11]) > COV_MIN_FLOAT) {
            Du.x = l1 - inverse[15];
            Du.y = inverse[11];
            Du.z = 0.0;
            Dv.x = l2 - inverse[15];
            Dv.y = inverse[11];
            Dv.z = 0.0;

            Du.Normalize();
            Dv.Normalize();

            Du = sqrt(l1)/(2.0*M_PI) * Du;
            Dv = sqrt(l2)/(2.0*M_PI) * Dv;
          } else {
            Du.x = sqrt(inverse[10])/(2.0*M_PI);
            Du.y = 0.0;
            Du.z = 0.0;
            Dv.x = 0.0;
            Dv.y = sqrt(inverse[15])/(2.0*M_PI);
            Dv.z = 0.0;
          }
      }


      /////////////////////
      // Spatial Filters //
      /////////////////////

      /*  Compute the spatial filter in primal space. This resumes to computing
       *  the inverse of the spatial submatrix from the inverse covariance
       *  matrix in frequency space.
       */
      void SpatialFilter(Float& sxx, Float& sxy, Float& syy) const {

         // Compute the inverse matrix. We need to add an epsilon to the
         // diagonal in order to ensure that the matrix can be inverted.
         Float inverse[16];
         inverse[ 0] = matrix[ 0] + COV_MIN_FLOAT;
         inverse[ 1] = matrix[ 1]; inverse[ 4] = matrix[ 1];
         inverse[ 5] = matrix[ 2] + COV_MIN_FLOAT;
         inverse[ 2] = matrix[ 3]; inverse[ 8] = matrix[ 3];
         inverse[ 6] = matrix[ 4]; inverse[ 9] = matrix[ 4];
         inverse[10] = matrix[ 5] + COV_MIN_FLOAT;
         inverse[ 3] = matrix[ 6]; inverse[12] = matrix[ 6];
         inverse[ 7] = matrix[ 7]; inverse[13] = matrix[ 7];
         inverse[11] = matrix[ 8]; inverse[14] = matrix[ 8];
         inverse[15] = matrix[ 9] + COV_MIN_FLOAT;

         if(!Inverse<Float>(inverse, 4)) { throw 1; }

         // The outgoing filter is the inverse submatrix of this inverse
         // matrix.
         Float det = (inverse[0]*inverse[5]-inverse[1]*inverse[1]) / pow(2.0*M_PI, 2);
         sxx =  inverse[5] / det;
         syy =  inverse[0] / det;
         sxy = -inverse[1] / det;
      }

      /*  Compute the spatial extent of the space related to the equivalent
       *  spatial filter to the covariance matrix. The extent is provided as
       * 'Dx' and 'Dy', the main axis of the filter's fooprint.
       *
       * 'Dx' and 'Dy' are express using the first two components:
       *        Dx = [x, y, 0]
       *  as they represent direction in the local frame x,y.
       *
       *  To extract Dx and Dy, we do an eigen-decomposition of the
       *  covariance matrix and use the normalized eigen-vectors as the
       *  axis of the extent and the eigen-values are the squared extent
       *  of the polygonal shape.
       */
      void SpatialExtent(Vector& Dx, Vector& Dy) const {
         // Compute the inverse matrix. We need to add an epsilon to the
         // diagonal in order to ensure that the matrix can be inverted.
         Float inverse[16];
         inverse[ 0] = matrix[ 0] + COV_MIN_FLOAT;
         inverse[ 1] = matrix[ 1]; inverse[ 4] = matrix[ 1];
         inverse[ 5] = matrix[ 2] + COV_MIN_FLOAT;
         inverse[ 2] = matrix[ 3]; inverse[ 8] = matrix[ 3];
         inverse[ 6] = matrix[ 4]; inverse[ 9] = matrix[ 4];
         inverse[10] = matrix[ 5] + COV_MIN_FLOAT;
         inverse[ 3] = matrix[ 6]; inverse[12] = matrix[ 6];
         inverse[ 7] = matrix[ 7]; inverse[13] = matrix[ 7];
         inverse[11] = matrix[ 8]; inverse[14] = matrix[ 8];
         inverse[15] = matrix[ 9] + COV_MIN_FLOAT;

         if(!Inverse<Float>(inverse, 4)) { throw 1; }

          // T = trace and D = det of the spatial submatrix
          Float T = inverse[0]+inverse[5];
          Float D = inverse[0]*inverse[5] - inverse[1]*inverse[1];

          // Solve the 2nd order polynomial roots of p(l) = l^2 - l T + D.
          // This gives us the eigen values.
          Float d  = 0.25*T*T - D;
          if(d < 0.0) { throw 1; } // No solution exists
          Float l1 = 0.5*T + sqrt(d);
          Float l2 = 0.5*T - sqrt(d);

          if(abs(inverse[1]) > COV_MIN_FLOAT) {
            Dx.x = l1 - inverse[5];
            Dx.y = inverse[1];
            Dx.z = 0.0;
            Dy.x = l2 - inverse[5];
            Dy.y = inverse[1];
            Dy.z = 0.0;

            Dx.Normalize();
            Dy.Normalize();

            Dx = sqrt(l1)/(2.0*M_PI) * Dx;
            Dy = sqrt(l2)/(2.0*M_PI) * Dy;
          } else {
            Dx.x = sqrt(inverse[0])/(2.0*M_PI);
            Dx.y = 0.0;
            Dx.z = 0.0;
            Dy.x = 0.0;
            Dy.y = sqrt(inverse[5])/(2.0*M_PI);
            Dy.z = 0.0;
          }
      }


      /////////////////////
      // Angular Filters //
      /////////////////////

      /* Compute the angular filter in primal space. The angular filter is
       * the 2D Gaussian parameters 'suu', 'suv', 'svv' such that the filter
       * is written as:
       *        f(u,v) = exp(- 0.5 (suu u^2 + 2 suv u v + svv v^2))
       * in the local tangent frame or directions..
       *
       * This resumes to computing the inverse of the spatial submatrix from
       * the inverse covariance matrix in frequency space.
       */
      void AngularFilter(Float& suu, Float& suv, Float& svv) const {

         // Compute the inverse matrix. We need to add an epsilon to the
         // diagonal in order to ensure that the matrix can be inverted.
         Float inverse[16];
         inverse[ 0] = matrix[ 0] + COV_MIN_FLOAT;
         inverse[ 1] = matrix[ 1]; inverse[ 4] = matrix[ 1];
         inverse[ 5] = matrix[ 2] + COV_MIN_FLOAT;
         inverse[ 2] = matrix[ 3]; inverse[ 8] = matrix[ 3];
         inverse[ 6] = matrix[ 4]; inverse[ 9] = matrix[ 4];
         inverse[10] = matrix[ 5] + COV_MIN_FLOAT;
         inverse[ 3] = matrix[ 6]; inverse[12] = matrix[ 6];
         inverse[ 7] = matrix[ 7]; inverse[13] = matrix[ 7];
         inverse[11] = matrix[ 8]; inverse[14] = matrix[ 8];
         inverse[15] = matrix[ 9] + COV_MIN_FLOAT;

         if(!Inverse<Float>(inverse, 4)) { throw 1; }

         // The outgoing filter is the inverse submatrix of this inverse
         // matrix.
         Float det = (matrix[5]*matrix[9]-matrix[8]*matrix[8]) / pow(2.0*M_PI,2);
         if(det > 0.0) {
            suu =  matrix[9] / det;
            svv =  matrix[5] / det;
            suv = -matrix[8] / det;
         } else {
            suu = COV_MAX_FLOAT;
            svv = COV_MAX_FLOAT;
            suv = 0.0;
         }
      }

      /*  Compute the angular extent of the angular component of the equivalent
       *  spatial filter to the covariance matrix. The extent is provided as
       * 'Du' and 'Dv', the main axis of the filter's fooprint.
       *
       * 'Du' and 'Dv' are express using the first two components:
       *        Du = [u, v, 0]
       *  as they represent direction in the local frame x,y.
       *
       *  To extract Du and Dv, we do an eigen-decomposition of the
       *  covariance matrix and use the normalized eigen-vectors as the
       *  axis of the extent and the eigen-values are the squared extent
       *  of the polygonal shape.
       */
      void AngularExtent(Vector& Du, Vector& Dv) const {
         // Compute the inverse matrix. We need to add an epsilon to the
         // diagonal in order to ensure that the matrix can be inverted.
         Float inverse[16];
         inverse[ 0] = matrix[ 0] + COV_MIN_FLOAT;
         inverse[ 1] = matrix[ 1]; inverse[ 4] = matrix[ 1];
         inverse[ 5] = matrix[ 2] + COV_MIN_FLOAT;
         inverse[ 2] = matrix[ 3]; inverse[ 8] = matrix[ 3];
         inverse[ 6] = matrix[ 4]; inverse[ 9] = matrix[ 4];
         inverse[10] = matrix[ 5] + COV_MIN_FLOAT;
         inverse[ 3] = matrix[ 6]; inverse[12] = matrix[ 6];
         inverse[ 7] = matrix[ 7]; inverse[13] = matrix[ 7];
         inverse[11] = matrix[ 8]; inverse[14] = matrix[ 8];
         inverse[15] = matrix[ 9] + COV_MIN_FLOAT;

         if(!Inverse<Float>(inverse, 4)) { throw 1; }

          // T = trace and D = det of the spatial submatrix
          Float T = inverse[10]+inverse[15];
          Float D = inverse[10]*inverse[15] - inverse[11]*inverse[11];

          // Solve the 2nd order polynomial roots of p(l) = l^2 - l T + D.
          // This gives us the eigen values.
          Float d  = 0.25*T*T - D;
          if(d < 0.0) { throw 1; } // No solution exists
          Float l1 = 0.5*T + sqrt(d);
          Float l2 = 0.5*T - sqrt(d);

          if(abs(inverse[11]) > COV_MIN_FLOAT) {
            Du.x = l1 - inverse[15];
            Du.y = inverse[11];
            Du.z = 0.0;
            Dv.x = l2 - inverse[15];
            Dv.y = inverse[11];
            Dv.z = 0.0;

            Du.Normalize();
            Dv.Normalize();

            Du = sqrt(l1)/(2.0*M_PI) * Du;
            Dv = sqrt(l2)/(2.0*M_PI) * Dv;
          } else {
            Du.x = sqrt(inverse[10])/(2.0*M_PI);
            Du.y = 0.0;
            Du.z = 0.0;
            Dv.x = 0.0;
            Dv.y = sqrt(inverse[15])/(2.0*M_PI);
            Dv.z = 0.0;
          }
      }

      /* Compute the volume (in frequency domain) spanned by the matrix.
       * Note: this volume should always be positive.
       */
      Float Volume() const {

         Float fmatrix[16];
         fmatrix[ 0] = matrix[ 0] + COV_MIN_FLOAT;
         fmatrix[ 1] = matrix[ 1]; fmatrix[ 4] = matrix[ 1];
         fmatrix[ 5] = matrix[ 2] + COV_MIN_FLOAT;
         fmatrix[ 2] = matrix[ 3]; fmatrix[ 8] = matrix[ 3];
         fmatrix[ 6] = matrix[ 4]; fmatrix[ 9] = matrix[ 4];
         fmatrix[10] = matrix[ 5] + COV_MIN_FLOAT;
         fmatrix[ 3] = matrix[ 6]; fmatrix[12] = matrix[ 6];
         fmatrix[ 7] = matrix[ 7]; fmatrix[13] = matrix[ 7];
         fmatrix[11] = matrix[ 8]; fmatrix[14] = matrix[ 8];
         fmatrix[15] = matrix[ 9] + COV_MIN_FLOAT;

         return Determinant<Float>(fmatrix, 4);
      }

      /////////////////////
      //   Constructors  //
      /////////////////////

      Covariance4D() {
         matrix = { 0.0f,
                    0.0f, 0.0f,
                    0.0f, 0.0f, 0.0f,
                    0.0f, 0.0f, 0.0f, 0.0f};
      }
      Covariance4D(Float sxx, Float syy, Float suu, Float svv) {
         matrix = {  sxx,
                    0.0f,  syy,
                    0.0f, 0.0f,  suu,
                    0.0f, 0.0f, 0.0f,  svv};
      }
      Covariance4D(std::array<Float, 10> matrix, const Vector& z) :
         matrix(matrix), z(z) {
         Vector::Frame(z, x, y);
      }
      Covariance4D(std::array<Float, 10> matrix,
                   const Vector& x,
                   const Vector& y,
                   const Vector& z) :
         matrix(matrix), x(x), y(y), z(z) {}
   };
}
