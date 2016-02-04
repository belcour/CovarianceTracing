#pragma once

// STL includes
#include <array>
#include <cmath>

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
   template<class Vector>
   struct Covariance4D {

      std::array<float, 10> matrix;
      Vector x, y, z;


      ////////////////////////
      //  Atomic operators  //
      ////////////////////////

      /* Travel operator
       *
       * 'd' distance of travel along the central ray
       */
      inline void Travel(float d) {
         ShearAngleSpace(d, d);
      }

      /* Curvature operator
       *
       * 'kx' curvature along the X direction
       * 'ky' curvature along the Y direction
       */
      inline void Curvature(float kx, float ky) {
         ShearSpaceAngle(kx, ky);
      }

      /* Cosine operator
       *
       * 'wz' The incident direction's elevation in the local frame
       */
      inline void Cosine(float wz) {
         const float theta = acos(wz);
         const float dist  = std::abs(0.5*M_PI-theta);
         const float frequ = 2.0 / M_PI;
         const float freqv = 1.0 / fmax(dist, 1.0E-10);
         matrix[5] += frequ*frequ;
         matrix[9] += freqv*freqv;
      }

      /* Reflection operator
       *
       * 'suu' the covariance of the BRDF along the X axis.
       * 'svv' the covariance of the BRDF along the Y axis.
       */
      inline void Reflection(float suu, float svv) {
         ProductUV(suu, svv);
      }

      /* Reflection operator
       */
      inline void Refraction() {
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
         if(cx == 0.0f) { return; }

         // Rotate the Frame to be aligned with plane.
         const float alpha = atan2(cx, cy);
         const float c = cos(alpha), s = -sin(alpha);
         Rotate(c, s);

         // Scale the componnent that project by the cosine of the ray direction
         // and the normal.
         const float cosine = std::abs(Vector::Dot(z, n));
         ScaleY(cosine);

         // Update direction vectors.
         const auto temp = c*x + s*y;
         z = n;
         y = c*y - s*x;
         x = temp;
      }

      //! Perform the projection of a lightfield defined on a surface to an
      //! outgoing direction. This function assumes that the lightfield main
      //! vector is the surface normal and that the outgoing vector is in the
      //! same direction.
      //!
      //! \param out_d the outgoing direction.
      inline void InverseProjection(const Vector& d) {

         const auto cx = Vector::Dot(x, d);
         const auto cy = Vector::Dot(y, d);
         if(cx == 0.0f) { return; }

         // Rotate the Frame to be aligned with plane.
         const float alpha = atan2(cx, cy);
         const float c = cos(alpha), s = -sin(alpha);
         Rotate(c, s); // Rotate of -alpha

         // Scale the componnent that project by the inverse cosine of the ray
         // direction and the normal.
         const float cosine = std::abs(Vector::Dot(z, d));
         ScaleY(1.0/cosine);

         // Update direction vectors.
         const auto temp = c*x + s*y;
         z = d;
         y = c*y - s*x;
         x = temp;
      }


      ////////////////////////
      //   Matrix scaling   //
      ////////////////////////

      inline void ScaleX(float alpha) {
         matrix[0] *= alpha*alpha;
         matrix[1] *= alpha;
         matrix[3] *= alpha;
         matrix[6] *= alpha;
      }

      inline void ScaleY(float alpha) {
         matrix[1] *= alpha;
         matrix[2] *= alpha*alpha;
         matrix[4] *= alpha;
         matrix[7] *= alpha;
      }

      inline void ScaleU(float alpha) {
         matrix[3] *= alpha;
         matrix[4] *= alpha;
         matrix[5] *= alpha*alpha;
         matrix[8] *= alpha;
      }

      inline void ScaleV(float alpha) {
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
      inline void ShearSpaceAngle(float cx, float cy) {
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
      inline void ShearAngleSpace(float cu, float cv) {
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
      inline void Rotate(float alpha) {
         const float c = cos(alpha);
         const float s = sin(alpha);
         Rotate(c, s);
      }

      // 'c' the cosine of the rotation angle
      // 's' the sine of the rotation angle
      inline void Rotate(float c, float s) {
         const float cs = c*s;
         const float c2 = c*c;
         const float s2 = s*s;

         const float cov_xx = matrix[0];
         const float cov_xy = matrix[1];
         const float cov_yy = matrix[2];
         const float cov_xu = matrix[3];
         const float cov_yu = matrix[4];
         const float cov_uu = matrix[5];
         const float cov_xv = matrix[6];
         const float cov_yv = matrix[7];
         const float cov_uv = matrix[8];
         const float cov_vv = matrix[9];

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
      //  Symmetry of space //
      ////////////////////////

      inline void Symmetry() {
         matrix[3] = -matrix[3];
         matrix[4] = -matrix[4];
         matrix[6] = -matrix[6];
         matrix[7] = -matrix[7];
      }


      ////////////////////////
      // Product of signals //
      ////////////////////////

      // Evaluate the covariance matrix of the product of the local lightfield
      // and a angularly varying only signal (like a BSDF).
      //
      // 'su' is the sigma u of the inverse angular signal's covariance matrix.
      // 'sv' the sigma v of the inverse angular signal's covariance matrix.
      inline void ProductUV(float su, float sv) {

         if(su == sv == 0.0) {
            return;
         }

         const float cov_xx = matrix[0];
         const float cov_xy = matrix[1];
         const float cov_yy = matrix[2];
         const float cov_xu = matrix[3];
         const float cov_yu = matrix[4];
         const float cov_uu = matrix[5];
         const float cov_xv = matrix[6];
         const float cov_yv = matrix[7];
         const float cov_uv = matrix[8];
         const float cov_vv = matrix[9];

         // The following is an application of the Woodbury matrix identity for
         // a rank 2 C matrix. See:
         //         http://en.wikipedia.org/wiki/Woodbury_matrix_identity
         const float sig_u = cov_uu+su, sig_v = cov_vv+sv;
         const float og = (sig_u*sig_v - cov_uv*cov_uv);
         const float  g = 1.0f / og;

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
      Covariance4D(float sxx, float syy, float suu, float svv) {
         matrix = {  sxx,
                    0.0f,  syy,
                    0.0f, 0.0f,  suu,
                    0.0f, 0.0f, 0.0f,  svv};
      }
      Covariance4D(std::array<float, 10> matrix,
                   const Vector& x,
                   const Vector& y,
                   const Vector& z) :
         matrix(matrix), x(x), y(y), z(z) {}
   };
}
