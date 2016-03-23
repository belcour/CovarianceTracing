// STL includes
#include <cstdlib>
#include <cstdio>
#include <random>
#include <utility>
#include <iostream>
#include <sstream>
#include <thread>

// Local includes
#include "common.hpp"
#include "tutorial2.hpp"

std::stringstream sout;

PosCov CovarianceFilter(const std::vector<Sphere>& spheres, const Ray &r,
                        const Cov4D& cov, int depth, int maxdepth,
                        std::stringstream& out) {
   double t;
   int id=0;
   if (!Intersect(spheres, r, t, id)) {
     return PosCov(Vector(), Cov4D());
   }
   const Sphere&   obj = spheres[id];
   const Material& mat = obj.mat;
   Vector x  = r.o+r.d*t,
          n  = (x-obj.c).Normalize(),
          nl = Vector::Dot(n,r.d) < 0 ? n:n*-1;
   const double k = 1.f/spheres[id].r;

   // Update the covariance with travel and project it onto the tangent plane
   // of the hit object.
   Cov4D cov2 = cov;
   cov2.Travel(t);
   out << "After travel of " << t << " meters" << std::endl;
   out << cov2 << std::endl;
   out << "Volume = " << cov2.Volume() << std::endl;
   out << std::endl;

   out << "After projection, cos=" << Vector::Dot(n, cov2.z) << std::endl;
   cov2.Projection(n);
   out << cov2 << std::endl;
   out << "Volume = " << cov2.Volume() << std::endl;
   out << std::endl;

   // if the max depth is reached
   if(depth >= maxdepth) {
      return PosCov(x, cov2);
   } else {
      // Sample a new direction
      auto wi = -r.d;
      auto wr = 2*Vector::Dot(wi, nl)*nl - wi;
      auto r2 = Ray(x, wr);

      cov2.Curvature(k, k);
      out << "After curvature" << std::endl;
      out << cov2 << std::endl;
      out << "Volume = " << cov2.Volume() << std::endl;
      out << std::endl;

      cov2.Cosine(1.0f);
      out << "After cosine multiplication" << std::endl;
      out << cov2 << std::endl;

      cov2.Symmetry();
      out << "After symmetry" << std::endl;
      out << cov2 << std::endl;
      out << "Volume = " << cov2.Volume() << std::endl;
      out << std::endl;

      const double rho = mat.exponent / (4*M_PI*M_PI);
      cov2.Reflection(rho, rho);
      out << "After BRDF convolution, sigma=" << rho << std::endl;
      out << cov2 << std::endl;
      out << "Volume = " << cov2.Volume() << std::endl;
      out << std::endl;

      cov2.Curvature(-k, -k);
      out << "After inverse curvature, k=" << -k << std::endl;
      out << cov2 << std::endl;
      out << "Volume = " << cov2.Volume() << std::endl;
      out << std::endl;

      cov2.InverseProjection(wr);
      out << "After inverse projection, cos=" << Vector::Dot(n, wr) << std::endl;
      out << cov2 << std::endl;
      out << "Volume = " << cov2.Volume() << std::endl;
      out << std::endl;

      return CovarianceFilter(spheres, r2, cov2, depth+1, maxdepth, out);
   }
}

// Should the covariance image display the Gaussian filter or the equivalent
// ray differential's footprint?
bool useCovFilter = true;

void CovarianceTexture(int x, int y) {
   // Generate a covariance matrix at the sampling position
   const auto t = (cx*((x+0.5)/double(width) - .5) + cy*((y+0.5)/double(height) - .5) + cam.d).Normalize();
   //*/
   const auto pixelCov = Cov4D({ 1.0E+5, 0.0, 1.0E+5, 0.0, 0.0, 1.0E+5, 0.0, 0.0, 0.0, 1.0E+5 }, t);
   /*/
   const auto pixelCov = Cov4D({ 1.0E-5, 0.0, 1.0E-5, 0.0, 0.0, 1.0E-5, 0.0, 0.0, 0.0, 1.0E-5 }, t);
   //*/
   sout.str("");
   const auto surfCov  = CovarianceFilter(spheres, Ray(cam.o, t), pixelCov, 0, 1, sout);
   sout << surfCov.second << std::endl;
   sout << "Volume = " << surfCov.second.Volume() << std::endl;
   sout << std::endl;

   double sxx = 0, syy = 0, sxy = 0;
   Vector Dx, Dy;
   try {
      surfCov.second.SpatialFilter(sxx, sxy, syy);
      surfCov.second.SpatialExtent(Dx, Dy);
      sout << "Spatial filter = [" << sxx << "," << sxy << "; " << sxy << ", " << syy << "]"<< std::endl;
      sout << "Extent = " << Dx << ", " << Dy << std::endl;
      sout << "|Dx| = " << Vector::Norm(Dx) << ", |Dy| = " << Vector::Norm(Dy) << std::endl;
   } catch (...) {
      std::cout << "Error: incorrect spatial filter" << std::endl;
      sout << surfCov.second << std::endl;
      return;
   }

   // Loop over the rows and columns of the image and evaluate radiance and
   // covariance per pixel using Monte-Carlo.
   #pragma omp parallel for schedule(dynamic, 1)
   for (int y=0; y<height; y++){
      for (int x=0; x<width; x++) {
         // Pixel index
         int i=(width-x-1)*height+y;

         // Generate the pixel direction
         Vector d = cx*( ( 0.5 + x)/width - .5) +
                    cy*( ( 0.5 + y)/height - .5) + cam.d;
         d.Normalize();

         Ray ray(cam.o, d);
         double t; int id;
         if(!Intersect(spheres, ray, t, id)){ continue; }
         Vector hitp = ray.o + t*ray.d;

         // Evaluate the covariance
         const Vector dx  = surfCov.first - hitp;
         const Vector dU = Vector(Vector::Dot(dx, surfCov.second.x), Vector::Dot(dx, surfCov.second.y), Vector::Dot(dx, surfCov.second.z));
         if(useCovFilter) {
            double bf = dU.x*dU.x*sxx + dU.y*dU.y*syy + 2*dU.x*dU.y*sxy;
            cov_img[i] = exp(-10.0*dU.z*dU.z) * exp(- 0.5* bf);
         } else {

            const double du  = Vector::Dot(dU, Dx) / Vector::Norm(Dx);
            const double dv  = Vector::Dot(dU, Dy) / Vector::Norm(Dy);
            cov_img[i] = (abs(du) < Vector::Norm(Dx) &&
                          abs(dv) < Vector::Norm(Dy)) ? exp(-10.0*dU.z*dU.z) : 0.0;
         }
      }
   }
}

void PrintHelp() {
   std::cout << "Covariance Tracing tutorial 2" << std::endl;
   std::cout << "Usage: ./CovFiltering" << std::endl;
}

int main(int argc, char** argv) {

   PrintHelp();

   // Variables
   int spp=1000, npasses=10;
   int x=410, y=175;

   //generateReference = true;

   // Display the covariance filter
   if(generateCovariance) {
      CovarianceTexture(x, y);
   }

   // Progressively display the background image and the brute force evaluation
   // of the indirect pixel filter.
   for(int i=0; i<std::max(spp, npasses); ++i) {
      if(generateBackground && i<spp)
         RadianceTexture();

      if(generateReference && i<npasses)
         BruteForceTexture(x, y);

      if(i % 10 == 9) { ExportImage(x, y); }
   }

   // Display result
   ExportImage(x, y);

   // Clean memory
   if(bcg_img) { delete[] bcg_img; }
   if(cov_img) { delete[] cov_img; }
   if(ref_img) { delete[] ref_img; }
   return EXIT_SUCCESS;
}
