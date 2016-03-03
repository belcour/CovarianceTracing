// STL includes
#include <cstdlib>
#include <cstdio>
#include <random>
#include <utility>

std::default_random_engine gen;
std::uniform_real_distribution<double> dist(0,1);

// Local includes
#include "common.hpp"

// Covariance Tracing includes
#include <Covariance/Covariance4D.hpp>
using namespace Covariance;
using Cov    = Covariance4D<Vector>;
using RadCov = std::pair<Vector, Cov>;


Material phong(Vector(), Vector(0,0,0), Vector(1,1,1)*.999, 100.0);

std::vector<Sphere> spheres = {
   Sphere(Vector( 1e5+1,40.8,81.6),  1e5,  Vector(),Vector(.75,.25,.25)),//Left
   Sphere(Vector(-1e5+99,40.8,81.6), 1e5,  Vector(),Vector(.25,.25,.75)),//Rght
   Sphere(Vector(50,40.8, 1e5),      1e5,  Vector(),Vector(.75,.75,.75)),//Back
   Sphere(Vector(50,40.8,-1e5+170),  1e5,  Vector(),Vector()           ),//Frnt
   Sphere(Vector(50, 1e5, 81.6),     1e5,  Vector(),Vector(.75,.75,.75)),//Botm
   Sphere(Vector(50,-1e5+81.6,81.6), 1e5,  Vector(),Vector(.75,.75,.75)),//Top
   Sphere(Vector(27,16.5,47),        16.5, phong),//Mirr
   Sphere(Vector(73,16.5,78),        16.5, Vector(),Vector(1,1,1)*.999),//Glas
   Sphere(Vector(50,681.6-.27,81.6), 600,  Vector(12,12,12),  Vector()) //Lite
};

RadCov radiance(const Ray &r, int depth, int maxdepth=1){
   double t;                               // distance to intersection
   int id=0;                               // id of intersected object
   if (!Intersect(spheres, r, t, id)) return RadCov(Vector(), Cov()); // if miss, return black
   const Sphere&   obj = spheres[id];      // the hit object
   const Material& mat = obj.mat;          // Its material
   
   Vector x  = r.o+r.d*t;
   Vector n  = (x-obj.c).Normalize();
   Vector nl = (Vector::Dot(n, r.d) < 0.f) ? n : (-1.f)*n;

   const double k = 1.f/spheres[id].r;

   /* Local Frame at the surface of the object */
   Vector w = nl;
   Vector u = Vector::Cross((fabs(w.x) > .1 ? Vector(0,1,0) : Vector(1,0,0)), w).Normalize();
   Vector v = Vector::Cross(w, u);

   // If the object is a source, return the its covariance. A source has a
   // constant angular emission but has bounded spatial extent. Since there
   // is no way here to infer the extent of the source, its frequency is
   // defined as 1.0E5.
   if(!mat.ke.IsNull()) {
      Cov cov({ 1.0E2, 0.0, 1.0E2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, u, v, w);
      cov.InverseProjection(-r.d);
      cov.Travel(t);
      return RadCov(mat.ke, cov) ;

   // Terminate the recursion after a finite number of call. Since this
   // implementation is recursive and passing covariance objects, the
   // number of max bounces is deterministic.
   } else if(depth > maxdepth) {
      Cov cov({ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, u, v, w);
      cov.InverseProjection(-r.d);
      return RadCov(Vector(), cov) ;

   // Main covariance computation. First this code generate a new direction
   // and query the covariance+radiance in that direction. Then, it computes
   // the covariance after the reflection/refraction.
   } else {
      /* Sampling a new direction + recursive call */
      double pdf = 0.f;
      const auto e  = Vector(dist(gen), dist(gen), dist(gen));
      const auto wo = -r.d;
      const auto wi = mat.Sample(wo, nl, e, pdf);
      if(Vector::Dot(wo, nl) <= 0.f || pdf <= 0.f) {
         Cov cov({ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, u, v, w);
         cov.InverseProjection(wo);
      	 return RadCov(Vector((pdf <= 0.f) ? 1.0 : 0.0,0.0,0.0), cov) ;
      }
      auto f = Vector::Dot(wi, nl)*mat.Reflectance(wi, wo, nl);
      const RadCov radcov = radiance(Ray(x, wi), depth+1);

      /* Covariance computation */
      Covariance4D<Vector> cov = radcov.second;
      cov.Projection(nl);
      cov.Curvature(k, k);
      cov.Cosine(1.0f);
      cov.Symmetry();
      const double rho = mat.exponent / (4*M_PI*M_PI);
      cov.Reflection(rho, rho); // TODO correct the formula
      cov.Curvature(-k, -k);
      cov.InverseProjection(-r.d);
      cov.Travel(t);
      return RadCov((1.f/pdf) * f.Multiply(radcov.first), cov);
   }
}

#include <xmmintrin.h>

int main(int argc, char** argv){
   int w=512, h=512, samps = argc==2 ? atoi(argv[1])/4 : 1; // # samples
   Ray cam(Vector(50,52,295.6), Vector(0,-0.042612,-1).Normalize()); // cam pos, dir
   cam.o = cam.o + 140.0*cam.d;
   double fovx = 1.2; // 0.5135;
   double fovy = 1.2; // 0.5135;
   Vector  cx  = Vector(fovx);
   Vector  cy  = Vector::Cross(cx, cam.d).Normalize()*fovy;
   Vector ncx  = cx; ncx.Normalize();
   Vector ncy  = cy; ncy.Normalize();
   Vector* img = new Vector[w*h];

   _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

   // Loop over the rows and columns of the image and evaluate radiance and
   // covariance per pixel using Monte-Carlo.
   #pragma omp parallel for schedule(dynamic, 1) private(gen)
   for (int y=0; y<h; y++){
      fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1));
      for (unsigned short x=0; x<w; x++) {

         // Sub pixel sampling
         for (int sy=0, i=(h-y-1)*w+x; sy<2; sy++) {
            for (int sx=0; sx<2; sx++){

               Vector _r;
               Covariance4D<Vector> _cov;

               for (int s=0; s<samps; s++){

                  // Generate a sub-pixel random position to perform super
                  // sampling.
                  double r1=2*dist(gen), dx=r1<1 ? sqrt(r1)-1: 1-sqrt(2-r1);
                  double r2=2*dist(gen), dy=r2<1 ? sqrt(r2)-1: 1-sqrt(2-r2);

                  // Generate the pixel direction
                  Vector d = ncx*fovx*(( (sx+.5 + dx)/2 + x)/w - .5) +
                             ncy*fovy*(( (sy+.5 + dy)/2 + y)/h - .5) + cam.d;
                  d.Normalize();

                  // Covariance tracing requires to know the pixel frame in order to
                  // align the orientation of the covariance matrix with respect to
                  // the image plane. (cx, cy, d) is not a proper frame and we need
                  // to correct it.
                  const Vector px = (ncx - Vector::Dot(d, ncx)*d).Normalize(),
                               py = (ncy - Vector::Dot(d, ncy)*d).Normalize();
                  const double scaleX = Vector::Norm(ncx) / double(w),
                              scaleY = Vector::Norm(ncy) / double(h);

                  // Evaluate the Covariance and Radiance at the pixel location
                  auto radcov = radiance(Ray(cam.o, d),0);
                  auto rad = radcov.first;
                  auto cov = radcov.second;

                  // Orient the covariance and scale it to be in pixel^{-2} and not
                  // in meter^{-2} or rad^{-2}.
                  double cr, sr;
                  cr = Vector::Dot(cov.x, px);
                  sr = Vector::Dot(cov.x, py);
                  cov.Rotate(cr, sr);
                  cov.ScaleU(scaleX);
                  cov.ScaleV(scaleY);

                  _cov.Add(cov, Vector::Norm(_r), Vector::Norm(rad));
                  _r = (_r*double(s) + rad)*(1.f/(s+1.f));
               }

               // What do you want to see? [Un]comment some of those line to
               // output a different part of aspect of frequency analysis.
               Vector c;

               //  1) The angular part of the covariance
               //c = Vector(std::fabs(_cov.matrix[5]),
               //           std::fabs(_cov.matrix[8]),
               //           std::fabs(_cov.matrix[9]));

               //  2) The spatial part of the covariance
               //c = Vector(std::fabs(_cov.matrix[0]),
               //           std::fabs(_cov.matrix[1]),
               //           std::fabs(_cov.matrix[2]));
               //c = Vector(std::fabs(cr), 0.0, std::fabs(sr));

               //  3) Predicted sampling density. This is what Belcour et al.
               //  [2013] used to generate the image space adaptive sampling.
               double den;
               den = _cov.matrix[0]*_cov.matrix[2]-pow(_cov.matrix[1], 2);
               den = sqrt(fmax(den, 0.0));
               c = Vector(den, den, den);

               //  4) Radiance
               c = _r;

               img[i] = img[i] + c*.25;
            }
         }
      }
   }

   // Output image
   const auto ret = SaveEXR(img, w, h, "image.exr");

   delete[] img;
   return ret;
}

