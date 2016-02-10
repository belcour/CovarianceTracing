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


Sphere spheres[] = {
   Sphere(Vector( 1e5+1,40.8,81.6),  1e5,  Vector(),Vector(.75,.25,.25)),//Left
   Sphere(Vector(-1e5+99,40.8,81.6), 1e5,  Vector(),Vector(.25,.25,.75)),//Rght
   Sphere(Vector(50,40.8, 1e5),      1e5,  Vector(),Vector(.75,.75,.75)),//Back
   Sphere(Vector(50,40.8,-1e5+170),  1e5,  Vector(),Vector()           ),//Frnt
   Sphere(Vector(50, 1e5, 81.6),     1e5,  Vector(),Vector(.75,.75,.75)),//Botm
   Sphere(Vector(50,-1e5+81.6,81.6), 1e5,  Vector(),Vector(.75,.75,.75)),//Top
   Sphere(Vector(27,16.5,47),        16.5, Vector(),Vector(1,1,1)*.999),//Mirr
   Sphere(Vector(73,16.5,78),        16.5, Vector(),Vector(1,1,1)*.999),//Glas
   Sphere(Vector(50,681.6-.27,81.6), 600,  Vector(12,12,12),  Vector()) //Lite
};


inline bool intersect(const Ray &r, double &t, int &id){
   double n=sizeof(spheres)/sizeof(Sphere), d, inf=t=1e20;
   for(int i=int(n);i--;) if((d=spheres[i].Intersect(r))&&d<t){t=d;id=i;}
   return t<inf;
}

RadCov radiance(const Ray &r, int depth){
   double t;                               // distance to intersection
   int id=0;                               // id of intersected object
   if (!intersect(r, t, id)) return RadCov(Vector(), Cov()); // if miss, return black
   const Sphere &obj = spheres[id];        // the hit object
   Vector x  = r.o+r.d*t,
          n  = (x-obj.c).Normalize(),
          nl = Vector::Dot(n,r.d) < 0 ? n:n*-1,
          f  = obj.kd;

   const float k = 1.f/spheres[id].r;

   /* Local Frame at the surface of the object */
   Vector w = nl;
   Vector u = Vector::Cross((fabs(w.x) > .1 ? Vector(0,1,0) : Vector(1,0,0)), w).Normalize();
   Vector v = Vector::Cross(w, u);

   // If the object is a source, return the its covariance. A source has a
   // constant angular emission but has bounded spatial extent. Since there
   // is no way here to infer the extent of the source, its frequency is
   // defined as 1.0E5.
   if(!obj.ke.IsNull()) {
      Cov cov({ 1.0E2, 0.0, 1.0E2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, u, v, w);
      cov.InverseProjection(-r.d);
      cov.Travel(t);
      return RadCov(obj.ke, cov) ;

   // Terminate the recursion after a finite number of call. Since this
   // implementation is recursive and passing covariance objects, the
   // number of max bounces is deterministic.
   } else if(depth>1) {
      Cov cov({ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, u, v, w);
      cov.InverseProjection(-r.d);
      return RadCov(Vector(), cov) ;

   // Main covariance computation. First this code generate a new direction
   // and query the covariance+radiance in that direction. Then, it computes
   // the covariance after the reflection/refraction.
   } else {
      /* Sampling a new direction + recursive call */
      double r1=2*M_PI*dist(gen), r2=dist(gen), r2s=sqrt(r2);
      Vector w = nl;
      Vector u = Vector::Cross((fabs(w.x) > .1 ? Vector(0,1,0) : Vector(1,0,0)), w).Normalize();
      Vector v = Vector::Cross(w, u);
      Vector d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).Normalize();
      const RadCov radcov = radiance(Ray(x,d), depth+1);

      /* Covariance computation */
      Covariance4D<Vector> cov = radcov.second;
      cov.Projection(nl);
      cov.Curvature(k, k);
      cov.Cosine(1.0f);
      cov.Symmetry();
      cov.Reflection(1.0, 1.0);
      cov.Curvature(-k, -k);
      cov.InverseProjection(-r.d);
      cov.Travel(t);
      return RadCov(f.Multiply(radcov.first), cov);
   }
}

#include <xmmintrin.h>

int main(int argc, char** argv){
   int w=1024, h=768, samps = argc==2 ? atoi(argv[1])/4 : 1; // # samples
   Ray cam(Vector(50,52,295.6), Vector(0,-0.042612,-1).Normalize()); // cam pos, dir
   Vector  cx  = Vector(w*.5135/h);
   Vector  cy  = Vector::Cross(cx, cam.d).Normalize()*.5135;
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
                  Vector d = cx*( ( (sx+.5 + dx)/2 + x)/w - .5) +
                             cy*( ( (sy+.5 + dy)/2 + y)/h - .5) + cam.d;
                  d.Normalize();

                  // Covariance tracing requires to know the pixel frame in order to
                  // align the orientation of the covariance matrix with respect to
                  // the image plane. (cx, cy, d) is not a proper frame and we need
                  // to correct it.
                  const Vector px = (cx - Vector::Dot(d, cx)*d).Normalize(),
                               py = (cy - Vector::Dot(d, cy)*d).Normalize();
                  const float scaleX = Vector::Norm(cx) / float(w),
                              scaleY = Vector::Norm(cy) / float(h);

                  // Evaluate the Covariance and Radiance at the pixel location
                  auto radcov = radiance(Ray(cam.o+d*140,d.Normalize()),0);
                  auto rad = radcov.first;
                  auto cov = radcov.second;

                  // Orient the covariance and scale it to be in pixel^{-2} and not
                  // in meter^{-2} or rad^{-2}.
                  float cr, sr;
                  cr = Vector::Dot(cov.x, px);
                  sr = Vector::Dot(cov.x, py);
                  cov.Rotate(cr, sr);
                  cov.ScaleU(scaleX);
                  cov.ScaleV(scaleY);

                  _cov.Add(cov, Vector::Norm(_r), Vector::Norm(rad));
                  _r = (_r*float(s) + rad)*(1.f/(s+1.f));
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
               float den;
               den = _cov.matrix[0]*_cov.matrix[2]-pow(_cov.matrix[1], 2);
               den = sqrt(fmax(den, 0.0));
               c = Vector(den, den, den);

               //  4) Radiance
               //c = _r;

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

