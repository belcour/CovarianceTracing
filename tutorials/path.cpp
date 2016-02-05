#include <cstdlib>
#include <cstdio>
#include <random>
#include <utility>

std::default_random_engine gen;
std::uniform_real_distribution<double> dist(0,1);

#include "common.hpp"

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

   Cov cov({ 0.0, 0.0, 0.0, 0.0, 0.0, 1.0E5, 0.0, 0.0, 0.0, 1.0E5 },
           u, v, w);

   cov.Curvature(k, k);
   cov.Symmetry();
   cov.Reflection(10.0, 10.0);
   cov.Curvature(-k, -k);
   cov.InverseProjection(-r.d);
   cov.Travel(t);
   
   return RadCov(f, cov);
   //return Vector(cov.matrix[0], cov.matrix[1], cov.matrix[2]);

   /*
   if (++depth>5) { return obj.ke; }

   double r1=2*M_PI*dist(gen), r2=dist(gen), r2s=sqrt(r2);
   Vector w = nl;
   Vector u = Vector::Cross((fabs(w.x) > .1 ? Vector(0,1,0) : Vector(1,0,0)), w).Normalize(),
   Vector v = Vector::Cross(w, u);
   Vector d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).Normalize();
   return obj.ke + f.Multiply(radiance(Ray(x,d),depth));
   */
 }

inline int toInt(double x){ return int(pow(Clamp(x),1/2.2)*255+.5); }

int main(int argc, char** argv){
   int w=1024, h=768, samps = argc==2 ? atoi(argv[1])/4 : 1; // # samples
   Ray cam(Vector(50,52,295.6), Vector(0,-0.042612,-1).Normalize()); // cam pos, dir
   Vector cx=Vector(w*.5135/h), cy=(Vector::Cross(cx, cam.d)).Normalize()*.5135, r, *c=new Vector[w*h];
   #pragma omp parallel for schedule(dynamic, 1) private(r)
   for (int y=0; y<h; y++){                       // Loop over image rows
      fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1));
      for (unsigned short x=0; x<w; x++) {  // Loop cols
         for (int sy=0, i=(h-y-1)*w+x; sy<2; sy++) {     // 2x2 subpixel rows
            for (int sx=0; sx<2; sx++, r=Vector()){        // 2x2 subpixel cols
               for (int s=0; s<samps; s++){
                  double r1=2*dist(gen), dx=r1<1 ? sqrt(r1)-1: 1-sqrt(2-r1);
                  double r2=2*dist(gen), dy=r2<1 ? sqrt(r2)-1: 1-sqrt(2-r2);
                  Vector d = cx*( ( (sx+.5 + dx)/2 + x)/w - .5) +
                             cy*( ( (sy+.5 + dy)/2 + y)/h - .5) + cam.d;

                  // Evaluate the Covariance and Radiance at the pixel location
                  auto radcov = radiance(Ray(cam.o+d*140,d.Normalize()),0);
                  auto cov = radcov.second;
                  
                  // Orient the covariance
                  float cr, sr;
                  cr = Vector::Dot(cov.x, cx);
                  sr = Vector::Dot(cov.x, cy);
                  cov.Rotate(cr, -sr);

                  // What do you want to see?
                  Vector rgb;

                  //  1) The angular part of the covariance
                  rgb = Vector(std::fabs(cov.matrix[5]),
                                          std::fabs(cov.matrix[8]),
                                          std::fabs(cov.matrix[9]));
                  //  2) The spatial part of the covariance
                  rgb = Vector(std::fabs(cov.matrix[0]),
                                          std::fabs(cov.matrix[1]),
                                          std::fabs(cov.matrix[2]));
                  //  3) Density
                  float den;
                  den = cov.matrix[0]*cov.matrix[2]-pow(cov.matrix[1], 2);
                  den = sqrt(den);
                  rgb = Vector(den, den, den);

                  r = r + rgb*(1.f/samps);
               } // Camera rays are pushed ^^^^^ forward to start in interior
               c[i] = c[i] + Vector(Clamp(r.x),Clamp(r.y),Clamp(r.z))*.25;
            }
         }
      }
   }
   FILE *f = fopen("image.ppm", "w");         // Write image to PPM file.
   fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
   for (int i=0; i<w*h; i++) {
      fprintf(f,"%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
   }
}

