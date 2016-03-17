#pragma once

// Include STL
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <utility>

double Clamp(double x) {
   return x<0 ? 0 : x>1 ? 1 : x;
}

struct Vector {
   double x, y, z;
   Vector(double x=0, double y=0, double z=0) : x(x), y(y), z(z) {}
   static double Dot(const Vector& w1, const Vector& w2) {
      return w1.x*w2.x + w1.y*w2.y + w1.z*w2.z;
   }
   static Vector Cross(const Vector& u, const Vector& v) {
      Vector r;
      r.x = u.y*v.z - u.z*v.y;
      r.y = u.z*v.x - u.x*v.z;
      r.z = u.x*v.y - u.y*v.x;
      return r;
   }
   Vector Normalize() {
      double norm = sqrt(Dot(*this, *this));
      x /= norm;
      y /= norm;
      z /= norm;
      return *this;
   }
   static double Norm(const Vector& a) {
       return sqrt(Vector::Dot(a, a));
   }
   bool IsNull() const {
    return x==0.0f && y==0.0f && z==0.0f;
   }
   Vector Reflect(const Vector& w, const Vector& n) {
      Vector r;
      const double dot = Dot(w, n);
      r = 2*dot*n - w;
      return r;
   }
   Vector Multiply(const Vector& m) {
      Vector r;
      r.x = this->x * m.x;
      r.y = this->y * m.y;
      r.z = this->z * m.z;
      return r;
   }
   friend Vector operator*(double a, const Vector& w) {
      Vector v;
      v.x = a*w.x;
      v.y = a*w.y;
      v.z = a*w.z;
      return v;
   }
   friend Vector operator*(const Vector& w, double a) {
      Vector v;
      v.x = a*w.x;
      v.y = a*w.y;
      v.z = a*w.z;
      return v;
   }
   friend Vector operator+(const Vector& a, const Vector& w) {
      Vector v;
      v.x = a.x+w.x;
      v.y = a.y+w.y;
      v.z = a.z+w.z;
      return v;
   }
   friend Vector operator-(const Vector& a, const Vector& w) {
      Vector v;
      v.x = a.x-w.x;
      v.y = a.y-w.y;
      v.z = a.z-w.z;
      return v;
   }
   friend Vector operator-(const Vector& w) {
      Vector v;
      v.x = -w.x;
      v.y = -w.y;
      v.z = -w.z;
      return v;
   }
   friend std::ostream& operator<<(std::ostream& out, const Vector& w) {
      out << "[" << w.x << ", " << w.y << ", " << w.z << "]";
      return out;
   }
   static void Frame(const Vector& z, Vector& x, Vector& y) {
      if(std::fabs(z.x) > 0.1f) {
         x = Vector::Cross(Vector(0,1,0), z).Normalize();
      } else {
         x = Vector::Cross(Vector(1,0,0), z).Normalize();
      }
      y = Vector::Cross(z, x);
   }
};

struct Ray {
   Vector o, d;
   Ray(const Vector& o, const Vector& d) : o(o), d(d) {}
};

struct Camera {
   Vector o, d;
   Vector cx, cy;
   double  fx, fy;

   Camera(const Vector& o, const Vector& d) : o(o), d(d), fx(1.0), fy(1.0) {
      Vector::Frame(d, cx, cy);
   }

   Vector PixelToDirection(const Vector uv) const {
      return (d + (uv.x-.5)*fx*cx + (uv.y-0.5)*fy*cy).Normalize();
   }

#ifdef TODO
   Vector DirectionToPixel(const Vector dir) const {
      double u,v;
      u = Vector::Dot(dir, cx);
      v = Vector::Dot(dir, cy);
   }
#endif
};

struct Material {
   // Diffuse and specular values
   Vector ke, kd, ks;
   double exponent;

   // Constructors
   Material(const Vector& ke, const Vector& kd, const Vector& ks, double exponent)
      : ke(ke), kd(kd), ks(ks), exponent(exponent) {}

   Material(const Vector& ke, const Vector& kd)
      : ke(ke), kd(kd), ks(), exponent(0.f) {}

   // Evaluate
   Vector Emission() const {
      return ke;
   }
   inline double PhongLobe(const Vector& wr, const Vector& wo) const {
      const auto dot = Vector::Dot(wr, wo);
      const auto f   = ((exponent+1.f)/(2.f*M_PI))*pow(dot, exponent);
      return f;
   }
   Vector Reflectance(const Vector& wi, const Vector& wo, const Vector& n) const {
      if(Vector::Dot(wi, n) <= 0.0f || Vector::Dot(wo, n) <= 0.0f) {
         return Vector();
      }
      const auto wr = 2.f*Vector::Dot(wi, n)*n - wi;
      const auto fs = PhongLobe(wr, wo);
      const auto fd = 1.f / (2.f * M_PI);
      return fd*kd + fs*ks;
   }

   // Sample phong
   Vector Sample(const Vector& wi, const Vector& n,
                 const Vector& e, double& pdf) const {
      Vector u, v;
      if(Vector::Dot(wi, n) <= 0.0) {
         pdf = 0.0f;
         return Vector();
      }

      // TODO here, should scale by the selection ratio
      if(e.z * Vector::Norm(kd+ks) <= Vector::Norm(ks)) {
         const auto cosT = pow(e.x, 1.f/(exponent+1.f));
         const auto sinT = sqrt(1.f - cosT*cosT);
         const auto phi  = 2*M_PI*e.y;
         const auto wr   = 2.f*Vector::Dot(wi, n)*n - wi;
         Vector::Frame(wr, u, v);
         const auto wo = cosT*wr + sinT*(cos(phi)*u + sin(phi)*v);
         pdf = PhongLobe(wr, wo);
         return wo;
      } else {
         Vector::Frame(n, u, v);
         const auto cosT = e.x;
         const auto sinT = sqrt(1.f - cosT*cosT);
         const auto phi  = 2*M_PI*e.y;
         const auto wo = cosT*n + sinT*(cos(phi)*u + sin(phi)*v);
         pdf = 1.f / (2.f * M_PI);
         return wo;
      }
   }
};

struct Sphere {
   // Geometric information for the sphere
   Vector c;
   double r;

   // Emission and reflectance color
   Material mat;

   Sphere(const Vector& c, double r, const Material& mat)
      : c(c), r(r), mat(mat) {}
   Sphere(const Vector& c, double r, const Vector& ke, const Vector& kd)
      : c(c), r(r), mat(ke, kd) {}

   // returns distance, -1 if nohit
   double Intersect(const Ray &ray) const {
      Vector op = c-ray.o; // Solve t^2*d.d + 2*t*(o-c).d + (o-c).(o-c)-R^2 = 0
      double t, eps = 1e-4;
      double b   = Vector::Dot(op, ray.d);
      double det = b*b - Vector::Dot(op,op) + r*r;
      if (det<0) { return 0; } else { det=sqrt(det); }
      return (t=b-det)>eps ? t : ((t=b+det)>eps ? t : 0);
   }
};

/* Random number generator using the STL random and unfiorm distribution
 */
struct Random {

   Random(unsigned int seed) : _gen(seed), _dist(0.0, 1.0) {}
   Random() : _gen(clock()), _dist(0.0, 1.0) {}

   double operator()() { return _dist(_gen); }

   std::default_random_engine _gen;
   std::uniform_real_distribution<double> _dist;
};

/* 'Intersect' routine return the intersection of a ray with a vector of spheres.
 * this method return true if a sphere is intersected and false if nothing is
 * hit. If a hit is found, 't' contain the distance to the hit point and 'id' the
 * index of the hit sphere.
 */
inline bool Intersect(const std::vector<Sphere>& spheres, const Ray &r, double &t, int &id){
   double d, inf=t=1e20;
   for(int i=spheres.size(); i--;) if((d=spheres[i].Intersect(r)) && d<t) { t=d; id=i; }
   return t<inf;
}

/* 'Radiance' evaluate the RGB throughput using path tracing with no explicit
 * connection to the light. This is a very crude way to do it, but it is the simplest
 * to rapidly prototype in practice. This code assumes that light sources do not
 * scatter light.
 */
Vector Radiance(const std::vector<Sphere>& spheres,
                const Ray &r,
                Random& rng,
                int depth,
                int maxdepth=1) {
   double t;                               // distance to intersection
   int id=0;                               // id of intersected object
   if (!Intersect(spheres, r, t, id)) return Vector(); // if miss, return black
   const Sphere&   obj = spheres[id];      // the hit object
   const Material& mat = obj.mat;          // Its material

   Vector x  = r.o+r.d*t;
   Vector n  = (x-obj.c).Normalize();
   Vector nl = (Vector::Dot(n, r.d) < 0.f) ? n : (-1.f)*n;

   // If the object is a source return radiance
   if(!mat.ke.IsNull()) {
      return mat.ke;

   // Terminate the recursion after a finite number of call
   } else if(depth > maxdepth) {
      return Vector();

   // Ray shooting
   } else {
      double pdf = 0.f;
      const auto e  = Vector(rng(), rng(), rng());
      const auto wo = -r.d;
      const auto wi = mat.Sample(wo, nl, e, pdf);
      if(Vector::Dot(wo, nl) <= 0.f || pdf <= 0.f) {
      	 return Vector();
      }
      auto f = Vector::Dot(wi, nl)*mat.Reflectance(wi, wo, nl);
      const Vector rad = Radiance(spheres, Ray(x, wi), rng, depth+1, maxdepth);

      return (1.f/pdf) * f.Multiply(rad);
   }
}

// Covariance Tracing includes
//*
#include <Covariance/Covariance4D.hpp>
using Cov4D   = Covariance::Covariance4D<Vector, double>;
/*/
#include <Covariance/InvCovariance4D.hpp>
using Cov4D  = Covariance::InvCovariance4D<Vector, double>;
//*/
using PosCov  = std::pair<Vector, Cov4D>;

std::ostream& operator<<(std::ostream& out, const Cov4D& cov) {
   out << std::scientific << std::showpos;
   out.precision(3);
   out << "[" << cov.matrix[0] << ",\t" << cov.matrix[1] << ",\t" << cov.matrix[3] << ",\t" << cov.matrix[6] << ";" << std::endl;
   out << " " << cov.matrix[1] << ",\t" << cov.matrix[2] << ",\t" << cov.matrix[4] << ",\t" << cov.matrix[7] << ";" << std::endl;
   out << " " << cov.matrix[3] << ",\t" << cov.matrix[4] << ",\t" << cov.matrix[5] << ",\t" << cov.matrix[8] << ";" << std::endl;
   out << " " << cov.matrix[6] << ",\t" << cov.matrix[7] << ",\t" << cov.matrix[8] << ",\t" << cov.matrix[9] << "]";
   return out;
}

PosCov CovarianceFilter(const std::vector<Sphere>& spheres, const Ray &r, const Cov4D& cov, int depth, int maxdepth=2) {
   double t;                               // distance to Intersection
   int id=0;                               // id of Intersected object
   if (!Intersect(spheres, r, t, id)) return PosCov(Vector(), Cov4D()); // if miss, return black
   const Sphere&   obj = spheres[id];      // the hit object
   const Material& mat = obj.mat;          // Its material
   Vector x  = r.o+r.d*t,
          n  = (x-obj.c).Normalize(),
          nl = Vector::Dot(n,r.d) < 0 ? n:n*-1;
   const double k = 1.f/spheres[id].r;

   // Update the covariance with travel and project it onto the tangent plane
   // of the hit object.
   Cov4D cov2 = cov;
   cov2.Travel(t);
   cov2.Projection(n);

   // if the max depth is reached
   if(depth >= maxdepth) {
      cov2.matrix[1] = - cov2.matrix[1];
      return PosCov(x, cov2);
   } else {

      // Sample a new direction
      auto wi = -r.d;
      auto wr = 2*Vector::Dot(wi, nl)*nl - wi;
      auto r2 = Ray(x, wr);

      cov2.Curvature(k, k);
      cov2.Cosine(1.0f);
      cov2.Symmetry();
      const double rho = mat.exponent / (4*M_PI*M_PI);
      cov2.Reflection(rho, rho);
      cov2.Curvature(-k, -k);
      cov2.InverseProjection(wr);
      return CovarianceFilter(spheres, r2, cov2, depth+1, maxdepth);
   }
}

// TinyEXR includes
#define TINYEXR_IMPLEMENTATION
#include <tinyexr/tinyexr.h>

int SaveEXR(const Vector* img, int w, int h, const std::string& filename) {
   EXRImage image;
   InitEXRImage(&image);
   image.num_channels  = 3;
   const char* names[] = {"B", "G", "R"};

    std::vector<float> images[3];
    images[0].resize(w * h);
    images[1].resize(w * h);
    images[2].resize(w * h);

    for (int i = 0; i < w * h; i++) {
      images[0][i] = img[i].x;
      images[1][i] = img[i].y;
      images[2][i] = img[i].z;
    }

    float* image_ptr[3];
    image_ptr[0] = &(images[2].at(0)); // B
    image_ptr[1] = &(images[1].at(0)); // G
    image_ptr[2] = &(images[0].at(0)); // R

    image.channel_names = names;
    image.images = (unsigned char**)image_ptr;
    image.width = w;
    image.height = h;
    image.compression = TINYEXR_COMPRESSIONTYPE_ZIP;

    image.pixel_types = (int *)malloc(sizeof(int) * image.num_channels);
    image.requested_pixel_types = (int *)malloc(sizeof(int) * image.num_channels);
    for (int i = 0; i < image.num_channels; i++) {
      image.pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT; // pixel type of input image
      image.requested_pixel_types[i] = TINYEXR_PIXELTYPE_HALF; // pixel type of output image to be stored in .EXR
    }

    const char* err;
    int ret = SaveMultiChannelEXRToFile(&image, filename.c_str(), &err);
    if (ret != 0) {
      fprintf(stderr, "Save EXR err: %s\n", err);
    }

    free(image.pixel_types);
    free(image.requested_pixel_types);
    return ret;
}
