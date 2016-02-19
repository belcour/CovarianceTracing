#pragma once

// Include STL
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

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

   // TODO: Finish implementation
   Vector DirectionToPixel(const Vector dir) const {
      double u,v;
      u = Vector::Dot(dir, cx);
      v = Vector::Dot(dir, cy);
   }
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
