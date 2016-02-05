#pragma once

#include <cmath>
#include <iostream>

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
      Vector r;
      double norm = sqrt(Dot(*this, *this));
      r.x = x/norm;
      r.y = y/norm;
      r.z = z/norm;
      return r;
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
      v.x = w.x;
      v.y = w.y;
      v.z = w.z;
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

struct Sphere {
   // Geometric information for the sphere
   Vector c;
   double r;

   // Emission and reflectance color
   Vector ke;
   Vector kd;

   Sphere(const Vector& c, double r, const Vector ke, const Vector& kd) 
      : c(c), r(r), ke(ke), kd(kd) {}

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
