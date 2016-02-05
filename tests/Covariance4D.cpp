// STL includes
#include <iostream>
#include <iomanip>
#include <cmath>

// Covariance includes
#include <Covariance/Covariance4D.hpp>

using namespace Covariance;

struct Vector {
   float x, y, z;
   Vector() {}
   Vector(float x, float y, float z) : x(x), y(y), z(z) {}
   static float Dot(const Vector& w1, const Vector& w2) {
      return w1.x*w2.x + w1.y*w2.y + w1.z*w2.z;
   }
   static Vector Cross(const Vector& u, const Vector& v) {
      Vector r;
      r.x = u.y*v.z - u.z*v.y;
      r.y = u.z*v.x - u.x*v.z;
      r.z = u.x*v.y - u.y*v.x;
      return r;
   }
   void Normalize() {
      float norm = sqrt(Dot(*this, *this));
      x /= norm;
      y /= norm;
      z /= norm;
   }
   friend Vector operator*(float a, const Vector& w) {
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
};

bool IsApprox(const Covariance4D<Vector>& A, const Covariance4D<Vector>& B, float Eps=1.0E-5) {
   bool isApprox = true;
   for(int i=0; i<10; ++i) {
      isApprox &= std::abs(A.matrix[i] - B.matrix[i]) < Eps;
   }
   return isApprox;
}

std::ostream& operator<<(std::ostream& out, const Covariance4D<Vector>& A) {
   for(int i=0; i<10; ++i) {
      out << A.matrix[i] << ", ";
   }
   return out;
}

int TestRotation() {
   int nb_fails = 0;

   Covariance4D<Vector> A, B;

   A = B = Covariance4D<Vector>(1.0, 0.0, 1.0, 0.0);
   A.Rotate(2.0*M_PI);
   if(!IsApprox(A, B)) {
      std::cerr << "Error: 2π rotation is not indempotent" << std::endl;
      std::cerr << A << std::endl;
      std::cerr << B << std::endl;
      ++nb_fails;
   }

   A = B;
   A.Rotate( 0.567);
   A.Rotate(-0.567);
   if(!IsApprox(A, B)) {
      std::cerr << "Error: Rotation and its inverse is not indempotent" << std::endl;
      std::cerr << "       cos = " << cos(0.567) << std::endl;
      std::cerr << "       sin = " << sin(0.567) << std::endl;
      A = B;
      std::cerr << A << std::endl;
      A.Rotate( 0.567);
      std::cerr << A << std::endl;
      A.Rotate(-0.567);
      std::cerr << A << std::endl;
      ++nb_fails;
   }

   return nb_fails;
}

int TestShear() {
   int nb_fails = 0;

   Covariance4D<Vector> A, B;
   float c = 123.456;
   float d = 765.432;

   A = B = Covariance4D<Vector>(1.0, 2.0, 3.0, 4.0);
   A.ShearAngleSpace( c,  d);
   A.ShearAngleSpace(-c, -d);

   if(!IsApprox(A, B)) {
      std::cerr << "Error: Shear Angle -> Space is not indempotent" << std::endl;
      std::cerr << A << std::endl;
      std::cerr << B << std::endl;
      ++nb_fails;
   }

   A = B = Covariance4D<Vector>(1.0, 2.0, 3.0, 4.0);
   A.ShearSpaceAngle( c,  d);
   A.ShearSpaceAngle(-c, -d);

   if(!IsApprox(A, B)) {
      std::cerr << "Error: Shear Space -> Angle is not indempotent" << std::endl;
      std::cerr << A << std::endl;
      std::cerr << B << std::endl;
      ++nb_fails;
   }

   A = Covariance4D<Vector>(1.0, 2.0, 0.0, 0.0);
   B = Covariance4D<Vector>(0.0, 0.0, 1.0, 2.0);
   A.Travel(1);
   A.Curvature(-1, -1);
   if(!IsApprox(A, B)) {
      std::cerr << "Error: Curvature + Travel incorrect" << std::endl;
      std::cerr << A << std::endl;
      std::cerr << B << std::endl;
      ++nb_fails;
   }

   // Mimicking a lens looking at the plane in focus
   A = Covariance4D<Vector>(1.0, 1.0, 0.0, 0.0);
   B = A;
   A.Travel(1);
   A.Curvature(-2, -2);
   A.Travel(1);
   if(!IsApprox(A, B)) {
      std::cerr << "Error: Lens operator incorrect" << std::endl;
      std::cerr << A << std::endl;
      std::cerr << B << std::endl;
      ++nb_fails;
   }

   return nb_fails;
}

int TestProjection() {
   std::array<float, 10> matrix = {1.0,
                                   0.1, 1.0,
                                   0.0, 0.0, 1.0,
                                   0.0, 0.0, 0.1, 1.0};
   int nb_fails = 0;

   Vector x(1,0,0), y(0,1,0), z(0,0,1), n(0,0,1);
   /*
   std::cout << Vector::Cross(x, y) << " : " << z << std::endl;
   std::cout << Vector::Cross(z, x) << " : " << y << std::endl;
   std::cout << Vector::Cross(y, z) << " : " << x << std::endl;
   */

   Covariance4D<Vector> A(matrix, x, y, z), B(matrix, x, y, z);

   A.Projection(z);
   if(!IsApprox(A, B)) {
      std::cerr << "Error: Project with n == z is not identity" << std::endl;
      std::cerr << A << std::endl;
      std::cerr << B << std::endl;
      ++nb_fails;
   }

   /* Perform a specular reflection on the identity covariance */
   {
   const float rho = std::numeric_limits<float>::max();
   const float k   = 0.0;
   x = Vector(0.0, 1.0,  0.0);
   y = Vector(0.5, 0.0,  0.5); y.Normalize();
   z = Vector(0.5, 0.0, -0.5); z.Normalize();
   Vector o = y;
   A = Covariance4D<Vector>(matrix, x, y, z);
   A.Travel(1);
   B = A; B.Symmetry(); // Reflection mirror along Y
   //std::cout << "A: " << A << std::endl;
   A.Projection(n);
   //std::cout << "After proj: " << A << std::endl;
   A.Curvature(k, k);
   //std::cout << "After curv: " << A << std::endl;
   A.Symmetry();
   //std::cout << "After symm: " << A << std::endl;
   A.Reflection(rho, rho);
   //std::cout << "After refl: " << A << std::endl;
   A.Curvature(-k, -k);
   //std::cout << "After icur: " << A << std::endl;
   A.InverseProjection(o);
   //std::cout << "After ipro: " << A << std::endl;
   if(!IsApprox(A, B)) {
      std::cerr << "Error: Specular reflection preserve content but negate Y-correlations" << std::endl;
      std::cerr << A << std::endl;
      std::cerr << B << std::endl;
      ++nb_fails;
   }
   }

   {
   const float rho = std::numeric_limits<float>::max();
   const float k   = 1.0;
   matrix[0] = 0.0;
   matrix[1] = 0.0;
   matrix[2] = 0.0;
   matrix[5] = 1.0;
   matrix[8] = 0.0;
   matrix[9] = 1.0;
   x = Vector( 1, 0, 0);
   y = Vector( 0,-1, 0);
   z = Vector( 0, 0,-1);
   n = Vector( 0, 0, 1);
   A = Covariance4D<Vector>(matrix, x, y, z);
   B = A;
   //std::cout << "A: " << A << std::endl;
   A.Projection(n);
   //std::cout << "After proj: " << A << std::endl;
   A.Curvature(k, k);
   //std::cout << "After curv: " << A << std::endl;
   A.Symmetry();
   //std::cout << "After symm: " << A << std::endl;
   A.Reflection(rho, rho);
   //std::cout << "After refl: " << A << std::endl;
   A.Curvature(-k, -k);
   //std::cout << "After icur: " << A << std::endl;
   A.InverseProjection(n);
   //std::cout << "After ipro: " << A << std::endl;
   if(A.matrix[5] != B.matrix[5] && A.matrix[9] != B.matrix[9] &&
      A.matrix[8] != B.matrix[8] && A.matrix[0] != 2.0*B.matrix[0] &&
      A.matrix[2] != 2*B.matrix[2] && A.matrix[1] != 2*B.matrix[1]) {
      std::cerr << "Error: Specular curved reflection preserve angular content" << std::endl;
      std::cerr << "       but increase spatial content." << std::endl;
      std::cerr << A << std::endl;
      std::cerr << B << std::endl;
      ++nb_fails;
   }}

   return nb_fails;
}


int TestReflection() {
   int nb_fails = 0;

   Covariance4D<Vector> A(0.0f, 0.0f, 1.0f, 1.0f);
   Covariance4D<Vector> B = A;
   Covariance4D<Vector> Z(0.0f, 0.0f, 0.0f, 0.0f);

   A.Reflection(0.0f, 0.0f);
   if(!IsApprox(A, Z)) {
      std::cerr << "Error: Diffuse reflection does not kill angular freqs" << std::endl;
      std::cerr << A << std::endl;
      std::cerr << Z << std::endl;
      ++nb_fails;
   }

   A = B;
   float rho = std::numeric_limits<float>::max();
   A.Reflection(rho, rho);
   if(!IsApprox(A, B)) {
      std::cerr << "Error: Specular reflection reduces angular freqs" << std::endl;
      std::cerr << A << std::endl;
      std::cerr << B << std::endl;
      ++nb_fails;
   }

   return nb_fails;
}

int main(int argc, char** argv) {
   int nb_fails = 0;
   std::cout << std::fixed << std::showpos << std::setprecision(2);

   nb_fails += TestRotation();
   nb_fails += TestShear();
   nb_fails += TestProjection();
   nb_fails += TestReflection();

   if(nb_fails > 0) {
      return EXIT_FAILURE;
   } else {
      return EXIT_SUCCESS;
   }
};
