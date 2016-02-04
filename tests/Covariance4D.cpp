// STL includes
#include <iostream>
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

   return nb_fails;
}

int TestProjection() {
   std::array<float, 10> matrix = {1.0,
                                   0.0, 1.0,
                                   0.0, 0.0, 1.0,
                                   0.0, 0.0, 0.0, 1.0};
   int nb_fails = 0;

   Vector x(1,0,0), y(0,1,0), z(0,0,1);
   Covariance4D<Vector> A(matrix, x, y, z), B(matrix, x, y, z);

   A.Projection(z);
   if(!IsApprox(A, B)) {
      std::cerr << "Error: Project with n == z is not identity" << std::endl;
      ++nb_fails;
   }
   return nb_fails;
}

int main(int argc, char** argv) {
   int nb_fails = 0;

   nb_fails += TestRotation();
   nb_fails += TestShear();
   nb_fails += TestProjection();

   if(nb_fails > 0) {
      return EXIT_FAILURE;
   } else {
      return EXIT_SUCCESS;
   }
};
