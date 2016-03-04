// OpenGL includes
#ifdef __APPLE__
#include <glut.h>
#include <gl.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#endif

// STL includes
#include <cstdlib>
#include <cstdio>
#include <random>
#include <utility>
#include <string>
#include <iostream>
#include <sstream>

std::default_random_engine gen;
std::uniform_real_distribution<double> dist(0,1);

// Local includes
#include "common.hpp"
#include "opengl.hpp"

// Scene description
Material phongH(Vector(), Vector(), Vector(1,1,1)*.999, 1.E3);
Material phongL(Vector(), Vector(), Vector(1,1,1)*.999, 1.E2);

std::vector<Sphere> spheres = {
   Sphere(Vector(27,16.5,47),        16.5, phongH),//RightSp
   Sphere(Vector(73,16.5,78),        16.5, phongL),//LeftSp
   Sphere(Vector( 1e5+1,40.8,81.6),  1e5,  Vector(),Vector(.75,.25,.25)),//Left
   Sphere(Vector(-1e5+99,40.8,81.6), 1e5,  Vector(),Vector(.25,.25,.75)),//Rght
   Sphere(Vector(50,40.8, 1e5),      1e5,  Vector(),Vector(.75,.75,.75)),//Back
   Sphere(Vector(50,40.8,-1e5+170),  1e5,  Vector(),Vector()           ),//Frnt
   Sphere(Vector(50, 1e5, 81.6),     1e5,  Vector(),Vector(.75,.75,.75)),//Botm
   Sphere(Vector(50,-1e5+81.6,81.6), 1e5,  Vector(),Vector(.75,.75,.75)),//Top
   Sphere(Vector(50,681.6-.27,81.6), 600,  Vector(12,12,12),  Vector()) //Lite
};

// Texture for the bcg_img + size
int  width = 512, height = 512;
bool generateBackground = true;
int  nPasses = 0;

ShaderProgram* program;

// Different buffer, the background image, the covariance filter and the brute
// force filter.
GLuint texs_id[3];
float* bcg_img = new float[width*height];
float* cov_img = new float[width*height];
float* ref_img = new float[width*height];

// Camera frame
Ray cam(Vector(50,52,295.6), Vector(0,-0.042612,-1).Normalize());
double fov  = 1.2;
Vector  cx  = Vector(width*fov/height);
Vector  cy  = Vector::Cross(cx, cam.d).Normalize()*fov;
Vector ncx  = cx;
Vector ncy  = cy;

std::stringstream sout;

void ExportImage() {
   Vector* img = new Vector[width*height];
   for(int i=0; i<width*height; ++i) {
      img[i].x = bcg_img[i] + cov_img[i];
      img[i].y = bcg_img[i];
      img[i].z = bcg_img[i];
   }

   int ret = SaveEXR(img, width, height, "output.exr");
   if(ret != 0) { std::cerr << "Unable to export image" << std::endl; }
}

void KeyboardKeys(unsigned char key, int x, int y) {
   if(key == 'b') {
      generateBackground = !generateBackground;
   } else if(key == '+') {
      Material phong(Vector(), Vector(), Vector(1,1,1)*.999, spheres[1].mat.exponent * 10);
      spheres[1].mat = phong;
      nPasses = 0;
   } else if(key == '-') {
      Material phong(Vector(), Vector(), Vector(1,1,1)*.999, fmax(spheres[1].mat.exponent / 10, 1.0));
      spheres[1].mat = phong;
      nPasses = 0;
   } else if(key == 'p') {
      ExportImage();
   } else if(key == 'd') {
      std::cout << sout.str() << std::endl;
   }
   glutPostRedisplay();
}

void RadianceTexture() {

   const float dx = dist(gen);
   const float dy = dist(gen);

   // Loop over the rows and columns of the image and evaluate radiance and
   // covariance per pixel using Monte-Carlo.
   #pragma omp parallel for schedule(dynamic, 1)
   for (int y=0; y<height; y++){
      for (int x=0; x<width; x++) {

         int i=(width-x-1)*height+y;
         float _r = 0.0f;
         // Generate the pixel direction
         Vector d = cx*((dx + x)/float(width)  - .5) +
            cy*((dy + y)/float(height) - .5) + cam.d;
         d.Normalize();

         Ray ray(cam.o, d);
         Vector radiance = Radiance(spheres, ray, 0, 1);

         bcg_img[i] = (float(nPasses)*bcg_img[i] + radiance.x) / float(nPasses+1);
      }
   }

   ++nPasses;
   glutPostRedisplay();
}

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

   cov2.Projection(n);
   out << "After projection" << std::endl;
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
      out << "After BRDF convolution of sigma=" << rho << std::endl;
      out << cov2 << std::endl;
      out << "Volume = " << cov2.Volume() << std::endl;
      out << std::endl;

      cov2.Curvature(-k, -k);
      out << "After inverse curvature" << std::endl;
      out << cov2 << std::endl;
      out << "Volume = " << cov2.Volume() << std::endl;
      out << std::endl;

      cov2.InverseProjection(wr);
      out << "After inverse projection" << std::endl;
      out << cov2 << std::endl;
      out << "Volume = " << cov2.Volume() << std::endl;
      out << std::endl;

      return CovarianceFilter(spheres, r2, cov2, depth+1, maxdepth, out);
   }
}

void CovarianceTexture() {
   // Generate a covariance matrix at the sampling position
   int x = width*mouse.X, y = height*mouse.Y;
   const auto t = (cx*((x+0.5)/double(width) - .5) + cy*((y+0.5)/double(height) - .5) + cam.d).Normalize();
   const auto u = (ncx - Vector::Dot(t, ncx)*t).Normalize();
   const auto v = (ncy - Vector::Dot(t, ncy)*t).Normalize();
   const auto pixelCov = Cov4D({ 1.0E5, 0.0, 1.0E5, 0.0, 0.0, 1.0E5, 0.0, 0.0, 0.0, 1.0E5 }, t);

   sout = std::stringstream();
   const auto surfCov  = CovarianceFilter(spheres, Ray(cam.o, t), pixelCov, 0, 1, sout);
   sout << surfCov.second << std::endl;
   sout << "Volume = " << surfCov.second.Volume() << std::endl;
   sout << std::endl;

   double sxx = 0, syy = 0, sxy = 0;
   try {
      surfCov.second.SpatialFilter(sxx, sxy, syy);
   } catch (...) {
      std::cout << "Error: incorrect spatial filter" << std::endl;
      sout << surfCov.second << std::endl;
      return;
   }
   //std::cout << "Filter = [" << sxx << ", " << sxy << ", " << syy << "]" << std::endl;

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
         const double du  = Vector::Dot(dx, surfCov.second.x);
         const double dv  = Vector::Dot(dx, surfCov.second.y);
         const double dt  = Vector::Dot(dx, surfCov.second.z);

         double bf = du*du*sxx + dv*dv*syy + 2*du*dv*sxy;
         cov_img[i] = exp(-10.0*dt*dt) * exp(- 0.5* bf);
      }
   }
}

void Draw() {

   if(generateBackground) {
      RadianceTexture();

      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_2D, texs_id[0]);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_LUMINANCE, GL_FLOAT, bcg_img);
   }

   CovarianceTexture();
   glActiveTexture(GL_TEXTURE1);
   glBindTexture(GL_TEXTURE_2D, texs_id[1]);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
   glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_LUMINANCE, GL_FLOAT, cov_img);

   program->use();

   glActiveTexture(GL_TEXTURE0);
   glBindTexture(GL_TEXTURE_2D, texs_id[0]);

   glActiveTexture(GL_TEXTURE1);
   glBindTexture(GL_TEXTURE_2D, texs_id[1]);

   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
   auto uniLocation = program->uniform("pointer");
   glUniform2f(uniLocation, mouse.Y, 1.0-mouse.X);

   glBegin(GL_QUADS);
   glVertex3f(-1.0f,-1.0f, 0.0f); glTexCoord2f(0, 0);
   glVertex3f( 1.0f,-1.0f, 0.0f); glTexCoord2f(1, 0);
   glVertex3f( 1.0f, 1.0f, 0.0f); glTexCoord2f(1, 1);
   glVertex3f(-1.0f, 1.0f, 0.0f); glTexCoord2f(0, 1);
   glEnd();
   program->disable();
   glutSwapBuffers();
}

// Create geometry and textures
void Init() {
   // Background color
   glClearColor(0.0f, 0.0f, 0.0f, 2.0f);

   // Create the shader programs
   program = new ShaderProgram(false);
   std::string vertShader =
      "void main(void) {"
      "   gl_TexCoord[0] = gl_MultiTexCoord0;"
      "   gl_Position    = vec4(gl_Vertex);"
      "}";
   std::string fragShader =
      "uniform sampler2D tex0;"
      "uniform sampler2D tex1;"
      "uniform sampler2D tex2;"
      "uniform vec2      pointer;"
      "uniform float     width;"
      "uniform float     height;"
      "void main(void) {"
      "  float fact = exp(- width*height * pow(length(gl_TexCoord[0].xy - pointer.xy), 2.0));"
      "  gl_FragColor = vec4(0,0,1,1)*fact + vec4(1,1,1,1)*texture2D(tex0, gl_TexCoord[0].st) + vec4(1,0,0,1)*texture2D(tex1, gl_TexCoord[0].st) + vec4(0,1,0,1)*texture2D(tex2, gl_TexCoord[0].st);"
      "}";
   program->initFromStrings(vertShader, fragShader);

   // Reserve textures on the GPU
   glGenTextures(3, texs_id);

   // Define the different uniform locations in the shader
   program->use();

   const auto t1Location = program->addUniform("tex0");
   glUniform1i(t1Location, 0);
   const auto t2Location = program->addUniform("tex1");
   glUniform1i(t2Location, 1);
   const auto t3Location = program->addUniform("tex2");
   glUniform1i(t3Location, 2);

   const auto uniWidth = program->addUniform("width");
   glUniform1f(uniWidth, float(width));
   const auto uniHeight = program->addUniform("height");
   glUniform1f(uniHeight, float(height));

   const auto uniLocation = program->addUniform("pointer");
   glUniform2f(uniLocation, width*mouse.X, height*mouse.Y);

   program->disable();

   // Clean buffer memory
   memset(ref_img, 0.0, width*height);
   memset(cov_img, 0.0, width*height);
}

void PrintHelp() {
   std::cout << "Covariance Tracing tutorial 2" << std::endl;
   std::cout << "----------------------------" << std::endl;
   std::cout << std::endl;
   std::cout << "This tutorial display the indirect pixel filter after one bounce for" << std::endl;
   std::cout << "non-diffuse surfaces. To display the filter, click on one of the two" << std::endl;
   std::cout << "shiny spheres." << std::endl;
   std::cout << std::endl;
   std::cout << " + 'b' stop/resume rendering the bcg_img image" << std::endl;
   std::cout << " + 'p' output image to EXR file" << std::endl;
   std::cout << " + 'd' print Covariance Tracing step by step" << std::endl;
   std::cout << " + '+' increase the BRDF exponent for the right sphere" << std::endl;
   std::cout << " + '-' decrease the BRDF exponent for the right sphere" << std::endl;
   std::cout << std::endl;
}

int main(int argc, char** argv) {
   cam.o = cam.o + 140.0*cam.d;
   ncx.Normalize();
   ncy.Normalize();

   PrintHelp();

   glutInit(&argc, argv);
   glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);

   glutInitWindowSize(width, height);
   glutCreateWindow("Covariance Tracing tutorial 2");

   Init();

   glutDisplayFunc(Draw);
   glutMouseFunc(MouseClicked);
   glutMotionFunc(MouseMoved);
   glutKeyboardFunc(KeyboardKeys);
   glutMainLoop();

   if(bcg_img) { delete[] bcg_img; }
   if(cov_img) { delete[] cov_img; }
   if(ref_img) { delete[] ref_img; }
   return EXIT_SUCCESS;
}
