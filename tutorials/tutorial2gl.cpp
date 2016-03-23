// STL includes
#include <cstdlib>
#include <cstdio>
#include <random>
#include <utility>
#include <string>
#include <iostream>
#include <sstream>
#include <thread>

// Local includes
#include "common.hpp"
#include "opengl.hpp"
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

void CovarianceTexture() {
   // Generate a covariance matrix at the sampling position
   int x = width*mouse.X, y = height*mouse.Y;
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

ShaderProgram* program;

// Different buffer, the background image, the covariance filter and the brute
// force filter.
GLuint texs_id[3];

void KeyboardKeys(unsigned char key, int x, int y) {
   if(key == 'c') {
      memset(cov_img, 0.0, width*height);
      generateCovariance = !generateCovariance;
   } else if(key == 'B') {
      generateReference = !generateReference;
   } else if(key == 'b') {
      generateBackground = !generateBackground;
   } else if(key == 'h') {
      displayBackground  = !displayBackground;
   } else if(key == 'f') {
       useCovFilter = !useCovFilter;
   } else if(key == '+') {
      Material phong(Vector(), Vector(), Vector(1,1,1)*.999, spheres[1].mat.exponent * 10);
      spheres[1].mat = phong;
      nPasses = 0;
   } else if(key == '-') {
      Material phong(Vector(), Vector(), Vector(1,1,1)*.999, fmax(spheres[1].mat.exponent / 10, 1.0));
      spheres[1].mat = phong;
      nPasses = 0;
   } else if(key == 'p') {
      ExportImage(width*mouse.X, height*mouse.Y);
   } else if(key == 'd') {
      std::cout << sout.str() << std::endl;
   }
   glutPostRedisplay();
}

void Draw() {

   if(generateBackground) {
      RadianceTexture();

      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_2D, texs_id[0]);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_LUMINANCE, GL_FLOAT, bcg_img);
   }

   if(generateCovariance) {
      CovarianceTexture();

      glActiveTexture(GL_TEXTURE1);
      glBindTexture(GL_TEXTURE_2D, texs_id[1]);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_LUMINANCE, GL_FLOAT, cov_img);
   }

   if(generateReference) {

      if(fabs(mouse.Dx) > 0.0f || fabs(mouse.Dy) > 0.0f) {
         nPassesFilter = 0;
         filterRadius  = 1.0f;
      }
      BruteForceTexture(width*mouse.X, height*mouse.Y);

      glActiveTexture(GL_TEXTURE2);
      glBindTexture(GL_TEXTURE_2D, texs_id[2]);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_LUMINANCE, GL_FLOAT, ref_img);
   }

   program->use();

   glActiveTexture(GL_TEXTURE0);
   glBindTexture(GL_TEXTURE_2D, texs_id[0]);

   glActiveTexture(GL_TEXTURE1);
   glBindTexture(GL_TEXTURE_2D, texs_id[1]);

   glActiveTexture(GL_TEXTURE2);
   glBindTexture(GL_TEXTURE_2D, texs_id[2]);

   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
   auto uniLocation = program->uniform("pointer");
   glUniform2f(uniLocation, mouse.Y, 1.0-mouse.X);

   // Update the scaling
   glUniform1f(program->uniform("tex0scale"), displayBackground  ? bcg_scale : 0.0f);
   glUniform1f(program->uniform("tex1scale"), generateCovariance ? cov_scale : 0.0f);
   glUniform1f(program->uniform("tex2scale"), generateReference  ? ref_scale : 0.0f);

   glBegin(GL_QUADS);
   glVertex3f(-1.0f,-1.0f, 0.0f); glTexCoord2f(0, 0);
   glVertex3f( 1.0f,-1.0f, 0.0f); glTexCoord2f(1, 0);
   glVertex3f( 1.0f, 1.0f, 0.0f); glTexCoord2f(1, 1);
   glVertex3f(-1.0f, 1.0f, 0.0f); glTexCoord2f(0, 1);
   glEnd();
   program->disable();
   glutSwapBuffers();

   if(generateReference || generateBackground || generateCovariance) {
      glutPostRedisplay();
   }
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
      "uniform sampler2D tex0; uniform float tex0scale;"
      "uniform sampler2D tex1; uniform float tex1scale;"
      "uniform sampler2D tex2; uniform float tex2scale;"
      "uniform vec2      pointer;"
      "uniform float     width;"
      "uniform float     height;"
      "void main(void) {"
      "  float fact = exp(- width*height * pow(length(gl_TexCoord[0].xy - pointer.xy), 2.0));"
      "  gl_FragColor = vec4(0,0,1,1)*fact + tex0scale*vec4(1,1,1,1)*texture2D(tex0, gl_TexCoord[0].st) + tex1scale*vec4(1,0,0,1)*texture2D(tex1, gl_TexCoord[0].st) + tex2scale*vec4(0,1,0,1)*texture2D(tex2, gl_TexCoord[0].st);"
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

   const auto t1sLocation = program->addUniform("tex0scale");
   glUniform1f(t1sLocation, bcg_scale);
   const auto t2sLocation = program->addUniform("tex1scale");
   glUniform1f(t2sLocation, cov_scale);
   const auto t3sLocation = program->addUniform("tex2scale");
   glUniform1f(t3sLocation, ref_scale);

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
   std::cout << " + 'b' stop/resume rendering the background image" << std::endl;
   std::cout << " + 'h' hide/show the background image" << std::endl;
   std::cout << " + 'c' stop/resume rendering the covariance filter" << std::endl;
   std::cout << " + 'B' brute-force indirect pixel filter (SLOW)" << std::endl;
   std::cout << " + 'p' output image to EXR file" << std::endl;
   std::cout << " + 'd' print Covariance Tracing step by step" << std::endl;
   std::cout << " + 'f' switch between displaying the Gaussian filter or the polygonal footprint" << std::endl;
   std::cout << " + '+' increase the BRDF exponent for the right sphere" << std::endl;
   std::cout << " + '-' decrease the BRDF exponent for the right sphere" << std::endl;
   std::cout << std::endl;
}

int main(int argc, char** argv) {

   mouse.X = 0.5;
   mouse.Y = 0.5;

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
