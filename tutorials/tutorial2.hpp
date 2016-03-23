#pragma once

// Local includes
#include "common.hpp"


/*****************************************************************************\

  Scene description:
  This scene is a classical 'Cornell box' with two glossy spheres inside. The
  scene defines also the camera position and the screen resolution.

   + 'spheres' contains the list of all objects with sphere geometry and phong
     materials.

   + 'width' and 'height' are the horizontal and vertical screen resolutions.

   + The camera is defined as a central ray 'cam' and a tangent plane 'cx',
     'cy'. Normalized versions of 'cx' and 'cy': 'ncx' and 'ncy' are provided.

\*****************************************************************************/

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

int width = 512, height = 512;

Ray cam(Vector(50.0f, 46.0f, 155.8f), Vector(0,0,-1));
double fov  = 1.2;
Vector  cx  = Vector(width*fov/height);
Vector  cy  = Vector::Cross(cx, cam.d).Normalize()*fov;
Vector ncx  = 1.0/Vector::Norm(cx) * cx;
Vector ncy  = 1.0/Vector::Norm(cy) * cy;



/*****************************************************************************\

  Predefined variables:

   + '*_img' are the buffers corresponding to the background image, the
     covariance filter image and the brute force pixel filter image. Each
     buffer is associated with a scaling factor to normalize the outputed
     image.

   + 'nPasses' is the number of passes used to render the background image
     using progressive rendering.

   + 'nPassesFilter' is the number of passes used to render the indirect pixel
     filter image using progressive rendering and 'filterRadius' is the radius
     used for density estimation.

\*****************************************************************************/

float* bcg_img = new float[width*height]; float bcg_scale = 1.0f;
float* cov_img = new float[width*height]; float cov_scale = 1.0f;
float* ref_img = new float[width*height]; float ref_scale = 1.0f;

int   nPasses       = 0;
int   nPassesFilter = 0;
float filterRadius  = 1.0f;

std::default_random_engine gen;
std::uniform_real_distribution<double> dist(0,1);

bool  displayBackground  = true;
bool  generateBackground = true;
bool  generateCovariance = true;
bool  generateReference  = false;



/*****************************************************************************\

   Render the indirect pixel filter after a specified number of bounces. To
   render the indirect pixel filter, I use density estimation.

   First the method 'indirect_filter' enables to compute the position in world
   space of a path tracing sample after 'maxdepth' bounces. Its return a pair
   of World space coordinate and RGB value corresponding to the importance of
   generate this point.

   Then, using the 'BruteForceTexture', a list of world space samples are used
   to perform density estimation for a complete frame. This method performs
   progressive density estimation to converge towards the correct filter image.

\*****************************************************************************/

using PosFilter = std::pair<Vector, Vector>;

PosFilter indirect_filter(const Ray &r, Random& rng, int depth, int maxdepth=2){
   double t;                               // distance to Intersection
   int id=0;                               // id of Intersected object
   if (!Intersect(spheres, r, t, id)) return PosFilter(Vector(), Vector()); // if miss, return black
   const Sphere&   obj = spheres[id];      // the hit object
   const Material& mat = obj.mat;          // Its material
   Vector x  = r.o+r.d*t,
          n  = (x-obj.c).Normalize(),
          nl = Vector::Dot(n,r.d) < 0 ? n:n*-1;

   // Once you reach the max depth, return the hit position and the filter's
   // value using the recursive form.
   if(depth >= maxdepth) {
      return PosFilter(x, Vector(1.0f, 1.0f, 1.0f));

   // Main covariance computation. First this code generate a new direction
   // and query the covariance+radiance in that direction. Then, it computes
   // the covariance after the reflection/refraction.
   } else {
      /* Sampling a new direction + recursive call */
      double pdf;
      const auto e   = Vector(rng(), rng(), rng());
      const auto wo  = -r.d;
      const auto wi  = mat.Sample(wo, nl, e, pdf);
            auto f   = Vector::Dot(wi, nl)*mat.Reflectance(wi, wo, nl);
      const auto res = indirect_filter(Ray(x, wi), rng, depth+1, maxdepth);
      return PosFilter(res.first, (1.f/pdf) * f.Multiply(res.second));
   }
}

void BruteForceTexture(int x, int y, int samps = 1000) {

   std::vector<PosFilter> _filter_elems;

   // Sub pixel sampling
   const int nthread = std::thread::hardware_concurrency();
   #pragma omp parallel for schedule(dynamic, 1)
   for(int t=0; t<nthread; ++t) {
      Random rng(t + nthread*clock());

      for(int s=0; s<samps/nthread; s++){

         // Create the RNG and get the sub-pixel sample
         float dx = rng();
         float dy = rng();

         // Generate the pixel direction
         Vector d = cx*((dx + x)/width  - .5) +
                    cy*((dy + y)/height - .5) + cam.d;
         d.Normalize();

         // Evaluate the Covariance and Radiance at the pixel location
         const auto filter = indirect_filter(Ray(cam.o, d), rng, 0, 1);
         #pragma omp critical
         {
            _filter_elems.push_back(filter);
         }
      }
   }

   // Loop over the rows and columns of the image and evaluate radiance and
   // covariance per pixel using Monte-Carlo.
   float max_ref = 0.0f;
   #pragma omp parallel for schedule(dynamic, 1), shared(max_ref)
   for (int y=0; y<height; y++){
      float max_temp = 0.0f;

      for (int x=0; x<width; x++) {
         int i=(width-x-1)*height+y;
         float _r = 0.0f;

         // Generate the pixel direction
         Vector d = cx*((0.5 + x)/width  - .5) +
                    cy*((0.5 + y)/height - .5) + cam.d;
         d.Normalize();

         Ray ray(cam.o, d);
         double t; int id;
         if(!Intersect(spheres, ray, t, id)){ continue; }
         Vector hitp = ray.o + t*ray.d;

         for(auto& elem : _filter_elems) {
            const auto& p = elem.first;
            const auto  x = Vector::Norm(hitp-p) / filterRadius;
            _r += (0.3989422804f/filterRadius) * exp(-0.5f * pow(x, 2)) * elem.second.x;
         }

         const auto scale = 20.f;
         const auto Nold  = nPassesFilter * samps;
         const auto Nnew  = (nPassesFilter+1) * samps;
         ref_img[i] = (ref_img[i]*Nold + scale*_r) / float(Nnew);
         max_temp   = std::max(ref_img[i], max_temp);
      }

      #pragma omp critical
      {
         max_ref = std::max(max_ref, max_temp);
      }
   }

   // Update the scaling
   ref_scale = 1.0f/max_ref;

   // Progressive refinement of the radius and number of passes
   filterRadius *= sqrt((nPassesFilter + 0.8) / (nPassesFilter + 1.0));
   ++nPassesFilter;
}



/*****************************************************************************\

      Render the background image using path tracing with implicit light
      connection.

\*****************************************************************************/

void RadianceTexture() {

   // Loop over the rows and columns of the image and evaluate radiance and
   // covariance per pixel using Monte-Carlo.
   #pragma omp parallel for schedule(dynamic, 1)
   for (int y=0; y<height; y++){
      Random rng(y + height*clock());

      for (int x=0; x<width; x++) {
         int i=(width-x-1)*height+y;

         // Create the RNG and get the sub-pixel sample
         float dx = rng();
         float dy = rng();

         // Generate the pixel direction
         Vector d = cx*((dx + x)/float(width)  - .5) +
                    cy*((dy + y)/float(height) - .5) + cam.d;
         d.Normalize();

         Ray ray(cam.o, d);
         Vector radiance = Radiance(spheres, ray, rng, 0, 1);

         bcg_img[i] = (float(nPasses)*bcg_img[i] + Vector::Dot(radiance, Vector(1,1,1))/3.0f) / float(nPasses+1);
      }
   }

   ++nPasses;
}



/*****************************************************************************\

               Export the image to EXR using TinyEXR library.

\*****************************************************************************/

void ExportImage(int x, int y) {
   Vector* img = new Vector[width*height];
   #pragma omp parallel for schedule(dynamic, 1)
   for(int i=0; i<width; ++i) {
      for(int j=0; j<height; ++j) {
         int id = i + j*width;
         int di = height-1-j + (width-1-i)*height;

         bool filter = ((abs(i-x)-4 < 0 && height-1-j==y) || (abs((height-1-j)-y)-4 < 0 && i==x)) &&
                       (generateReference || generateCovariance);

         // Background img
         img[id] = (displayBackground) ? bcg_img[di]*Vector(1,1,1) : Vector(0,0,0);

         // Adding the filters
         img[id].x += generateCovariance ? cov_scale*cov_img[di] : 0.0f;
         img[id].y += generateReference  ? ref_scale*ref_img[di] : 0.0f;
         img[id].z += (filter ? 1.0f : 0.0f);
      }
   }

   int ret = SaveEXR(img, width, height, "output.exr");
   if(ret != 0) { std::cerr << "Unable to export image" << std::endl; }
}
