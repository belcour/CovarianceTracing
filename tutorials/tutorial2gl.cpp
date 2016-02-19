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

std::default_random_engine gen;
std::uniform_real_distribution<double> dist(0,1);

// Local includes
#include "common.hpp"
#include "opengl.hpp"

// Texture for the background + size
int width = 512, height = 512;
float* background = nullptr;

ShaderProgram* program;

float gl_time = 0.0f;

void Draw() {

   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
   program->use();
   auto uniLocation = program->uniform("pointer");
   glUniform2f(uniLocation, mouse.Y, 1.0-mouse.X);

   glBegin(GL_QUADS);              // Drawing Using Triangles
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

   // Load textures
   background = new float[width*height];

   // Create the shader programs
   program = new ShaderProgram();
   std::string vertShader =
   "void main(void) {"
   "   gl_TexCoord[0] = gl_MultiTexCoord0;"
   "   gl_Position    = vec4(gl_Vertex);"
   "}";
   std::string fragShader =
   "uniform sampler2D texture;"
   "uniform vec2      pointer;"
   "void main(void) {"
   //"  gl_FragColor = vec4(1,0,0,1);"
   //"  gl_FragColor = vec4(gl_FragCoord.xy, 0.0, 1.0);"
   //"  gl_FragColor = vec4(gl_TexCoord[0].st - pointer, 0.0, 1.0);"
   "  float fact = exp(- pow(length(gl_TexCoord[0].xy - pointer.xy), 2.0));"
   "  gl_FragColor = fact*texture2D(texture, gl_TexCoord[0].st);"
   "}";
   program->initFromStrings(vertShader, fragShader);

   int w, h;
   float* img; const char* err;
   int ret = LoadEXR(&img, &w, &h, "ref.exr", &err);
   if(width != w || height != h) exit(EXIT_FAILURE);
   float* pixels = new float[width*height];
   for(int y=0; y<h; ++y)
      for(int x=0; x<w; ++x) {
         int i=(h-y-1)*w+x;
         int j=(w-x-1)*h+y;
         pixels[j] = img[i*4];
      }
   GLuint texture_id;
   glGenTextures(1, &texture_id);
   glBindTexture(GL_TEXTURE_2D, texture_id);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
   glTexImage2D(
      GL_TEXTURE_2D, // target
		0,             // level, 0 = base, no minimap,
		GL_RGBA,       // internal format
		width,         // width
		height,        // height
		0,             // border, always 0 in OpenGL ES
		GL_LUMINANCE,  // format
		GL_FLOAT,      // type
		pixels);       // ptr to data
   glActiveTexture(GL_TEXTURE0);
   glBindTexture(GL_TEXTURE_2D, texture_id);
   delete[] img;
   delete[] pixels;

   auto uniLocation = program->addUniform("pointer");
   glUniform2f(uniLocation, 0.1, 0.1);
}

int main(int argc, char** argv) {
   glutInit(&argc, argv);
   glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);

   glutInitWindowSize(width, height);
   glutCreateWindow("Tutorial number 2");

   Init();

   glutDisplayFunc(Draw);
   glutMouseFunc(MouseClicked);
   glutMainLoop();

   if(background) {
      delete[] background;
   }
   return EXIT_SUCCESS;
}
