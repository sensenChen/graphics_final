#ifndef __ARG_PARSER_H__
#define __ARG_PARSER_H__

#include <cassert>
#include <string>
#include <glm/glm.hpp>
#include "mtrand.h"
#include <math.h>

// ================================================================================
// ================================================================================

inline void separatePathAndFile(const std::string &input, std::string &path, std::string &file) {
  // we need to separate the filename from the path
  // (we assume the vertex & fragment shaders are in the same directory)
  // first, locate the last '/' in the filename
  size_t last = std::string::npos;
  while (1) {
    int next = input.find('/',last+1);
    if (next != (int)std::string::npos) {
      last = next;
      continue;
    }
    next = input.find('\\',last+1);
    if (next != (int)std::string::npos) {
      last = next;
      continue;
    }
    break;
  }
  if (last == std::string::npos) {
    // if there is no directory in the filename
    file = input;
    path = ".";
  } else {
    // separate filename & path
    file = input.substr(last+1,input.size()-last-1);
    path = input.substr(0,last);
  }
}

// ====================================================================================

class ArgParser {

public:

  ArgParser() { DefaultValues(); }

  ArgParser(int argc, char *argv[]) {
    DefaultValues();

    for (int i = 1; i < argc; i++) {
      if (argv[i] == std::string("-cloth")) {
        i++; assert (i < argc);
        separatePathAndFile(argv[i],path,cloth_file);
      } else if (argv[i] == std::string("-fluid")) {
	i++; assert (i < argc);
        separatePathAndFile(argv[i],path,fluid_file);
      } else if (argv[i] == std::string("-size")) {
        i++; assert (i < argc);
	width = height = atoi(argv[i]);
      } else if (argv[i] == std::string("-timestep")) {
	i++; assert (i < argc);
	timestep = atof(argv[i]);
        assert (timestep > 0);
      } else {
	printf ("whoops error with command line argument %d: '%s'\n",i,argv[i]);
	assert(0);
      }
    }
  }

  // ===================================
  // ===================================

  void DefaultValues() {
    pipex = 0;
    pipey = 0;
    pipez = 0;
    theta = 0*M_PI;

    //rotation matrix
    // r1 = glm::vec4(1,0,0,0);
    // r2 = glm::vec4(0,cos(theta),-1*sin(theta),0);
    // r3 = glm::vec4(0,sin(theta),cos(theta),0);

    //transformation matrix
    r1 = glm::vec4(1,0,0,pipex);
    r2 = glm::vec4(0,1,0,pipey);
    r3 = glm::vec4(0,0,1,pipez);


    width = 500;
    height = 500;

    timestep = 0.01;
    animate = false;

    particles = true;
    velocity = true;
    force = true;

    face_velocity = 0;
    dense_velocity = 0;

    surface = false;
    isosurface = 0.7;

    wireframe = false;
    bounding_box = true;
    cubes = false;
    pressure = false;

    gravity = glm::vec3(0,-9.8,0);

    // uncomment for deterministic randomness
    // mtrand = MTRand(37);

  }

  // ===================================
  // ===================================
  // REPRESENTATION
  // all public! (no accessors)

  glm::vec3 mult(glm::vec3& pos) {
    glm::vec3 ans = glm::vec3(
      r1[0] * pos.x + r1[1] * pos.y + r1[2] * pos.z + r1[3],
      r2[0] * pos.x + r2[1] * pos.y + r2[2] * pos.z + r2[3],
      r3[0] * pos.x + r3[1] * pos.y + r3[2] * pos.z + r3[3]);

    return ans;
  }

  std::string cloth_file;
  std::string fluid_file;
  std::string path;
  int width;
  int height;

  double pipex;
  double pipey;
  double pipez;
  double theta;


  glm::vec4 r1;
  glm::vec4 r2;
  glm::vec4 r3;
  // glm::vec4 r4;

  // animation control
  double timestep;
  bool animate;
  glm::vec3 gravity;

  // display option toggles
  // (used by both)
  bool particles;
  bool velocity;
  bool surface;
  bool bounding_box;

  // used by cloth
  bool force;
  bool wireframe;

  // used by fluid
  int face_velocity;
  int dense_velocity;
  double isosurface;
  bool cubes;
  bool pressure;

  // default initialization
  MTRand mtrand;

};

// ================================================================================

#endif
