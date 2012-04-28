// Source file for the particle system



// Include files

#include "R2/R2.h"
#include "R3/R3.h"
#include "R3Scene.h"
#include "particle.h"
using namespace std;
#ifdef _WIN32
#   include <windows.h>
#else
#   include <sys/time.h>
#endif




////////////////////////////////////////////////////////////
// Random Number Generator
////////////////////////////////////////////////////////////

double RandomNumber(void)
{
#ifdef _WIN32
  // Seed random number generator
  static int first = 1;
  if (first) {
    srand(GetTickCount());
    first = 0;
  }

  // Return random number
  int r1 = rand();
  double r2 = ((double) rand()) / ((double) RAND_MAX);
  return (r1 + r2) / ((double) RAND_MAX);
#else 
  // Seed random number generator
  static int first = 1;
  if (first) {
    struct timeval timevalue;
    gettimeofday(&timevalue, 0);
    srand48(timevalue.tv_usec);
    first = 0;
  }

  // Return random number
  return drand48();
#endif
}



////////////////////////////////////////////////////////////
// Generating Particles
////////////////////////////////////////////////////////////

void GenerateParticles(R3Scene *scene, double current_time, double delta_time)
{
  // Generate new particles for every source

  // FILL IN CODE HERE
}



////////////////////////////////////////////////////////////
// Updating Particles
////////////////////////////////////////////////////////////

void UpdateParticles(R3Scene *scene, double current_time, double delta_time, int integration_type)
{
  // Update position for every particle

  // FILL IN CODE HERE
}



////////////////////////////////////////////////////////////
// Rendering Particles
////////////////////////////////////////////////////////////

void RenderParticles(R3Scene *scene, double current_time, double delta_time)
{
  // Draw every particle

  // REPLACE CODE HERE
  glDisable(GL_LIGHTING);
  glPointSize(5);
  glBegin(GL_POINTS);
  for (int i = 0; i < scene->NParticles(); i++) {
    R3Particle *particle = scene->Particle(i);
    glColor3d(particle->material->kd[0], particle->material->kd[1], particle->material->kd[2]);
    const R3Point& position = particle->position;
    glVertex3d(position[0], position[1], position[2]);
  }   
  glEnd();
}



