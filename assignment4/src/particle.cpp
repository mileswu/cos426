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
#include "raytrace.h"




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

  for (int i = 0; i < scene->NParticleSources(); i++) {
    R3ParticleSource *source = scene->ParticleSource(i);
    R3Shape *shape = source->shape;
    
    double idealnumtogen = source->rate*delta_time;
    int numtogen = 0;
    numtogen += idealnumtogen;
    
    if(RandomNumber() < idealnumtogen-numtogen)
      numtogen++;
    
    for(int j=0; j<numtogen; j++) {
      R3Particle *particle = new R3Particle();
      
      if(shape->type == R3_SPHERE_SHAPE) {
        R3Sphere *sphere = shape->sphere;
        R3Point center = sphere->Center();
        double r = sphere->Radius();
        
        double z,phi,d;
        z = rand()%2 == 0 ? -RandomNumber()*r : RandomNumber()*r;
        phi = RandomNumber()*2.0*M_PI;
        d = sqrt(r*r - z*z);
        particle->position = R3Point(center[0] + d*cos(phi), center[1] + d*sin(phi), center[2] + z);
      
        R3Vector normal = particle->position - center;
        normal.Normalize();
        R3Vector tangentplanevector = normal;
        tangentplanevector[2] += 1;
        tangentplanevector[0] += 1;
        tangentplanevector[1] += 1;
        tangentplanevector.Cross(normal);
        tangentplanevector.Normalize();
      
        double t1, t2;
        t1 = RandomNumber()*2.0*M_PI;
        t2 = RandomNumber()*sin(source->angle_cutoff);
      
        R3Vector direction = tangentplanevector;
        direction.Rotate(normal, t1);
        R3Vector vcrossn = direction;
        vcrossn.Cross(normal);
        direction.Rotate(vcrossn, acos(t2));
        direction.Normalize();
        particle->velocity = source->velocity*direction;
      }
      else if(shape->type == R3_CIRCLE_SHAPE) {
        R3Circle *circle = shape->circle;
        R3Point center = circle->Center();
        double rmax = circle->Radius();
        
        R3Vector normal = circle->Normal();
        normal.Normalize();
        R3Vector tangentplanevector1 = normal;
        tangentplanevector1[2] += 1;
        tangentplanevector1[0] += 1;
        tangentplanevector1[1] += 1;
        tangentplanevector1.Cross(normal);
        tangentplanevector1.Normalize();
        R3Vector tangentplanevector2 = tangentplanevector1;
        tangentplanevector2.Cross(normal);
        
        double theta = RandomNumber()*2.0*M_PI;
        double r = RandomNumber()*rmax;
        
        particle->position = center + sqrt(r)*cos(theta)*tangentplanevector1 + sqrt(r)*sin(theta)*tangentplanevector2;
      
        double t1, t2;
        t1 = RandomNumber()*2.0*M_PI;
        t2 = RandomNumber()*sin(source->angle_cutoff);
      
        R3Vector direction = tangentplanevector1;
        direction.Rotate(normal, t1);
        R3Vector vcrossn = direction;
        vcrossn.Cross(normal);
        direction.Rotate(vcrossn, acos(t2));
        direction.Normalize();
        if(rand() % 2 == 0)
          particle->velocity = source->velocity*direction;
        else
          particle->velocity = -source->velocity*direction;
      }
      else {
        delete particle;
        continue;
      }
      
      particle->material = source->material;
      particle->mass = source->mass;
      particle->drag = source->drag;
      particle->elasticity = source->elasticity;
      particle->fixed = source->fixed;
      if(particle->fixed) {
        particle->velocity = R3Vector(0,0,0);
      }
      if(source->lifetime > 0) {
        particle->lifetime = source->lifetime + current_time;
      } else {
        particle->lifetime = source->lifetime;
      }
    
      scene->particles.push_back(particle);
    }
  }
}



////////////////////////////////////////////////////////////
// Updating Particles
////////////////////////////////////////////////////////////

void UpdateParticles(R3Scene *scene, double current_time, double delta_time, int integration_type, int mutualattraction)
{
  // Lifetime deletion
  for (int i = 0; i < scene->NParticles(); i++) {
    R3Particle *particle = scene->Particle(i);
    if(particle->lifetime > 0 && particle->lifetime < current_time) {
      delete particle;
      scene->particles.erase(scene->particles.begin() + i);
      i--;
      break;
    }
  }
  
  
  std::vector<R3Vector> accels; 
  // Compute forces
  for (int i = 0; i < scene->NParticles(); i++) {
    R3Particle *particle = scene->Particle(i);
    if(particle->fixed) {
      accels.push_back(R3Vector(0,0,0));
      continue;
    }
    
    // Drag + Gravity
    R3Vector drag = -particle->drag*particle->velocity;
    R3Vector gravity = scene->gravity*particle->mass;
    
    // Sinks
    R3Vector sinks(0,0,0);
    for(int j=0; j < scene->NParticleSinks(); j++) {
      R3ParticleSink* sink = scene->ParticleSink(j);
      R3Shape *shape = sink->shape;
      if(shape->type == R3_SPHERE_SHAPE) {
        R3Sphere *sphere = shape->sphere;
        double d = R3Distance(particle->position, sphere->Center()) - sphere->Radius();
        double magnitude = sink->intensity / (sink->constant_attenuation + sink->linear_attenuation*d + sink->quadratic_attenuation*d*d);
        
        R3Vector v = sphere->Center() - particle->position;
        v.Normalize();
        sinks += v*magnitude;
      }
    }
    
    // Springs
    R3Vector springs(0,0,0);
    for(unsigned int j=0; j < particle->springs.size(); j++) {
      R3ParticleSpring *spring = particle->springs[j];
      R3Particle *particle2;
      if(particle == spring->particles[0])
        particle2 = spring->particles[1];
      else
        particle2 = spring->particles[0];
      
      double d = R3Distance(particle2->position, particle->position);
      R3Vector D = particle2->position - particle->position;
      D.Normalize();
      
      double hooke = spring->ks*(d - spring->rest_length);
      double damphooke = spring->kd*(particle2->velocity - particle->velocity).Dot(D);
      
      springs += (hooke + damphooke) * D;
    }
    
    R3Vector attraction(0,0,0);
    if(mutualattraction == 1) {
      for(int j=0; j < scene->NParticles(); j++) {
        R3Particle *particle2 = scene->Particle(j);
        if(particle == particle2) continue;
      
        double r = R3Distance(particle->position, particle2->position);
        R3Vector rvec = particle2->position - particle->position;
        rvec.Normalize();
        double newton = 6.67E-11*particle->mass*particle2->mass/r/r;
        attraction = newton*rvec;
      }
    }
    
    accels.push_back(drag+gravity+sinks+springs+attraction);
  }
    
  // Update velocity for every particle
  for (int i = 0; i < scene->NParticles(); i++) {
    R3Particle *particle = scene->Particle(i);
    if(particle->fixed) continue;
    particle->velocity += accels[i]*delta_time/particle->mass;
  }
  
  // Update position for every particle
  
  for (int i = 0; i < scene->NParticles(); i++) {
    R3Particle *particle = scene->Particle(i);
    if(particle->fixed) continue;
    
    // Intersection
    double time_left = delta_time;
    for(;;) {
      R3Ray ray(particle->position, particle->velocity);
      int intersects;
      R3Point intersectionpoint;
      R3Node *intersectingnode;
      R3Vector intersectionnormal;
      double t;
      intersects = IntersectScene(scene, ray, &intersectionpoint, &intersectionnormal, &t, &intersectingnode, NULL);
      
      if(intersects == 1) {        
        t /= particle->velocity.Length();
        if(t < time_left) {
          R3Vector perpvector1 = intersectionnormal;
          perpvector1[2] += 1;
          perpvector1.Cross(intersectionnormal);
          perpvector1.Normalize();

          R3Vector perpvector2 = intersectionnormal;
          perpvector2.Cross(perpvector1);

          double vperp1, vperp2, vpar;
          vperp1 = particle->velocity.Dot(perpvector1);
          vperp2 = particle->velocity.Dot(perpvector2);
          vpar = particle->velocity.Dot(intersectionnormal);

          R3Vector v = vperp1*perpvector1 + vperp2*perpvector2 - particle->elasticity*vpar*intersectionnormal;
          time_left -= t;
          particle->position = intersectionpoint;
          particle->velocity = v;
        }
        else {
          particle->position += particle->velocity*time_left;
          break;
        }
      }
      else {
        particle->position += particle->velocity*time_left;
        break;
      }
    }
    
    // Sink Deletion
    for(int j=0; j < scene->NParticleSinks(); j++) {
      R3ParticleSink* sink = scene->ParticleSink(j);
      R3Shape *shape = sink->shape;
      if(shape->type == R3_SPHERE_SHAPE) {
        R3Sphere *sphere = shape->sphere;
        double d = R3Distance(particle->position, sphere->Center());
        if(d < sphere->Radius()) {
          delete particle;
          scene->particles.erase(scene->particles.begin() + i);
          i--;
          break;
        }
      }
    }
  }
}



////////////////////////////////////////////////////////////
// Rendering Particles
////////////////////////////////////////////////////////////

void RenderParticles(R3Scene *scene, double current_time, double delta_time, int dynamic)
{
  // Draw every particle

  glDisable(GL_LIGHTING);
  glPointSize(5);
  glBegin(GL_POINTS);
  for (int i = 0; i < scene->NParticles(); i++) {
    R3Particle *particle = scene->Particle(i);
    
    double scalef = 1;
    if(dynamic == 1 && particle->lifetime > 0 && particle->lifetime - current_time < 10) {
      scalef = (particle->lifetime - current_time)/10;
    }
    glColor3d(scalef*particle->material->kd[0], scalef*particle->material->kd[1], scalef*particle->material->kd[2]);
    const R3Point& position = particle->position;
    glVertex3d(position[0], position[1], position[2]);
  }   
  glEnd();

}



