// Include file for the R3 scene stuff

#define R3Rgb R2Pixel



// Constant definitions

typedef enum {
  R3_BOX_SHAPE,
  R3_SPHERE_SHAPE,
  R3_CYLINDER_SHAPE,
  R3_CONE_SHAPE,
  R3_MESH_SHAPE,
  R3_SEGMENT_SHAPE,
  R3_CIRCLE_SHAPE,
  R3_NUM_SHAPE_TYPES
} R3ShapeType;

typedef enum {
  R3_DIRECTIONAL_LIGHT,
  R3_POINT_LIGHT,
  R3_SPOT_LIGHT,
  R3_AREA_LIGHT,
  R3_NUM_LIGHT_TYPES
} R3LightType;



// Scene element definitions

struct R3Shape {
  R3ShapeType type;
  R3Box *box;
  R3Sphere *sphere;
  R3Cylinder *cylinder;
  R3Cone *cone;
  R3Mesh *mesh;
  R3Segment *segment;
  R3Circle *circle;
};  

struct R3Material {
  R3Rgb ka;
  R3Rgb kd;
  R3Rgb ks;
  R3Rgb kt;
  R3Rgb emission;
  double shininess;
  double indexofrefraction;
  R2Image *texture;
  int texture_index;
  int id;
};

struct R3Light {
  R3LightType type;
  R3Point position;
  R3Vector direction;
  double radius;
  R3Rgb color;
  double constant_attenuation;
  double linear_attenuation;
  double quadratic_attenuation;
  double angle_attenuation;
  double angle_cutoff;
};

struct R3Camera {
  R3Point eye;
  R3Vector towards;
  R3Vector right;
  R3Vector up;
  double xfov, yfov;
  double neardist, fardist;
};

struct R3Node {
  struct R3Node *parent;
  vector<struct R3Node *> children;
  R3Shape *shape;
  R3Matrix transformation;
  R3Material *material;
  R3Box bbox;
};



// Particle system definitions

struct R3Particle {
  R3Point position;
  R3Vector velocity;
  double mass;
  bool fixed;
  double drag;
  double elasticity;
  double lifetime;
  R3Material *material;
  vector<struct R3ParticleSpring *> springs;
};

struct R3ParticleSource {
  R3Shape *shape;
  double rate;
  double velocity;
  double angle_cutoff;
  double mass;
  bool fixed;
  double drag;
  double elasticity;
  double lifetime;
  R3Material *material;
};

struct R3ParticleSink {
  R3Shape *shape;
  double intensity;
  double constant_attenuation;
  double linear_attenuation;
  double quadratic_attenuation;
};

struct R3ParticleSpring {
  R3Particle *particles[2];
  double rest_length;
  double ks;
  double kd;
};



// Scene graph definition

struct R3Scene {
 public:
  // Constructor functions
  R3Scene(void);

  // Access functions
  R3Node *Root(void) const;
  int NLights(void) const;
  R3Light *Light(int k) const;
  R3Camera& Camera(void);
  R3Box& BBox(void);

  // Particle stuff
  int NParticleSources(void) const;
  R3ParticleSource *ParticleSource(int k) const;
  int NParticleSinks(void) const;
  R3ParticleSink *ParticleSink(int k) const;
  int NParticles(void) const;
  R3Particle *Particle(int k) const;

  // I/O functions
  int Read(const char *filename, R3Node *root = NULL);

 public:
  R3Node *root;
  vector<R3Particle *> particles;
  vector<R3ParticleSource *> particle_sources;
  vector<R3ParticleSink *> particle_sinks;
  vector<R3ParticleSpring *> particle_springs;
  vector<R3Light *> lights;
  R3Vector gravity;
  R3Camera camera;
  R3Box bbox;
  R3Rgb background;
  R3Rgb ambient;
};



// Inline functions 

inline R3Node *R3Scene::
Root(void) const
{
  // Return root node
  return root;
}



inline int R3Scene::
NLights(void) const
{
  // Return number of lights
  return lights.size();
}



inline R3Light *R3Scene::
Light(int k) const
{
  // Return kth light
  return lights[k];
}



inline R3Camera& R3Scene::
Camera(void) 
{
  // Return camera
  return camera;
}



inline R3Box& R3Scene::
BBox(void) 
{
  // Return bounding box 
  return bbox;
}



inline int R3Scene::
NParticleSources(void) const
{
  // Return number of particle sources
  return particle_sources.size();
}



inline R3ParticleSource *R3Scene::
ParticleSource(int k) const
{
  // Return kth particle source
  return particle_sources[k];
}



inline int R3Scene::
NParticleSinks(void) const
{
  // Return number of particle sinks
  return particle_sinks.size();
}



inline R3ParticleSink *R3Scene::
ParticleSink(int k) const
{
  // Return kth particle sink
  return particle_sinks[k];
}



inline int R3Scene::
NParticles(void) const
{
  // Return number of particles
  return particles.size();
}



inline R3Particle *R3Scene::
Particle(int k) const
{
  // Return kth particle 
  return particles[k];
}




