// Source file for raytracing code


// Include files

#include "R2/R2.h"
#include "R3/R3.h"
#include "R3Scene.h"
#include "raytrace.h"
#include "float.h"

using namespace std;

R3Ray ConstructRay(R3Camera c, int x, int y, int width, int height) {
//	c.xfov = 1;
//	c.yfov = 1;
	R3Point pcenter = c.eye + c.towards * c.neardist;
	R3Vector vright = c.right * c.neardist * tan(c.xfov);
	R3Vector vup = c.up * c.neardist * tan(c.yfov);

	double halfwidth = (double) width / 2.0;
	double halfheight = (double) height / 2.0;
	R3Point p = pcenter + ((double)x - halfwidth)/halfwidth * vright + ((double)y - halfheight)/halfheight * vup;
	
	R3Vector v = p - c.eye;
	v.Normalize();
	return R3Ray(c.eye, v);
}

int IntersectPlane(R3Plane *p, R3Ray r, R3Point *position, R3Vector *normal, double *t) {
	double denom = r.Vector().Dot(p->Normal());
	if(denom == 0) {
		return 0;
	} else {
		*t = - (r.Start().Vector().Dot(p->Normal()) + p->D()) / denom;
		*position = r.Start() + *t * r.Vector();
		*normal = p->Normal();
		return 1;
	}
}

int IntersectCylinder(R3Cylinder *c, R3Ray r, R3Point *position, R3Vector *normal, double *t) {
	R3Vector v = r.Vector();
	R3Point s = r.Start();
	*t = DBL_MAX;

	// check cylinder sides on infinite
	double eqa = v[0]*v[0] + v[2]*v[2];
	double eqb = 2.0*v[0]*(s[0] - c->Center()[0]) + 2.0*v[2]*(s[2] - c->Center()[2]);
	double eqc = pow(s[0] - c->Center()[0], 2.0) + pow(s[2] - c->Center()[2], 2.0) - c->Radius()*c->Radius();

	double t1, t2;
	t1 = (-eqb + sqrt(eqb*eqb - 4.0*eqa*eqc))/2.0/eqa;
	t2 = (-eqb - sqrt(eqb*eqb - 4.0*eqa*eqc))/2.0/eqa;

	R3Point possible1, possible2;
	possible1 = r.Start() + r.Vector()*t1;
	possible2 = r.Start() + r.Vector()*t2;

	// Check if in y range
	if(possible1[1] < c->Center()[1] - c->Height()/2.0 || possible1[1] > c->Center()[1] + c->Height()/2.0) {
		t1 = -1;
	}
	if(possible2[1] < c->Center()[1] - c->Height()/2.0 || possible2[1] > c->Center()[1] + c->Height()/2.0) {
		t2 = -1;
	}
	if(t1 > 0 && t2 > 0) {
		if(t1 > t2) {
			*t = t2;
			*position = possible2;
			*normal = R3Vector(possible2[0] - c->Center()[0], 0, possible2[2] - c->Center()[2]);
			normal->Normalize();
		}
		else {
			*t = t1;
			*position = possible1;
			*normal = R3Vector(possible1[0] - c->Center()[0], 0, possible1[2] - c->Center()[2]);
			normal->Normalize();
		}
	}
	else if(t1 > 0) {
		*t = t1;
		*position = possible1;
		*normal = R3Vector(possible1[0] - c->Center()[0], 0, possible1[2] - c->Center()[2]);
		normal->Normalize();
	}
	else if(t2 > 0) {
		*t = t2;
		*position = possible2;
		*normal = R3Vector(possible2[0] - c->Center()[0], 0, possible2[2] - c->Center()[2]);
		normal->Normalize();
	}

	//Check endcaps
	for(int i=0; i<2; i++) {
		R3Vector norm_y(0, i == 0 ? 1 : -1, 0);
		double sign = (i == 0 ? 1.0 : -1.0);
		R3Point planepoint(c->Center()[0], c->Center()[1] + sign*c->Height()/2.0, c->Center()[2]);
		R3Plane plane(planepoint, norm_y);
		
		R3Point endcap_p;
		double endcap_t;
		if(IntersectPlane(&plane, r, &endcap_p, &norm_y, &endcap_t) == 1) {
			if(pow(endcap_p[0] - c->Center()[0], 2.0) + pow(endcap_p[2] - c->Center()[2], 2.0) < c->Radius()*c->Radius() && endcap_t > 0) {
				if(endcap_t < *t) {
					*t = endcap_t;
					*position = endcap_p;
					*normal = norm_y;
				}
			}
		}
	}

	if(*t < DBL_MAX) {
		return 1;
	}
	else {
		return 0;
	}
}

int IntersectCone(R3Cone *c, R3Ray r, R3Point *position, R3Vector *normal, double *t) {
	R3Vector v = r.Vector();
	R3Point s = r.Start();
	*t = DBL_MAX;

	double k = pow(c->Radius() / c->Height(), 2.0);
	double y_apex = c->Center()[1] + c->Height() / 2.0;
	// check cone side
	double eqa = v[0]*v[0] + v[2]*v[2] - k*v[1]*v[1];
	double eqb = 2.0*v[0]*(s[0] - c->Center()[0]) + 2.0*v[2]*(s[2] - c->Center()[2]) - 2.0*k*v[1]*(s[1] - y_apex);
	double eqc = pow(s[0] - c->Center()[0], 2.0) + pow(s[2] - c->Center()[2], 2.0) - k*pow(s[1] - y_apex, 2.0);

	double t1, t2;
	t1 = (-eqb + sqrt(eqb*eqb - 4.0*eqa*eqc))/2.0/eqa;
	t2 = (-eqb - sqrt(eqb*eqb - 4.0*eqa*eqc))/2.0/eqa;

	R3Point possible1, possible2;
	possible1 = r.Start() + r.Vector()*t1;
	possible2 = r.Start() + r.Vector()*t2;

	// Check if in y range
	if(possible1[1] < c->Center()[1] - c->Height()/2.0 || possible1[1] > c->Center()[1] + c->Height()/2.0) {
		t1 = -1;
	}
	if(possible2[1] < c->Center()[1] - c->Height()/2.0 || possible2[1] > c->Center()[1] + c->Height()/2.0) {
		t2 = -1;
	}
	if(t1 > 0 && t2 > 0) {
		if(t1 > t2) {
			*t = t2;
			*position = possible2;
			*normal = R3Vector(possible2[0] - c->Center()[0], 0, possible2[2] - c->Center()[2]);
			normal->Normalize();
			(*normal)[1] = k;
		}
		else {
			*t = t1;
			*position = possible1;
			*normal = R3Vector(possible1[0] - c->Center()[0], 0, possible1[2] - c->Center()[2]);
			normal->Normalize();
			(*normal)[1] = k;
		}
	}
	else if(t1 > 0) {
		*t = t1;
		*position = possible1;
		*normal = R3Vector(possible1[0] - c->Center()[0], 0, possible1[2] - c->Center()[2]);
		normal->Normalize();
		(*normal)[1] = k;
	}
	else if(t2 > 0) {
		*t = t2;
		*position = possible2;
		*normal = R3Vector(possible2[0] - c->Center()[0], 0, possible2[2] - c->Center()[2]);
		normal->Normalize();
		(*normal)[1] = k;
	}

	//Check endcaps
	R3Vector norm_y(0, -1, 0);
	R3Point planepoint(c->Center()[0], c->Center()[1] - c->Height()/2.0, c->Center()[2]);
	R3Plane plane(planepoint, norm_y);
	
	R3Point endcap_p;
	double endcap_t;
	if(IntersectPlane(&plane, r, &endcap_p, &norm_y, &endcap_t) == 1) {
		if(pow(endcap_p[0] - c->Center()[0], 2.0) + pow(endcap_p[2] - c->Center()[2], 2.0) < c->Radius()*c->Radius() && endcap_t > 0) {
			if(endcap_t < *t) {
				*t = endcap_t;
				*position = endcap_p;
				*normal = norm_y;
			}
		}
	}

	if(*t < DBL_MAX) {
		return 1;
	}
	else {
		return 0;
	}
}

int IntersectSphere(R3Sphere *s, R3Ray r, R3Point *position, R3Vector *normal, double *t) {
	R3Vector l = s->Center() - r.Start();
	double tca = l.Dot(r.Vector());
	if(tca < 0) return 0;

	double d2 = l.Dot(l) - tca*tca;
	if(d2 > s->Radius() * s->Radius()) return 0;

	double thc = sqrt(s->Radius() * s->Radius() - d2);
	*t = thc > 0 ? tca - thc : tca + thc;
	if(*t < 0) {
		*t = thc > 0 ? tca + thc : tca - thc;
	}
	if(*t < 0) {
		return 0;
	}
	
	*position = r.Start() + *t * r.Vector();
	*normal = *position - s->Center();
	normal->Normalize();

	return 1;
}


int IntersectBox(R3Box *b, R3Ray r, R3Point *position, R3Vector *normal, double *t) {
	*t = DBL_MAX;
	for(int i=0; i<3; i++) { // each direction
		for(int j=0; j<2; j++) { 
			R3Vector facenormal(0,0,0);
			facenormal.SetCoord(i, j == 0 ? -1.0 : 1.0);
			R3Point pointonface;
			if(i == 0)
				pointonface = b->Corner(j == 0 ? 0 : 1, 0, 0);
			else if(i == 1)
				pointonface = b->Corner(0, j == 0 ? 0 : 1, 0);
			else if(i == 2)
				pointonface = b->Corner(0, 0, j == 0 ? 0 : 1);

			R3Point intersectposition;
			double t_intersect;
			R3Plane plane(pointonface, facenormal);
			if(IntersectPlane(&plane, r, &intersectposition, &facenormal, &t_intersect) != 0) {
				// check inside rectangle
				int withinface = 1;
				for(int i2=0; i2<3; i2++) {
					if(i2 == i) continue; //we know it's on the plane
					double min = b->Min()[i2] > b->Max()[i2] ? b->Max()[i2] : b->Min()[i2];
					double max = b->Min()[i2] > b->Max()[i2] ? b->Min()[i2] : b->Max()[i2];
					if(intersectposition[i2] < min || intersectposition[i2] > max) {
						withinface = 0;
						break;
					}
				}
				
				if(withinface == 1 && t_intersect > 0 && t_intersect < *t) {
					*t = t_intersect;
					*position = intersectposition;
					*normal = facenormal;
				}
			}
		}
	}
	if(*t < DBL_MAX) {
		return 1;
	}
	else {
		return 0;
	}
}

int IntersectMesh(R3Mesh *m, R3Ray r, R3Point *position, R3Vector *normal, double *t) {
	*t = DBL_MAX;
	for(int i=0; i < m->NFaces(); i++) {
		R3MeshFace *f = m->Face(i);
		if(f->vertices.size() != 3) continue;
		
		R3Vector trianglenormal = (f->vertices[1]->position - f->vertices[0]->position);
		trianglenormal.Cross(f->vertices[2]->position - f->vertices[0]->position);
		trianglenormal.Normalize();

		R3Plane triangleplane(f->vertices[0]->position, trianglenormal);

		R3Point intersectionpoint;
		double t_intersection;
		if(IntersectPlane(&triangleplane, r, &intersectionpoint, &trianglenormal, &t_intersection) !=0 ) {
			// check inside triangle
			int withintriangle = 1;
			for(int j=0; j<3; j++) {
				R3Vector v1 = f->vertices[j%3]->position - r.Start();
				R3Vector v2 = f->vertices[(j+1)%3]->position - r.Start();
				R3Vector n1 = v2;
				n1.Cross(v1);
				n1.Normalize();
				R3Plane p(r.Start(), n1);
				if(R3SignedDistance(p, intersectionpoint) < 0) {
					withintriangle = 0;
					break;
				}
			}
			
			if(withintriangle == 1 && t_intersection > 0 && t_intersection < *t) {
				*t = t_intersection;
				*position = intersectionpoint;
				*normal = trianglenormal;
			}
		}

	}

	if(*t < DBL_MAX) {
		return 1;
	}
	else {
		return 0;
	}
}

int IntersectNode(R3Node *node, R3Ray r, R3Point *position, R3Vector *normal, double *t, R3Node **intersectingnode, R3Node *excludenode) {
	*t = DBL_MAX;
	R3Ray orig_r = r;
	R3Point intersectionpoint;
	R3Vector intersectionnormal;
	double t_intersection;

	R3Matrix tmatrix = node->transformation;
	R3Matrix tinvmatrix = tmatrix.Inverse();
	r.Transform(tinvmatrix);

	if(node->shape != NULL && excludenode != node) {
		R3Shape *shape = node->shape;
		int intersects = -1;
		if(shape->type == R3_SPHERE_SHAPE) {
			R3Sphere *s = shape->sphere;
			intersects = IntersectSphere(s, r, &intersectionpoint, &intersectionnormal, &t_intersection);
		}
		else if(shape->type == R3_BOX_SHAPE) {
			R3Box *b = shape->box;
			intersects = IntersectBox(b, r, &intersectionpoint, &intersectionnormal, &t_intersection);
		}
		else if(shape->type == R3_MESH_SHAPE) {
			R3Mesh *m = shape->mesh;
			intersects = IntersectMesh(m, r, &intersectionpoint, &intersectionnormal, &t_intersection);
		}
		else if(shape->type == R3_CYLINDER_SHAPE) {
			R3Cylinder *c = shape->cylinder;
			intersects = IntersectCylinder(c, r, &intersectionpoint, &intersectionnormal, &t_intersection);
		}
		else if(shape->type == R3_CONE_SHAPE) {
			R3Cone *c = shape->cone;
			intersects = IntersectCone(c, r, &intersectionpoint, &intersectionnormal, &t_intersection);
		}

		if(intersects == 1) {
			if(t_intersection > 0 && t_intersection < *t) {
				*t = t_intersection;
				*position = intersectionpoint;
				*normal = intersectionnormal;
				*intersectingnode = node;
			}
		}
	}

	for(unsigned int i=0; i<node->children.size(); i++) {
		R3Node *child = node->children[i];
		R3Node *temp_intersectingnode;
		
		if(IntersectBox(&(child->bbox), r, &intersectionpoint, &intersectionnormal, &t_intersection) == 1) {
      if(t_intersection > *t) continue;
		}
		else {
      continue;
		}
		
		if(IntersectNode(child, r, &intersectionpoint, &intersectionnormal, &t_intersection, &temp_intersectingnode, excludenode) == 1) {
			if(t_intersection > 0 && t_intersection < *t) {
				*t = t_intersection;
				*position = intersectionpoint;
				*normal = intersectionnormal;
				*intersectingnode = temp_intersectingnode;
			}
		}
	}

	if(*t < DBL_MAX) {
		normal->Transform(tmatrix);
		normal->Normalize();
		position->Transform(tmatrix);

		if(orig_r.Vector()[0] != 0)
			*t = (*position - orig_r.Start())[0] / orig_r.Vector()[0];
		else if(orig_r.Vector()[1] != 0)
			*t = (*position - orig_r.Start())[1] / orig_r.Vector()[1];
		else
			*t = (*position - orig_r.Start())[2] / orig_r.Vector()[2];

		return 1;
	}
	else {
		return 0;
	}
}

int IntersectScene(R3Scene *scene, R3Ray r, R3Point *position, R3Vector *normal, double *t, R3Node **intersectingnode, R3Node *excludenode) {
	return(IntersectNode(scene->Root(), r, position, normal, t, intersectingnode, excludenode));
}

R3Rgb ComputeRadiance(R3Scene *scene, R3Ray r, int max_depth, R3Node *excludenode, int hardshadows_enabled) {
	R3Point intersectionpoint;
	R3Vector intersectionnormal;
	double t;
	R3Node *n;

	if(IntersectScene(scene, r, &intersectionpoint, &intersectionnormal, &t, &n, excludenode) != 0) {

		R3Rgb emission = n->material->emission;
		R3Rgb ambient = n->material->ka * scene->ambient;
		R3Rgb diffuse(0,0,0,1);
		R3Rgb specular(0,0,0,1);
		R3Rgb reflection(0,0,0,1);
		R3Rgb transmission(0,0,0,1);

		for(unsigned int i=0; i<scene->lights.size(); i++) {
			R3Rgb light = scene->lights[i]->color;

			R3Vector l;
			if(scene->lights[i]->type == R3_DIRECTIONAL_LIGHT) {
				l = -scene->lights[i]->direction;
			}
			else if(scene->lights[i]->type == R3_POINT_LIGHT) {
				l = scene->lights[i]->position - intersectionpoint;
				double dist = sqrt(l.Dot(l));
				l.Normalize();
				light *= 1.0 / (scene->lights[i]->constant_attenuation + scene->lights[i]->linear_attenuation * dist + scene->lights[i]->quadratic_attenuation * dist * dist);
				//printf("%f\n", dist);
			} else if(scene->lights[i]->type == R3_SPOT_LIGHT) {
				l = scene->lights[i]->position - intersectionpoint;
				double dist = sqrt(l.Dot(l));
				l.Normalize();

				double theta = acos((-l).Dot(scene->lights[i]->direction));
				if(theta <= scene->lights[i]->angle_cutoff) {
					light *= 1.0 / (scene->lights[i]->constant_attenuation + scene->lights[i]->linear_attenuation * dist + scene->lights[i]->quadratic_attenuation * dist * dist);
					light *= pow(cos(theta), scene->lights[i]->angle_attenuation);
				}
				else {
					l *= 0;
				}
			} else if(scene->lights[i]->type == R3_AREA_LIGHT) {
				l = scene->lights[i]->position - intersectionpoint;
				double dist = sqrt(l.Dot(l));
				
        l = -scene->lights[i]->direction;
				
				light *= 1.0 / (scene->lights[i]->constant_attenuation + scene->lights[i]->linear_attenuation * dist + scene->lights[i]->quadratic_attenuation * dist * dist);
			}

			R3Plane surfplane(intersectionpoint, intersectionnormal);
			R3Vector reflected = l;
			reflected.Mirror(surfplane);
			reflected = -reflected;

			R3Vector v = r.Start() - intersectionpoint;
			v.Normalize();


			if(hardshadows_enabled == 1) {
				R3Ray shadowray(intersectionpoint, l);
				int blocked;
				R3Point shadowipoint; R3Vector shadowinormal; double shadow_t; R3Node *shadownode;
				blocked = IntersectScene(scene, shadowray, &shadowipoint, &shadowinormal, &shadow_t, &shadownode, n);
				if(blocked == 1 && shadow_t > 0 ) {
					if(shadow_t < shadowray.T(scene->lights[i]->position));
						continue;
				}
			}

			if(intersectionnormal.Dot(l) > 0) {
				diffuse += n->material->kd * light * intersectionnormal.Dot(l);
				if(v.Dot(reflected) > 0) {
					specular += n->material->ks * light * pow(v.Dot(reflected), n->material->shininess);
				}
			}

		}
		if(max_depth > 0) {
			//Reflection
			R3Vector v = r.Start() - intersectionpoint;
			R3Plane surfplane(intersectionpoint, intersectionnormal);
			R3Vector reflected = v;
			reflected.Mirror(surfplane);
			reflected.Normalize();
			reflected = -reflected;
			R3Ray r_reflection(intersectionpoint, reflected);
			
			reflection += n->material->ks * ComputeRadiance(scene, r_reflection, max_depth - 1, n, hardshadows_enabled);


			//Transmission
			R3Rgb blank(0,0,0,1);
			if(n->material->kt != blank) {
				R3Ray r_transmission(intersectionpoint, r.Vector());
				transmission += n->material->kt * ComputeRadiance(scene, r_transmission, max_depth - 1, n, hardshadows_enabled);
			}
		}

		return (emission + ambient + diffuse + specular + reflection + transmission);
	} else {
		return scene->background;
	}
}


////////////////////////////////////////////////////////////////////////
// Create image from scene
//
// This is the main ray tracing function called from raypro
// 
// "width" and "height" indicate the size of the ray traced image
//   (keep these small during debugging to speed up your code development cycle)
//
// "max_depth" indicates the maximum number of secondary reflections/transmissions to trace for any ray
//   (i.e., stop tracing a ray if it has already been reflected max_depth times -- 
//   0 means direct illumination, 1 means one bounce, etc.)
//
// "num_primary_rays_per_pixel" indicates the number of random rays to generate within 
//   each pixel during antialiasing.  This argument can be ignored if antialiasing is not implemented.
//
// "num_distributed_rays_per_intersection" indicates the number of secondary rays to generate
//   for each surface intersection if distributed ray tracing is implemented.  
//   It can be ignored otherwise.
// 
////////////////////////////////////////////////////////////////////////

R2Image *RenderImage(R3Scene *scene, int width, int height, int max_depth,
  int num_primary_rays_per_pixel, int num_distributed_rays_per_intersection, int hardshadows_enabled, int progress)
{
  // Allocate  image
  R2Image *image = new R2Image(width, height);
  if (!image) {
    fprintf(stderr, "Unable to allocate image\n");
    return NULL;
  }

	int counter = 0;
	for(int i=0; i<width; i++) {
		for(int j=0; j<height; j++) {
			counter++;
			if(progress == 1) 
				printf("Progress: %f percent\n", (double)counter*100.0/((double)width*(double)height));
			R3Ray ray = ConstructRay(scene->camera, i, j, width, height);
			R3Rgb radiance = ComputeRadiance(scene, ray, max_depth, NULL, hardshadows_enabled);
			image->SetPixel(i, j, radiance);
		}
	}

  // Return image
  return image;
}
