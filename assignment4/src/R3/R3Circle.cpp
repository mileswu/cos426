// Source file for the R3 circle class 



// Include files 

#include "R3.h"



// Public variables 

const R3Circle R3null_circle(R3Point(0.0, 0.0, 0.0), -1.0, R3Vector(0.0, 0.0, 0.0));
const R3Circle R3zero_circle(R3Point(0.0, 0.0, 0.0), 0.0, R3Vector(0.0, 0.0, 1.0));
const R3Circle R3unit_circle(R3Point(0.0, 0.0, 0.0), 1.0, R3Vector(0.0, 0.0, 1.0));



R3Circle::
R3Circle(void)
{
}



R3Circle::
R3Circle(const R3Circle& circle)
  : center(circle.center),
    radius(circle.radius),
    plane(circle.plane)
{
}



R3Circle::
R3Circle(const R3Point& center, double radius, const R3Vector& normal)
  : center(center),
    radius(radius),
    plane(center, normal)
{
}



const double R3Circle::
Area(void) const
{
  // Return area of circle
  return (M_PI * radius * radius);
}



const R3Box R3Circle::
BBox(void) const
{
  // Return bounding box of circle ???
  return R3Box(center.X() - radius, center.Y() - radius, center.Z() - radius,
               center.X() + radius, center.Y() + radius, center.Z() + radius);
}



void R3Circle::
Flip (void) 
{
  // Flip orientation of plane
  plane.Flip();
}



void R3Circle::
Reposition (const R3Point& center) 
{
  // Reposition circle's center
  this->center = center;
  plane.Reposition(center);
}



void R3Circle::
Translate (const R3Vector& offset) 
{
  // Translate circle
  center.Translate(offset);
  plane.Translate(offset);
}



void R3Circle::
Align (const R3Vector& normal) 
{
  // Rotate around center to align with normal
  plane.Reset(center, normal);
}



void R3Circle::
Resize (double radius)
{
  // Set radius
  this->radius = radius;
}



void R3Circle::
Draw(void) const
{
  // Draw circle
  int nsteps = 16;
  R3Vector axis1 = R3zero_vector;
  axis1[plane.Normal().MinDimension()] = 1;
  axis1.Cross(plane.Normal());
  axis1.Normalize();
  R3Vector axis2 = plane.Normal() % axis1;
  axis2.Normalize();
  glNormal3d(plane[0], plane[1], plane[2]);
  glBegin(GL_POLYGON);
  for (int i = 0; i < 16; i++) {
    double angle = i * 2 * M_PI / nsteps;
    double x = radius * cos(angle);
    double y = radius * sin(angle);
    R3Point p = center + x*axis1 + y*axis2;
    glVertex3d(p[0], p[1], p[2]);
  }
  glEnd();
}



void R3Circle::
Outline(void) const
{
  // Outline circle
  int nsteps = 16;
  R3Vector axis1 = R3zero_vector;
  axis1[plane.Normal().MinDimension()] = 1;
  axis1.Cross(plane.Normal());
  axis1.Normalize();
  R3Vector axis2 = plane.Normal() % axis1;
  axis2.Normalize();
  glBegin(GL_LINE_LOOP);
  for (int i = 0; i < 16; i++) {
    double angle = i * 2 * M_PI / nsteps;
    double x = radius * cos(angle);
    double y = radius * sin(angle);
    R3Point p = center + x*axis1 + y*axis2;
    glVertex3d(p[0], p[1], p[2]);
  }
  glEnd();
}



void R3Circle::
Print(FILE *fp) const
{
  // Print parameters
  fprintf(fp, "(%g %g %g) (%g %g %g) %g", center[0], center[1], center[2], plane[0], plane[1], plane[2], radius);
}


