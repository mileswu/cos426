// Include file for the R3 circle class 



// Class definition 

class R3Circle {
public:
  // Constructor functions
  R3Circle(void);
  R3Circle(const R3Circle& circle);
  R3Circle(const R3Point& center, double radius, const R3Vector& normal);

  // Circle property functions/operators
  const R3Point& Center(void) const;
  const double Radius(void) const;
  const R3Plane& Plane(void) const;
  const R3Vector& Normal(void) const;
  const bool IsEmpty(void) const;

  // Shape property functions/operators
  virtual const double Area(void) const;
  virtual const R3Box BBox(void) const;

  // Manipulation functions/operators
  virtual void Flip(void);
  virtual void Reposition(const R3Point& center);
  virtual void Translate(const R3Vector& offset);
  virtual void Align(const R3Vector& normal);
  virtual void Resize(double radius);

  // Draw functions/operators
  virtual void Draw(void) const;
  virtual void Outline(void) const;
  void Print(FILE *fp = stdout) const;

private:
  R3Point center;
  double radius;
  R3Plane plane;
};



// Public variables 

extern const R3Circle R3null_circle;
extern const R3Circle R3zero_circle;
extern const R3Circle R3unit_circle;



// Inline functions 

inline const R3Point& R3Circle::
Center(void) const
{
  // Return circle center
  return center;
}



inline const double R3Circle::
Radius(void) const
{
  // Return circle radius
  return radius;
}



inline const R3Plane& R3Circle::
Plane(void) const
{
  // Return plane containing circle
  return plane;
}



inline const R3Vector& R3Circle::
Normal(void) const
{
  // Return circle normal
  return plane.Normal();
}



inline const bool R3Circle::
IsEmpty(void) const
{
  // Return whether the sphere is null
  return (radius < 0.0);
}



