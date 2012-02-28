// COS 426, Spring 2007, Thomas Funkhouser
// Assignment 2: Mesh Processing



// Include files
#ifdef _WIN32
#include <windows.h>
#endif

#include "R2/R2.h"
#include "R3/R3.h"
#include "R3Mesh.h"



// Program arguments

static char options[] =
"  -help\n"
"  -translate <real:dx> <real:dy> <real:dz>\n"
"  -scale <real:sx> <real:sy> <real:sz>\n"
"  -rotate <real:angle> <real:px> <real:py> <real:pz> <real:vx> <real:vy> <real:vz>\n"
"  -noise <real:factor>\n"
"  -inflate <real:factor>\n"
"  -fun\n"
"  -smooth\n"
"  -sharpen\n"
"  -smoothbilateral\n"
"  -truncate <real:t>\n"
"  -bevel <real:t>\n"
"  -splitfaces\n"
"  -starfaces <real:factor>\n"
"  -splitlongedges <real:max_edge_length>\n (max_edge_length is relative to mesh radius)"
"  -collapesshortedges <real:min_edge_length>\n (min_edge_length is relative to mesh radius)"
"  -vertexcluster <real:grid_cell_size> (grid_cell_size is relative to mesh radius)\n"
"  -bezier <file:controlmesh> <int:m> <int:n>\n"
"  -bspline <file:controlmesh> <int:m> <int:n>\n"
"  -loop\n"
"  -catmullclark\n"
"  -revolution <file:profile> <real:px> <real:py> <real:pz> <real:vx> <real:vy> <real:vz> <real:angle_step>\n"
"  -sweep <file:crosssection> <file:centerline>\n"
"  -fixholes\n"
"  -fixcracks <real:epsilon> (epsilon is relative to mesh radius)\n"
"  -fixintersections\n"
"  -intersect <file:mesh>\n"
"  -subtract <file:mesh>\n"
"  -union <file:mesh>\n"
"  -crop <real:a> <real:b> <real:c> <real:d>\n";



static void 
ShowUsage(void)
{
  // Print usage message and exit
  fprintf(stderr, "Usage: meshpro input_mesh output_mesh [  -option [arg ...] ...]\n");
  fprintf(stderr, "%s", options);
  exit(EXIT_FAILURE);
}



static void 
CheckOption(char *option, int argc, int minargc)
{
  // Check if there are enough remaining arguments for option
  if (argc < minargc)  {
    fprintf(stderr, "Too few arguments for %s\n", option);
    ShowUsage();
    exit(-1);
  }
}



int 
main(int argc, char **argv)
{
  // Look for help
  for (int i = 0; i < argc; i++) {
    if (!strcmp(argv[i], "-help")) {
      ShowUsage();
    }
  }

  // Read input and output mesh filenames
  if (argc < 3)  ShowUsage();
  argv++, argc--; // First argument is program name
  char *input_mesh_name = *argv; argv++, argc--; 
  char *output_mesh_name = *argv; argv++, argc--; 

  // Allocate mesh
  R3Mesh *mesh = new R3Mesh();
  if (!mesh) {
    fprintf(stderr, "Unable to allocate mesh\n");
    exit(-1);
  }

  // Read input mesh
  if (!mesh->Read(input_mesh_name)) {
    fprintf(stderr, "Unable to read mesh from %s\n", input_mesh_name);
    exit(-1);
  }

  // Parse arguments and perform operations 
  while (argc > 0) {
    if (!strcmp(*argv, "-translate")) {
      CheckOption(*argv, argc, 4);
      double dx = atof(argv[1]);
      double dy = atof(argv[2]);
      double dz = atof(argv[3]);
      argv += 4, argc -= 4;
      mesh->Translate(dx, dy, dz);
    }
    else if (!strcmp(*argv, "-scale")) {
      CheckOption(*argv, argc, 4);
      double sx = atof(argv[1]);
      double sy = atof(argv[2]);
      double sz = atof(argv[3]);
      argv += 4, argc -= 4;
      mesh->Scale(sx, sy, sz);
    }
    else if (!strcmp(*argv, "-rotate")) {
      CheckOption(*argv, argc, 8);
      double angle = atof(argv[1]);
      double px = atof(argv[2]);
      double py = atof(argv[3]);
      double pz = atof(argv[4]);
      double vx = atof(argv[5]);
      double vy = atof(argv[6]);
      double vz = atof(argv[7]);
      argv += 8, argc -= 8;
      R3Point p(px, py, pz);
      R3Vector v(vx, vy, vz);
      R3Line axis(p, v);
      mesh->Rotate(angle, axis);
    }
    else if (!strcmp(*argv, "-noise")) {
      CheckOption(*argv, argc, 2);
      double factor = atof(argv[1]);
      argv += 2, argc -= 2;
      mesh->RandomNoise(factor);
    }
    else if (!strcmp(*argv, "-inflate")) {
      CheckOption(*argv, argc, 2);
      double factor = atof(argv[1]);
      argv += 2, argc -= 2;
      mesh->Inflate(factor);
    }
    else if (!strcmp(*argv, "-fun")) {
      CheckOption(*argv, argc, 1);
      argv += 1, argc -= 1;
      mesh->Fun();
    }
    else if (!strcmp(*argv, "-smooth")) {
      CheckOption(*argv, argc, 1);
      argv += 1, argc -= 1;
      mesh->Smooth();
    }
    else if (!strcmp(*argv, "-sharpen")) {
      CheckOption(*argv, argc, 1);
      argv += 1, argc -= 1;
      mesh->Sharpen();
    }
    else if (!strcmp(*argv, "-smoothbilateral")) {
      CheckOption(*argv, argc, 1);
      argv += 1, argc -= 1;
      mesh->SmoothBilateral();
    }
    else if (!strcmp(*argv, "-truncate")) {
      CheckOption(*argv, argc, 2);
      double t = atof(argv[1]);
      argv += 2, argc -= 2;
      mesh->Truncate(t);
    }
    else if (!strcmp(*argv, "-bevel")) {
      CheckOption(*argv, argc, 2);
      double t = atof(argv[1]);
      argv += 2, argc -= 2;
      mesh->Bevel(t);
    }
    else if (!strcmp(*argv, "-splitfaces")) {
      CheckOption(*argv, argc, 1);
      argv += 1, argc -= 1;
      mesh->SplitFaces();
    }
    else if (!strcmp(*argv, "-starfaces")) {
      CheckOption(*argv, argc, 1);
      double factor = atof(argv[1]);
      argv += 2, argc -= 2;
      mesh->StarFaces(factor);
    }
    else if (!strcmp(*argv, "-splitlongedges")) {
      CheckOption(*argv, argc, 2);
      double max_edge_length = mesh->Radius() * atof(argv[1]);
      argv += 2, argc -= 2;
      mesh->SplitLongEdges(max_edge_length);
    }
    else if (!strcmp(*argv, "-collapseshortedges")) {
      CheckOption(*argv, argc, 2);
      double min_edge_length = mesh->Radius() * atof(argv[1]);
      argv += 2, argc -= 2;
      mesh->CollapseShortEdges(min_edge_length);
    }
    else if (!strcmp(*argv, "-vertexcluster")) {
      CheckOption(*argv, argc, 2);
      double grid_cell_size = mesh->Radius() * atof(argv[1]);
      argv += 2, argc -= 2;
      mesh->ClusterVertices(grid_cell_size);
    }
    else if (!strcmp(*argv, "-bezier")) {
      CheckOption(*argv, argc, 4);
      R3Mesh control_mesh;
      if (!control_mesh.Read(argv[1])) exit(-1);
      int m = atoi(argv[2]);
      int n = atoi(argv[3]);
      argv += 4, argc -= 4;
      mesh->Bezier(control_mesh, m, n);
    }
    else if (!strcmp(*argv, "-bspline")) {
      CheckOption(*argv, argc, 4);
      R3Mesh control_mesh;
      if (!control_mesh.Read(argv[1])) exit(-1);
      int m = atoi(argv[2]);
      int n = atoi(argv[3]);
      argv += 4, argc -= 4;
      mesh->BSpline(control_mesh, m, n);
    }
    else if (!strcmp(*argv, "-loop")) {
      CheckOption(*argv, argc, 1);
      argv += 1, argc -= 1;
      mesh->SubdivideLoop();
    }
    else if (!strcmp(*argv, "-catmullclark")) {
      CheckOption(*argv, argc, 1);
      argv += 1, argc -= 1;
      mesh->SubdivideCatmullClark();
    }
    else if (!strcmp(*argv, "-revolution")) {
      CheckOption(*argv, argc, 9);
      R3Mesh profile;
      if (!profile.Read(argv[1])) exit(-1);
      double px = atof(argv[2]);
      double py = atof(argv[3]);
      double pz = atof(argv[4]);
      double vx = atof(argv[5]);
      double vy = atof(argv[6]);
      double vz = atof(argv[7]);
      double angle_step = atof(argv[8]);
      argv += 9, argc -= 9;
      R3Point p(px, py, pz);
      R3Vector v(vx, vy, vz);
      R3Line axis(p, v);
      mesh->SurfaceOfRevolution(profile, axis, angle_step);
    }
    else if (!strcmp(*argv, "-sweep")) {
      CheckOption(*argv, argc, 3);
      R3Mesh crosssection, centerline;
      if (!crosssection.Read(argv[1])) exit(-1);
      if (!centerline.Read(argv[2])) exit(-1);
      argv += 3, argc -= 3;
      mesh->SurfaceSweep(crosssection, centerline);
    }
    else if (!strcmp(*argv, "-fixholes")) {
      CheckOption(*argv, argc, 1);
      argv += 1, argc -= 1;
      mesh->FixHoles();
    }
    else if (!strcmp(*argv, "-fixcracks")) {
      CheckOption(*argv, argc, 2);
      double epsilon = mesh->Radius() * atof(argv[1]);
      argv += 2, argc -= 2;
      mesh->FixCracks(epsilon);
    }
    else if (!strcmp(*argv, "-fixintersections")) {
      CheckOption(*argv, argc, 1);
      argv += 1, argc -= 1;
      mesh->FixIntersections();
    }
    else if (!strcmp(*argv, "-intersect")) {
      CheckOption(*argv, argc, 2);
      R3Mesh mask;
      if (!mask.Read(argv[1])) exit(-1);
      argv += 2, argc -= 2;
      mesh->Intersect(mask);
    }
    else if (!strcmp(*argv, "-subtract")) {
      CheckOption(*argv, argc, 2);
      R3Mesh mask;
      if (!mask.Read(argv[1])) exit(-1);
      argv += 2, argc -= 2;
      mesh->Subtract(mask);
    }
    else if (!strcmp(*argv, "-union")) {
      CheckOption(*argv, argc, 2);
      R3Mesh mask;
      if (!mask.Read(argv[1])) exit(-1);
      argv += 2, argc -= 2;
      mesh->Union(mask);
    }
    else if (!strcmp(*argv, "-crop")) {
      CheckOption(*argv, argc, 2);
      double a = atof(argv[1]);
      double b = atof(argv[2]);
      double c = atof(argv[3]);
      double d = atof(argv[4]);
      argv += 5, argc -= 5;
      R3Plane plane(a, b, c, d);
      mesh->Crop(plane);
    }
    else {
      // Unrecognized program argument
      fprintf(stderr, "meshpro: invalid option: %s\n", *argv);
      ShowUsage();
    }
  }

  // Write output mesh
  if (!mesh->Write(output_mesh_name)) {
    fprintf(stderr, "Unable to write mesh to %s\n", output_mesh_name);
    exit(-1);
  }

  // Delete mesh
  delete mesh;

  // Return success
  return EXIT_SUCCESS;
}



