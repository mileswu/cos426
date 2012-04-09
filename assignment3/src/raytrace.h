// Include file for ray tracing code

R2Image *RenderImage(R3Scene *scene, int width, int height, int max_depth, 
  int num_primary_rays_per_pixel, int num_distributed_rays_per_intersection, int hardshadows_eanbled);
