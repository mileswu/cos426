#
# Makefile for producing Assignment #2 output
#

EXE=src/meshpro

IMGS=output/ico_rotate_0.5.jpg \
     output/ico_scale_1_2_3.jpg \
     output/cow.jpg \
     output/cheetah.jpg \
     output/cow_noise_0.1.jpg \
     output/cow_noise_0.4.jpg \
     output/cow_noise_0.8.jpg \
     output/cow_inflate_0.2.jpg \
     output/cow_inflate_0.8.jpg \
     output/cow_inflate_1.6.jpg \
     output/cow_smooth_1.jpg \
     output/cow_smooth_2.jpg \
     output/cow_smooth_3.jpg \
     output/cow_sharpen_1.jpg \
     output/cow_sharpen_2.jpg \
     output/cow_sharpen_3.jpg \
     output/ico.jpg \
     output/ico_truncate_0.1.jpg \
     output/ico_truncate_0.2.jpg \
     output/ico_truncate_0.4.jpg \
     output/sphere.jpg \
     output/sphere_star_0.1.jpg \
     output/sphere_star_0.2.jpg \
     output/sphere_star_0.4.jpg \
     output/cow_loop.jpg \
     output/cheetah_loop.jpg \
     output/cow_collapse_0.01.jpg \
     output/cow_collapse_0.02.jpg \
     output/cow_collapse_0.04.jpg \
     output/wineglass_surface_0.2.jpg \
     output/wineglass_surface_0.1.jpg \
     output/wineglass_surface_0.01.jpg \
     output/ico_edges.jpg \
     output/ico_splitfaces.jpg \
     output/ico_normals.jpg \
     output/sphere_normals.jpg \
     output/cow_normals.jpg \
     output/movie.webm \
     art/movie.mp4


all: $(IMGS) 

output/ico_rotate_0.5.jpg: input/ico.off $(EXE)
	$(EXE) input/ico.off output/ico_rotate_0.5.off -rotate  0.5  0 0 0  0 0 1
	src/meshview output/ico_rotate_0.5.off -output_image output/ico_rotate_0.5.jpg -exit_immediately

output/ico_scale_1_2_3.jpg: input/ico.off $(EXE)
	$(EXE) input/ico.off output/ico_scale_1_2_3.off -scale 1 2 3
	src/meshview output/ico_scale_1_2_3.off -output_image output/ico_scale_1_2_3.jpg -exit_immediately

output/cow.jpg: input/ico.off $(EXE)
	src/meshview input/cow.off -output_image output/cow.jpg -exit_immediately

output/ico.jpg: input/ico.off $(EXE)
	src/meshview input/ico.off -output_image output/ico.jpg -exit_immediately

output/ico_edges.jpg: input/ico.off $(EXE)
	$(EXE) input/ico.off output/ico_edges.off -rotate -0.1 0 0 0 1 0 0
	src/meshview output/ico_edges.off -output_image output/ico_edges.jpg -exit_immediately -edges

output/sphere.jpg: input/ico.off $(EXE)
	src/meshview input/sphere.off -output_image output/sphere.jpg -exit_immediately

output/cheetah.jpg: input/ico.off $(EXE)
	src/meshview input/cheetah.off -output_image output/cheetah.jpg -exit_immediately

output/cow_noise_0.1.jpg: input/ico.off $(EXE)
	$(EXE) input/cow.off output/cow_noise_0.1.off -noise 0.1
	src/meshview output/cow_noise_0.1.off -output_image output/cow_noise_0.1.jpg -exit_immediately

output/cow_noise_0.4.jpg: input/ico.off $(EXE)
	$(EXE) input/cow.off output/cow_noise_0.4.off -noise 0.4
	src/meshview output/cow_noise_0.4.off -output_image output/cow_noise_0.4.jpg -exit_immediately

output/cow_noise_0.8.jpg: input/ico.off $(EXE)
	$(EXE) input/cow.off output/cow_noise_0.8.off -noise 0.8
	src/meshview output/cow_noise_0.8.off -output_image output/cow_noise_0.8.jpg -exit_immediately

output/cow_inflate_0.2.jpg: input/ico.off $(EXE)
	$(EXE) input/cow.off output/cow_inflate_0.2.off -inflate 0.2
	src/meshview output/cow_inflate_0.2.off -output_image output/cow_inflate_0.2.jpg -exit_immediately

output/cow_inflate_0.8.jpg: input/ico.off $(EXE)
	$(EXE) input/cow.off output/cow_inflate_0.8.off -inflate 0.8
	src/meshview output/cow_inflate_0.8.off -output_image output/cow_inflate_0.8.jpg -exit_immediately

output/cow_inflate_1.6.jpg: input/ico.off $(EXE)
	$(EXE) input/cow.off output/cow_inflate_1.6.off -inflate 1.6
	src/meshview output/cow_inflate_1.6.off -output_image output/cow_inflate_1.6.jpg -exit_immediately

output/cow_smooth_1.jpg: input/ico.off $(EXE)
	$(EXE) input/cow.off output/cow_smooth_1.off -smooth
	src/meshview output/cow_smooth_1.off -output_image output/cow_smooth_1.jpg -exit_immediately

output/cow_smooth_2.jpg: input/ico.off $(EXE)
	$(EXE) input/cow.off output/cow_smooth_2.off -smooth -smooth
	src/meshview output/cow_smooth_2.off -output_image output/cow_smooth_2.jpg -exit_immediately

output/cow_smooth_3.jpg: input/ico.off $(EXE)
	$(EXE) input/cow.off output/cow_smooth_3.off -smooth -smooth -smooth
	src/meshview output/cow_smooth_3.off -output_image output/cow_smooth_3.jpg -exit_immediately

output/cow_sharpen_1.jpg: input/ico.off $(EXE)
	$(EXE) input/cow.off output/cow_sharpen_1.off -sharpen
	src/meshview output/cow_sharpen_1.off -output_image output/cow_sharpen_1.jpg -exit_immediately

output/cow_sharpen_2.jpg: input/ico.off $(EXE)
	$(EXE) input/cow.off output/cow_sharpen_2.off -sharpen -sharpen
	src/meshview output/cow_sharpen_2.off -output_image output/cow_sharpen_2.jpg -exit_immediately

output/cow_sharpen_3.jpg: input/ico.off $(EXE)
	$(EXE) input/cow.off output/cow_sharpen_3.off -sharpen -sharpen -sharpen
	src/meshview output/cow_sharpen_3.off -output_image output/cow_sharpen_3.jpg -exit_immediately

output/ico_truncate_0.1.jpg: input/ico.off $(EXE)
	$(EXE) input/ico.off output/ico_truncate_0.1.off -truncate 0.1
	src/meshview output/ico_truncate_0.1.off -output_image output/ico_truncate_0.1.jpg -exit_immediately

output/ico_truncate_0.2.jpg: input/ico.off $(EXE)
	$(EXE) input/ico.off output/ico_truncate_0.2.off -truncate 0.2
	src/meshview output/ico_truncate_0.2.off -output_image output/ico_truncate_0.2.jpg -exit_immediately

output/ico_truncate_0.4.jpg: input/ico.off $(EXE)
	$(EXE) input/ico.off output/ico_truncate_0.4.off -truncate 0.4
	src/meshview output/ico_truncate_0.4.off -output_image output/ico_truncate_0.4.jpg -exit_immediately

output/sphere_star_0.1.jpg: input/sphere.off $(EXE)
	$(EXE) input/sphere.off output/sphere_star_0.1.off -starfaces 0.1
	src/meshview output/sphere_star_0.1.off -output_image output/sphere_star_0.1.jpg -exit_immediately

output/sphere_star_0.2.jpg: input/sphere.off $(EXE)
	$(EXE) input/sphere.off output/sphere_star_0.2.off -starfaces 0.2
	src/meshview output/sphere_star_0.2.off -output_image output/sphere_star_0.2.jpg -exit_immediately

output/sphere_star_0.4.jpg: input/sphere.off $(EXE)
	$(EXE) input/sphere.off output/sphere_star_0.4.off -starfaces 0.4
	src/meshview output/sphere_star_0.4.off -output_image output/sphere_star_0.4.jpg -exit_immediately

output/cow_loop.jpg: input/sphere.off $(EXE)
	$(EXE) input/cow.off output/cow_loop.off -loop
	src/meshview output/cow_loop.off -output_image output/cow_loop.jpg -exit_immediately

output/cheetah_loop.jpg: input/sphere.off $(EXE)
	$(EXE) input/cheetah.off output/cheetah_loop.off -loop
	src/meshview output/cheetah_loop.off -output_image output/cheetah_loop.jpg -exit_immediately

output/cow_collapse_0.01.jpg: input/cow.off $(EXE)
	$(EXE) input/cow.off output/cow_collapse_0.01.off -collapseshortedges 0.01
	src/meshview output/cow_collapse_0.01.off -output_image output/cow_collapse_0.01.jpg -exit_immediately

output/cow_collapse_0.02.jpg: input/cow.off $(EXE)
	$(EXE) input/cow.off output/cow_collapse_0.02.off -collapseshortedges 0.02
	src/meshview output/cow_collapse_0.02.off -output_image output/cow_collapse_0.02.jpg -exit_immediately

output/cow_collapse_0.04.jpg: input/cow.off $(EXE)
	$(EXE) input/cow.off output/cow_collapse_0.04.off -collapseshortedges 0.04
	src/meshview output/cow_collapse_0.04.off -output_image output/cow_collapse_0.04.jpg -exit_immediately

output/wineglass_surface_0.2.jpg: input/cow.off $(EXE)
	$(EXE) input/empty.off output/wineglass_surface_0.2.off -revolution input/wineglass_curve.off 0 0 0 0 0 1 0.2 -rotate -1.6 0 0 0 1 0 0
	src/meshview output/wineglass_surface_0.2.off -output_image output/wineglass_surface_0.2.jpg -exit_immediately

output/wineglass_surface_0.1.jpg: input/cow.off $(EXE)
	$(EXE) input/empty.off output/wineglass_surface_0.1.off -revolution input/wineglass_curve.off 0 0 0 0 0 1 0.1 -rotate -1.6 0 0 0 1 0 0
	src/meshview output/wineglass_surface_0.1.off -output_image output/wineglass_surface_0.1.jpg -exit_immediately

output/wineglass_surface_0.01.jpg: input/cow.off $(EXE)
	$(EXE) input/empty.off output/wineglass_surface_0.01.off -revolution input/wineglass_curve.off 0 0 0 0 0 1 0.01 -rotate -1.6 0 0 0 1 0 0
	src/meshview output/wineglass_surface_0.01.off -output_image output/wineglass_surface_0.01.jpg -exit_immediately

output/ico_splitfaces.jpg: input/cow.off $(EXE)
	$(EXE) input/ico.off output/ico_splitfaces.off -splitfaces -rotate -0.1 0 0 0 1 0 0
	src/meshview output/ico_splitfaces.off -output_image output/ico_splitfaces.jpg -exit_immediately -edges

output/ico_normals.jpg: input/cow.off $(EXE)
	src/meshview input/ico.off -output_image output/ico_normals.jpg -exit_immediately -normals 

output/sphere_normals.jpg: input/cow.off $(EXE)
	src/meshview input/sphere.off -output_image output/sphere_normals.jpg -exit_immediately -normals 

output/cow_normals.jpg: input/cow.off $(EXE)
	src/meshview input/cow.off -output_image output/cow_normals.jpg -exit_immediately -normals 

art/%.mp4 output/%.webm: input/cow.off $(EXE)
	./make_video.sh

clean:
	rm -f $(IMGS)
	rm -f $(IMGS:.jpg=.off)
