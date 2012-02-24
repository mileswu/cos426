#
# Makefile for producing Assignment #0 output
#

EXE=src/imgpro

IMGS=output/globos_blackandwhite.jpg \
		 output/globos_gamma2.jpg output/globos_gamma0.5.jpg output/globos_gamma1.jpg \
		 output/globos_saturation2.5.jpg output/globos_saturation0.5.jpg output/globos_saturation1.jpg \
		 output/globos_contrast0.5.jpg output/globos_contrast1.5.jpg output/globos_contrast1.jpg\
		 output/globos_extractred.jpg output/globos_extractgreen.jpg output/globos_extractblue.jpg\
		 output/globos_gaussian3.jpg output/globos_gaussian6.jpg output/globos_gaussian0.jpg  \
		 output/globos_sharpen.jpg \
		 output/globos_edge.jpg output/globos_motionblur.jpg output/globos_sobel.jpg \
		 output/globos_median_5.jpg output/globos_median_2.jpg output/globos_median_0.jpg \
		 output/globos_bilateral_2_0.1.jpg output/globos_bilateral_2_0.3.jpg output/globos_bilateral_6_0.3.jpg output/globos_bilateral_6_0.1.jpg \
		 output/globos_scale_gaussian0.7.jpg output/globos_scale_nearest0.7.jpg output/globos_scale_bilinear0.7.jpg \
		 output/globos_scale_gaussian1.5.jpg output/globos_scale_nearest1.5.jpg output/globos_scale_bilinear1.5.jpg \
		 output/globos_rotate_nearest_5.3.jpg output/globos_rotate_bilinear_5.3.jpg output/globos_rotate_gaussian_5.3.jpg \
		 output/globos_rotate_nearest_0.3.jpg output/globos_rotate_bilinear_0.3.jpg output/globos_rotate_gaussian_0.3.jpg \
		 output/globos_quantize_4.jpg output/globos_quantize_3.jpg output/globos_quantize_2.jpg output/globos_quantize_1.jpg \
		 output/globos_randomdither_4.jpg output/globos_randomdither_3.jpg output/globos_randomdither_2.jpg output/globos_randomdither_1.jpg \
		 output/globos_ordereddither_4.jpg output/globos_ordereddither_3.jpg output/globos_ordereddither_2.jpg output/globos_ordereddither_1.jpg \
		 output/globos_floyddither_4.jpg output/globos_floyddither_3.jpg output/globos_floyddither_2.jpg output/globos_floyddither_1.jpg \
		 output/globos_fun.jpg \
		 output/globos_composite.jpg \
		 output/morph_0.jpg output/morph_2.jpg output/morph_4.jpg output/morph_6.jpg output/morph_8.jpg output/morph_10.jpg \
		 output/morph.webm


all: $(IMGS)

output/globos_blackandwhite.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -blackandwhite

output/globos_gamma2.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -gamma 2
output/globos_gamma0.5.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -gamma 0.5
output/globos_gamma1.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -gamma 1

output/globos_saturation2.5.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -saturation 2.5
output/globos_saturation0.5.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -saturation 0.5
output/globos_saturation1.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -saturation 1

output/globos_contrast0.5.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -contrast 0.5
output/globos_contrast1.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -contrast 1.0
output/globos_contrast1.5.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -contrast 1.5

output/globos_extractred.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -extract 0
output/globos_extractblue.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -extract 2
output/globos_extractgreen.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -extract 1

output/globos_gaussian3.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -blur 3
output/globos_gaussian0.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -blur 0
output/globos_gaussian6.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -blur 6

output/globos_sharpen.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sharpen
output/globos_edge.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -edge
output/globos_sobel.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sobel
output/globos_motionblur.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -motionblur 5

output/globos_median_0.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -median 0
output/globos_median_2.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -median 2
output/globos_median_5.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -median 5

output/globos_bilateral_6_0.3.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -bilateral 6 0.3
output/globos_bilateral_6_0.1.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -bilateral 6 0.1
output/globos_bilateral_2_0.3.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -bilateral 2 0.3
output/globos_bilateral_2_0.1.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -bilateral 2 0.1

output/globos_scale_gaussian0.7.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sampling 2 -scale 0.7 0.7
output/globos_scale_bilinear0.7.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sampling 1 -scale 0.7 0.7
output/globos_scale_nearest0.7.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sampling 0 -scale 0.7 0.7
output/globos_scale_gaussian1.5.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sampling 2 -scale 1.5 1.5
output/globos_scale_bilinear1.5.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sampling 1 -scale 1.5 1.5
output/globos_scale_nearest1.5.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sampling 0 -scale 1.5 1.5


output/globos_rotate_nearest_5.3.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sampling 0 -rotate 5.3
output/globos_rotate_bilinear_5.3.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sampling 1 -rotate 5.3
output/globos_rotate_gaussian_5.3.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sampling 2 -rotate 5.3
output/globos_rotate_nearest_0.3.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sampling 0 -rotate 0.3
output/globos_rotate_bilinear_0.3.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sampling 1 -rotate 0.3
output/globos_rotate_gaussian_0.3.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sampling 2 -rotate 0.3
	
output/globos_quantize_4.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -quantize 4
output/globos_quantize_3.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -quantize 3
output/globos_quantize_2.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -quantize 2
output/globos_quantize_1.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -quantize 1

output/globos_randomdither_4.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -dither 0 4
output/globos_randomdither_3.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -dither 0 3
output/globos_randomdither_2.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -dither 0 2
output/globos_randomdither_1.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -dither 0 1

output/globos_ordereddither_4.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -dither 1 4
output/globos_ordereddither_3.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -dither 1 3
output/globos_ordereddither_2.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -dither 1 2
output/globos_ordereddither_1.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -dither 1 1

output/globos_floyddither_4.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -dither 2 4
output/globos_floyddither_3.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -dither 2 3
output/globos_floyddither_2.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -dither 2 2
output/globos_floyddither_1.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -dither 2 1

output/globos_fun.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -fun

output/globos_composite.jpg: input/globos_de_colores_resized.jpg
	$(EXE) $^ $@ -composite input/islas_del_rosario_resized.jpg input/bottommask.jpg input/topmask.jpg 0

output/morph_0.jpg: input/morph1.jpg
	$(EXE) $^ $@ -sampling 1 -morph input/morph2.jpg input/morph_correspondences.txt 0
output/morph_2.jpg: input/morph1.jpg
	$(EXE) $^ $@ -sampling 1 -morph input/morph2.jpg input/morph_correspondences.txt 0.2
output/morph_4.jpg: input/morph1.jpg
	$(EXE) $^ $@ -sampling 1 -morph input/morph2.jpg input/morph_correspondences.txt 0.4
output/morph_6.jpg: input/morph1.jpg
	$(EXE) $^ $@ -sampling 1 -morph input/morph2.jpg input/morph_correspondences.txt 0.6
output/morph_8.jpg: input/morph1.jpg
	$(EXE) $^ $@ -sampling 1 -morph input/morph2.jpg input/morph_correspondences.txt 0.8
output/morph_10.jpg: input/morph1.jpg
	$(EXE) $^ $@ -sampling 1 -morph input/morph2.jpg input/morph_correspondences.txt 1

output/morph.webm:
	./make_video.sh

clean:
	rm -f $(IMGS)
