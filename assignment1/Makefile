#
# Makefile for producing Assignment #0 output
#

EXE=src/imgpro

IMGS=output/globos_blackandwhite.jpg output/globos_gamma2.jpg output/globos_gamma0.5.jpg \
		 output/globos_saturation2.5.jpg output/globos_saturation0.5.jpg \
		 output/globos_contrast0.5.jpg output/globos_contrast1.5.jpg \
		 output/globos_extractred.jpg \
		 output/globos_gaussian.jpg output/globos_sharpen.jpg \
		 output/globos_edge.jpg output/globos_motionblur.jpg \
		 output/globos_median_5.jpg output/globos_median_1.jpg \
		 output/globos_bilateral_6_0.15.jpg \
		 output/globos_scale_gaussian0.7.jpg output/globos_scale_nearest0.7.jpg output/globos_scale_bilinear0.7.jpg \
		 output/globos_scale_gaussian1.5.jpg output/globos_scale_nearest1.5.jpg output/globos_scale_bilinear1.5.jpg
		 


all: $(IMGS)

output/globos_blackandwhite.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -blackandwhite

output/globos_gamma2.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -gamma 2
output/globos_gamma0.5.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -gamma 0.5

output/globos_saturation2.5.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -saturation 2.5
output/globos_saturation0.5.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -saturation 0.5

output/globos_contrast0.5.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -contrast 0.5
output/globos_contrast1.5.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -contrast 1.5

output/globos_extractred.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -extract 0

output/globos_gaussian.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -blur 3
output/globos_sharpen.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -sharpen
output/globos_edge.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -edge
output/globos_motionblur.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -motionblur 5

output/globos_median_1.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -median 1
output/globos_median_5.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -median 5

output/globos_bilateral_6_0.15.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -bilateral 6 0.15

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

clean:
	rm -f $(IMGS)
