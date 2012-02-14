#
# Makefile for producing Assignment #0 output
#

EXE=src/imgpro

IMGS=output/globos_blackandwhite.jpg output/globos_gamma2.jpg output/globos_gamma0.5.jpg \
		 output/globos_saturation2.5.jpg output/globos_saturation0.5.jpg \
		 output/globos_contrast0.5.jpg output/globos_contrast1.5.jpg

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

clean:
	rm -f $(IMGS)
