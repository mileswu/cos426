#
# Makefile for producing Assignment #0 output
#

EXE=src/imgpro

IMGS=output/globos_blackandwhite.jpg output/globos_gamma2.jpg output/globos_gamma0.5.jpg

all: $(IMGS)

output/globos_blackandwhite.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -blackandwhite

output/globos_gamma2.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -gamma 2
output/globos_gamma0.5.jpg: input/globos_de_colores.jpg
	$(EXE) $^ $@ -gamma 0.5

clean:
	rm -f $(IMGS)
