#
# Makefile for producing Assignment #2 output
#

EXE=src/meshpro

IMGS=output/ico_rotate_0.5.jpg \
     output/ico_scale_1_2_3.jpg 

all: $(IMGS) 

output/ico_rotate_0.5.jpg: input/ico.off $(EXE)
	$(EXE) input/ico.off output/ico_rotate_0.5.off -rotate  0.5  0 0 0  0 0 1
	src/meshview output/ico_rotate_0.5.off -output_image output/ico_rotate_0.5.jpg -exit_immediately

output/ico_scale_1_2_3.jpg: input/ico.off $(EXE)
	$(EXE) input/ico.off output/ico_scale_1_2_3.off -scale 1 2 3
	src/meshview output/ico_scale_1_2_3.off -output_image output/ico_scale_1_2_3.jpg -exit_immediately

clean:
	rm -f $(IMGS)
	rm -f $(IMGS:.jpg=.off)
