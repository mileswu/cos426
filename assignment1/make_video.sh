#!/bin/bash

echo "You need ffmpeg"

for i in `seq 0 100`; do
	T=`echo "scale=4;${i}.0/100.0" | bc`
	echo "$i/100"
	./src/imgpro input/morph1.jpg morphvideo_$i.jpg -sampling 1 -morph input/morph2.jpg input/morph_correspondences.txt $T
done

ffmpeg -i morphvideo_%d.jpg art/morph.mp4
ffmpeg -i morphvideo_%d.jpg output/morph.webm

rm morphvideo*.jpg
