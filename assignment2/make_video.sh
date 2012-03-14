#!/bin/bash

echo "You need ffmpeg"

for i in `seq 0 100`; do
	T=`echo "scale=4;${i}.0/100.0" | bc`
	echo "$i/100"
	./src/meshpro input/sphere.off movie.off -starfaces $T -starfaces $T -starfaces $T
	./src/meshview movie.off -output_image movie_$i.jpg -exit_immediately
done

ffmpeg -i movie_%d.jpg art/movie.mp4
ffmpeg -i movie_%d.jpg output/movie.webm

rm movie.off
rm movie_*.jpg
