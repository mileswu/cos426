#!/bin/bash

# First Step
./src/particleview ./input/particlerace.scn -recordandquit 50
convert -delay 3 -loop 0 video-frames/frame*.jpg output/particlerace.gif
rm video-frames/*

./src/particleview ./input/gravityrace.scn -recordandquit 50
convert -delay 3 -loop 0 video-frames/frame*.jpg output/gravityrace.gif
rm video-frames/*

./src/particleview ./input/dragrace.scn -recordandquit 50
convert -delay 3 -loop 0 video-frames/frame*.jpg output/dragrace.gif
rm video-frames/*

# Sources/sinks
./src/particleview ./input/sources.scn -recordandquit 1000 -videoskip 10
convert -delay 3 -loop 0 video-frames/frame*.jpg output/sources.gif
rm video-frames/*

./src/particleview ./input/sinks.scn -recordandquit 14000 -videoskip 160 -videofps 400 -nomutualattraction
convert -delay 3 -loop 0 video-frames/frame*.jpg output/sinks.gif
rm video-frames/*

# Lifetime
./src/particleview ./input/particlerace-lifetime.scn -nodynamic -recordandquit 200
convert -delay 3 -loop 0 video-frames/frame*.jpg output/particlerace-lifetime.gif
rm video-frames/*

# Collisions
./src/particleview ./input/cannonball-collisions.scn -recordandquit 200
convert -delay 3 -loop 0 video-frames/frame*.jpg output/cannonball-collisions.gif
rm video-frames/*

# Interactions
./src/particleview ./input/opposites.scn -recordandquit 100
convert -delay 3 -loop 0 video-frames/frame*.jpg output/opposites.gif
rm video-frames/*

./src/particleview ./input/rope.scn -recordandquit 1000 -videofps 50 -videoskip 5
convert -delay 3 -loop 0 video-frames/frame*.jpg output/rope.gif
rm video-frames/*

# Rendering
./src/particleview ./input/sinks-lifetime.scn -recordandquit 1000 -videoskip 5
convert -delay 3 -loop 0 video-frames/frame*.jpg output/sinks-lifetime.gif
rm video-frames/*

# 

