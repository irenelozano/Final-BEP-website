#!/bin/bash
# Convert PPM frames from SPARTA dump image to video
# Usage: ./make_video.sh [framerate]
# Requires: ffmpeg

FRAMERATE=${1:-30}

echo "Converting image.*.ppm to particle_animation.mp4 at ${FRAMERATE} fps..."

ffmpeg -y -framerate ${FRAMERATE} -pattern_type glob -i 'image.*.ppm' \
    -c:v libx264 -preset medium -crf 18 -pix_fmt yuv420p \
    particle_animation.mp4

echo "Done! Output: particle_animation.mp4"
