#!/bin/bash

i=000

for file in `ls *.pdf -r` 
do 
	pdftoppm $file $i -png -rx 300 -ry 300
	printf -v i "%03d" $((10#$i+1))
done

ffmpeg -y\
  -framerate 20 \
  -pattern_type glob \
  -i '*.png' \
  -r 10 \
  -vf scale=512:-1 \
  out.gif \
;

rm *.png