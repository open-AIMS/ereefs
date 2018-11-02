#!/bin/sh

for i in $(seq 100001 100029)
do 
   convert -resize 1000x1000 true*${i}.png true_${i}.gif
   convert -resize 1000x1000 fit*${i}.png fit_${i}.gif
   convert +append true*${i}.gif fit*${i}.gif both_${i}.gif
done
convert -delay 25 -loop 0 both*gif Yasi.gif

