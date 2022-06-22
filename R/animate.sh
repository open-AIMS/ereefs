#!/bin/sh
 for i in $(seq 100001 100090)
do 
   convert -resize 1500x1500 ToAnimate/tul*${i}.png ToAnimate/tul_${i}.gif
   convert -resize 1500x1500 ToAnimate/her*${i}.png ToAnimate/her_${i}.gif
#   convert +append true*${i}.gif fit*${i}.gif both_${i}.gif
done
convert -delay 15 -loop 0 ToAnimate/tul*gif Tully.gif
convert -delay 15 -loop 0 ToAnimate/her*gif Herbert.gif
