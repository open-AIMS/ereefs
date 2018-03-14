
for i in $(seq 100001 100032)
do
    convert +append {DIN,PhyL_N}_$i.png ab_out_$i.png
done

convert -delay 10 -loop 0 ab_out*.png combined.gif
