#!/bin/bash

for i in `seq 0 18`;
    do
		AMP=$(echo "scale=9; -0.05*$i+0.1" | bc)
		echo "amp = "$AMP
		./sphaleron4 $AMP
		./pi $i > data/pi_$i.out
		if [ "$?" = "0" ]; then
			./pi3 $i > data/pi3_$i.out
		else
			echo pi failed, value returned is $?
		fi
    done
