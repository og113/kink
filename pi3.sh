#!/bin/bash

FILE="results/15.01.15_Tb_0.79_AMP.txt"
echo "output to" $FILE

#AMP=0.065
#echo "amp = "$AMP > $FILE
#./sphaleron4 $AMP

echo "N = 160" > $FILE
echo "Nb = 160" >> $FILE
echo "Tb = 0.79" >> $FILE
echo "" >> $FILE

for i in `seq 0 20`;
    do
    	echo "-------------------------------------------------------------------------------------------------------" >> $FILE
		AMP=$(echo "scale=9; -0.1*$i+1.0" | bc)
		echo "amp = "$AMP >> $FILE
		echo "" >> $FILE
		./sphaleron4 -t1 0.79 -amp $AMP
		#let Nb=80+$i*20
		#echo "Tb = "$Tb >> $FILE
		#./changeInputs Tb $Tb
		echo "pi output:" >> $FILE
		echo "" >> $FILE
		./pi $i >> $FILE
		if [ "$?" = "0" ]; then
			echo "pi3 output:" >> $FILE
			echo "" >> $FILE
			./pi3 $i >> $FILE
			if [ "$?" = "0" ]; then
				echo "success: solution tunnelled" >> $FILE
				echo "" >> $FILE
			else
				echo "solution didn't tunnel" >> $FILE
				echo "" >> $FILE
			fi
			if [ "$i" -gt "0" ]; then 
				echo "compare vector with previous loop output:" >> $FILE
				echo "" >> $FILE
				let im=$i-1
				./compareVector -fileA data/"$i"pip_0.dat -fileB data/"$im"pip_0.dat -c 1 -colA 4 -colB 4 >> $FILE
				echo "" >> $FILE
			fi
			#echo "compare vector with static sphaleron:" >> $FILE
			#echo "" >> $FILE
			#let im=$i-1
			#./compareVector -fileA data/"$i"pip_0.dat -fileB data/staticSphaleron.dat -c 1 -colA 4 -colB 4 >> $FILE
			#echo "" >> $FILE
		else
			echo pi failed, value returned is $? >> $FILE
		fi
		echo "" >> $FILE
    done
