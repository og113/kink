#!/bin/bash

FILE="results/16.01.15_Tb_0.74_0.7_AMP_1.0_0.0.txt"
SUMMARY="results/16.01.15_summary.txt"
echo "output to" $FILE
echo "summary to" $SUMMARY

#AMP=0.065
#echo "amp = "$AMP > $FILE
#./sphaleron4 $AMP
N=160
Nb=160

echo "N =" $N > $FILE
echo "Nb =" $Nb >> $FILE
#printf '%-10s %-10s %-10s %-10s ' N Nb Tb AMP "T/N/F" >> $SUMMARY

for j in `seq 0 4`
	do
	echo "-------------------------------------------------------------------------------------------------------" >> $FILE
	echo "-------------------------------------------------------------------------------------------------------" >> $FILE
	Tb=$(echo "scale=9; -0.01*$j+0.74" | bc)
	echo "Tb =" $Tb >> $FILE
	./changeInputs Tb $Tb
	echo "" >> $FILE

	for i in `seq 0 10`;
		do
			echo "-------------------------------------------------------------------------------------------------------" >> $FILE
			AMP=$(echo "scale=9; -0.1*$i+1.0" | bc)
			echo "Tb =" $Tb >> $FILE
			echo "amp = "$AMP >> $FILE
			echo "" >> $FILE
			./sphaleron4 -t1 $Tb -amp $AMP
			#let Nb=80+$i*20
			#echo "Tb = "$Tb >> $FILE
			#./changeInputs Tb $Tb
			echo "pi output:" >> $FILE
			echo "" >> $FILE
			TIMENUMBER="$j"_"$i"
			./pi $TIMENUMBER >> $FILE
			if [ "$?" = "0" ]; then
				echo "pi3 output:" >> $FILE
				echo "" >> $FILE
				./pi3 $TIMENUMBER >> $FILE
				if [ "$?" = "0" ]; then
					echo "success: solution tunnelled" >> $FILE
					echo "" >> $FILE
					printf '%-10s %-10s %-10s %-10s ' $N $Nb $Tb $AMP "T" >> $SUMMARY
				else
					echo "solution didn't tunnel" >> $FILE
					echo "" >> $FILE
					printf '%-10s %-10s %-10s %-10s ' $N $Nb $Tb $AMP "N" >> $SUMMARY
				fi
				if [ "$i" -gt "0" ]; then 
					echo "compare vector with previous loop output:" >> $FILE
					echo "" >> $FILE
					let im=$i-1
					TIMENUMBER_MINUS_ONE="$j"_"$im"
					./compareVector -fileA data/"$TIMENUMBER"pip_0.dat -fileB data/"$TIMENUMBER_MINUS_ONE"pip_0.dat -c 1 -colA 4 -colB 4 >> $FILE
					echo "" >> $FILE
				fi
				#echo "compare vector with static sphaleron:" >> $FILE
				#echo "" >> $FILE
				#let im=$i-1
				#./compareVector -fileA data/"$i"pip_0.dat -fileB data/stable/staticSphaleron.dat -c 1 -colA 4 -colB 4 >> $FILE
				#echo "" >> $FILE
			else
				echo pi failed, value returned is $? >> $FILE
				echo "" >> $FILE
				printf '%-10s %-10s %-10s %-10s ' $N $Nb $Tb $AMP "F" >> $SUMMARY
			fi
			echo "" >> $FILE
		done
		echo "" >> $SUMMARY
	done
