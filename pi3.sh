#!/bin/bash

FILE="results/21.01.15_Tb_0.8_AMP_0.4.txt"
SUMMARY="results/21.01.15_summary.txt"
echo "output to" $FILE
echo "summary to" $SUMMARY

#AMP=0.065
#echo "amp = "$AMP > $FILE
#./sphaleron4 $AMP
N=160
Nb=160
L=10

echo "N =" $N > $FILE
echo "Nb =" $Nb >> $FILE
echo "L =" $L >> $FILE
echo "" >> $FILE
echo "summary from pi3.sh"> $SUMMARY
printf '%-10s %-10s %-10s %-10s %-10s \n' N Nb Tb AMP "T/N/F" >> $SUMMARY

LoR=$(echo "scale=1; $L/10.0" | bc)
./changeInputs LoR $LoR
./changeInputs N $N
./changeInputs Nb $Nb
./sphaleron -r1 $L >> $FILE
cp data/sphaleron.dat data/stable/sphaleron.dat
cp data/D1.dat ../mpi/data/D1.dat
cp data/D2.dat ../mpi/data/D2.dat
./mx "[D1,D2,diff,maximum] = compareDDS;"
./mx "[V,D] = eigs(D2,2,-20);"
./mx "D"
./mx "V0 = V(:,1);"
./mx "printVector(V0,'../kink/data/stable/sphaleronEigVec.dat');"

for j in `seq 0 4`
	do
	echo "-------------------------------------------------------------------------------------------------------" >> $FILE
	echo "-------------------------------------------------------------------------------------------------------" >> $FILE
	Tb=$(echo "scale=9; -0.01*$j+0.80" | bc)
	echo "Tb =" $Tb >> $FILE
	./changeInputs Tb $Tb
	echo "" >> $FILE

	for i in `seq 3 7`;
		do
			echo "-------------------------------------------------------------------------------------------------------" >> $FILE
			AMP=$(echo "scale=9; 0.1*$i-1.0" | bc)
			echo "Tb =" $Tb >> $FILE
			echo "amp = "$AMP >> $FILE
			echo "" >> $FILE
			./sphaleron4 -r1 $L -t1 $Tb -amp $AMP >> $FILE
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
				./pi3 -tn $TIMENUMBER >> $FILE
				if [ "$?" = "0" ]; then
					echo "success: solution tunnelled" >> $FILE
					echo "" >> $FILE
					printf '%-10s %-10s %-10s %-10s %-10s \n' $N $Nb $Tb $AMP "T" >> $SUMMARY
				else
					echo "solution didn't tunnel" >> $FILE
					echo "" >> $FILE
					printf '%-10s %-10s %-10s %-10s %-10s \n' $N $Nb $Tb $AMP "N" >> $SUMMARY
				fi
				if [ "$i" -gt "0" ]; then 
					echo "compare vector with previous loop output:" >> $FILE
					echo "" >> $FILE
					let im=$i-1
					TIMENUMBER_MINUS_ONE="$j"_"$im"
					./compareVector -fileA data/"$TIMENUMBER"pip_0.dat -fileB data/"$TIMENUMBER_MINUS_ONE"pip_0.dat -c 1 -colA 4 -colB 4 >> $FILE
					echo "" >> $FILE
				fi
				echo "compare vector with static sphaleron:" >> $FILE
				echo "" >> $FILE
				let im=$i-1
				./compareVector -fileA data/"$i"pip_0.dat -fileB data/stable/staticSphaleron.dat -c 1 -colA 4 -colB 4 >> $FILE
				echo "" >> $FILE
			else
				echo pi failed, value returned is $? >> $FILE
				echo "" >> $FILE
				printf '%-10s %-10s %-10s %-10s %-10s \n' $N $Nb $Tb $AMP "F" >> $SUMMARY
			fi
			echo "" >> $FILE
		done
		echo "" >> $SUMMARY
	done
