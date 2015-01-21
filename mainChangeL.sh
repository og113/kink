#!/bin/bash

#tmux new -s matlab "matlab -nodesktop -nojvm"

FILE="results/20.01.15_L_5.0_10.0_Tb_0.8.txt"
SUMMARY="results/20.01.15_summary.txt"
echo "output to" $FILE
echo "summary to" $SUMMARY
echo "output from mainChangeL.sh" > $FILE
echo "" >> $FILE
echo "output from mainChangeL.sh" > $SUMMARY
echo "" >> $SUMMARY
printf '%-10s %-10s \n' "L" "S/F" >> $SUMMARY

for j in `seq 0 10`
	do
	echo "-------------------------------------------------------------------------------------------------------" >> $FILE
	L=$(echo "scale=1; -0.5*$j+10.0" | bc)
	LoR=$(echo "scale=1; $L/10.0" | bc)
	echo "L =" $L >> $FILE
	echo "" >> $FILE
	./changeInputs LoR $LoR
	echo "./sphaleron" >> $FILE
	echo "" >> $FILE
	./sphaleron -r1 $L >> $FILE
	cp data/sphaleron.dat data/stable/sphaleron.dat
	cp data/D1.dat ../mpi/data/D1.dat
	cp data/D2.dat ../mpi/data/D2.dat
	./mx "[D1,D2,diff,maximum] = compareDDS;"
	./mx "[D1,D2,diff,maximum] = compareDDS;"
	./mx "[V,D] = eigs(D2,2,-20);"
	./mx "D"
	./mx "V0 = V(:,1);"
	./mx "printVector(V0,'../kink/data/stable/sphaleronEigVec.dat');"
	AMP="0.3"
	echo "amp = "$AMP >> $FILE
	echo "./sphaleron4" >> $FILE
	echo "" >> $FILE
	./sphaleron4 -r1 $L -amp $AMP
	echo "./pi" >> $FILE
	echo "" >> $FILE
	TIMENUMBER=$j
	./pi $TIMENUMBER >> $FILE
	if [ "$?" = "0" ]; then
		echo "pi3 output:" >> $FILE
		echo "" >> $FILE
		./pi3 -tn $TIMENUMBER >> $FILE
		if [ "$?" = "0" ]; then
			echo "success: solution tunnelled" >> $FILE
			echo "" >> $FILE
			./changeInputs -f mainInputs -n minfileNo -v "$TIMENUMBER"
			echo "#################################################################################################" >> $FILE
			./main >> $FILE
			echo "#################################################################################################" >> $FILE
			if [ "$?" = "0" ]; then
				printf '%-10s %-10s \n' $L "S" >> $SUMMARY
			else
				printf '%-10s %-10s \n' $L "FM" >> $SUMMARY
			fi
		else
			echo "solution didn't tunnel" >> $FILE
			echo "" >> $FILE
			printf '%-10s %-10s ' $L "FT" >> $SUMMARY
		fi
	else
		echo pi failed, value returned is $? >> $FILE
		echo "" >> $FILE
		printf '%-10s %-10s \n' $L "FP" >> $SUMMARY
	fi
	echo "" >> $FILE
	done
