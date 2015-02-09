#!/bin/bash

#tmux new -s matlab "matlab -nodesktop -nojvm"

FILE="results/6.02.15_NL_output.txt"
SUMMARY="results/6.02.15_summary.txt"
echo "output to" $FILE
echo "summary to" $SUMMARY
echo "output from mainChangeNL.sh" >> $FILE
echo "" >> $FILE
#echo "output from mainChangeNL.sh" > $SUMMARY
#echo "" >> $SUMMARY
#printf '%-10s%-10s%-10s%-10s%-10s%-10s%-10s\n' "N" "Na" "Nb" "Nc" "L" "Ta" "S/F" >> $SUMMARY

function getEigenvectors {
	echo "./sphaleron" >> $FILE
	echo "" >> $FILE
	./sphaleron -r1 $1 >> $FILE #1st argument is L
	cp data/sphaleron.dat data/stable/sphaleron.dat
	cp data/D1.dat ../mpi/data/D1.dat
	cp data/D2.dat ../mpi/data/D2.dat
	./mx "[D1,D2,diff,maximum] = compareDDS;"
	./mx "[V,D] = eigs(D2,2,-20);"
	./mx "D"
	./mx "V0 = V(:,1);"
	./mx "printVector(V0,'../kink/data/stable/sphaleronEigVec.dat');"
}

function changeParameter {
	./changeInputs $1 $2
	echo $1 " = " $2 >> $FILE
}

Tb=0.80
changeParameter Tb $Tb
loops=0

for j in `seq 0 $loops`
	do
	echo "-------------------------------------------------------------------------------------------------------" >> $FILE
	let N=150
	let Nb=80
	let Na=300
	let Nc=4
	Ta=$(echo "scale=3; $Na*$Tb/($Nb-1.0)" | bc)
	L=$(echo "scale=3; 5" | bc)
	LoR=$(echo "scale=3; $L/10.0" | bc)
	changeParameter "N" $N
	changeParameter "Na" $Na
	changeParameter "Nb" $Nb
	changeParameter "Nc" $Nc
	./changeInputs LoR $LoR
	echo "L =" $L >> $FILE
	echo "" >> $FILE
	getEigenvectors $L
	AMP=0.5
	echo "./sphaleron4" >> $FILE
	echo "" >> $FILE
	./sphaleron4 -t1 $Tb -r1 $L -amp $AMP  >> $FILE
	echo "./pi" >> $FILE
	echo "" >> $FILE
	TIMENUMBER=$j
	./pi $TIMENUMBER >> $FILE
	if [ "$?" = "0" ]; then
		echo "pi3 output:" >> $FILE
		echo "" >> $FILE
		./pi3 -tn $TIMENUMBER -test 1 >> $FILE
		if [ "$?" = "0" ]; then
			echo "success: solution tunnelled" >> $FILE
			echo "" >> $FILE
			./pi3 -tn $TIMENUMBER -lin 1 >> $FILE
			echo "" >> $FILE
			./pi3 -tn $TIMENUMBER -test 0 >> $FILE
			echo "" >> $FILE
			./changeInputs -f mainInputs -n minFileNo -v $TIMENUMBER
			./changeInputs -f mainInputs -n maxFileNo -v $TIMENUMBER
			echo "#################################################################################################" >> $FILE
			./main >> $FILE
			if [ "$?" = "0" ]; then
				printf '%-10s%-10s%-10s%-10s%-10s%-10s%-10s\n' $N $Na $Nb $Nc $L $Ta "S" >> $SUMMARY
			else
				printf '%-10s%-10s%-10s%-10s%-10s%-10s%-10s\n' $N $Na $Nb $Nc $L $Ta "FM" >> $SUMMARY
			fi
			echo "#################################################################################################" >> $FILE
		else
			echo "solution didn't tunnel" >> $FILE
			echo "" >> $FILE
			printf '%-10s%-10s%-10s%-10s%-10s%-10s%-10s\n' $N $Na $Nb $Nc $L $Ta "FT" >> $SUMMARY
		fi
	else
		echo pi failed, value returned is $? >> $FILE
		echo "" >> $FILE
		printf '%-10s%-10s%-10s%-10s%-10s%-10s%-10s\n' $N $Na $Nb $Nc $L $Ta "FP" >> $SUMMARY
	fi
	echo "" >> $FILE
	done
