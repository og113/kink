#!/bin/bash

#tmux new -s matlab "matlab -nodesktop -nojvm"

#note that this has been written wrongly as only main needs to be repeated for each value of cutoff, not pi or pi3 - should rewrite if reusing

DATE="22.02.15"
FILE="results/"$DATE"_cutoff.txt"
echo "output to" $FILE
echo "output from mainChangeCutoff.sh" >> $FILE
echo "" >> $FILE

./changeInputs -f mainInputs -n inF -v p

N=130
Na=320
Nb=80
Nc=2
Tb=0.8
L=$(echo "scale=3; 5" | bc)
LoR=$(echo "scale=3; $L/10.0" | bc)
./changeInputs LoR $LoR
./changeInputs Tb $Tb
./changeInputs N $N
./changeInputs Na $Na
./changeInputs Nb $Nb
./changeInputs Nc $Nc
echo "Tb =" $Tb >> $FILE
echo "L =" $L >> $FILE
echo "N =" $N >> $FILE
echo "Na =" $Na >> $FILE
echo "Nb =" $Nb >> $FILE
echo "Nc =" $Nc >> $FILE
echo "" >> $FILE

./sphaleron -r1 $L >> $FILE
cp data/sphaleron.dat data/stable/sphaleron.dat
cp data/D1.dat ../mpi/data/D1.dat
cp data/D2.dat ../mpi/data/D2.dat
./mx "[D1,D2,diff,maximum] = compareDDS;"
./mx "[V,D] = eigs(D2,2,-20);"
./mx "D"
./mx "V0 = V(:,1);"
./mx "printVector(V0,'../kink/data/stable/sphaleronEigVec.dat');"

AMP=0.5
echo "./sphaleron4" >> $FILE
echo "" >> $FILE
./sphaleron4 -t1 $Tb -r1 $L -amp $AMP  >> $FILE

loops=4

for j in `seq 0 $loops`
	do
	echo "-------------------------------------------------------------------------------------------------------" >> $FILE
	let CUTOFF=18-j*2
	echo "CUTOFF = $CUTOFF" >> $FILE
	./changeInputs cutoff $CUTOFF
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
			./pi3 -tn $TIMENUMBER -test 0 >> $FILE
			echo "" >> $FILE
			./changeInputs -f mainInputs -n minFileNo -v $TIMENUMBER
			./changeInputs -f mainInputs -n maxFileNo -v $TIMENUMBER
			echo "#################################################################################################" >> $FILE
			./main >> $FILE
			echo "#################################################################################################" >> $FILE
		else
			echo "solution didn't tunnel" >> $FILE
			echo "" >> $FILE
		fi
	else
		echo pi failed, value returned is $? >> $FILE
		echo "" >> $FILE
	fi
	echo "" >> $FILE
	done
