#!/bin/bash

#tmux new -s matlab "matlab -nodesktop -nojvm"

FILE="results/24.01.15_Bcs.txt"
SUMMARY="results/24.01.15_summary.txt"
echo "output to" $FILE
echo "summary to" $SUMMARY
echo "output from mainChangeBcs.sh" > $FILE
echo "" >> $FILE
echo "------------------------------------------------------------------------------------------------------------------" >> SUMMARY
echo "output from mainChangeBcs.sh" >> $SUMMARY
echo "------------------------------------------------------------------------------------------------------------------" >> SUMMARY
echo "" >> $SUMMARY
printf '%-10s %-10s \n' "zmt" "S/F" >> $SUMMARY

N=120
Na=320
Nb=100
Nc=2
Tb=0.8
L=$(echo "scale=3; 4.5" | bc)
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
echo "./pi" >> $FILE
echo "" >> $FILE
TIMENUMBER=0
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
	else
		echo "solution didn't tunnel" >> $FILE
		echo "" >> $FILE
		printf '%-10s %-10s \n' $zmt "FT" >> $SUMMARY
	fi
else
	echo pi failed, value returned is $? >> $FILE
	echo "" >> $FILE
	printf '%-10s %-10s \n' $zmt "FP" >> $SUMMARY
fi
		
./changeInputs -f mainInputs -n minFileNo -v $TIMENUMBER
./changeInputs -f mainInputs -n maxFileNo -v $TIMENUMBER

loops=4

for j in `seq 0 $loops`
	do
	echo "-------------------------------------------------------------------------------------------------------" >> $FILE
	let lines=2+j*j*j*j
	zmt=nA$lines
	./changeInputs -f mainInputs -n zmt -v $zmt
	echo "zmt =" $zmt >> $FILE
	echo "" >> $FILE
	echo "#################################################################################################" >> $FILE
	./main >> $FILE
	if [ "$?" = "0" ]; then
		printf '%-10s %-10s \n' $zmt "S" >> $SUMMARY
	else
		printf '%-10s %-10s \n' $zmt "FM" >> $SUMMARY
	fi
	echo "#################################################################################################" >> $FILE
	echo "" >> $FILE
	done
