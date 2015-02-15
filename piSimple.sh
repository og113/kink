#!/bin/bash

#tmux new -s matlab "matlab -nodesktop -nojvm"

DATE="13.02.15"

FILE="results/"$DATE"_pi_output.txt"
SUMMARY="results/"$DATE"_summary.txt"
echo "output to" $FILE
echo "summary to" $SUMMARY
#echo "output from piSimple.sh on "$DATE > $FILE
#echo "" >> $FILE
#echo "output from piSimple.sh on "$DATE > $SUMMARY
#echo "" >> $SUMMARY
#printf '%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s\n' "N" "Na" "Nb" "Nc" "L" "Tb" "S/F" "E" "Num" "im(S)" >> $SUMMARY

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

#let N=130
#let Nb=80
#let Na=300
#let Nc=2
#L=$(echo "scale=3; 5" | bc)
#LoR=$(echo "scale=3; $L/10.0" | bc)
#changeParameter "N" $N
#changeParameter "Na" $Na
#changeParameter "Nb" $Nb
#changeParameter "Nc" $Nc
#./changeInputs LoR $LoR

loops=12

for j in `seq 0 $loops`
	do
	echo "-------------------------------------------------------------------------------------------------------" >> $FILE
	#AMP=0.5
	#echo "./sphaleron4" >> $FILE
	#echo "" >> $FILE
	#./sphaleron4 -t1 $Tb -r1 $L -amp $AMP  >> $FILE
	#echo "./pi" >> $FILE
	#echo "" >> $FILE
	TIMENUMBER="150213133251"
	let LOOP=28+$j
	Tb=$(echo "scale=5; 0.655-0.005*$j" | bc)
	changeParameter "Tb" $Tb
	echo "pi3 output:" >> $FILE
	echo "" >> $FILE
	./pi3 -tn $TIMENUMBER -loop $LOOP -test 1 >> $FILE
	if [ "$?" = "0" ]; then
		echo "success: solution tunnelled" >> $FILE
		echo "" >> $FILE
		./pi3 -tn $TIMENUMBER -loop $LOOP -lin 1 -changeNa 0 -N 300 -approxOmega 0 >> $FILE
		echo "" >> $FILE
		printf '%-10s%-10s%-10s%-10s%-10s%-10s%-10s\n' $N $Na $Nb $Nc $L $Tb "S" >> $SUMMARY
		echo "-------------------------------------------------------------------------------------------------------" >> $FILE
	else
		echo "solution didn't tunnel" >> $FILE
		echo "" >> $FILE
		printf '%-10s%-10s%-10s%-10s%-10s%-10s%-10s\n' $N $Na $Nb $Nc $L $Tb "FT" >> $SUMMARY
		echo "-------------------------------------------------------------------------------------------------------" >> $FILE
	fi
	echo "" >> $FILE
	done
