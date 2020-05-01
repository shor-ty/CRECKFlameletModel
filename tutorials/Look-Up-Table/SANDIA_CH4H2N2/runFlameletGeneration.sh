#!/bin/bash
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Autor: Tobias Holzmann												#
# Datum: 26.11.2012													#
#															#
# Description:														#
# This script builds non-adiabatic flamelet libraries for OpenFOAM and fluent used by the binaries from Alberto 	#
# Cuoci. The defects has to be defined like below:									#
#															#
# fEd[0]	 	>>> adiabat flame										#
# fEd[1]   -> fEd[x]	>>> positiv enthalpy defects ( 10 20 30 100 200 ...)						#
# fEd[x+1] -> fEd[end]  >>> negativ enthalpy defects (-10 -20 -30 -100 -200 ...)					#
#															#
# its important to set the right order of the values									#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


#
#		EDIT HERE - EDIT HERE - EDIT HERE - EDIT HERE
#			
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# delete Fluent Libs yes = true, No = false										#
dFl=true														#
															#
# used Kinetics:													#
#kinetics="PolimiH2CO"													#
kinetics="PolimiC1C3"													#
															#
# aiabate flame														#
fEd[0]=-800														#
															#
# defects														#
#fEd[1]=-2														#
#fEd[2]=-5														#
#fEd[3]=-10														#
#fEd[3]=-50														#
#fEd[5]=-50														#
#fEd[6]=-80														#
#fEd[7]=-100														#
#fEd[8]=-150														#
#fEd[9]=-250														#
#fEd[4]=-500														#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



# SCRIPT START
#
#
#
# Calculate Flamelets all Flamelets and all enthalpy defects  				
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
./runCleanCase.sh
mkdir Logs 
#mkdir $kinetics

k=0
j=1
for i in ${fEd[@]} 
do
	echo "Calculation of flamelet lib with enthalpy defect ${fEd[$k]}"	
	sed "s/0   kJ\/kg/${fEd[$k]}   kJ\/kg/g" Data/Data.inp > Data/DataEnthalpyDefect${fEd[$k]}.inp 
	./OpenSMOKE_SteadyStateFlamelet.sh --kinetics "$kinetics" --input "Data/DataEnthalpyDefect${fEd[$k]}.inp" --schedule "Operations.inp"  > Logs/Log-Enthalpydefect_${fEd[$k]}
	echo "Finished..."  

	mv Flamelet_$kinetics/Profiles Flamelet_${fEd[$k]} 
	rm -rf Flamelet_$kinetics 
	k=$((k+1))
done
echo "Flamelets generated..."
echo "LookUpTable-Ensemble"

#Ensemble Folder
mkdir PDF-Library
mkdir PDF-Library/Fluent
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


#
# 
# Generate Look-Up-Tables using the Beta-Probability-Function  				
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# insert Flameletlist into Input.inp and rename the output path
for i in ${fEd[*]}
do
	echo "Look-Up-Table for Flamelets with enthalpy defect $i"
	cp Input/Input.inp Input/Input_$i.inp

	sed -i "s/#FlameletLibrary.PDF-0.fla/#FlameletLibrary        PDF-Library\/Fluent\/PDF$i.fla/g" Input/Input_$i.inp
	sed -i "s/#LookUpTable.*lib/#LookUpTable            PDF-Library\/Fluent\/PDF$i.lib/g" Input/Input_$i.inp
	sed -i "s/#LookUpTableOpenFOAM.*/#LookUpTableOpenFOAM    PDF-Library\/PDF$i/g" Input/Input_$i.inp

	#FlameletArray with list of Flamelets
	fA=(`ls Flamelet_$i`)
	k=0
	len=${#fA[*]}

	#Insert Flameletpath to the Input.inp
	while [ $k -lt $len ]; 
	do
		sed -i "8aFlamelet_$i/${fA[$k]}" Input/Input_$i.inp
		let k++
	done

	#Using beta-PDF on flamelet to get L-U-T
	./OpenSMOKE_LookUpTable.sh --input "Input/Input_$i.inp" --kinetics "$kinetics" -build-library > Logs/Log-LookUpTable$i
	echo "Look-Up-Table for enthalpy defect $i generated"
done
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


#
# 
# Prepare the PDF-Library/LookUpTable.out file for OpenFOAM
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


#Create and insert first lines
echo "EnthalpyDefects" > PDF-Library/LookUpTable.out
echo "nphi ${#fEd[*]}" >> PDF-Library/LookUpTable.out
echo "" >> PDF-Library/LookUpTable.out

#search positiv enthalpy values
j=0
k=0

for i in ${fEd[*]}
do
	if [ ${fEd[$j]} -gt 0 ]
	then
		#temp positiv Enthalpy Defects Array
		tmpp[$k]=${fEd[$j]}
		let k++
	fi
	let j++
done

j=0
k=0
for i in ${tmpp[*]}
do
	k=${#tmpp[*]}
	k=$((k-j-1))

	#reverse Array
	sEd[$j]=${tmpp[$k]}
	let j++
done

#adiabat flamelet
j=0
k=${#tmpp[*]}

for i in ${fEd[*]}
do
	if [ ${fEd[$j]} -eq 0 ]
	then
		#add adiabat flame to array
		sEd[$k]=${fEd[$j]}
	fi
	let j++
done

#search for negativ enthalpy values
j=0
k=0

for i in ${fEd[*]}
do
	if [ ${fEd[$j]} -lt 0 ]
	then
		#temp negativ Enthalpy Defects Array
		tmpn[$k]=${fEd[$j]}
		let k++
	fi
	let j++
done

#add negativ enthaply values after adiabat 
j=${#tmpp[*]}
j=$((j+1))
k=0
for i in ${tmpn[*]}
do
	#noreverse of Array
	sEd[$j]=${tmpn[$k]}
	let j++
	let k++
done

#insert sorted array into LookUpTable.out
j=0
for i in ${fEd[*]}
do
	if [ ${sEd[$j]} -eq 0 ] 
	then
		echo "${sEd[$j]}.   PDF${sEd[$j]}" >> PDF-Library/LookUpTable.out
	else
		echo "${sEd[$j]}.e3   PDF${sEd[$j]}" >> PDF-Library/LookUpTable.out
	fi
	let j++
done

#remove fluent libs

if [ $dFl == true ]
then
	echo "Remove Fluent libs"
	rm -rf PDF-Library/Fluent
fi
echo "Finished..."
