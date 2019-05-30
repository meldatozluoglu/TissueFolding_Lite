#!/bin/bash
numThreads=4 #number of threads to use in sparse matrix solving.
mainPath=pathOfClone-please change this accordingly: eg. /home/{user}/Documents/TissueFolding/
simulationOption=5

# Simulation to run options:
#	1: morphogenesis - from scratch
#	2: pipette aspiration for stiffness measurement - from scratch
#	3: morphogenesis - continue from save
#	4: pipette aspiration for stiffness measurement - continue from save
#	5: morphogenesis - one time step only

displayOption=1
# Simulation can be run with or without the user interface: 
#	0: no display
#	1: with display

case $simulationOption in
	1)
		currSimToTun=$mainPath/ExampleSimulationsToRun/morphogenesis/
		currMode=SimulationOnTheGo
		currModelInputFile=modelinput
		;;
	2)
		currSimToTun=$mainPath/ExampleSimulationsToRun/pipette/
		currMode=SimulationOnTheGo
		currModelInputFile=modelinput
		;;			
	3)
		currSimToTun=$mainPath/ExampleSimulationsToRun/morphogenesis/
		currMode=ContinueFromSave
		currModelInputFile=modelinput
		;;
	4)
		currSimToTun=$mainPath/ExampleSimulationsToRun/pipette/
		currMode=ContinueFromSave
		;;	
	5)
		currSimToTun=$mainPath/ExampleSimulationsToRun/morphogenesis/
		currMode=SimulationOnTheGo
		currModelInputFile=modelinput-onestep	
		;;	
	*)
		echo "invalid simulation option, should be 1 - 5"
		exit 1
		;;	
esac

case $displayOption in
	0)
		executable=$mainPath/TissueFolding/Release/TissueFolding
		;;
	1)
		executable=$mainPath/UserInterface/Release/TissueFoldingUI
		;;
	*)
		echo "invalid display mode option, should be 0 or 1"
		exit 1
		;;			
esac

mkdir $currSimToTun/ScreenShots/
rm $currSimToTun/ScreenShots/frame*
export OMP_NUM_THREADS=$numThreads

if [ "$currMode" = SimulationOnTheGo ] 
then
	cd $currSimToTun
	$executable -mode $currMode -i $currSimToTun/$currModelInputFile -od $currSimToTun > $currSimToTun/tmpRunOutput
fi

if [ "$currMode" = ContinueFromSave ] 
then
	cd $currSimToTun
	$executable -mode $currMode -i $currSimToTun/$currModelInputFile -dInput $currSimToTun/savedData/ -od $currSimToTun > $currSimToTun/tmpRunOutput
fi

