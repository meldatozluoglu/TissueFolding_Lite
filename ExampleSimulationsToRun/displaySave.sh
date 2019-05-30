#!/bin/bash
#mainPath=pathOfClone-please change this accordingly: eg. /home/{user}/Documents/TissueFolding/

simulationOption=1
# Simulation to display options:
#	1: morphogenesis - from scratch
#	2: pipette aspiration for stiffness measurement - from scratch
#	3: morphogenesis - continued from save
#	4: pipette aspiration for stiffness measurement - continued from save
#	5: morphogenesis - one time step only

case $simulationOption in
	1)
		currInputToDisplay=$mainPath/ExampleSimulationsToRun/morphogenesis
		;;
	2)
		currInputToDisplay=$mainPath/ExampleSimulationsToRun/pipette
		;;	
	3)
		currInputToDisplay=$mainPath/ExampleSimulationsToRun/morphogenesis/
		;;
	4)
		currInputToDisplay=$mainPath/ExampleSimulationsToRun/pipette/
		;;	
	5)
		currInputToDisplay=$mainPath/ExampleSimulationsToRun/morphogenesis/	
		;;
	*)
		echo "invalid simulation to display option, should be 1 - 5"
		exit 1
		;;			
esac

mkdir  $currInputToDisplay/ScreenShots
rm $currInputToDisplay/ScreenShots/frame*

$mainPath/UserInterface/Release/TissueFoldingUI  -mode DisplaySave -dInput $currInputToDisplay/  -i $mainPath/ExampleSimulations/modelinputImageSave -od $currInputToDisplay/ > $currInputToDisplay/tmpDisplayOutput
