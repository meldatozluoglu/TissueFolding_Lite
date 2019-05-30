#!/bin/bash
#mainPath=pathOfClone-please change this accordingly: eg. /home/{user}/Documents/TissueFolding/

simulationOption=2
# Simulation to display options:
#	1: morphogenesis
#	2: pipette aspiration for stiffness measurement

case $simulationOption in
	1)
		currInputToDisplay=$mainPath/ExampleSimulations/morphogenesis
		
		#unzipping the large save files if necessary
		echo unzipping the save files if have not been unzipped yet. If the display fails please check Save_Frame and Save_Growth are unzipped.
		zippedSaveFrameFile=$mainPath/ExampleSimulations/morphogenesis/Save_Frame.zip
		zippedSaveGrowthFile=$mainPath/ExampleSimulations/morphogenesis/Save_Growth.zip
		unzippedSaveFrameFile=$mainPath/ExampleSimulations/morphogenesis/Save_Frame
		unzippedSaveGrowthFile=$mainPath/ExampleSimulations/morphogenesis/Save_Growth
		if [ -e "$zippedSaveFrameFile" -a ! -e "$unzippedSaveFrameFile" ]
		then
			unzip $zippedSaveFrameFile -d $mainPath/ExampleSimulations/morphogenesis/
		else
			echo Save_Frame is already unzipped
		fi
		if [ -e "$zippedSaveGrowthFile" -a ! -e "$unzippedSaveGrowthFile" ]
		then
			unzip $zippedSaveGrowthFile -d $mainPath/ExampleSimulations/morphogenesis/
		else
			echo Save_Growth is already unzipped
		fi
		;;
	2)
		currInputToDisplay=$mainPath/ExampleSimulations/pipette
		;;	
	*)
		echo "invalid simulation to display option, should be 1 - 2"
		exit 1
		;;			
esac

mkdir  $currInputToDisplay/ScreenShots
rm $currInputToDisplay/ScreenShots/frame*

$mainPath/UserInterface/Release/TissueFoldingUI  -mode DisplaySave -dInput $currInputToDisplay/  -i $mainPath/ExampleSimulations/modelinputImageSave -od $currInputToDisplay/ > $currInputToDisplay/tmpDisplayOutput
