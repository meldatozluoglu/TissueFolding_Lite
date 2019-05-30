/*
 * ModelInputObject.cpp
 *
 *  Created on: 5 Aug 2014
 *      Author: melda
 */

#include "ModelInputObject.h"
#include "Simulation.h"
#include "GrowthFunctionBase.h"
#include "GrowthFunctionTypes.h"

#include <vector>
#include <iostream>
#include <fstream>

ModelInputObject::ModelInputObject(){
}

ModelInputObject::~ModelInputObject(){
}

bool ModelInputObject::readParameters(){
	/**
	 *  This function will read all available model inputs from the file ModelInputObject#parameterFileName. \n
	 *  It will start by opening the model input file, after each attempt to open a file, there will be a health check to ensure the file
	 *  could be opened. In case there are issues with the file (most common one being the file is not opened due to a path error),
	 *  the function will throw an error with corresponding explanatory error message, and quit the simulation.
	 */
	bool Success = true;
    std::ifstream parametersFile;
    parametersFile.open(parameterFileName, std::ifstream::in);
	Success = checkFileStatus(parametersFile,parameterFileName);
	if (!Success){
		return Success;
	}
	/**
	 *  After successfully opening the input file, the function will read it until it reaches to the end of the file.
	 *  This will involve a series of private functions, which are thoroughly documented in source code, while
	 *  the processed documentation may or may not be available in this user interface documentation structure.
	 */
	while (!parametersFile.eof()){
        std::string currline;
		getline(parametersFile,currline);
		if (!parametersFile.eof()){
			if (currline.empty()){
				//skipping empty line
				continue;
			}
            std::istringstream currSStrem(currline);
            std::string currParameterHeader;
			currSStrem >> currParameterHeader;
			if (currParameterHeader == ""){
				//just in case an empty line stores two consecutive line ends (/n/n),
				//and thus, currline is not technically empty, while the line is still empty
				continue;
			}
			/**
			 * Depending on the header the function it encounters it will read:
			 */
			if(currParameterHeader == "InputMeshParameters:"){
				/**
				 * Mesh geometry related parameters through the private function ModelInputObject#readMeshParameters
				 */
				Success  = readMeshParameters(parametersFile);
			}
			else if(currParameterHeader == "PeripodialMembraneParameters:"){
				/**
				 * Peripodial membrane structure related parameters through the private function ModelInputObject#readPeripodialMembraneParameters
				 */
				Success  = readPeripodialMembraneParameters(parametersFile);
			}
			else if(currParameterHeader == "LinkerZoneParameters:"){
				/**
				 * Linker zone physical properties through the private function ModelInputObject#readLinkerZoneParameters
				 */
				Success  = readLinkerZoneParameters(parametersFile);
			}
			else if(currParameterHeader == "NodeFixingOptions:"){
				/**
				 * Inputs relating to fixing the nodes of the tissue through the private function ModelInputObject#readNodeFixingParameters
				 */
				Success  = readNodeFixingParameters(parametersFile);
			}
			else if(currParameterHeader == "NodeBindingOptions:"){
				/**
				 * Inputs relating to binding the displacement the nodes of the tissue through the private function ModelInputObject#readNodeBindingParameters
				 */
				Success  = readNodeBindingParameters(parametersFile);
			}
			else if(currParameterHeader == "ExternalViscositySetup:"){
				/**
				 * Inputs relating to the external viscosity felt by the tissue through private function ModelInputObject#readExternalViscosityParameters
				 */
				Success  = readExternalViscosityParameters(parametersFile);
			}
			else if(currParameterHeader == "TimeParameters:"){
				/**
				 * Simulation time and time step related parameters through the private function ModelInputObject#readTimeParameters
				 */
				Success  = readTimeParameters(parametersFile);
			}
			else if(currParameterHeader == "PysicalProperties:"){
				/**
				 * Physical parameters of the tissue  through the private function ModelInputObject#readPysicalProperties
				 */
				Success  = readPysicalProperties(parametersFile);
			}
			else if (currParameterHeader == "SaveOptions:"){
				/**
				 * Save options of the simulation through the private function ModelInputObject#readSaveOptions
				 */
				Success  = readSaveOptions(parametersFile);
			}
			else if (currParameterHeader == "GrowthOptions:"){
				/**
				 * Growth functions and related parameters of the simulation through the private function ModelInputObject#readGrowthOptions
				 */
				Success  = readGrowthOptions(parametersFile);
			}
			else if(currParameterHeader == "ShapeChangeOptions:"){
				/**
				 * Shape change functions and related parameters of the simulation through the private function ModelInputObject#readShapeChangeOptions
				 */
				Success  = readShapeChangeOptions(parametersFile);
			}
			else if (currParameterHeader == "PlasticDeformationOptions:"){
				/**
				* Plastic deformation options through the private function ModelInputObject#readPlasticDeformationOptions
				*/
				Success  = readPlasticDeformationOptions(parametersFile);
			}
			else if(currParameterHeader == "Stretcher:"){
				/**
				 * Stretcher experimental setup parameters of the simulation through the private function ModelInputObject#readStretcherSetup
				 */
				Success  = readStretcherSetup(parametersFile);
			}
			else if(currParameterHeader == "Pipette_Aspiration:"){
				/**
				 * Pipette aspiration experimental setup parameters of the simulation through the private function ModelInputObject#readPipetteSetup
				 */
				Success  = readPipetteSetup(parametersFile);
			}
			else if(currParameterHeader == "ECM_Perturbation:"){
				/**
				 * Setting perturbations to ECM, current setup includes softening of apical or basal ECM at a given time point.
				 */
				Success  = readECMPerturbation(parametersFile);
			}
			else if(currParameterHeader == "Stiffness_Perturbation:"){
				/**
				 * Setting perturbations to tissue stiffness, current setup includes changing Young's modulus for apical, basal or whole tissue at a given time point.
				 */
				Success  = readStiffnessPerturbation(parametersFile);
			}
			else if(currParameterHeader == "Marker_Ellipses:"){
				/**
				 * Setting perturbations to ECM, current setup includes softening of apical or basal ECM at a given time point.
				 */
				Success  = readMarkerEllipseBandOptions(parametersFile);
			}
			else if(currParameterHeader == "ExplicitECMOptions:"){
				/**
				 * Setting perturbations to ECM, current setup includes softening of apical or basal ECM at a given time point.
				 */
				Success  = readExplicitECMOptions(parametersFile);
			}
			else if(currParameterHeader == "AdhesionOptions:"){
				/**
				 * Setting the flag to activate adhesion in the model
				 */
				Success  = readAdhesionOptions(parametersFile);
			}
			else if(currParameterHeader == "NodeCollapseOptions:"){
				/**
				 * Setting the flag to activate node collapsing for elements dangerously close to flipping
				 */
				Success  = readNodeCollapseOptions(parametersFile);
			}
			else if(currParameterHeader == "ExplicitActinOptions:"){
				/**
				 * Setting explicit actin options, The actin layer will not grow in z.
				 */
				Success  = readExplicitActinOptions(parametersFile);
			}
			else if(currParameterHeader == "ColumnViseVolumeConservationOptions:"){
				/**
				 * Setting volume redistributions options, The volume will be redistributed within a column of elements, starting from the apical surface.
				 */
				Success  =readColumnViseVolumeConservationOptions(parametersFile);
			}
			else if(currParameterHeader == "MutationOptions:"){
				/**
				 * Setting mutation options.
				 */
				Success  = readMutationOptions(parametersFile);
			}
			else if(currParameterHeader == "zShellOptions:"){
				/**
				 * Setting z enclosement options.
				 */
				Success  = readEnclosementOptions(parametersFile);
			}
			else if (currParameterHeader == "ArtificialRelaxationOptions:"){
				Success  = readartificialRelaxationOptions(parametersFile);;
			}
			else {
				/**
				 * In the case that the function, or any of the above listed parameter reading functions, encounters an
				 * unexpected line in the model input file, it will throw an error with a corresponding explanatory message,
				 * and quit the simulation.
				 */
                std::cerr<<"Unidentified parameter input line: "<<std::endl;
                std::cerr<<"		"<<currParameterHeader<<std::endl;
				return false;
			}
			if (!Success){
				return false;
			}
		}
	}
	return Success;
}

bool ModelInputObject::checkFileStatus(std::ifstream& file, std::string fileName){
	if (!file.is_open()) {
        std::cerr<<"Cannot open parameter input file, "<<fileName<<std::endl;
		return false;
	}
	if (!file.good()) {
        std::cerr<<"File does not exist, "<<fileName<<std::endl;
		return false;
	}
	return true;
}


bool ModelInputObject::readGrowthOptions(std::ifstream& file){
    std::string currHeader;
    file >> currHeader;
	int n;
    std::cout<<"entered read growth options, current header: "<<currHeader<<std::endl;
	if(currHeader == "NumberofGrowthFunctions(int):"){
		file >> n;
		Sim->nGrowthFunctions = n;
	}
	else{
        printErrorMessage(currHeader,"Growth options","NumberofGrowthFunctions(int):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "GridGrowthsPinnedOnInitialMesh(bool):"){
		bool tmpBool;
		file >> tmpBool;
		Sim->GridGrowthsPinnedOnInitialMesh = tmpBool;
	}
	else{
        printErrorMessage(currHeader,"Growth options","GridGrowthsPinnedOnInitialMesh(bool):");
        return false;
	}
	file >> currHeader;
	if(currHeader == "PinningUpdateTimes(number-times(sec)):"){
		file >> Sim->nGrowthPinning;
		if ( Sim->nGrowthPinning > 0){
            for (int j=0; j<Sim->nGrowthPinning; ++j){
                int currentUpdateTime;
                file >> currentUpdateTime;
                Sim->growthPinUpdateTime.push_back(currentUpdateTime);
                Sim->growthPinUpdateBools.push_back(false);
            }
		}
	}
	else{
        printErrorMessage(currHeader,"Growth options","PinningUpdateTimes(number-times(sec):");
        return false;
	}
	file >> currHeader;
	if (currHeader == "GridGrowthsInterpolationType(0=step,1=linear):"){
		file >> Sim->gridGrowthsInterpolationType;
	}
	else{
        printErrorMessage(currHeader,"Growth options","GridGrowthsInterpolationType(0=step,1=linear):");
		return false;
	}
	for (int i = 0; i<n; ++i){
		file >> currHeader;
		int type;
		cout<<"inside the loop, read growth options, current header: "<<currHeader<<endl;
		if(currHeader == "GrowthFunctionType(int-seeDocumentation):"){
			file >> type;
		}
		else{
            printErrorMessage(currHeader,"Growth options","GrowthFunctionType(int-seeDocumentation):");
            return false;
		}
		bool Success = true;
		if (type == 1){
			Success = readGrowthType1(file);
		}
		else if (type == 2){
			Success = readGrowthType2(file);
		}
		else if (type == 3){
			Success = readGrowthType3(file);
		}
		else{
            printErrorMessage(std::to_string(type),"Growth options","a valid type: {1, 2, 3}");
            return false;
		}
		if (!Success){
			return false;
		}
	}
	return true;
}

bool ModelInputObject::readGrowthType1(std::ifstream& file){
    std::string currHeader;
	file >> currHeader;
	float initialtime;
	float finaltime;
	bool applyToColumnarLayer = false;
	bool applyToPeripodialMembrane = false;
	bool applyToBasalECM = false;
	bool applyToLateralECM = false;
	float DVRate;
	float APRate;
	float ABRate;
	float angle;
	if(currHeader == "InitialTime(sec):"){
		file >> initialtime;
	}
	else{
        printErrorMessage(currHeader,"Growth type 1 options","InitialTime(sec)");
		return false;
	}
	file >> currHeader;
	if(currHeader == "FinalTime(sec):"){
		file >> finaltime;
	}
	else{
        printErrorMessage(currHeader,"Growth type 1 options","FinalTime(sec)");
        return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToColumnarLayer(bool):"){
			file >> applyToColumnarLayer;
	}
	else{
        printErrorMessage(currHeader,"Growth type 1 options","ApplyToColumnarLayer(bool):");
        return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToPeripodialMembrane(bool):"){
			file >> applyToPeripodialMembrane;
	}
	else{
        printErrorMessage(currHeader,"Growth type 1 options","ApplyToPeripodialMembrane(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToBasalECM(bool):"){
			file >> applyToBasalECM;
	}
	else{
        printErrorMessage(currHeader,"Growth type 1 options","ApplyToBasalECM(bool):");
        return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToLateralECM(bool):"){
			file >> applyToLateralECM;
	}
	else{
        printErrorMessage(currHeader,"Growth type 1 options","ApplyToLateralECM(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "MaxValue(fractionPerHour-DV,AP,AB):"){
		double timeMultiplier = 1.0 / 3600.0;
		file >> DVRate;
		DVRate *= timeMultiplier;
		file >> APRate;
		APRate *= timeMultiplier;
		file >> ABRate;
		ABRate *= timeMultiplier;
	}
	else{
        printErrorMessage(currHeader,"Growth type 1 options","MaxValue(fractionPerHour-DV,AP,AB):");
        return false;
	}
	file >> currHeader;
	if(currHeader == "Angle(degrees):"){
		file >> angle;
		angle *= M_PI/180.0; 	 // converting to radians
	}
    else{
        printErrorMessage(currHeader,"Growth type 1 options","Angle(degrees):");
		return false;
	}   
    int Id = Sim->GrowthFunctions.size();
    //type is 1
    std::unique_ptr<GrowthFunctionBase> GSBp = std::make_unique<UniformGrowthFunction>(Id, 1, initialtime, finaltime, applyToColumnarLayer, applyToPeripodialMembrane, applyToBasalECM, applyToLateralECM, DVRate, APRate,  ABRate, angle);
    Sim->GrowthFunctions.push_back(std::move(GSBp));
    return true;
}

bool ModelInputObject::readGrowthType2(std::ifstream& file){
    std::string currHeader;
	file >> currHeader;
    std::cout<<"entered read growth type 2, current header: "<<currHeader<<std::endl;
	float initialtime;
	float finaltime;
	bool applyToColumnarLayer = false;
	bool applyToPeripodialMembrane = false;
	bool applyToBasalECM = false;
	bool applyToLateralECM = false;
	float CentreX,CentreY;
	float innerR, outerR;
	float DVRate;
	float APRate;
	float ABRate;
	if(currHeader == "InitialTime(sec):"){
		file >> initialtime;
	}
	else{
        printErrorMessage(currHeader,"Growth type 2 options","InitialTime(sec):");
        return false;
	}
	file >> currHeader;
	if(currHeader == "FinalTime(sec):"){
		file >> finaltime;
	}
	else{
        printErrorMessage(currHeader,"Growth type 2 options","FinalTime(sec):");
        return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToColumnarLayer(bool):"){
			file >> applyToColumnarLayer;
	}
	else{
        printErrorMessage(currHeader,"Growth type 2 options","ApplyToColumnarLayer(bool):");
        return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToPeripodialMembrane(bool):"){
			file >> applyToPeripodialMembrane;
	}
	else{
        printErrorMessage(currHeader,"Growth type 2 options","ApplyToPeripodialMembrane(sec):");
        return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToBasalECM(bool):"){
			file >> applyToBasalECM;
	}
	else{
        printErrorMessage(currHeader,"Growth type 2 options","ApplyToBasalECM(sec):");
        return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToLateralECM(bool):"){
			file >> applyToLateralECM;
	}
	else{
        printErrorMessage(currHeader,"Growth type 2 options","ApplyToLateralECM(sec):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "Centre:"){
		file >> CentreX;
		file >> CentreY;
	}
	else{
        printErrorMessage(currHeader,"Growth type 2 options","Centre:");
		return false;
	}
	file >> currHeader;
	if(currHeader == "InnerRadius:"){
		file >> innerR;
	}
	else{
        printErrorMessage(currHeader,"Growth type 2 options","InnerRadius:");
        return false;
	}
	file >> currHeader;
	if(currHeader == "OuterRadius:"){
		file >> outerR;
	}
	else{
        printErrorMessage(currHeader,"Growth type 2 options","OuterRadius:");
        return false;
	}
	file >> currHeader;
	if(currHeader == "MaxValue(fractionPerHour-DV,AP,AB):"){
		double timeMultiplier = 1.0 / 3600.0;
		file >> DVRate;
		DVRate *= timeMultiplier;
		file >> APRate;
		APRate *= timeMultiplier;
		file >> ABRate;
		ABRate *= timeMultiplier;
	}
	else{
        printErrorMessage(currHeader,"Growth type 2 options","MaxValue(fractionPerHour-DV,AP,AB):");
		return false;
	}
	double angle;
	file >> currHeader;
	if(currHeader == "Angle(degrees):"){
		file >> angle;
		angle *= M_PI/180.0; 	 // converting to radians
	}
	else{
        printErrorMessage(currHeader,"Growth type 2 options","Angle(degrees):");
		return false;
	}
    int Id = Sim->GrowthFunctions.size();
    //type is 2
    std::unique_ptr<GrowthFunctionBase> GSBp = std::make_unique<RingGrowthFunction>(Id, 2, initialtime, finaltime, applyToColumnarLayer, applyToPeripodialMembrane, applyToBasalECM, applyToLateralECM, CentreX, CentreY, innerR,  outerR, DVRate, APRate,  ABRate, angle);
    Sim->GrowthFunctions.push_back(std::move(GSBp));
    return true;
}

bool ModelInputObject::readGrowthType3(std::ifstream& file){
    std::string currHeader;
	file >> currHeader;
    std::cout<<"entered read growth type 3, current header: "<<currHeader<<std::endl;
	float initialtime;
	float finaltime;
	bool applyToColumnarLayer = false;
	bool applyToPeripodialMembrane = false;
	bool applyToLateralECM = false;
	bool applyToBasalECM = false;
    size_t gridX, gridY;
    std::vector<std::vector<std::array<double,3>>> GrowthMatrix; //[grid_i][grid_j][x,y,z]
    std::vector<std::vector<double>>  AngleMatrix; //[grid_i][grid_j][tetha]
	if(currHeader == "InitialTime(sec):"){
		file >> initialtime;
	}
	else{
        printErrorMessage(currHeader,"Growth type 2 options","InitialTime(sec):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "FinalTime(sec):"){
		file >> finaltime;
	}
	else{
        printErrorMessage(currHeader,"Growth type 2 options","FinalTime(sec):");
        return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToColumnarLayer(bool):"){
			file >> applyToColumnarLayer;
	}
	else{
        printErrorMessage(currHeader,"Growth type 2 options","ApplyToColumnarLayer(bool):");
        return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToPeripodialMembrane(bool):"){
			file >> applyToPeripodialMembrane;
	}
	else{
        printErrorMessage(currHeader,"Growth type 2 options","ApplyToPeripodialMembrane(bool):");
        return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToBasalECM(bool):"){
			file >> applyToBasalECM;
	}
	else{
        printErrorMessage(currHeader,"Growth type 2 options","ApplyToBasalECM(bool):");
        return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToLateralECM(bool):"){
			file >> applyToLateralECM;
	}
	else{
        printErrorMessage(currHeader,"Growth type 2 options","ApplyToLateralECM(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "Filename(full-path):"){
        std::string filepath;
		file >> filepath;
        std::cerr<<" filename is: "<<filepath<<std::endl;
		const char* name_growthRates = filepath.c_str();
        std::ifstream GrowthRateFile;
        GrowthRateFile.open(name_growthRates, std::ifstream::in);
		if (!(GrowthRateFile.good() && GrowthRateFile.is_open())){
            std::cerr<<"could not open growth rate file: "<<name_growthRates<<std::endl;
			return false;
		}
		//adding the indice of the growth matrix
        //std::cout<<"reading from growth file"<<std::endl;
		GrowthRateFile >> gridX;
		GrowthRateFile >> gridY;
		float rate;
        //cout<<"initiating growth and angle matrices"<<endl;

        for (size_t i= 0; i<gridX; ++i){
            std::vector<std::array<double,3>> tmpGridMatrixY(gridY,std::array<double,3>{0.0});
            GrowthMatrix.push_back(tmpGridMatrixY);
            std::vector<double> tmpShearY(gridY,0.0);
            AngleMatrix.push_back(tmpShearY);
		}
		double timeMultiplier = 1.0 / 3600.0; // converting rate per hour to rate per second
        //cout<<"reading growth matrix"<<endl;
        for (int j=gridY-1; j>-1; --j){
            for (size_t i=0; i<gridX; ++i){
                for (size_t k=0; k<3; ++k){
					//cout<<"i :"<<i<<" j: "<<j<<" k: "<<k<<" ";
					GrowthRateFile >> rate;
					//cout<<"rate: "<<rate<<" ";
					//GrowthMatrix[i][j][k] = rate*timeMultiplier;
					GrowthMatrix[i][j][k] = rate;
					//cout<<"matrix value: "<<GrowthMatrix[i][j][k]<<endl;
				}
			}
		}
		double angle;
        for (int j=gridY-1; j>-1; --j){
            for (size_t i=0; i<gridX; ++i){
				GrowthRateFile >> angle;
				AngleMatrix[i][j] = angle; // angles in degrees!
			}
		}
		GrowthRateFile.close();
        for (auto &grid_i : GrowthMatrix){
            for (auto &grid_j : grid_i){
                for (auto &currCoord : grid_j){
                    currCoord *= timeMultiplier;
				}
			}
		}
	}
	else{
        printErrorMessage(currHeader,"Growth options","Filename(full-path):");
		return false;
	}
	double zMin, zMax;
	file >> currHeader;
	if(currHeader == "zRange:"){
		file >> zMin;
		file >> zMax;
	}
	else{
        printErrorMessage(currHeader,"Growth type 2 options","zRange:");
        return false;
	}
    int Id = Sim->GrowthFunctions.size();
    //type is 3
    std::unique_ptr<GrowthFunctionBase> GSBp = std::make_unique<GridBasedGrowthFunction>(Id, 3, initialtime, finaltime, applyToColumnarLayer, applyToPeripodialMembrane, applyToBasalECM, applyToLateralECM, gridX, gridY, GrowthMatrix, AngleMatrix);
    GSBp->zMin = zMin;
    GSBp->zMax = zMax;
    Sim->GrowthFunctions.push_back(std::move(GSBp));
    return true;
}

bool ModelInputObject::readMeshParameters(std::ifstream& file){
    std::string currHeader;
	file >> currHeader;
	if(currHeader == "MeshInputMode(int-seeDocumentation):"){
		file >> Sim->MeshType;
	}
	else{
        printErrorMessage(currHeader,"reading mesh options","MeshInputMode(int-seeDocumentation):");
		return false;
	}
	if ( Sim->MeshType == 2 ){
		bool Success  = readMeshType2(file);
		if (!Success){
			return false;
		}
	}
	else if ( Sim->MeshType == 4){
		bool Success  = readMeshType4(file);
		if (!Success){
			return false;
		}
	}
	file >> currHeader;
	if(currHeader == "symmetricInX(bool):"){
		file >> Sim->symmetricX;
	}
	else{
        printErrorMessage(currHeader,"reading mesh options","symmetricInX(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "symmetricInY(bool):"){
		file >> Sim->symmetricY;
	}
	else{
        printErrorMessage(currHeader,"reading mesh options","symmetricInY(bool):");
		return false;
	}
	return true;
}

bool ModelInputObject::readExternalViscosityParameters(std::ifstream& file){
    std::string currHeader;
	file >> currHeader;
	if(currHeader == "ExtendToWholeTissue:"){
		file >>Sim->extendExternalViscosityToInnerTissue;
	}
	else{
        printErrorMessage(currHeader,"reading external viscosity options","ExtendToWholeTissue:");
		return false;
	}
	file >> currHeader;
	if(currHeader == "DiscProperApicalExternalViscosity:"){
		file >>Sim->externalViscosityDPApical;
	}
	else{
        printErrorMessage(currHeader,"reading external viscosity  options","DiscProperApicalExternalViscosity:");
        return false;
	}
	file >> currHeader;
	if(currHeader == "DiscProperBasalExternalViscosity:"){
		file >>Sim->externalViscosityDPBasal;
	}
	else{
        printErrorMessage(currHeader,"reading external viscosity  options","DiscProperBasalExternalViscosity:");
        return false;
	}
	file >> currHeader;
	if(currHeader == "PeripodialMembraneApicalExternalViscosity:"){
		file >>Sim->externalViscosityPMApical;
	}
	else{
        printErrorMessage(currHeader,"reading external viscosity  options","PeripodialMembraneApicalExternalViscosity:");
		return false;
	}
	file >> currHeader;
	if(currHeader == "PeripodialMembraneBasalExternalViscosity:"){
		file >>Sim->externalViscosityPMBasal;
	}
	else{
        printErrorMessage(currHeader,"reading external viscosity  options","PeripodialMembraneBasalExternalViscosity:");
		return false;
	}
	file >> currHeader;
	if(currHeader == "LinkerZoneApicalExternalViscosity:"){
		file >>Sim->externalViscosityLZApical;
	}
	else{
        printErrorMessage(currHeader,"reading external viscosity  options","LinkerZoneApicalExternalViscosity:");
		return false;
	}
	file >> currHeader;
	if(currHeader == "LinkerZoneBasalExternalViscosity:"){
		file >>Sim->externalViscosityLZBasal;
	}
	else{
        printErrorMessage(currHeader,"reading external viscosity  options","LinkerZoneBasalExternalViscosity:");
		return false;
	}
	return true;
}

void ModelInputObject::printErrorMessage(std::string currentInput, std::string sourceFuction, std::string expectedInput){
    std::cerr<<"Error in reading "<<sourceFuction<<" current input: "<<currentInput<<", should have been: "<<expectedInput<<endl;
}

bool ModelInputObject::readNodeBindingParameters(std::ifstream& file){
    std::string currHeader;
	file >> currHeader;
	if(currHeader == "bindCircumferenceXYToBasal(bool):"){
		file >>Sim->thereIsCircumferenceXYBinding;
	}
	else{
		printErrorMessage(currHeader,"NodeBinding options","bindCircumferenceXYToBasal(bool):");
		return false;
	}
	return true;
}



bool ModelInputObject::readNodeFixingParameters(std::ifstream& file){
    std::string currHeader;
	file >> currHeader;
	if(currHeader == "FixingViscosity(x,y,z):"){
		file >>Sim->fixingExternalViscosity[0];
		file >>Sim->fixingExternalViscosity[1];
		file >>Sim->fixingExternalViscosity[2];
	}
	else{
		printErrorMessage(currHeader,"Fixing options","FixingViscosity(x,y,z):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApicSurfaceFix(bool-x,y,z):"){
		file >>Sim->ApicalNodeFix[0];
		file >>Sim->ApicalNodeFix[1];
		file >>Sim->ApicalNodeFix[2];
	}
	else{
		printErrorMessage(currHeader,"Fixing options","ApicSurfaceFix(bool-x,y,z):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "FixApicalExtVisc(bool):"){
		file >>Sim->ApicalNodeFixWithExternalViscosity;
	}
	else{
		printErrorMessage(currHeader,"Fixing options","FixApicalExtVisc(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "BasalSurfaceFix(bool-x,y,z):"){
		file >>Sim->BasalNodeFix[0];
		file >>Sim->BasalNodeFix[1];
		file >>Sim->BasalNodeFix[2];
	}
	else{
		printErrorMessage(currHeader,"Fixing options","BasalSurfaceFix(bool-x,y,z):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "FixBasalExtVisc(bool):"){
		file >>Sim->BasalNodeFixWithExternalViscosity;
	}
	else{
		printErrorMessage(currHeader,"Fixing options","FixBasalExtVisc(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "CircumferenceFix(bool-x,y,z):"){
		file >>Sim->CircumferentialNodeFix[4][0];
		file >>Sim->CircumferentialNodeFix[4][1];
		file >>Sim->CircumferentialNodeFix[4][2];
	}
	else{
		printErrorMessage(currHeader,"Fixing options","CircumferenceFix(bool-x,y,z):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "FixCircWithExtVisc(bool):"){
		file >>Sim->CircumferentialNodeFixWithHighExternalViscosity[4];
	}
	else{
		printErrorMessage(currHeader,"Fixing options","FixCircWithExtVisc(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApicCircumFix(bool-x,y,z):"){
		file >>Sim->CircumferentialNodeFix[0][0];
		file >>Sim->CircumferentialNodeFix[0][1];
		file >>Sim->CircumferentialNodeFix[0][2];
	}
	else{
		printErrorMessage(currHeader,"Fixing options","ApicCircumFix(bool-x,y,z):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "FixApicCircWithExtVisc(bool):"){
		file >>Sim->CircumferentialNodeFixWithHighExternalViscosity[0];
	}
	else{
		printErrorMessage(currHeader,"Fixing options","FixApicCircWithExtVisc(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "BasalCircumFix(bool-x,y,z):"){
		file >>Sim->CircumferentialNodeFix[1][0];
		file >>Sim->CircumferentialNodeFix[1][1];
		file >>Sim->CircumferentialNodeFix[1][2];
	}
	else{
		printErrorMessage(currHeader,"Fixing options","BasalCircumFix(bool-x,y,z):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "FixBasalCircWithExtVisc(bool):"){
		file >>Sim->CircumferentialNodeFixWithHighExternalViscosity[1];
	}
	else{
		printErrorMessage(currHeader,"Fixing options","FixBasalCircWithExtVisc(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "LinkerApicCircumFix(bool-x,y,z):"){
		file >>Sim->CircumferentialNodeFix[2][0];
		file >>Sim->CircumferentialNodeFix[2][1];
		file >>Sim->CircumferentialNodeFix[2][2];
	}
	else{
		printErrorMessage(currHeader,"Fixing options","LinkerApicCircumFix(bool-x,y,z):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "FixLinkerApicCircWithExtVisc(bool):"){
		file >>Sim->CircumferentialNodeFixWithHighExternalViscosity[2];
	}
	else{
		printErrorMessage(currHeader,"Fixing options","FixLinkerApicCircWithExtVisc(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "LinkerBasalCircumFix(bool-x,y,z):"){
		file >>Sim->CircumferentialNodeFix[3][0];
		file >>Sim->CircumferentialNodeFix[3][1];
		file >>Sim->CircumferentialNodeFix[3][2];
	}
	else{
		printErrorMessage(currHeader,"Fixing options","LinkerBasalCircumFix(bool-x,y,z):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "FixLinkerBasalCircWithExtVisc(bool):"){
		file >>Sim->CircumferentialNodeFixWithHighExternalViscosity[3];
	}
	else{
		printErrorMessage(currHeader,"Fixing options","FixLinkerBasalCircWithExtVisc(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "NotumFix(bool-x,y,z,double-xFracMin,xFracMax):"){
		file >>Sim->NotumNodeFix[0];
		file >>Sim->NotumNodeFix[1];
		file >>Sim->NotumNodeFix[2];
		file >>Sim->notumFixingRange[0];
		file >>Sim->notumFixingRange[1];
	}
	else{
		printErrorMessage(currHeader,"Fixing options","NotumFix(bool-x,y,z,double-xFracMin,xFracMax):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "FixNotumExtVisc(bool):"){
		file >>Sim->NotumNodeFixWithExternalViscosity;
	}
	else{
		printErrorMessage(currHeader,"Fixing options","FixNotumExtVisc(bool):");
		return false;
	}
	return true;
}

bool ModelInputObject::readMeshType4(std::ifstream& file){
    std::string currHeader;
	file >> currHeader;
	if(currHeader == "MeshFile(full-path):"){
		file >> Sim->inputMeshFileName;
	}
	else{
        printErrorMessage(currHeader,"reading mesh path","MeshFile(full-path):");
        return false;
	}
	return true;
}


bool ModelInputObject::readMeshType2(std::ifstream& file){
    std::string currHeader;
	file >> currHeader;
	if(currHeader == "MeshRow(int):"){
		file >> Sim->Row;
	}
	else{
        printErrorMessage(currHeader,"reading mesh type 2","MeshRow(int):");
        return false;
	}

	file >> currHeader;
	if(currHeader == "MeshColumn(int):"){
		file >> Sim->Column;
	}
	else{
        printErrorMessage(currHeader,"reading mesh type 2","MeshColumn(int):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "SideLength:"){
		file >> Sim->SideLength;
	}
	else{
        printErrorMessage(currHeader,"reading mesh type 2","SideLength:");
		return false;
	}
	file >> currHeader;
	if(currHeader == "zHeight:"){
		file >> Sim->zHeight;
	}
	else{
        printErrorMessage(currHeader,"reading mesh type 2","zHeight:");
        return false;
	}
	//checking consistency:
	if (Sim->Column>Sim->Row-2){
		Sim->Column = Sim->Row-2;
		Sim->outputFile<<"Too few rows vs. columns, column count cannot be higher than Row-2"<<endl;
		Sim->outputFile<<"Updated to have a mesh: Row: "<<Sim->Row<<" Column: "<<Sim->Column<<endl;
	}
	float aspectratio = Sim->zHeight/Sim->SideLength;
	if ( aspectratio > 10 || aspectratio < 0.01 ){
		Sim->outputFile<<"Warning: The aspect ratio of the shapes are too high or low (aspectratio (z/side): "<<aspectratio<<endl;
	}
	return true;
}

bool ModelInputObject::readPeripodialMembraneParameters(std::ifstream& file){
    std::string currHeader;
	file >> currHeader;
	if(currHeader == "PeripodialMembraneYoungsModulus:"){
		file >>Sim->PeripodialElasticity;
	}
	else{
        printErrorMessage(currHeader,"reading peripodial options","PeripodialMembraneYoungsModulus:");
        return false;
	}
	file >> currHeader;
	if(currHeader == "PeripodialMembraneApicalViscosity:"){
		file >>Sim->peripodialApicalViscosity;
	}
	else{
        printErrorMessage(currHeader,"reading peripodial options","PeripodialMembraneApicalViscosity:");
		return false;
	}
	file >> currHeader;
	if(currHeader == "PeripodialMembraneBasalViscosity:"){
		file >>Sim->peripodialBasalViscosity;
	}
	else{
        printErrorMessage(currHeader,"reading peripodial options","PeripodialMembraneBasalViscosity:");
		return false;
	}
	file >> currHeader;
	if(currHeader == "PeripodialMembraneMidlineViscosity:"){
		file >>Sim->peripodialMidlineViscosity;
	}
	else{
        printErrorMessage(currHeader,"reading peripodial options","PeripodialMembraneMidlineViscosity:");
		return false;
	}
	return true;
}

bool ModelInputObject::readLinkerZoneParameters(std::ifstream& file){
    std::string currHeader;
	file >> currHeader;
	if(currHeader == "BaseOnPeripodialness(bool):"){
		file >>Sim->BaseLinkerZoneParametersOnPeripodialness;
	}
	else{
        printErrorMessage(currHeader,"reading linker zone options","BaseOnPeripodialness(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "LinkerZoneApicalYoungsModulus:"){
		file >>Sim->LinkerZoneApicalElasticity;
	}
	else{
        printErrorMessage(currHeader,"reading linker zone options","LinkerZoneApicalYoungsModulus:");
        return false;
	}
	file >> currHeader;
	if(currHeader == "LinkerZoneBasalYoungsModulus:"){
		file >>Sim->LinkerZoneBasalYoungsModulus;
	}
	else{
        printErrorMessage(currHeader,"reading linker zone options","LinkerZoneBasalYoungsModulus:");
		return false;
	}
	file >> currHeader;
	if(currHeader == "LinkerZoneApicalViscosity:"){
		file >>Sim->linkerZoneApicalViscosity;
	}
	else{
        printErrorMessage(currHeader,"reading linker zone options","LinkerZoneApicalViscosity:");
		return false;
	}
	file >> currHeader;
	if(currHeader == "LinkerZoneBasalViscosity:"){
		file >>Sim->linkerZoneBasalViscosity;
	}
	else{
        printErrorMessage(currHeader,"reading linker zone options","LinkerZoneBasalViscosity:");
        return false;
	}
	file >> currHeader;
	if(currHeader == "LinkerZoneMidlineViscosity:"){
		file >>Sim->linkerZoneMidlineViscosity;
	}
	else{
        printErrorMessage(currHeader,"reading linker zone options","LinkerZoneMidlineViscosity:");
		return false;
	}
	return true;
}

bool ModelInputObject::readTimeParameters(std::ifstream& file){
    std::string currHeader;

	file >> currHeader;
	if(currHeader == "TimeStep(sec):"){
		file >> Sim->dt;
	}
	else{
        printErrorMessage(currHeader,"reading simulation time options","TimeStep(sec):");
        return false;
	}
	file >> currHeader;
	if(currHeader == "SimulationLength(sec):"){
		file >> Sim->SimLength;
	}
	else{
        printErrorMessage(currHeader,"reading simulation time options","SimulationLength(sec):");
        return false;
	}
    std::cout<<"Simulation time step	: "<<Sim->dt<<std::endl;
    std::cout<<"Simulation Length	: "<<Sim->SimLength<<std::endl;
	return true;
}

bool ModelInputObject::readPysicalProperties(std::ifstream& file){
    std::string currHeader;
	file >> currHeader;
	if(currHeader == "YoungsModulusApical:"){
		file >> Sim->EApical;
	}
	else{
        printErrorMessage(currHeader,"reading physical property options","YoungsModulusApical:");
		return false;
	}
	file >> currHeader;
	if(currHeader == "YoungsModulusBasal:"){
		file >> Sim->EBasal;
	}
	else{
        printErrorMessage(currHeader,"reading physical property options","YoungsModulusBasal:");
        return false;
	}
	file >> currHeader;
	if(currHeader == "YoungsModulusMid:"){
		file >> Sim->EMid;
	}
	else{
        printErrorMessage(currHeader,"reading physical property options","YoungsModulusMid:");
        return false;
	}
	file >> currHeader;
	if(currHeader == "Noise(%-int):"){
		file >> Sim->noiseOnPysProp[0];
	}
	else{
        printErrorMessage(currHeader,"reading physical property options","Noise(%-int):");
        return false;
	}

	file >> currHeader;
	if(currHeader == "PoissonsRatio:"){
		file >> Sim->poisson;
	}
	else{
        printErrorMessage(currHeader,"reading physical property options","PoissonsRatio:");
		return false;
	}
	file >> currHeader;
	if(currHeader == "Noise(%-int):"){
		file >> Sim->noiseOnPysProp[1];
	}
	else{
        printErrorMessage(currHeader,"reading physical property options","Noise(%-int):");
        return false;
	}

	file >> currHeader;
	if(currHeader == "ApicalViscosity:"){
		file >> Sim->discProperApicalViscosity;
	}
	else{
        printErrorMessage(currHeader,"reading physical property options","ApicalViscosity:");
        return false;
	}
	file >> currHeader;
	if(currHeader == "Noise(%-int):"){
		file >> Sim->noiseOnPysProp[2];
	}
	else{
        printErrorMessage(currHeader,"reading physical property options","Noise(%-int):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "BasalViscosity:"){
		file >> Sim->discProperBasalViscosity;
	}
	else{
        printErrorMessage(currHeader,"reading physical property options","BasalViscosity:");
		return false;
	}
	file >> currHeader;
	if(currHeader == "MidLineViscosity:"){
		file >> Sim->discProperMidlineViscosity;
	}
	else{
        printErrorMessage(currHeader,"reading physical property options","MidLineViscosity:");
        return false;
	}
	return true;
}

bool ModelInputObject::readSaveOptions(std::ifstream& file){
    std::string currHeader;
	file >> currHeader;
	if(currHeader == "SaveImages(bool):"){
		file >> Sim->saveImages;
	}
	else{
        printErrorMessage(currHeader,"reading save options","SaveImages(bool):");
		return false;
	}
    file >> currHeader;
    if(currHeader == "ImagesSaveDirectory(relative_path_to_od):"){
        file >> Sim->saveScreenshotsDirectory;
        /** The input directory for saving the screenshots is relative, such that the user can re-use the same modelinput file
         * with ease, by changing  the exetutable input output directory. The path will be completed during image saving.
         * in MainWindow#takeScreenshot function.
         */
    }
    else{
        printErrorMessage(currHeader,"reading save options","ImagesSaveDirectory(relative_path_to_od):");
        return false;
    }
	file >> currHeader;
    if(currHeader == "SaveData(bool):"){
		file >> Sim->saveData;
	}
	else{
        printErrorMessage(currHeader,"reading save options","SaveData(bool):");
        return false;
	}

	file >> currHeader;
	if(currHeader == "ImageSaveInterval(sec):"){
		float timeInSec;
		file >> timeInSec;
		Sim->imageSaveInterval = timeInSec/Sim->dt;
        //the image save interval cannot be smaller than time step!, if this is the case,
        //dataSaveInterval will be zero! Check and correct if necessary
        if (Sim->imageSaveInterval < 1 ){
          Sim->imageSaveInterval =1;
        }
	}
	else{
        printErrorMessage(currHeader,"reading save options","ImageSaveInterval(sec):");
		return false;
	}

	file >> currHeader;
	if(currHeader == "DataSaveInterval(sec):"){
		float timeInSec;
		file >> timeInSec;
		Sim->dataSaveInterval = timeInSec/Sim->dt;
        //the data save interval cannot be smaller than time step!, if this is the case,
        //dataSaveInterval will be zero! Check and correct if necessary
        if (Sim->dataSaveInterval < 1 ){
          Sim->dataSaveInterval =1;
        }
	}
	else{
        printErrorMessage(currHeader,"reading save options","DataSaveInterval(sec):");
		return false;
	}
	return true;
}

bool ModelInputObject::readShapeChangeOptions(std::ifstream& file){
    std::string currHeader;
	file >> currHeader;
	int n;
	if(currHeader == "NumberofShapeChangeFunctions(int):"){
		file >> n;
		Sim->nShapeChangeFunctions = n;
	}
	else{
        printErrorMessage(currHeader,"reading shape change options","NumberofShapeChangeFunctions(int):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "ShapeChangeStartsBelowECMLevel(fraction):"){
			file >> Sim->shapeChangeECMLimit;
		}
		else{
			printErrorMessage(currHeader,"read shape change options","ShapeChangeStartsBelowECMLevel(fraction)");
			return false;
		}
	for (int i = 0; i<n; ++i){
		file >> currHeader;
		int type;
		if(currHeader == "ShapeChangeFunctionType(int-seeDocumentation):"){
			file >> type;
		}
		else{
            printErrorMessage(currHeader,"read shape change options","ShapeChangeFunctionType(int-seeDocumentation)");
            return false;
		}
		if (type == 1){
			bool success = readShapeChangeType1(file);
			if (!success){
				return false;
			}
		}
		else if (type ==2){
			bool success = readShapeChangeType2(file);
			if (!success){
				return false;
			}
		}
		else{
            printErrorMessage(std::to_string(type),"read shape change options","please enter a valid type: {1},{2}");
            return false;
		}
	}
	return true;
}

bool ModelInputObject::readPlasticDeformationOptions(std::ifstream& file){
    std::string currHeader;
	file >> currHeader;
	if(currHeader == "ThereIsPlasticDeformation(bool):"){
		file >> Sim->thereIsPlasticDeformation;
	}
	else{
        printErrorMessage(currHeader,"read splastic deformation options","ThereIsPlasticDeformation(bool):");
        return false;
	}


	file >> currHeader;
        if(currHeader == "ApplyToColumnarLayer(bool):"){
				file >> Sim->plasticDeformationAppliedToColumnar;
		}
		else{
            printErrorMessage(currHeader,"read splastic deformation options","ApplyToColumnarLayer(bool):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "ApplyToPeripodialMembrane(bool):"){
				file >> Sim->plasticDeformationAppliedToPeripodial;;
		}
		else{
            printErrorMessage(currHeader,"read splastic deformation options","ApplyToPeripodialMembrane(bool):");
			return false;
		}

	file >> currHeader;
	if(currHeader == "VolumeConserved(bool):"){
		file >> Sim->volumeConservedInPlasticDeformation;
	}
	else{
        printErrorMessage(currHeader,"read splastic deformation options","VolumeConserved(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "DeformationHalfLife(hour):"){
		file >> Sim->plasticDeformationHalfLife;
		Sim->plasticDeformationHalfLife *= 3600; //converting to seconds.
	}
	else{
        printErrorMessage(currHeader,"read splastic deformation options","DeformationHalfLife(hour):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "zDeformationLimits(lowerFraction-upperFraction):"){
		file >> Sim->zRemodellingLowerThreshold;
		file >> Sim->zRemodellingUpperThreshold;
	}
	else{
        printErrorMessage(currHeader,"read splastic deformation options","zDeformationLimits(lowerFraction-upperFraction):");
        return false;
	}
	return true;
}

bool ModelInputObject::readShapeChangeType1(std::ifstream& file){
    std::string currHeader;
	file >> currHeader;
	float initialtime;
	float finaltime;
	bool applyToColumnarLayer = false;
	bool applyToPeripodialMembrane = false;
	float Rate;
	if(currHeader == "InitialTime(sec):"){
		file >> initialtime;
	}
	else{
        printErrorMessage(currHeader,"read shape change type 1 options","InitialTime(sec):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "FinalTime(sec):"){
		file >> finaltime;
	}
	else{
        printErrorMessage(currHeader,"read shape change type 1 options","FinalTime(sec):");
        return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToColumnarLayer(bool):"){
			file >> applyToColumnarLayer;
	}
	else{
        printErrorMessage(currHeader,"read shape change type 1 options","ApplyToColumnarLayer(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToPeripodialMembrane(bool):"){
			file >> applyToPeripodialMembrane;
	}
	else{
        printErrorMessage(currHeader,"read shape change type 1 options","ApplyToPeripodialMembrane(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "MaxValue(fractionPerHour):"){
		double timeMultiplier = 1.0 / 3600.0;
		file >> Rate;
		Rate  *= timeMultiplier;
	}
	else{
        printErrorMessage(currHeader,"read shape change type 1 options","MaxValue(fractionPerHour):");
		return false;
	}
    int Id = Sim->ShapeChangeFunctions.size();
    //type is 1
    unique_ptr<GrowthFunctionBase> GSBp = std::make_unique<UniformShapeChangeFunction>(Id, 1, initialtime, finaltime, applyToColumnarLayer, applyToPeripodialMembrane, false /*applyToBasalECM*/, false /*applyToLateralECM*/, 1, Rate);
    Sim->ShapeChangeFunctions.push_back(std::move(GSBp));
	return true;
}


bool ModelInputObject::readShapeChangeType2(std::ifstream& file){
    std::string currHeader;
	file >> currHeader;
	float initialtime;
	float finaltime;
	bool applyToBasalECM{false};
	bool applyToLateralECM{false};
	bool applyTissueApical{false};
	bool applyTissueBasal{false};
	bool applyTissueMidline{false};
	bool conserveVolume{false};
    std::vector <int> markerEllipses;
	double ShapeChangeFractionPerHr;
	if(currHeader == "InitialTime(sec):"){
		file >> initialtime;
	}
	else{
		printErrorMessage(currHeader,"shape change type 2","InitialTime(sec):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "FinalTime(sec):"){
		file >> finaltime;
	}
	else{
		printErrorMessage(currHeader,"shape change type 2","FinalTime(sec):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyTissueApical(bool):"){
			file >> applyTissueApical;
	}
	else{
		printErrorMessage(currHeader,"shape change type 2","ApplyTissueApical(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyTissueBasal(bool):"){
			file >> applyTissueBasal;
	}
	else{
		printErrorMessage(currHeader,"shape change type 2","ApplyTissueBasal(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyTissueMidline(bool):"){
			file >> applyTissueMidline;
	}
	else{
		printErrorMessage(currHeader,"shape change type 2","ApplyTissueMidline(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToBasalECM(bool):"){
			file >> applyToBasalECM;
	}
	else{
		printErrorMessage(currHeader,"shape change type 2","ApplyToBasalECM(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToLateralECM(bool):"){
			file >> applyToLateralECM;
	}
	else{
		printErrorMessage(currHeader,"shape change type 2","ApplyToLateralECM(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "ShapeChangeAppliedToEllipses(number,[ellipseId][ellipseId]):"){
		int n;
		file >> n;
		for (int i=0; i<n; ++i){
			int ellipseId;
			file >> ellipseId;
			markerEllipses.push_back(ellipseId);
		}
	}
	else{
		printErrorMessage(currHeader,"shape change type 2","ShapeChangeAppliedToEllipses(number,[ellipseId][ellipseId]):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "xyShapeChange(fractionPerHour):"){
		file >> ShapeChangeFractionPerHr;
	}
	else{
		printErrorMessage(currHeader,"shape change type 2","xyShapeChange(fractionPerHour)");
		return false;
	}
	file >> currHeader;
	if(currHeader == "ConserveVolume(bool):"){
		file >> conserveVolume;
	}
	else{
		printErrorMessage(currHeader,"shape change type 2","ConserveVolume(bool)");
		return false;
	}
    //GrowthFunctionBase* GSBp;
    //int Id = Sim->ShapeChangeFunctions.size();
	//type is 1
    //GSBp = new 	markerEllipseBasedShapeChangeFunction(Id, 2, initialtime, finaltime, applyTissueApical, applyTissueBasal, applyTissueMidline, applyToBasalECM, applyToLateralECM, 2, ShapeChangeFractionPerHr, markerEllipses, conserveVolume);
    //Sim->ShapeChangeFunctions.push_back(GSBp);

    int Id = Sim->ShapeChangeFunctions.size();
    //type is 2
    std::unique_ptr<GrowthFunctionBase> GSBp = std::make_unique<markerEllipseBasedShapeChangeFunction>(Id, 2, initialtime, finaltime, applyTissueApical, applyTissueBasal, applyTissueMidline, applyToBasalECM, applyToLateralECM, 2, ShapeChangeFractionPerHr, markerEllipses, conserveVolume);
    Sim->ShapeChangeFunctions.push_back(std::move(GSBp));
	return true;
}

bool ModelInputObject::readStretcherSetup(std::ifstream& file){
    std::string currHeader;
	file >> currHeader;
	if(currHeader == "StretcherAttached(bool):"){
		bool stretcherAttached;
		file >> stretcherAttached;
		Sim->stretcherAttached = stretcherAttached;
	}
	else{
        printErrorMessage(currHeader,"read stretcher options","StretcherAttached(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "ClampedOnDV(bool):"){
		file >> Sim->DVClamp;
	}
	else{
        printErrorMessage(currHeader,"read stretcher options","ClampedOnDV(bool):");
        return false;
	}
	file >> currHeader;
	if(currHeader == "InitialTime(sec):"){
		double inittime;
		file >> inittime;
		Sim->StretchInitialTime = inittime;
	}
	else{
        printErrorMessage(currHeader,"read stretcher options","InitialTime(sec):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "FinalTime(sec):"){
		double endtime;
		file >> endtime;
		Sim->StretchEndTime = endtime;
	}
	else{
        printErrorMessage(currHeader,"read stretcher options","FinalTime(sec):");
        return false;
	}
	file >> currHeader;
	if(currHeader == "DVClampMin:"){
		double ClampPos;
		file >> ClampPos;
		Sim->StretchMin = ClampPos;
	}
	else{
        printErrorMessage(currHeader,"read stretcher options","DVClampMin:");
		return false;
	}
	file >> currHeader;
	if(currHeader == "DVClampMax:"){
		double ClampPos;
		file >> ClampPos;
		Sim->StretchMax = ClampPos;
	}
	else{
        printErrorMessage(currHeader,"read stretcher options","DVClampMax:");
		return false;
	}
	file >> currHeader;
	if(currHeader == "MaxStrain:"){
		double MaxStrain;
		file >> MaxStrain;
		Sim->StretchStrain = MaxStrain;
	}
	else{
        printErrorMessage(currHeader,"read stretcher options","MaxStrain:");
		return false;
	}
	return true;
}


bool ModelInputObject::readMarkerEllipseBandOptions(std::ifstream& file){
    std::string currHeader;
	file >> currHeader;
	if(currHeader == "numberOfMarkerEllipses(int):"){
		file >> Sim->nMarkerEllipseRanges;
	}
	else{
        printErrorMessage(currHeader,"read marker ellipse options","numberOfMarkerEllipses(int):");
		return false;
	}
	
	file >> currHeader;
	if(currHeader == "MarkerEllipseXCenters(fractionOfTissueSize):"){
		double dummy;
		for (int i=0; i<Sim->nMarkerEllipseRanges; ++i){
			file >> dummy;	
			Sim->markerEllipseBandXCentres.push_back(dummy);
		}
	}
	else{
        printErrorMessage(currHeader,"read marker ellipse options","MarkerEllipseXCenters(fractionOfTissueSize):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "MarkerEllipseBandR1Ranges(fractionOfTissueSize-x1Low-x1High):"){
		double dummy;
		for (int i=0; i<Sim->nMarkerEllipseRanges; ++i){
			file >> dummy;	
			Sim->markerEllipseBandR1Ranges.push_back(dummy);
			file >> dummy;	
			Sim->markerEllipseBandR1Ranges.push_back(dummy);
		}
	}
	else{
        printErrorMessage(currHeader,"read marker ellipse options","MarkerEllipseBandR1Ranges(fractionOfTissueSize-x1Low-x1High):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "MarkerEllipseBandR2Ranges(fractionOfTissueSize-y2Low-y2High):"){
		double dummy;
		for (int i=0; i<Sim->nMarkerEllipseRanges; ++i){
			file >> dummy;	
			Sim->markerEllipseBandR2Ranges.push_back(dummy);
			file >> dummy;	
			Sim->markerEllipseBandR2Ranges.push_back(dummy);
		}
	}
	else{
        printErrorMessage(currHeader,"read marker ellipse options","MarkerEllipseBandR2Ranges(fractionOfTissueSize-y2Low-y2High):");
        return false;
	}
	return true;
}

bool ModelInputObject::readStiffnessPerturbation(std::ifstream& file){
    std::string currHeader;
	file >> currHeader;
	if(currHeader == "ThereIsStiffnessPerturbation(bool):"){
		file >> Sim->ThereIsStiffnessPerturbation;
	}
	else{
        printErrorMessage(currHeader,"stiffness perturbations","ThereIsStiffnessPerturbation(bool):");
		return false;
	}
	int nStiffnessFunctions;
	file >> currHeader;
	if(currHeader == "NumberOfStiffnessPerturbations(int):"){
		file >>nStiffnessFunctions;
	}
	else{
		printErrorMessage(currHeader,"stiffness perturbations","NumberOfStiffnessPerturbations(int):");
		return false;
	}

	for (int i=0; i<nStiffnessFunctions; ++i){
		file >> currHeader;
		if(currHeader == "ApplyToApically(bool):"){
			bool ThereIsApicalStiffnessPerturbation;
			file >>ThereIsApicalStiffnessPerturbation;
			Sim->ThereIsApicalStiffnessPerturbation.push_back(ThereIsApicalStiffnessPerturbation);
		}
		else{
			printErrorMessage(currHeader,"stiffness perturbations","ApplyToApically(bool):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "ApplyBasally(bool):"){
			bool ThereIsBasalStiffnessPerturbation;
			file >>ThereIsBasalStiffnessPerturbation;
			Sim->ThereIsBasalStiffnessPerturbation.push_back(ThereIsBasalStiffnessPerturbation);
		}
		else{
			printErrorMessage(currHeader,"stiffness perturbations","ApplyBasally(bool):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "ApplyToWholeTissue(bool):"){
			bool ThereIsWholeTissueStiffnessPerturbation;
			file >>ThereIsWholeTissueStiffnessPerturbation;
			Sim->ThereIsWholeTissueStiffnessPerturbation.push_back(ThereIsWholeTissueStiffnessPerturbation);
		}
		else{
			printErrorMessage(currHeader,"stiffness perturbations","ApplyToWholeTissue(bool):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "Basolateral(bool):"){
			bool basolateral;
			file >>basolateral;
			Sim->ThereIsBasolateralStiffnessPerturbation.push_back(basolateral);
		}
		else{
			printErrorMessage(currHeader,"stiffness perturbations","BasolateralWithApicalRelaxation");
			return false;
		}
		file >> currHeader;
		if(currHeader == "BasolateralWithApicalRelaxation(bool):"){
			bool basolateral;
			file >>basolateral;
			Sim->ThereIsBasolateralWithApicalRelaxationStiffnessPerturbation.push_back(basolateral);
		}
		else{
			printErrorMessage(currHeader,"stiffness perturbations","BasolateralWithApicalRelaxation");
			return false;
		}
		file >> currHeader;
		if(currHeader == "timeOfStiffeningPerturbation(hr):"){
			double timeInHr;
			file >> timeInHr;
			Sim->stiffnessPerturbationBeginTimeInSec.push_back(timeInHr*3600);
			file >> timeInHr;
			Sim->stiffnessPerturbationEndTimeInSec.push_back(timeInHr*3600);
		}
		else{
			printErrorMessage(currHeader,"stiffness perturbations","timeOfStiffeningPerturbation(hr):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "stiffnessPerturbationAppliedToEllipses(number,[ellipseId][ellipseId]):"){
			int numberOfStiffnessPerturbationAppliesEllipseBands;
			file >>numberOfStiffnessPerturbationAppliesEllipseBands;
			Sim->numberOfStiffnessPerturbationAppliesEllipseBands.push_back(numberOfStiffnessPerturbationAppliesEllipseBands);
			double ellipseBandId;
			Sim->stiffnessPerturbationEllipseBandIds.push_back(vector<int>(0));
			for (int aa=0; aa<Sim->numberOfStiffnessPerturbationAppliesEllipseBands[i]; ++aa){
				file >>ellipseBandId;
				Sim->stiffnessPerturbationEllipseBandIds[i].push_back(ellipseBandId);
				if (ellipseBandId >= 100){
					Sim->thereIsEmergentEllipseMarking = true;
				}
			}
		}
		else{
			printErrorMessage(currHeader,"stiffness perturbations","stiffnessPerturbationAppliedToEllipses(number,[ellipseId][ellipseId]):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "stiffnessChangedToFractionOfOriginal(double):"){
			double fraction;
			file >> fraction;
			if (fraction <=0.0) {
				fraction = 0.0001;
			}
			Sim->stiffnessChangedToFractionOfOriginal.push_back(fraction);
		}
		else{
			printErrorMessage(currHeader,"stiffness perturbations","stiffnessChangedToFractionOfOriginal(double):");
			return false;
		}
		Sim->startedStiffnessPerturbation.push_back(false);
	}
	return true;
}

bool ModelInputObject::readECMPerturbation(std::ifstream& file){
    std::string currHeader;
	file >> currHeader;
	if(currHeader == "ThereIsECMStiffnessChange(bool):"){
		file >> Sim->thereIsECMChange;
	}
    else{
        printErrorMessage(currHeader,"ECM perturbations","ThereIsECMStiffnessChange(bool):");
        return false;
	}

	int nECMFunctions;
	file >> currHeader;
	if(currHeader == "NumberOfECMPerturbations(int):"){
		file >>nECMFunctions;
	}
	else{
		printErrorMessage(currHeader,"ECM perturbations","NumberOfECMPerturbations(int):");
		return false;
	}
	for (int i=0; i<nECMFunctions; ++i){
		Sim->changedECM.push_back(false);
		file >> currHeader;
		if(currHeader == "ApplyToApicalECM(bool):"){
			bool changeApicalECM;
			file >> changeApicalECM;
			Sim->changeApicalECM.push_back(changeApicalECM);
		}
		else{
			printErrorMessage(currHeader,"ECM perturbations","ApplyToApicalECM(bool):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "ApplyToBasalECM(bool):"){
			bool changeBasalECM;
			file >> changeBasalECM;
			Sim->changeBasalECM.push_back(changeBasalECM);
		}
		else{
			printErrorMessage(currHeader,"ECM perturbations","ApplyToBasalECM(bool):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "AppliedElementsAreEmergent(bool):"){
			bool emergentApplication;
			file >> emergentApplication;
			Sim->ECMChangeTypeIsEmergent.push_back(emergentApplication);
			if (emergentApplication){
				Sim->thereIsEmergentEllipseMarking = true;
			}
		}
		else{
			printErrorMessage(currHeader,"ECM perturbations","AppliedElementsAreEmergent(bool):");
			return false;
		}

		file >> currHeader;
		if(currHeader == "timeOfStiffnessChange(hr):"){
			double timeInHr;
			file >> timeInHr;
			Sim->ECMChangeBeginTimeInSec.push_back(timeInHr*3600);
			file >> timeInHr;
			Sim->ECMChangeEndTimeInSec.push_back(timeInHr*3600);
		}
		else{
			printErrorMessage(currHeader,"ECM perturbations","timeOfStiffnessChange(hr):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "stiffnessChangeAppliedToEllipses(number,[ellipseId][ellipseId]):"){
			int numberOfECMChangeEllipseBands;
			file >> numberOfECMChangeEllipseBands;
			Sim->numberOfECMChangeEllipseBands.push_back(numberOfECMChangeEllipseBands);
			int ellipseBandId;
			Sim->ECMChangeEllipseBandIds.push_back(vector<int>(0));
			for (int aa=0; aa<Sim->numberOfECMChangeEllipseBands[i]; ++aa){
				file >>ellipseBandId;
				Sim->ECMChangeEllipseBandIds[i].push_back(ellipseBandId);
				if (ellipseBandId >= 100){
					Sim->thereIsEmergentEllipseMarking = true;
				}
			}
		}
		else{
			printErrorMessage(currHeader,"ECM perturbations","stiffnessChangeAppliedToEllipses(number,[ellipseId][ellipseId]):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "stiffnessChangeFraction(double(0-1.0)):"){
			double fraction;
			file >> fraction;
			if (fraction <=0.0) {
				fraction = 0.0001;
			}
			Sim->ECMStiffnessChangeFraction.push_back(fraction);
		}
		else{
			printErrorMessage(currHeader,"ECM perturbations","stiffnessChangeFraction(double(0-1.0)):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "ECMRenewalHalfLifeTargetFraction(double(0-1.0)):"){
			double fraction;
			file >> fraction;
			if (fraction <=0.0) {
				fraction = 0.0001;
			}
			Sim->ECMRenewalHalfLifeTargetFraction.push_back(fraction);
		}
		else{
			printErrorMessage(currHeader,"ECM perturbations","ECMRenewalHalfLifeTargetFraction(double(0-1.0)):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "ECMViscosityChangeFraction(double):"){
			double fraction;
			file >> fraction;
			if (fraction <=0.0) {
				fraction = 0.0001;
			}
			Sim->ECMViscosityChangeFraction.push_back(fraction);
		}
		else{
			printErrorMessage(currHeader,"ECM perturbations","ECMViscosityChangeFraction(double):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "changeNotumECM(time,fraction):"){
			file >> Sim->notumECMChangeInitTime;
			Sim->notumECMChangeInitTime*=3600; //convet to hrs
			file >> Sim->notumECMChangeEndTime;
			Sim->notumECMChangeEndTime*=3600; //convet to hrs
			file >> Sim->notumECMChangeFraction;
		}
		else{
			printErrorMessage(currHeader,"ECM perturbations","changeNotumECM(time,fraction):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "changeHingeECM(time,fraction):"){
			file >> Sim->hingeECMChangeInitTime;
			Sim->hingeECMChangeInitTime*=3600; //convet to hrs
			file >> Sim->hingeECMChangeEndTime;
			Sim->hingeECMChangeEndTime*=3600; //convet to hrs
			file >> Sim->hingeECMChangeFraction;
		}
		else{
			printErrorMessage(currHeader,"ECM perturbations","changeNotumECM(time,fraction):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "changePouchECM(time,fraction):"){
			file >> Sim->pouchECMChangeInitTime;
			Sim->pouchECMChangeInitTime*=3600; //convet to hrs
			file >> Sim->pouchECMChangeEndTime;
			Sim->pouchECMChangeEndTime*=3600; //convet to hrs
			file >> Sim->pouchECMChangeFraction;
		}
		else{
			printErrorMessage(currHeader,"ECM perturbations","changeNotumECM(time,fraction):");
			return false;
		}
	}
	return true;
}


bool ModelInputObject::readExplicitActinOptions(std::ifstream& file){
    std::string currHeader;
	file >> currHeader;
	if(currHeader == "ThereIsExplicitActin(bool):"){
		file >> Sim->thereIsExplicitActin;
	}
	else{
        printErrorMessage(currHeader,"explicit actin options","ThereIsExplicitActin(bool):");
		return false;
	}
	return true;
}


bool ModelInputObject::readColumnViseVolumeConservationOptions(std::ifstream& file){
    std::string currHeader;
	file >> currHeader;
	if(currHeader == "ThereIsColumnViseVolumeConservation(bool):"){
		file >> Sim->conservingColumnVolumes;
	}
	else{
        printErrorMessage(currHeader,"column-vise volume conservation options","ThereIsColumnViseVolumeConservation(bool):");
		return false;
	}
	return true;
}

bool ModelInputObject::readartificialRelaxationOptions(std::ifstream& file){
    std::string currHeader;
	file >> currHeader;
	if(currHeader == "ThereIsArtificaialRelaxation(bool):"){
		file >> Sim->thereIsArtificaialRelaxation;
	}
	else{
		printErrorMessage(currHeader,"artificial relaxation options","ThereIsArtificaialRelaxation(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "ArtificialRelaxationTime(sec):"){
		file >> Sim->artificialRelaxationTime;
	}
	else{
		printErrorMessage(currHeader,"artificial relaxation options","ArtificialRelaxationTime(sec):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "relaxECM(bool):"){
		file >> Sim->relaxECMInArtificialRelaxation;
	}
	else{
		printErrorMessage(currHeader,"artificial relaxation options","relaxECM(bool):");
		return false;
	}
	return true;
}

bool ModelInputObject::readEnclosementOptions(std::ifstream& file){
    std::string currHeader;
	file >> currHeader;
	if(currHeader == "thereIsEnclosementOfTheTissue(bool):"){
		file >> Sim->encloseTissueBetweenSurfaces;
	}
	else{
        printErrorMessage(currHeader,"enclosement options","thereIsEnclosementOfTheTissue(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "initialLimits(lowerBound,upperBound):"){
		file >> Sim->initialZEnclosementBoundaries[0];
		file >> Sim->initialZEnclosementBoundaries[1];
	}
	else{
        printErrorMessage(currHeader,"enclosement options","initialLimits(lowerBound,upperBound):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "finalLimits(lowerBound,upperBound):"){
		file >> Sim->finalZEnclosementBoundaries[0];
		file >> Sim->finalZEnclosementBoundaries[1];
	}
	else{
        printErrorMessage(currHeader,"enclosement options","finalLimits(lowerBound,upperBound):");
		return false;
	}

	file >> currHeader;
	if(currHeader == "initialTime(sec):"){
		file >> Sim->initialTimeToEncloseTissueBetweenSurfacesSec;
	}
	else{
        printErrorMessage(currHeader,"enclosement options","initialTime(sec):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "finalTime(sec):"){
		file >> Sim->finalTimeToEncloseTissueBetweenSurfacesSec;
	}
	else{
        printErrorMessage(currHeader,"enclosement options","finalTime(sec):");
		return false;
	}
	return true;
}

bool ModelInputObject::readMutationOptions(std::ifstream& file){
    std::string currHeader;
	file >> currHeader;
	if(currHeader == "numberOfClones(int):"){
		file >> Sim->numberOfClones;
	}
	else{
		printErrorMessage(currHeader,"Mutation options","numberOfClones(int):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "cloneInformation(double-relativeX,relativeY,micronRadius,usingAbsoluteGrowth(bool),growthRatePerHour_OR_growthFoldIncrease):"){
		double relativeX, relativeY, micronRadius, growthRateORFold,useAbsoluteGrowthRate;
        for (size_t i=0; i<Sim->numberOfClones; ++i){
			file >> relativeX;
			file >> relativeY;
			file >> micronRadius;
			file >> useAbsoluteGrowthRate;
			file >> growthRateORFold;
			Sim->cloneInformationX.push_back(relativeX);
			Sim->cloneInformationY.push_back(relativeY);
			Sim->cloneInformationR.push_back(micronRadius);
			Sim->cloneInformationUsingAbsolueGrowth.push_back(useAbsoluteGrowthRate);
			Sim->cloneInformationGrowth.push_back(growthRateORFold);
		}
	}
	else{
		printErrorMessage(currHeader,"Mutation options","cloneInformation(double-relativeX,relativeY,micronRadius,usingAbsoluteGrowth(bool),growthRatePerHour_OR_growthFoldIncrease):");
		return false;
	}
	return true;
}

bool ModelInputObject::readAdhesionOptions(std::ifstream& file){
    std::string currHeader;
	file >> currHeader;
	if(currHeader == "ThereIsAdhesion(bool):"){
		file >> Sim->thereIsAdhesion;
	}
	else{
		printErrorMessage(currHeader,"adhesion options","ThereIsAdhesion(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "CollapseNodesOnAdhesion(bool):"){
		file >> Sim->collapseNodesOnAdhesion;
	}
	else{
		printErrorMessage(currHeader,"adhesion options","CollapseNodesOnAdhesion(bool):");
		return false;
	}
	return true;
}

bool ModelInputObject::readNodeCollapseOptions(std::ifstream& file){
    std::string currHeader;
	file >> currHeader;
	if(currHeader == "ThereIsNodeCollapse(bool):"){
		file >> Sim->thereNodeCollapsing;
	}
	else{
		printErrorMessage(currHeader,"node collapse options","ThereIsNodeCollapse(bool):");
		return false;
	}
	return true;
}


bool ModelInputObject::readExplicitECMOptions(std::ifstream& file){
    std::string currHeader;
	file >> currHeader;
	if(currHeader == "ThereIsExplicitECM(bool):"){
		file >> Sim->thereIsExplicitECM;
	}
	else{
		printErrorMessage(currHeader,"explicit ECM options","ThereIsExplicitECM(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "ECMRemodellingHalfLife(hour):"){
		file >> Sim->ECMRenawalHalfLife;
		Sim->ECMRenawalHalfLife *= 3600; //converting to seconds.
	}
	else{
		printErrorMessage(currHeader,"explicit ECM options","ECMRemodellingHalfLife(hour):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "ECMColumnarYoungsModulus:"){
		file >> Sim->EColumnarECM;
	}
	else{
		printErrorMessage(currHeader,"explicit ECM options","ECMColumnarYoungsModulus:");
		return false;
	}
	file >> currHeader;
	if(currHeader == "ECMPeripodialYoungsModulus:"){
		file >> Sim->EPeripodialECM;
	}
	else{
		printErrorMessage(currHeader,"explicit ECM options","ECMPeripodialYoungsModulus:");
		return false;
	}
	return true;
}

bool ModelInputObject::readPipetteSetup(std::ifstream& file){
    std::string currHeader;
	file >> currHeader;
	if(currHeader == "PipetteAspitarionActive(bool):"){
		bool PipetteSuction;
		file >> PipetteSuction;
		Sim->PipetteSuction = PipetteSuction;
	}
	else{
        printErrorMessage(currHeader," pipette aspiration setup","PipetteAspitarionActive(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "NumberOfPressureStages(int):"){
		file >> Sim->nPipetteSuctionSteps;
	}
	else{
        printErrorMessage(currHeader," pipette aspiration setup","NumberOfPressureStages(int):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "InitiationTimes(sec):"){
		for (int i=0;i<Sim->nPipetteSuctionSteps;++i){
			double pressureInitiationTime;
			file >> pressureInitiationTime;
			Sim->pipetteSuctionTimes.push_back(pressureInitiationTime);
		}
	}
	else{
        printErrorMessage(currHeader," pipette aspiration setup","InitiationTimes(sec):");
		return false;
	}
	file >> currHeader;
    if(currHeader == "Pressures(Pa):"){
		for (int i=0;i<Sim->nPipetteSuctionSteps;++i){
			double pipetePressure;
			file >> pipetePressure;
			Sim->pipetteSuctionPressures.push_back(pipetePressure);
		}
	}
	else{
        printErrorMessage(currHeader," pipette aspiration setup","Pressures(Pa):");
        return false;
	}
	if (Sim->nPipetteSuctionSteps>0){
		Sim->PipetteInitialStep = Sim->pipetteSuctionTimes[0]/Sim->dt;
	}
	file >> currHeader;
	if(currHeader == "ApicalSuction(bool-will_set_up_basal_suction_if_false):"){
		file >> Sim->ApicalSuction;
	}
	else{
        printErrorMessage(currHeader," pipette aspiration setup","ApicalSuction(bool-will_set_up_basal_suction_if_false):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "TissueStuck(bool-will_fix_the_opposite_surface_in_z):"){
		file >> Sim->TissueStuckOnGlassDuringPipetteAspiration;
	}
	else{
        printErrorMessage(currHeader," pipette aspiration setup","TissueStuck(bool-will_fix_the_opposite_surface_in_z):");
		return false;
	}

	file >> currHeader;
	if(currHeader == "Centre_Position(x,y,z):"){
		double dummy;
		for (int j=0;j<3;++j){
			file >> dummy;
			Sim->pipetteCentre[j] = dummy;
		}
	}
	else{
        printErrorMessage(currHeader," pipette aspiration setup","Centre_Position(x,y,z):");
		return false;
	}
	file >> currHeader;
    if(currHeader == "Pipette_InnerRadius(micron):"){
		file >> Sim->pipetteInnerRadius;
	}
	else{
        printErrorMessage(currHeader," pipette aspiration setup","Pipette_InnerRadius(micron):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "Pipette_OuterRadius(micron):"){
		double pippetOuterRad;
		file >> pippetOuterRad;
		Sim->pipetteThickness = pippetOuterRad - Sim->pipetteInnerRadius;
	}
	else{
        printErrorMessage(currHeader," pipette aspiration setup","Pipette_OuterRadius(micron):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "Pipette_Effect_Depth(micron):"){
		double dummy;
		file >> dummy;
		Sim->pipetteDepth = dummy;
	}
	else{
        printErrorMessage(currHeader," pipette aspiration setup","Pipette_Effect_Depth(micron):");
		return false;
	}
	return true;
}
