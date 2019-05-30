#include "Simulation.h"
#include "Prism.h"
#include "RandomGenerator.h"
#include <string.h>
#include <vector>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

Simulation::Simulation(){
    savedStepwise = true;
	currElementId = 0;
    ModInp = std::make_unique<ModelInputObject>();
	SystemCentre[0]=0.0; SystemCentre[1]=0.0; SystemCentre[2]=0.0;
	TissueHeight = 0.0;
	TissueHeightDiscretisationLayers= 1;
	timestep = 0;
	currSimTimeSec = 0;
	reachedEndOfSaveFile = false;
	thereIsPeripodialMembrane = false;
	needPeripodialforInputConsistency = false;
	conservingColumnVolumes = false;
	boundingBoxSize[0]=1000.0; boundingBoxSize[1]=1000.0; boundingBoxSize[2]=1000.0;
	ContinueFromSave = false;
    growthRotationUpdateFrequency = 60.0/dt;
    nElements = 0;
    nNodes = 0;
    if (growthRotationUpdateFrequency<1) {growthRotationUpdateFrequency =1;}
	setDefaultParameters();
	implicitPacking = true;
	thereIsAdhesion = false;
	collapseNodesOnAdhesion = false;
	thereNodeCollapsing = false;
	thereIsEmergentEllipseMarking = false;
	thereIsArtificaialRelaxation = false;
	artificialRelaxationTime = -1;
	relaxECMInArtificialRelaxation = false;
	checkedForCollapsedNodesOnFoldingOnce = false;
}

Simulation::~Simulation(){
    std::cout<<"destructor for simulation called"<<std::endl;
}

void Simulation::setDefaultParameters(){
	dt = 0.01;				//sec
    SimLength = dt; 		//10 sec of simulation
	saveImages = false;		//do not save simulation images by default
	saveData = false;		//do not save simulation data by default
	imageSaveInterval = 60.0/dt;	//save images every minute
    dataSaveInterval  = 60.0/dt;	//save data every minute
	saveDirectory = "Not-Set";	//the directory to save the images and data points
    saveScreenshotsDirectory = "Not-Set";	//the directory to save the images and data points
    saveDirectoryToDisplayString  = "Not-Set"; //the file whcih will be read and displayed - no simulation
    EApical = 250.0;
    EBasal = 250.0;
    EMid = 250.0;
    PeripodialElasticity = 10000.0;
    EColumnarECM = 10000.0;
    EPeripodialECM = 10000.0;
	poisson = 0.3;
    discProperApicalViscosity = 10000.0;
    discProperBasalViscosity = 10000.0;
	for (int i=0; i<3; ++i){
		zeroExternalViscosity[i] = true;
	}
    noiseOnPysProp.fill(0);
	// The default input is a calculated mesh of width 4 elements, each element being 2.0 unit high
	// and having 1.0 unit sides of the triangles.
	MeshType = 2;
	Row = 4;
	Column = Row-2;
	SideLength=1.0;
	zHeight = 2.0;
	ApicalNodeFixWithExternalViscosity = false;
	BasalNodeFixWithExternalViscosity = false;
	NotumNodeFixWithExternalViscosity = false;
	for (int i=0; i<5; ++i){
		CircumferentialNodeFixWithHighExternalViscosity[i] = false;
		for (int j=0; j<3; j++){
			CircumferentialNodeFix[i][j] = false;
			ApicalNodeFix[j] = false;
			BasalNodeFix[j] = false;
			NotumNodeFix[j] = false;
			fixingExternalViscosity[j] = 0;
		}
	}
	notumFixingRange[0] = 1.1;
	notumFixingRange[1] = -0.1;
	nGrowthFunctions = 0;
	GridGrowthsPinnedOnInitialMesh = false;
	nGrowthPinning = 0;
	gridGrowthsInterpolationType = 1;
	nShapeChangeFunctions = 0;
	TensionCompressionSaved = true;
    GrowthSaved = true;
    GrowthRateSaved = true;
	ForcesSaved = true;
	PackingSaved = true;
	growthRedistributionSaved = true;
	nodeBindingSaved = true;
	physicalPropertiesSaved = true;
	collapseAndAdhesionSaved = true;
	PeripodialElasticity = 0.0;
	peripodialApicalViscosity = discProperApicalViscosity;
	peripodialBasalViscosity  = discProperBasalViscosity;
	dorsalTipIndex = 0;
	ventralTipIndex = 1;
	anteriorTipIndex = 0;
	posteriorTipIndex = 0;
	stretcherAttached = false;
	distanceIndex = false;
	StretchDistanceStep = 0.0;
	DVClamp = true;
	StretchInitialTime = -100;
	StretchEndTime = -100;
	PipetteSuction = false;
	PipetteInitialStep= -100;
	nPipetteSuctionSteps = 0;
	ApicalSuction = true;
	TissueStuckOnGlassDuringPipetteAspiration = true;
	pipetteCentre[0] = 0.0;
	pipetteCentre[1] = 0.0;
	pipetteCentre[2] = 0.0;
	pipetteDepth = 0.0;
	pipetteInnerRadius =0.0;
	SuctionPressure[0] = 0.0;
	SuctionPressure[1] = 0.0;
	SuctionPressure[2] = 0.0;
	pipetteInnerRadiusSq = pipetteInnerRadius*pipetteInnerRadius;
    pipetteThickness = 11.7;
    effectLimitsInZ[0] = pipetteCentre[2] - pipetteDepth;
    effectLimitsInZ[1] = pipetteCentre[2] + pipetteDepth;

	boundingBox[0][0] =  1000.0;	//left x
	boundingBox[0][1] =  1000.0;	//low y
	boundingBox[0][2] =  1000.0;	//bottom z
	boundingBox[1][0] = -1000.0;	//right x
	boundingBox[1][1] = -1000.0;	//high y
	boundingBox[1][2] = -1000.0;	//top z

	symmetricY = false;
	symmetricX = false;
	addingRandomForces = false;
	randomForceMean = 0.0;
	randomForceVar = 0.0;
	packingDetectionThreshold = 5.5;
	packingThreshold = 6.0;

	softPeriphery = false;
	softDepth = 0.0;
	softnessFraction = 0.0;
	softPeripheryBooleans[0] = false; //apical
	softPeripheryBooleans[1] = false; //basal
	softPeripheryBooleans[2] = false; //columnar
	softPeripheryBooleans[3] = false; //peripodial

	thereIsPlasticDeformation = false;
	plasticDeformationAppliedToPeripodial = false;
	plasticDeformationAppliedToColumnar = false;
	volumeConservedInPlasticDeformation = false;
	plasticDeformationHalfLife = 0.0;
	zRemodellingLowerThreshold = 0.5;
	zRemodellingUpperThreshold = 2.0;

	BaseLinkerZoneParametersOnPeripodialness = true;
	LinkerZoneApicalElasticity = 0.0;
	LinkerZoneBasalYoungsModulus = 0.0;
	linkerZoneApicalViscosity = discProperApicalViscosity;
	linkerZoneBasalViscosity = discProperBasalViscosity;

	extendExternalViscosityToInnerTissue = false;
	externalViscosityDPApical = 0.0;
	externalViscosityDPBasal  = 0.0;
	externalViscosityPMApical = 0.0;
	externalViscosityPMBasal  = 0.0;
	externalViscosityLZApical = 0.0;
	externalViscosityLZBasal  = 0.0;

	thereIsECMChange = false;

	thereIsExplicitECM = false;
	lateralECMThickness = 0.0;
	ECMRenawalHalfLife = 0.0;
	thereIsExplicitActin = false;

	nMarkerEllipseRanges = 0;
	ThereIsStiffnessPerturbation = false;

	encloseTissueBetweenSurfaces =  false;
	zEnclosementBoundaries[0] = -1000;
	zEnclosementBoundaries[1] = 1000;
	initialZEnclosementBoundaries[0] = -1000;
	initialZEnclosementBoundaries[1] = 1000;
	finalZEnclosementBoundaries[0] = -1000;
	finalZEnclosementBoundaries[1] = 1000;

	thereIsCircumferenceXYBinding= false;

	notumECMChangeInitTime = 10000000.0;
	notumECMChangeEndTime = 0.0;
	notumECMChangeFraction =1.0;
	hingeECMChangeInitTime = 10000000.0;
	hingeECMChangeEndTime = 0.0;
	hingeECMChangeFraction =1.0;
	pouchECMChangeInitTime = 10000000.0;
	pouchECMChangeEndTime = 0.0;
	pouchECMChangeFraction =1.0;

	boundLateralElements = false;
}

bool Simulation::readExecutableInputs(int argc, char **argv){
	int i = 1;
	bool Success = true;
	while(i<argc){
		const char *inptype = argv[i];
        if (std::string(inptype) == "-mode"){
            /**
             * The mode of simulation can be "DisplaySave","SimulationOnTheGo" or "ContinueFromSave". As the
             *nemas suggest, "DisplaySave" option will dsplay a simulation from saved files, without running the
             *simulation. "SimulationOnTheGo" will start a fresh simulation and "ContinueFromSave" will continue simulating
             *from a save file.\n
             *
              */
			Success = readModeOfSim(i, argc, argv);
		}
        else if (std::string(inptype) == "-i"){
            /**
             * The tag "-i" defines the input file, this should be a modelinput file
             * detailing the parameters of the simulation.\n
             */
			Success = readParameters(i, argc, argv);
		}
        else if (std::string(inptype) == "-od"){
            /**
             * If the input file requires saving, then an output direcotry should be specified with the "-od" tag.
             */
			Success = readOutputDirectory(i, argc, argv);
		}
        else if (std::string(inptype) == "-dInput"){
            /**
             * In case of simulations contining from save or when the tool is called to display an existing simulation,
             * the input directory should be specified, with the tag "-dInput".
             */
			Success = readSaveDirectoryToDisplay(i, argc, argv);
		}
		else {
            std::cerr<<"Please enter a valid option key: {-mode,-i, -od, -dInput}, current string: "<<inptype<<std::endl;
			return false;
		}
		i++;
		if (!Success){
			return Success;
		}
	}
	Success = checkInputConsistency();
	return Success;
}

bool Simulation::readModeOfSim(int& i, int argc, char **argv){
	i++;
	if (i >= argc){
        std::cerr<<" input the mode of simulation: {DisplaySave, SimulationOnTheGo, ContinueFromSave, Default}"<<std::endl;
		return false;
	}
	const char* inpstring = argv[i];
    if (std::string(inpstring) == "DisplaySave"){
		DisplaySave = true;
		return true;
	}
    else if (std::string(inpstring) == "SimulationOnTheGo" || string(inpstring) == "Default"){
		DisplaySave = false;
		return true;
	}
    else if (std::string(inpstring) == "ContinueFromSave"){
		ContinueFromSave = true;
		DisplaySave = false;
		return true;
	}
	else{
        std::cerr<<"Please provide input mode: -mode {DisplaySave, SimulationOnTheGo, ContinueFromSave, Default}";
		return false;
	}
}

bool Simulation::readParameters(int& i, int argc, char **argv){
	i++;
	if (i >= argc){
        std::cerr<<" input the model input file"<<std::endl;
		return false;
	}
    ModInp->Sim=shared_from_this();
    ModInp->parameterFileName =  argv[i];
	bool Success = ModInp->readParameters();
	if (!Success){
		return Success;
	}
	return true;
}

bool Simulation::readSaveDirectoryToDisplay(int& i, int argc, char **argv){
	i++;
	if (i >= argc){
        std::cerr<<" input the save directory, contents of which will be displayed"<<std::endl;
		return false;
	}
    const char* inpstring = argv[i];
    saveDirectoryToDisplayString = string(inpstring);
	return true;
}

bool Simulation::readOutputDirectory(int& i, int argc, char **argv){
    /**
     * This function will read in the save directory. The boolean for saving files will
     * not be toggles. If your model input file states no saving, then the error and output files
     * will be directed into this directory, but the frame saving must be toggled independently.
     */
    i++;
	if (i >= argc){
        std::cerr<<" input the save directory"<<std::endl;
		return false;
	}
    const char* inpstring = argv[i];
    saveDirectory= string(inpstring);
	return true;
}

bool Simulation::readFinalSimulationStep(){
	bool success  = openFilesToDisplay();
	if (!success){
		return false;
	}
    /**
     * The data save interval and time step from the model input file
     * are backed up, then the properties of the saved system are read.
     */
	double timeStepCurrentSim = dt;
	int dataSaveIntervalCurrentSim = dataSaveInterval;
    //reading system properties:
	success  = readSystemSummaryFromSave();
	if (!success){
		return false;
	}
	string currline;
	while(reachedEndOfSaveFile == false){
        /**
         * Each frame is then read until the end of Simulation is reached. For each fraom first
         * line is discarded as it is the header.
         */
		//skipping the header:
		getline(saveFileToDisplayMesh,currline);
        //std::cerr<<" currline in read last step: "<<currline<<std::endl;
		if(saveFileToDisplayMesh.eof()){
			reachedEndOfSaveFile = true;
			break;
		}
        /**
         * If end of file is not reached, the node data is read.
         */
		success = readNodeDataToContinueFromSave();
        std::cout<<"dt after  readNodeDataToContinueFromSave "<<dt<<" timeStepCurrentSim: "<<timeStepCurrentSim<<" dataSaveInterval: "<<dataSaveInterval<<" dataSaveIntervalCurrentSim: "<<dataSaveIntervalCurrentSim<<std::endl;
		if (!success){
			return false;
		}
        /**
         * If node data is successfully read, then element data is read.
         */
		readElementDataToContinueFromSave();
        /**
         * Depending on which properties are saved, all the trmainder physocal porperies and
         * states of the tissue are read.
         */
		if (TensionCompressionSaved){
			readTensionCompressionToContinueFromSave();
		}
        if (GrowthSaved){
            readGrowthToContinueFromSave();
        }
        if (GrowthRateSaved){
            readGrowthRateToContinueFromSave();
        }
        if (physicalPropertiesSaved){
        	readPhysicalPropToContinueFromSave();
        }
    	if (growthRedistributionSaved){
    		readGrowthRedistributionToContinueFromSave();
    	}
        //std::cout<<" reading node binding"<<std::endl;
    	if (nodeBindingSaved){
    		readNodeBindingToContinueFromSave();
    	}
        /**
         * The time step are iterated with the save time step read from the save summary.
         */
        timestep = timestep + dataSaveInterval;
        currSimTimeSec += dt*dataSaveInterval;
        std::cout<<"current time step "<<timestep<<" currSimTimeSec: "<<currSimTimeSec<<std::endl;
        /**
         * Finally the footer is skipped.
         */
		//skipping the footer:
		getline(saveFileToDisplayMesh,currline);
		while (currline.empty() && !saveFileToDisplayMesh.eof()){
			//skipping empty line
			getline(saveFileToDisplayMesh,currline);
		}
        /**
         * Now with all this information, we will save the file again, as such, we can keep all the Simulation
         * in a continuous file and dont work on appending folders later. This reduces significant amount of book keeping,
         * specifically for the simulations that are run on high throughput servers at muliple steps.
         */
		//Now I will save it again:
        std::cout<<"saving"<<std::endl;
		saveStep();
	}
    /**
     * One the last step is reached, the element positions, volumes, physical properties, nodal masses, and other geometry dependent properties such as
     * the shape function derivatives, and bounding box are updated.
     */
    std::cout<<"carry out updates"<<std::endl;
	updateElementVolumesAndTissuePlacements();
	updateElasticPropertiesForAllNodes();
	clearNodeMassLists();
	assignNodeMasses();
	assignConnectedElementsAndWeightsToNodes();
    calculateShapeFunctionDerivatives();
	updateElementPositions();
	calculateBoundingBox();
    /**
     * Finally, the data save interval and time step from the model input of the current simulation are
     * instated back.
     */
	dataSaveInterval = dataSaveIntervalCurrentSim;
	dt = timeStepCurrentSim;

	return true;
}

void Simulation::updateMasterSlaveNodesInBinding(){
    for (const auto& itNode : Nodes){
		for (int dim = 0; dim<3; ++dim){
            if (itNode->slaveTo[dim] > -1){
                int slaveDof = itNode->Id*3+dim;
                int masterDof = itNode->slaveTo[dim]*3+dim;
                std::vector <int> fix;
				fix.push_back(slaveDof);
				fix.push_back(masterDof);
				NRSolver->slaveMasterList.push_back(fix);
				NRSolver->boundNodesWithSlaveMasterDefinition = true;
			}
		}
	}
}

bool Simulation::checkInputConsistency(){
	if (ContinueFromSave &&  saveDirectoryToDisplayString == "Not-Set"){
        std::cerr <<"The mode is set to continue from saved simulation, please provide an input directory containing the saved profile, using -dInput"<<std::endl;
		return false;
	}
	if (saveData || saveImages){
		if (saveDirectory == "Not-Set"){
            std::cerr <<"Modelinput file requires saving, please provide output directory, using -od tag"<<std::endl;
			return false;
		}
	}
	if (DisplaySave){
		if (saveDirectoryToDisplayString == "Not-Set"){
            std::cerr <<"The mode is set to display from save, please provide an input directory, using -dInput"<<std::endl;
			return false;
		}
	}
    for (int i=0; i<nGrowthFunctions; ++i){
        if(GrowthFunctions[i]->applyToPeripodialMembrane){
            std::cerr<<"There is no peripodial membrane, while growth function "<<i<<" is applicable to peropodial membrane, further checks needed"<<std::endl;
            needPeripodialforInputConsistency = true;
        }
    }
	return true;
}

bool Simulation::initiateSystem(){
	bool Success = openFiles();
	if (!Success){
		return Success;
	}
    if (MeshType == 2){
        /**
         * If the mesh type is 2, then a mesh with selected number of rows and columns will be initiated.
         *This is the default option is no input file is provided to the system.
         */
		Success = initiateMesh(MeshType, Row, Column,  SideLength,  zHeight);
	}
	else if(MeshType == 4){
        /**
         * Most commonly used input mesh type requires an input mesh file.
         */
		Success = initiateMesh(MeshType);
	}
	if (!Success){
		return Success;
	}
	if (symmetricY || symmetricX){
        /**
         * If there is symetricity in the system, a checkpoint for circumferential node definitions
         * is in place. The symmetricity boundary should not be considered as at the outer circumference.
         */
		 clearCircumferenceDataFromSymmetricityLine();
	}
    /**
     * Then the existance of the peripodial membrane in the system is checked.
     */
	Success = checkIfThereIsPeripodialMembrane();
    /**
     * The tissue height is calculated and the tissue bounding box is obtained. Then the relative z positions
     * of the elements are calculated.
     */
	Success = calculateTissueHeight(); //Calculating how many layers the columnar layer has, and what the actual height is.
	calculateBoundingBox();
	assignInitialZPositions();
	if (!Success){
		return Success;
	}
    /**
     * At this point consistency is checked again to see if the model inputs necessiate a peripodial membane, and if they do, if
     * it has been implemented.
     */
	if (needPeripodialforInputConsistency){
		if (!thereIsPeripodialMembrane){
            std::cerr<<"There is no peripodial membrane but at least one growth function desires one"<<std::endl;
			Success = false;
		}
	}
    /**
     * If the model inputs enable an explicit definition of extracellular matrix, this is labelled on the
     * mesh here via Simulation#setUpECMMimicingElements.
     */
	if (thereIsExplicitECM){
		setUpECMMimicingElements();
	}
    /**
     * If the model inputs enable an explicit definition of an actin rich top layer, this is labelled on the
     * mesh here via Simulation#setUpActinMimicingElements.
     */
	if (thereIsExplicitActin){
		setUpActinMimicingElements();
	}
    /**
     * Then element neighbourhoods are filled in.
     */
	setBasalNeighboursForApicalElements();
	fillInNodeNeighbourhood();
	fillInElementColumnLists();
    /**
     * The simulation can fix degrees of freedom for certain nodes, as boundary condition. These
     * node fixing options are set up via function Simulation#checkForNodeFixing.
     */
	checkForNodeFixing();
    /**
     * The tip nodes are assigned via Simulation#assignTips.
     */
	assignTips();
	if (!Success){
		return Success;
	}
    /**
     * If all the initiation setup conpleted without errors up to this point, then the mesh is characterised with flags, and now
     * the actual physical parametrs are set up. First the system forces matrix is set with the sytem node size, via function
     * Simulation#initiateSystemForces. Then the system center is calcualted, and physical parameters are assigned to the elements via
     * Simulation#assignPhysicalParameters. \n
     */
	initiateSystemForces();
	calculateSystemCentre();
	assignPhysicalParameters();
    /**
     * Once the viscoelastic properties are set, thje system is checked against zero external viscosity. If there is no external
     * viscosity in the system, then additional boundry conditions should e implemented to reach a unique solution in the NR iteration to
     * solve for nodal displacements. Oncce viscosities are updated, the node massess, and the surfaces that are expoed the the external
     * friction are set in Simulation#assignNodeMasses, Simulation#assignElementalSurfaceAreaIndices, and
     * Simulation#assignConnectedElementsAndWeightsToNodes functions.
     */
	checkForZeroExternalViscosity();
	calculateShapeFunctionDerivatives();
	assignNodeMasses();
	assignElementalSurfaceAreaIndices();
	assignConnectedElementsAndWeightsToNodes();
    alignTissueDVToXPositive();
    /**
     * Then the boiunding box of the tissue is calculated and the relative positions are set.
     */
    calculateBoundingBox();
    calculateDVDistance();
    for (auto &itElement : Elements){
        itElement->calculateRelativePosInBoundingBox(boundingBox[0][0],boundingBox[0][1],boundingBoxSize[0],boundingBoxSize[1]);
    }
    updateRelativePositionsToApicalPositioning();
    for(const auto& itElement : Elements){
        itElement->setInitialRelativePosInBoundingBox();
    }
    /**
     * If there are any additional manipulations, they are implemented next. Thee include growth mutant induction, experiemtal stretcher attachement,
     * experimental pippete aspiration setup, and setting up symettic axes of the tissue to have fixed boundaries, and assigningn labelling by
     * ellipse bands for further manipulation if these are implemented.
     */
    induceClones();
	if (stretcherAttached){
		setStretch();
	}
    std::cout<<"setting the pipette"<<std::endl;
    //setUpAFM();
	if(PipetteSuction){
		setupPipetteExperiment();
	}

	if (symmetricY){
		setupYsymmetricity();
	}
	if (symmetricX){
		setupXsymmetricity();
	}
	assignIfElementsAreInsideEllipseBands();
    /**
     * If the simultion is continuing from an already saved setup (Simualtion#ContinueFromSave = true), the last
     * frame of the save is read via Simulation#readFinalSimulationStep.
     */
	if (ContinueFromSave){
        cout<<"Reading Final SimulationStep: "<<std::endl;
		Success = readFinalSimulationStep();
		if (!Success){
			return Success;
		}
	}
    /**
     * If data is being saved, then the simulation summary will be written next.
     */
	if (saveData){
        cout<<"writing the summary current simulation parameters"<<std::endl;
		writeSimulationSummary();
		writeSpecificNodeTypes();
	}
    /**
     * Fianlly, the Newton Raphson solver object is initiated, and node DoF bonding as boundary condition is checked.
     */
	nNodes = Nodes.size();
	nElements = Elements.size();
    NRSolver = make_unique<NewtonRaphsonSolver>(Nodes[0]->nDim,nNodes);
    if (ContinueFromSave){
    	updateMasterSlaveNodesInBinding();
    }
    checkForNodeBinding();
    std::cout<<" system initiated"<<std::endl;
	return Success;
}

void Simulation::checkForNodeBinding(){
    /**
     * If the model inputs require node binding to eliminate rotation at tissue boundary (Simulation#thereIsCircumferenceXYBinding = true),
     * then this is implemented in here. The nodes are updated via Simulation#bindCircumferenceXY, and the NR solver data is updated.
     */
	if (thereIsCircumferenceXYBinding){
		bool thereIsBinding = bindCircumferenceXY();
		if (thereIsBinding){
			NRSolver->boundNodesWithSlaveMasterDefinition = true;
		}
	}
}

bool Simulation::areNodesToCollapseOnLateralECM(int slaveNodeId, int masterNodeId){
	if(Nodes[slaveNodeId]->hasLateralElementOwner || Nodes[masterNodeId]->hasLateralElementOwner){
        cout<<"binding nodes: "<<slaveNodeId<<" "<<masterNodeId<<" on single element collapse, but will not collapse nodes, as they are too close on lateral element"<<std::endl;
		return true;
	}
	int nElement = Nodes[slaveNodeId]->connectedElementIds.size();
	for (int elementCounter = 0; elementCounter<nElement; elementCounter++){
		int elementId = Nodes[slaveNodeId]->connectedElementIds[elementCounter];
		if (Elements[elementId]->isECMMimimcingAtCircumference){
            cout<<"binding nodes: "<<slaveNodeId<<" "<<masterNodeId<<" on single element collapse, but will not collapse nodes, as they are too close on circumferential element - slave"<<std::endl;
			return true;
		}
	}
	nElement = Nodes[masterNodeId]->connectedElementIds.size();
	for (int elementCounter = 0; elementCounter<nElement; elementCounter++){
		int elementId = Nodes[masterNodeId]->connectedElementIds[elementCounter];
		if (Elements[elementId]->isECMMimimcingAtCircumference){
            cout<<"binding nodes: "<<slaveNodeId<<" "<<masterNodeId<<" on single element collapse, but will not collapse nodes, as they are too close on circumferential element - master"<<std::endl;
			return true;
		}
	}
	return false;
}

bool Simulation::checkEdgeLenghtsForBindingPotentiallyUnstableElements(){
    /**
     * This function checks the side lengths of elements, and if they are shrinking below the collapse threshold, the
     * nodes of the element are collapsed to avoid flipping.
     */
	bool thereIsBinding = false;
    size_t dim = 3;
	vector<int> masterIdsBulk,slaveIdsBulk;
	if (boundLateralElements == false){
        /**
          * During the check, the corner ECM elements are treated specially. These elements are very thin, they are prone to flips,
          * therefore their apical and basal surfaces are bound to convert them effectively into sheet elements.
          * They are bound without collapse as they are very thin. This is carried out once at
          * the beginning of the simulation, and the operation is flagged with Simulation#boundLateralElements boolean.
          */
        for (auto &itElement : Elements){
			boundLateralElements = true;
            //cout<<"checking element "<<(*itElement)->Id<<std::endl;
            /**
             * The length is checked with ShapeBase#checkEdgeLenghtsForBinding.
             */
            itElement->checkEdgeLenghtsForBinding(masterIdsBulk,slaveIdsBulk);
			int selectedPair[2] = {0,0};
            if(itElement->isECMMimimcingAtCircumference && itElement->tissuePlacement == 0){
                /** Then in the following loop the corner ECM elements are listed to be bound.
                 */
                const std::vector<int>& nodeIds = itElement->getNodeIds();
                std::array<double,3> vec {0.0};
				for (int idim =0; idim<3; idim++){
					vec[idim] = Nodes[nodeIds[0]]->Position[idim]-Nodes[nodeIds[1]]->Position[idim];
				}
                double L = itElement->normaliseVector3D(vec);
				if (L < 0.5){
					selectedPair[0] = 0;
					selectedPair[1] = 1;
				}
				else{
                    for (size_t idim =0; idim<3; idim++){
						vec[idim] = Nodes[nodeIds[0]]->Position[idim]-Nodes[nodeIds[2]]->Position[idim];
					}
                    double L = itElement->normaliseVector3D(vec);
					if (L < 0.5){
						selectedPair[0] = 0;
						selectedPair[1] = 2;
					}
					else{
						selectedPair[0] = 1;
						selectedPair[1] = 2;
					}
				}
				masterIdsBulk.push_back(nodeIds[selectedPair[0]]);
				slaveIdsBulk.push_back(nodeIds[selectedPair[1]]);
				masterIdsBulk.push_back(nodeIds[selectedPair[0]+3]);
				slaveIdsBulk.push_back(nodeIds[selectedPair[1]+3]);
                for (size_t i=1; i < TissueHeightDiscretisationLayers;++i){
                    int idOfElementOnSameColumn = itElement->elementsIdsOnSameColumn[i];
                    const std::vector<int>&	nodeIdsOfElementOnSameColumn = Elements[idOfElementOnSameColumn]->getNodeIds();
                    masterIdsBulk.push_back(nodeIdsOfElementOnSameColumn[selectedPair[0]]);
                    slaveIdsBulk.push_back(nodeIdsOfElementOnSameColumn[selectedPair[1]]);
                    masterIdsBulk.push_back(nodeIdsOfElementOnSameColumn[selectedPair[0]+3]);
                    slaveIdsBulk.push_back(nodeIdsOfElementOnSameColumn[selectedPair[1]+3]);
				}
			}
		}
	}
    /**
     * Once the potential collapse list is obtained, the duplicates are cleaned.
     */
    vector<int> masterIds,slaveIds;
    size_t n = masterIdsBulk.size();
    for (size_t i=0;i<n;++i){
		bool duplicate = false;
        size_t nMaster = masterIds.size();
        for (size_t j=0; j<nMaster; ++j){
			if (masterIdsBulk[i] == masterIds[j] && slaveIdsBulk[i] == slaveIds[j]){
				duplicate = true;
			}
			if (masterIdsBulk[i] == slaveIds[j] && slaveIdsBulk[i] == masterIds[j]){
				duplicate = true;
			}
		}
		if (!duplicate){
			masterIds.push_back(masterIdsBulk[i]);
			slaveIds.push_back(slaveIdsBulk[i]);
		}
	}
    /**
     * In the clean list, the node degree of freedom couples are checked to see if they are already bound, or already flagged as masters or slaves
     * to other degrees of freedoms. The collapse is managed by function Node#collapseOnNode.
     */
	//Now go through the list, check if the binding is feasible:
	n = masterIds.size();
    for (size_t idMasterSlaveCouple=0;idMasterSlaveCouple<n;++idMasterSlaveCouple){
		int masterNodeId = masterIds[idMasterSlaveCouple];
		int slaveNodeId = slaveIds[idMasterSlaveCouple];
		if (binary_search(Nodes[masterNodeId]->collapsedWith.begin(), Nodes[masterNodeId]->collapsedWith.end(),slaveNodeId)){
			//the couple is already collapsed;
			continue;
		}
		bool  nodesHaveLateralOwners = areNodesToCollapseOnLateralECM(slaveNodeId,masterNodeId);
		if(!nodesHaveLateralOwners){
            cout<<"the nodes do not have lateral owners, continue to collapse"<<std::endl;
			Nodes[slaveNodeId]->collapseOnNode(Nodes, masterNodeId);
		}
        for(size_t i=0; i<dim; ++i){
            /** The fixed node degree of freedoms of the master node are reflected on the slave.
             */
			if (Nodes[masterNodeId]->FixedPos[i]){
				Nodes[slaveNodeId]->FixedPos[i]=true;
			}
			//not using an else, as the slave could be fixed in given dimension independent of the master
			if (!Nodes[slaveNodeId]->FixedPos[i]){
				int dofmaster = masterNodeId*dim+i;
				int dofslave  = slaveNodeId*dim+i;
                /**
                 * If the slave is alrady slave to another node (Node#slaveTo is not set to -1, the default non-slave indice), then the master of this processed couple is set to
                 * be the slave of the existing master of the slave, and all three degree of freedoms are bound to each other.
                 */
				if (Nodes[slaveNodeId]->slaveTo[i] > -1){
                    /**
                     * One possiblity is that the potential slave is already a slave to the potential master, then nothing is needed to be done.
                     */
					if(Nodes[slaveNodeId]->slaveTo[i] == masterNodeId || Nodes[slaveNodeId]->slaveTo[i] == Nodes[masterNodeId]->slaveTo[i] ){
						continue;
					}
					dofslave = Nodes[slaveNodeId]->slaveTo[i]*dim+i;
					slaveNodeId = Nodes[slaveNodeId]->slaveTo[i];
				}
                /**
                 * If the potential master node is already a slave to another node, then the potential master is moved to
                 * the original master of the potential master DoF in function NewtonRaphsonSolver#checkMasterUpdate.
                 */
				NRSolver->checkMasterUpdate(dofmaster,masterNodeId);
                /**
                 * In this last check point, if the original couple was flipped such that the potential master vas the slave of the potential slave,
                 * then now both master and slave nodes are equal to the initial potential slave id ( as it has been the master prior to this step).
                 * Then again, nothing needs to be implemented, slave master coupling is already in place.
                 */
				if (dofmaster != dofslave){
                    /**
                     * If this is not the case and the implementation shoulf continue, the next check point is for duplicate couples, to check
                     * if this couple already exists, via  NewtonRaphsonSolver#checkIfCombinationExists.
                     */
					bool continueAddition =  NRSolver->checkIfCombinationExists(dofslave,dofmaster);
					if (continueAddition){
                        /**
                         * If the couple is not implemented, then the check for the status of the slave is carried out, if the slave is already master of
                         * other nodes, the coupling is checked and corrected in function NewtonRaphsonSolver#checkIfSlaveIsAlreadyMasterOfOthers. If the
                         * function updates the node couple, then the potential slave was a master, and now is a slave to the potential master.
                         * All slaves of the potential master are moved on to the potential master. \n
                         */
						bool madeChange = NRSolver->checkIfSlaveIsAlreadyMasterOfOthers(dofslave,dofmaster);
						if (madeChange){
                            for (size_t nodeIt = 0 ; nodeIt<nNodes; ++nodeIt){
								if(Nodes[nodeIt]->slaveTo[i]==slaveNodeId){
									Nodes[nodeIt]->slaveTo[i]=masterNodeId;
								}
							}
						}
						vector <int> fixDOF;
						fixDOF.push_back(dofslave);
						fixDOF.push_back(dofmaster);
						NRSolver->slaveMasterList.push_back(fixDOF);
						Nodes[slaveNodeId]->slaveTo[i] = masterNodeId;
						Nodes[masterNodeId]->isMaster[i] = true;
						thereIsBinding = true;
					}
				}
			}
		}
	}
	if(thereIsBinding){
        /**
         * After finalisation of all node binding checks, if there has been binding, then node positions are updated to bring the positions of
         * the collapsed nodes together.\m
         * Please note there is no hierarchy between master and slaves, the choice siply defines which array of the
         * sparse matrix will be utilised in solving the cumulative forces of all master-slave groups.
         */
		updateElementPositions();
	}
	return thereIsBinding;
}

bool Simulation::bindCircumferenceXY(){
    /**
     * This is a boundary condition setup function, where a special case of node binding is implemented. For each column of nodes art tissue
     * circumference, the x and y degrees of freedom of each node is bound. The most basal (bottom - min z) node is declared as the master of the
     * remaining nodes, and they slaves. Please note there is no hierarchy between master and slaves, the choice siply defines which array of the
     * sparse matrix will be utilised in solving the cumulative forces of all master-slave groups.
     */
	bool thereIsBinding = false;
    size_t dim = 3;
	vector <int> nodeIds;
    for (size_t i=0; i<nNodes; ++i){
		if (Nodes[i]->tissuePlacement == 0){
			//basal node
			if (Nodes[i]->atCircumference){
				//at basal circumference:
				nodeIds.push_back(i);
			}
		}
	}
    size_t n = nodeIds.size();
    for (size_t i=0; i<n; ++i){
		int masterNodeId = nodeIds[i];
		int dofXmaster = masterNodeId*dim;
		int dofYmaster = masterNodeId*dim+1;
		int currNodeId = masterNodeId;
		bool reachedTop = false;
		int a = 0;
        while (!reachedTop && a<5){
			a++;
			int nConnectedElements = Nodes[currNodeId]->connectedElementIds.size();
			for (int j=0;j<nConnectedElements;++j){
				int elementId = Nodes[currNodeId]->connectedElementIds[j];
				bool IsBasalOwner = Elements[elementId]->IsThisNodeMyBasal(currNodeId);
				if (IsBasalOwner){
					int slaveNodeId = Elements[elementId]->getCorrecpondingApical(currNodeId);
					if (Nodes[masterNodeId]->FixedPos[0]){
						//master node is fixed in X, slave needs to be fixed as well, and no need for binding calculations.
						Nodes[slaveNodeId]->FixedPos[0]=true;
					}
					if (Nodes[masterNodeId]->FixedPos[1]){
						Nodes[slaveNodeId]->FixedPos[1]=true;
					}
					//If the DoF is fixed rigidly by node fixing functions, then I
					//will not make it a slave of any node. Making it a slave leads to displacement
					//I would need to check for fixed nodes and correct my forces/jacobian twice.
					//First fix jacobian and forces, then bind nodes, then fix jacobian again to
					//clear up all the fixed nodes that were slaves and bound to other nodes.
					//If I do the binding first,without the first fixing I will carry the load of the fixed node onto the master node unnecessarily
					//Simpler and cleaner to not make them slaves at all.
					if (!Nodes[slaveNodeId]->FixedPos[0]){
						int dofXslave  = slaveNodeId*dim;
						vector <int> fixX;
						fixX.push_back(dofXslave);
						fixX.push_back(dofXmaster);
						NRSolver->slaveMasterList.push_back(fixX);
						Nodes[slaveNodeId]->slaveTo[0] = masterNodeId;
						Nodes[masterNodeId]->isMaster[0] = true;
						thereIsBinding = true;
					}
					if (!Nodes[slaveNodeId]->FixedPos[1]){
						int dofYslave  = slaveNodeId*dim+1;
						vector <int> fixY;
						fixY.push_back(dofYslave);
						fixY.push_back(dofYmaster);
						NRSolver->slaveMasterList.push_back(fixY);
						Nodes[slaveNodeId]->slaveTo[1] = masterNodeId;
						Nodes[masterNodeId]->isMaster[1] = true;
						thereIsBinding = true;
					}
					currNodeId = slaveNodeId;
					break;
				}
			}
			if (!thereIsPeripodialMembrane && Nodes[currNodeId]->tissuePlacement == 1){
				//there is no peripodial and I have reached apical nodes
				reachedTop = true;
			}
			if (thereIsPeripodialMembrane && Nodes[currNodeId]->tissuePlacement == 0){
				//there is peripodial and I have reached basal nodes again
				reachedTop = true;
			}
		}
	}
	return thereIsBinding;
}

void Simulation::assignNodeMasses(){
    /**
     * Masses to all nodes are asigned through thier owner elements via ShapeBase#assignVolumesToNodes.
     */
    for(auto& itElement : Elements){
        itElement->assignVolumesToNodes(Nodes);
	}
}

void Simulation::assignElementalSurfaceAreaIndices(){
    /**
     * Exposed surface indices to all elements are asigned via ShapeBase#assignExposedSurfaceAreaIndices.
     */
    //not using range based loops here to ensure openMP comaptibility
    #pragma omp parallel for
    for(std::vector<std::unique_ptr<ShapeBase>>::iterator itElement = Elements.begin(); itElement < Elements.end(); ++itElement){
        (*itElement)->assignExposedSurfaceAreaIndices(Nodes);
    }

}

void Simulation::assignConnectedElementsAndWeightsToNodes(){
    /**
     * All connected elements are asigned nodes via ShapeBase#assignElementToConnectedNodes.
     */
    for(auto& itElement : Elements){
        itElement->assignElementToConnectedNodes(Nodes);
	}
}

bool Simulation::openFiles(){
	bool Success;
	if (saveData){
		//Mesh information at each step:
		string saveFileString = saveDirectory +"/Save_Frame";
		const char* name_saveFileMesh = saveFileString.c_str();
		saveFileMesh.open(name_saveFileMesh, ofstream::out);
		if (saveFileMesh.good() && saveFileMesh.is_open()){
			Success = true;
		}
		else{
            std::cerr<<"could not open file: "<<name_saveFileMesh<<std::endl;
			Success = false;
		}

		saveFileString = saveDirectory +"Save_Summary";
		const char* name_saveFileSimulationSummary = saveFileString.c_str();
		saveFileSimulationSummary.open(name_saveFileSimulationSummary, ofstream::out);
		if (saveFileSimulationSummary.good() && saveFileSimulationSummary.is_open()){
			Success = true;
		}
		else{
            std::cerr<<"could not open file: "<<name_saveFileSimulationSummary<<std::endl;
			Success = false;
		}

		//tension compression information at each step:
		saveFileString = saveDirectory +"/Save_TensionCompression";
		const char* name_saveFileTenComp = saveFileString.c_str();
        cout<<"opening the file" <<name_saveFileTenComp<<std::endl;
		saveFileTensionCompression.open(name_saveFileTenComp, ofstream::binary);
		if (saveFileTensionCompression.good() && saveFileTensionCompression.is_open()){
			Success = true;
		}
		else{
            std::cerr<<"could not open file: "<<name_saveFileTenComp<<std::endl;
			Success = false;
		}

        //Growth information at each step (Fg):
        saveFileString = saveDirectory +"/Save_Growth";
        const char* name_saveFileGrowth = saveFileString.c_str();
        cout<<"opening the file" <<name_saveFileGrowth<<std::endl;
        saveFileGrowth.open(name_saveFileGrowth, ofstream::binary);
        if (saveFileGrowth.good() && saveFileGrowth.is_open()){
            Success = true;
        }
        else{
            std::cerr<<"could not open file: "<<name_saveFileGrowth<<std::endl;
            Success = false;
        }

        //Growth rate information at each step (rx,ry,rz, (as in exp(rx*dt) )for display purposes only):
        saveFileString = saveDirectory +"/Save_GrowthRate";
        const char* name_saveFileGrowthRate = saveFileString.c_str();
        cout<<"opening the file" <<name_saveFileGrowthRate<<std::endl;
        saveFileGrowthRate.open(name_saveFileGrowthRate, ofstream::binary);
        if (saveFileGrowthRate.good() && saveFileGrowthRate.is_open()){
            Success = true;
        }
        else{
            std::cerr<<"could not open file: "<<name_saveFileGrowthRate<<std::endl;
            Success = false;
        }
		//opening the physical property information file:
		saveFileString = saveDirectory +"/Save_PhysicalProp";
		const char* name_saveFilePhysicalProp = saveFileString.c_str();
		saveFilePhysicalProp.open(name_saveFilePhysicalProp, ofstream::binary);
		if (saveFilePhysicalProp.good() && saveFilePhysicalProp.is_open()){
			Success = true;
		}
		else{
            std::cerr<<"could not open file: "<<name_saveFilePhysicalProp<<std::endl;
			Success = false;
		}

		//opening the specific Node/Element Type file:
		saveFileString = saveDirectory +"/Save_SpecificElementAndNodeTypes";
		const char* name_saveFileSpecificType = saveFileString.c_str();
		saveFileSpecificType.open(name_saveFileSpecificType, ofstream::binary);
		if (saveFileSpecificType.good() && saveFileSpecificType.is_open()){
			Success = true;
		}
		else{
            std::cerr<<"could not open file: "<<name_saveFileSpecificType<<std::endl;
			Success = false;
		}

		//opeining the growth redistribution information file:
		saveFileString = saveDirectory +"/Save_GrowthRedistribution";
		const char* name_saveFileGrowRedist = saveFileString.c_str();
        cout<<"opening the file" <<name_saveFileGrowRedist<<std::endl;
		saveFileGrowthRedistribution.open(name_saveFileGrowRedist, ofstream::binary);
		if (saveFileGrowthRedistribution.good() && saveFileGrowthRedistribution.is_open()){
			Success = true;
		}
		else{
            std::cerr<<"could not open file: "<<name_saveFileGrowRedist<<std::endl;
			Success = false;
		}

		//node binding information at each step:
		saveFileString = saveDirectory +"/Save_NodeBinding";
		const char* name_saveFileNodeBind = saveFileString.c_str();
        cout<<"opening the file" <<name_saveFileNodeBind<<std::endl;
		//saveFileNodeBinding.open(name_saveFileNodeBind, ofstream::binary);
		saveFileNodeBinding.open(name_saveFileNodeBind, ofstream::out);
		if (saveFileNodeBinding.good() && saveFileNodeBinding.is_open()){
			Success = true;
		}
		else{
            std::cerr<<"could not open file: "<<name_saveFileNodeBind<<std::endl;
			Success = false;
		}
	}
	if (saveDirectory == "Not-Set"){
        std::cerr<<"Output directory is not set, outputting on Out file in current directory"<<std::endl;
		outputFileString = "./Out";
	}
	else {
		outputFileString = saveDirectory +"/Out";
	}
	const char* name_outputFile = outputFileString.c_str();
	outputFile.open(name_outputFile, ofstream::out);
	if (outputFile.good() && outputFile.is_open()){
		Success = true;
	}
	else{
        std::cerr<<"could not open file: "<<name_outputFile<<std::endl;
		Success = false;
	}
	return Success;
}

void Simulation::writeSimulationSummary(){
	saveFileSimulationSummary<<"TimeStep(sec):  ";
    saveFileSimulationSummary<<dt<<std::endl;
	saveFileSimulationSummary<<"DataSaveInterval(sec):  ";
    saveFileSimulationSummary<<dataSaveInterval*dt<<std::endl;
	saveFileSimulationSummary<<"ModelinputName:  ";
    saveFileSimulationSummary<<ModInp->parameterFileName<<std::endl<<std::endl;
	saveFileSimulationSummary<<"Mesh_Type:  ";
    saveFileSimulationSummary<<MeshType<<std::endl;
    saveFileSimulationSummary<<"	Symmetricity-x: "<<symmetricX<<" Symmetricity-y: "<<symmetricY<<std::endl;
	writeMeshFileSummary();
	writeGrowthRatesSummary();
	writeECMSummary();
	writeActinSummary();
	writeExperimentalSummary();
}

void Simulation::writeSpecificNodeTypes(){
    //Counting the number of elements and nodes for, actin elements, ECM mimicing elements
    size_t counterForActinMimicingElements = 0;
    size_t counterForECMMimicingElements = 0;

    for(const auto& itElement : Elements){
        if (itElement->isActinMimicing){
			counterForActinMimicingElements++;
		}
        if (itElement->isECMMimicing){
			counterForECMMimicingElements++;
		}
	}
	//Writing explicit actin layer:
	saveFileSpecificType.write((char*) &counterForActinMimicingElements, sizeof counterForActinMimicingElements);
    for (const auto& itElement : Elements){
        if (itElement->isActinMimicing){
            saveFileSpecificType.write((char*) & itElement->Id, sizeof itElement->Id);
		}
	}
	//Writing explicit ECM layer:
	saveFileSpecificType.write((char*) &counterForECMMimicingElements, sizeof counterForECMMimicingElements);
    for (const auto& itElement : Elements){
        if (itElement->isECMMimicing){
            saveFileSpecificType.write((char*) & itElement->Id, sizeof itElement->Id);
		}
	}
    saveFileSpecificType.close();
}

void Simulation::writeMeshFileSummary(){
	if ( MeshType == 2){
		saveFileSimulationSummary<<"	Row: ";
		saveFileSimulationSummary<<Row;
		saveFileSimulationSummary<<"	Column: ";
		saveFileSimulationSummary<<Column;
		saveFileSimulationSummary<<"	SideLength: ";
		saveFileSimulationSummary<<SideLength;
		saveFileSimulationSummary<<"	zHeight: ";
		saveFileSimulationSummary<<zHeight;
        saveFileSimulationSummary<<std::endl;
	}
	if ( MeshType == 4){
		saveFileSimulationSummary<<"	MeshFileName: ";
		saveFileSimulationSummary<<inputMeshFileName;
        saveFileSimulationSummary<<std::endl;
	}
}

void Simulation::writeGrowthRatesSummary(){
    saveFileSimulationSummary<<std::endl;
	for (int i=0; i<nGrowthFunctions; ++i){
		if(GrowthFunctions[i]->Type == 3){ //grid based function does not need dt in summary
			GrowthFunctions[i]->writeSummary(saveFileSimulationSummary);
		}
		else{
			GrowthFunctions[i]->writeSummary(saveFileSimulationSummary,dt);
		}
	}
}

void Simulation::writeECMSummary(){
    saveFileSimulationSummary<<std::endl;
    saveFileSimulationSummary<<"There is explicit ECM:	"<<thereIsExplicitECM<<std::endl;
	if (thereIsExplicitECM){
		writeECMProperties();
	}
}

void Simulation::writeECMProperties(){
    saveFileSimulationSummary<<"	Columnar ECM Stiffness(Pa): "<<EColumnarECM<<std::endl;
    saveFileSimulationSummary<<"	Peripodial ECM Stiffness(Pa): "<<EPeripodialECM<<std::endl;
    saveFileSimulationSummary<<"	ECM remodelling half life (hour): "<<ECMRenawalHalfLife/3600.0<<std::endl;
    saveFileSimulationSummary<<"	Is there perturbation to the ECM: "<<thereIsECMChange<<std::endl;
	if (thereIsECMChange){
		int n = ECMChangeBeginTimeInSec.size();
        saveFileSimulationSummary<<"there are "<<n<<" ECM perturbations"<<std::endl;
		for (int ECMperturbationIndex =0;ECMperturbationIndex<n; ++ECMperturbationIndex ){
            saveFileSimulationSummary<<"		stiffness alteration begins at  "<<ECMChangeBeginTimeInSec[ECMperturbationIndex]/3600.0<<" hrs and ends at "<<ECMChangeEndTimeInSec[ECMperturbationIndex]/3600.0<<std::endl;
            saveFileSimulationSummary<<"		final fraction of ECM stiffness  "<<ECMStiffnessChangeFraction[ECMperturbationIndex]<<" times original values."<<std::endl;
            saveFileSimulationSummary<<"		final fraction of ECM renewal time  "<<ECMRenewalHalfLifeTargetFraction[ECMperturbationIndex]<<" times original values."<<std::endl;
            saveFileSimulationSummary<<"		final fraction of ECM viscosity  "<<ECMViscosityChangeFraction[ECMperturbationIndex]<<" times original values."<<std::endl;
            saveFileSimulationSummary<<"		stiffness alteration applied to apical ECM: "	<<changeApicalECM[ECMperturbationIndex]<<std::endl;
            saveFileSimulationSummary<<"		stiffness alteration applied to basal  ECM: "	<<changeBasalECM[ECMperturbationIndex]<<std::endl;
			saveFileSimulationSummary<<"		stiffness alteration applied to ellipse bands: ";
			for (int i=0; i<numberOfECMChangeEllipseBands[ECMperturbationIndex]; ++i){
				saveFileSimulationSummary<<" "<<ECMChangeEllipseBandIds[ECMperturbationIndex][i];
			}
            saveFileSimulationSummary<<std::endl;
		}
	}
}

void Simulation::writeActinSummary(){
    saveFileSimulationSummary<<std::endl;
    saveFileSimulationSummary<<"There is explicit actin:  "<<thereIsExplicitActin<<std::endl;
}

void Simulation::writeExperimentalSummary(){
	writePipetteSumary();
}

void Simulation::writePipetteSumary(){
    saveFileSimulationSummary<<std::endl;
    saveFileSimulationSummary<<"There is pipette aspiration:  "<<PipetteSuction<<std::endl;
	if (PipetteSuction){
        saveFileSimulationSummary<<"	is the tissue stuck on the glass: "<<TissueStuckOnGlassDuringPipetteAspiration<<std::endl;
        saveFileSimulationSummary<<"	suction on apical side: "<<ApicalSuction<<std::endl;
        saveFileSimulationSummary<<"	pipette position [centre(x,y,z), effective pippette suction depth, microns]: "<<pipetteCentre[0]<<" "<<pipetteCentre[1]<<" "<<pipetteCentre[2]<<" "<<pipetteDepth<<std::endl;
        saveFileSimulationSummary<<"	pipette radia [inner & outer(microns)]: "<<pipetteInnerRadius<<" "<<pipetteInnerRadius+pipetteThickness<<std::endl;
        saveFileSimulationSummary<<"	Steps for pipette suction "<<nPipetteSuctionSteps<<std::endl;
		saveFileSimulationSummary<<"	suction times (sec): ";
		for (int i=0; i<nPipetteSuctionSteps; ++i){
			saveFileSimulationSummary<<pipetteSuctionTimes[i]<<" ";
		}
        saveFileSimulationSummary<<std::endl;
		saveFileSimulationSummary<<"	suction pressures (Pa): ";
		for (int i=0; i<nPipetteSuctionSteps; ++i){
			saveFileSimulationSummary<<pipetteSuctionPressures[i]<<" ";
		}
	}

}

bool Simulation::openFilesToDisplay(){
	string saveFileString = saveDirectoryToDisplayString +"/Save_Frame";
	const char* name_saveFileToDisplayMesh = saveFileString.c_str();
	saveFileToDisplayMesh.open(name_saveFileToDisplayMesh, ifstream::in);
	if (!(saveFileToDisplayMesh.good() && saveFileToDisplayMesh.is_open())){
        std::cerr<<"Cannot open the save file to display: "<<name_saveFileToDisplayMesh<<std::endl;
		return false;
	}
	saveFileString = saveDirectoryToDisplayString +"/Save_Summary";
	const char* name_saveFileToDisplaySimSum = saveFileString.c_str();;
	saveFileToDisplaySimSum.open(name_saveFileToDisplaySimSum, ifstream::in);
	if (!(saveFileToDisplaySimSum.good() && saveFileToDisplaySimSum.is_open())){
        std::cerr<<"Cannot open the simulation summary: "<<name_saveFileToDisplaySimSum<<std::endl;
		return false;
	}
	saveFileString = saveDirectoryToDisplayString +"/Save_SpecificElementAndNodeTypes";
	const char* name_saveSpecificElementAndNodeTypes = saveFileString.c_str();;
	saveFileToDisplaySpecificNodeTypes.open(name_saveSpecificElementAndNodeTypes, ifstream::in);
	specificElementTypesRecorded = true;
	if (!(saveFileToDisplaySpecificNodeTypes.good() && saveFileToDisplaySpecificNodeTypes.is_open())){
        std::cerr<<"Cannot open the specific node types: "<<name_saveSpecificElementAndNodeTypes<<std::endl;
		specificElementTypesRecorded = false;
		//return false;
	}
	saveFileString = saveDirectoryToDisplayString +"/Save_TensionCompression";
	const char* name_saveFileToDisplayTenComp = saveFileString.c_str();
	saveFileToDisplayTenComp.open(name_saveFileToDisplayTenComp, ifstream::in);
	if (!(saveFileToDisplayTenComp.good() && saveFileToDisplayTenComp.is_open())){
        std::cerr<<"Cannot open the save file to display: "<<name_saveFileToDisplayTenComp<<std::endl;
		TensionCompressionSaved = false;
	}

    saveFileString = saveDirectoryToDisplayString +"/Save_Growth";
    const char* name_saveFileToDisplayGrowth = saveFileString.c_str();
    saveFileToDisplayGrowth.open(name_saveFileToDisplayGrowth, ifstream::in);
    if (!(saveFileToDisplayGrowth.good() && saveFileToDisplayGrowth.is_open())){
        std::cerr<<"Cannot open the save file to display: "<<name_saveFileToDisplayGrowth<<std::endl;
        GrowthSaved = false;
    }

    saveFileString = saveDirectoryToDisplayString +"/Save_GrowthRate";
	const char* name_saveFileToDisplayGrowthRate = saveFileString.c_str();
	saveFileToDisplayGrowthRate.open(name_saveFileToDisplayGrowthRate, ifstream::in);
	if (!(saveFileToDisplayGrowthRate.good() && saveFileToDisplayGrowthRate.is_open())){
        std::cerr<<"Cannot open the save file to display: "<<name_saveFileToDisplayGrowthRate<<std::endl;
		GrowthRateSaved = false;
	}
	saveFileString = saveDirectoryToDisplayString +"/Save_PhysicalProp";
	const char* name_saveFileToDisplayPhysicalProp = saveFileString.c_str();;
	saveFileToDisplayPhysicalProp.open(name_saveFileToDisplayPhysicalProp, ifstream::in);
	if (!(saveFileToDisplayPhysicalProp.good() && saveFileToDisplayPhysicalProp.is_open())){
        std::cerr<<"Cannot open the save file to display: "<<name_saveFileToDisplayPhysicalProp<<std::endl;
		physicalPropertiesSaved = false;
	}
	saveFileString = saveDirectoryToDisplayString +"/Save_GrowthRedistribution";
	const char* name_saveFileToDisplayGrowRedist = saveFileString.c_str();
	saveFileToDisplayGrowthRedistribution.open(name_saveFileToDisplayGrowRedist, ifstream::in);
	if (!(saveFileToDisplayTenComp.good() && saveFileToDisplayTenComp.is_open())){
        std::cerr<<"Cannot open the save file to display: "<<name_saveFileToDisplayGrowRedist<<std::endl;
		growthRedistributionSaved = false;
	}
	saveFileString = saveDirectoryToDisplayString +"/Save_NodeBinding";
	const char* name_saveFileToDisplayNodeBinding = saveFileString.c_str();
	saveFileToDisplayNodeBinding.open(name_saveFileToDisplayNodeBinding, ifstream::in);
	if (!(saveFileToDisplayNodeBinding.good() && saveFileToDisplayNodeBinding.is_open())){
        std::cerr<<"Cannot open the save file to display: "<<name_saveFileToDisplayNodeBinding<<std::endl;
		nodeBindingSaved = false;
	}
	return true;
}

bool Simulation::initiateSavedSystem(){
    /**
     * This function follows the logic of the Simulation#initiateSystem function, please refer to its documentation for details.
     * The element and node properties are read from saveed files where available.
     */
	bool success  = openFilesToDisplay();
	if (!success){
		return false;
	}
	//reading system properties:
	success  = readSystemSummaryFromSave();
	if (!success){
		return false;
	}
	string currline;
	//skipping the header:
	getline(saveFileToDisplayMesh,currline);
	initiateNodesFromSave();
	initiateElementsFromSave();
	fillInNodeNeighbourhood();
	assignPhysicalParameters();
	initiateSystemForces();
	if (specificElementTypesRecorded){
		success  = readSpecificNodeTypesFromSave();
	}
	if (!success){
		return false;
	}
	if (TensionCompressionSaved){
		updateTensionCompressionFromSave();
	}

    if (GrowthSaved){
        updateGrowthFromSave();
    }
    if (GrowthRateSaved){
        updateGrowthRateFromSave();
    }
    cout<<" reading forces"<<std::endl;
	if (growthRedistributionSaved){
		updateGrowthRedistributionFromSave();
	}
	if (nodeBindingSaved){
		updateNodeBindingFromSave();
	}
	updateElementVolumesAndTissuePlacements();
	clearNodeMassLists();
	assignNodeMasses();
	assignConnectedElementsAndWeightsToNodes();
    calculateShapeFunctionDerivatives();
	updateElementPositions();
	//skipping the footer:
	getline(saveFileToDisplayMesh,currline);
	while (currline.empty() && !saveFileToDisplayMesh.eof()){
		//skipping empty line
		getline(saveFileToDisplayMesh,currline);
	}
	if(saveFileToDisplayMesh.eof()){
		reachedEndOfSaveFile = true;
		return true;
	}
    return true;
}

bool Simulation::readSystemSummaryFromSave(){
	string dummystring;
	if(saveFileToDisplaySimSum.eof()){
        std::cerr<<"reached the end of summary file, expecting: TimeStep(sec):"<<std::endl;
		return false;
	}
	saveFileToDisplaySimSum >> dummystring; //reading "TimeStep(sec):"
	if(saveFileToDisplaySimSum.eof()){
        std::cerr<<"reached the end of summary file, expecting: dt value"<<std::endl;
		return false;
	}
	saveFileToDisplaySimSum >> dt;

	if(saveFileToDisplaySimSum.eof()){
        std::cerr<<"reached the end of summary file, expecting: DataSaveInterval(sec):"<<std::endl;
		return false;
	}
	saveFileToDisplaySimSum >> dummystring; //reading "DataSaveInterval(sec):"
	double dummydouble;
	if(saveFileToDisplaySimSum.eof()){
        std::cerr<<"reached the end of summary file, expecting: save interval value"<<std::endl;
		return false;
	}
	saveFileToDisplaySimSum >> dummydouble;
	dataSaveInterval =  (int) ceil(dummydouble/dt);


	if(saveFileToDisplaySimSum.eof()){
        std::cerr<<"reached the end of summary file, expecting: ModelinputName:"<<std::endl;
		return false;
	}
	saveFileToDisplaySimSum >> dummystring; //reading "ModelinputName: "
	if(saveFileToDisplaySimSum.eof()){
        std::cerr<<"reached the end of summary file, expecting name of the model input file"<<std::endl;
		return false;
	}
	saveFileToDisplaySimSum >> dummystring; //reading the model input file
	if(saveFileToDisplaySimSum.eof()){
        std::cerr<<"reached the end of summary file, expecting: Mesh_Type:"<<std::endl;
		return false;
	}
	saveFileToDisplaySimSum >> dummystring; //reading "Mesh_Type: "
	if(saveFileToDisplaySimSum.eof()){
        std::cerr<<"reached the end of summary file, expecting: the mesh type(int)"<<std::endl;
		return false;
	}
	int dummyint;
	saveFileToDisplaySimSum >> dummyint; //reading the mesh type
	if(saveFileToDisplaySimSum.eof()){
        std::cerr<<"reached the end of summary file, expecting: Symmetricity-x:"<<std::endl;
		return false;
	}
	saveFileToDisplaySimSum >> dummystring; //reading "Symmetricity-x: "
	if(saveFileToDisplaySimSum.eof()){
        std::cerr<<"reached the end of summary file, expecting: symmetricitX boolean:"<<std::endl;
		return false;
	}
	saveFileToDisplaySimSum >> symmetricX;
	if(saveFileToDisplaySimSum.eof()){
        std::cerr<<"reached the end of summary file, expecting: Symmetricity-y:"<<std::endl;
		return false;
	}
	saveFileToDisplaySimSum >> dummystring; //reading "Symmetricity-y: "
	if(saveFileToDisplaySimSum.eof()){
        std::cerr<<"reached the end of summary file, expecting: symmetricitY boolean:"<<std::endl;
		return false;
	}
	saveFileToDisplaySimSum >> symmetricY;
    cout<<"read the summary, symmetricity data: "<<symmetricX<<" "<<symmetricY<<std::endl;
	return true;
}

bool Simulation::readSpecificNodeTypesFromSave(){
    int currElementId;
    size_t counterForActinMimicingElements;
	saveFileToDisplaySpecificNodeTypes.read((char*) &counterForActinMimicingElements, sizeof counterForActinMimicingElements);
	if (counterForActinMimicingElements>0){
		thereIsExplicitActin = true;
	}
    for (size_t i=0; i<counterForActinMimicingElements; i++){
		saveFileToDisplaySpecificNodeTypes.read((char*) &currElementId, sizeof currElementId);
        int currIndice = getElementArrayIndexFromId(currElementId);
		if (currIndice < 0){
            std::cout<<" error in reading actin mimicking elements "<<std::endl;
            std::cerr<<" error in reading actin mimicking elements "<<std::endl;
			return false;
		}
		Elements[currIndice]->isActinMimicing = true;
	}
	//read in ECM layer to display:
    size_t counterForECMMimicingElements;
	saveFileToDisplaySpecificNodeTypes.read((char*) &counterForECMMimicingElements, sizeof counterForECMMimicingElements);
	if (counterForECMMimicingElements>0){
		thereIsExplicitECM = true;
	}
    for (size_t i=0; i<counterForECMMimicingElements; i++){
		saveFileToDisplaySpecificNodeTypes.read((char*) &currElementId, sizeof currElementId);
		int currIndice = getElementArrayIndexFromId(currElementId);
		if (currIndice < 0){
            std::cout<<" error in reading ECM mimicking elements "<<std::endl;
            std::cerr<<" error in reading ECM mimicking elements "<<std::endl;
			return false;
		}
		Elements[currIndice]->isECMMimicing = true;
	}
	assigneElementsAtTheBorderOfECM();
	assigneElementsAtTheBorderOfActin();
	return true;
}

int Simulation::getElementArrayIndexFromId(int currId){
	int currIndice = 0;
    for (const auto& itElement : Elements){
        if (itElement->getId() == currId){
			return currIndice;
		}
		currIndice++;
	}
    std::cerr<<"Error in getElementArrayIndexFromId, required element Id : "<<currId<<" returnin -10, inducing a crash"<<std::endl;
	return -10;
}

int Simulation::getNodeArrayIndexFromId(int currId){
	int currIndice = 0;
    for (const auto& itNode : Nodes){
        if (itNode->getId() == currId){
			return currIndice;
		}
		currIndice++;
	}
    std::cerr<<"Error in getElementArrayIndexFromId, required node Id : "<<currId<<" returnin -10, inducing a crash"<<std::endl;
	return -10;
}

bool Simulation::readNodeDataToContinueFromSave(){
    /**
     * The order of the node data is : \n
     * [number of nodes] \n
     * [pos_x] [pos_y] [pos_z] [Node#tissuePlacement] [Node#tissueType] [Node#atCircumference]
     */
    size_t n;
	saveFileToDisplayMesh >> n;
    std::cout<<"number of nodes: "<<n<<std::endl;
	if(nNodes != n){
        std::cerr<<"The node number from save file("<<n<<") and model input("<<Nodes.size()<<") are not consistent - cannot continue simulation from save"<<std::endl;
		return false;
	}
    for (size_t i=0; i<n; ++i){
		int tissuePlacement, tissueType;
		bool atCircumference;
		saveFileToDisplayMesh >> Nodes[i]->Position[0];
		saveFileToDisplayMesh >> Nodes[i]->Position[1];
		saveFileToDisplayMesh >> Nodes[i]->Position[2];
		saveFileToDisplayMesh >> tissuePlacement;
		saveFileToDisplayMesh >> tissueType;
		saveFileToDisplayMesh >> atCircumference;
		if (Nodes[i]->tissuePlacement != tissuePlacement || Nodes[i]->atCircumference != atCircumference || Nodes[i]->tissueType != tissueType ){
            std::cerr<<"Node "<<i<<" properties are not consistent  - cannot continue simulation from save "<<std::endl;
			return false;
		}
	}
	return true;
}

void Simulation::initiateNodesFromSave(){
	saveFileToDisplayMesh >> nNodes;
    for (size_t i=0; i<nNodes; ++i){
        std::array<double,3> pos = {0.0, 0.0, 0.0};
        int tissuePlacement, tissueType;
		bool atCircumference;
		saveFileToDisplayMesh >> pos[0];
		saveFileToDisplayMesh >> pos[1];
		saveFileToDisplayMesh >> pos[2];
		saveFileToDisplayMesh >> tissuePlacement;
		saveFileToDisplayMesh >> tissueType;
		saveFileToDisplayMesh >> atCircumference;
        std::unique_ptr<Node> tmp_nd = make_unique<Node>(i, 3, pos,tissuePlacement, tissueType);
        tmp_nd-> atCircumference = atCircumference;
        Nodes.push_back(std::move(tmp_nd));
	}
}

void Simulation::initiateNodesFromMeshInput(){
    size_t n;
	inputMeshFile >> n;
    for (size_t i=0; i<n; ++i){
        std::array<double,3> pos = {0.0, 0.0, 0.0};
		int tissuePos = -2;
		int tissueType = -2;
		int atCircumference;
		inputMeshFile >> pos[0];
		inputMeshFile >> pos[1];
		inputMeshFile >> pos[2];
		inputMeshFile >> tissuePos;
		inputMeshFile >> tissueType;
		inputMeshFile >> atCircumference;
        std::unique_ptr<Node> tmp_nd = make_unique<Node>(i, 3, pos,tissuePos, tissueType);
        tmp_nd-> atCircumference = atCircumference;
        Nodes.push_back(std::move(tmp_nd));
		nNodes = Nodes.size();
	}
}

void Simulation::initiateElementsFromMeshInput(){
    size_t n;
	inputMeshFile >> n;
    for (size_t i=0; i<n; ++i){
		int shapeType;
		inputMeshFile >> shapeType;
		if (shapeType == 1){
			initiatePrismFromMeshInput();
		}
		else{
            std::cerr<<"Element "<<i<<" of "<<n<<": Error in shape type, corrupt save file! - currShapeType: "<<shapeType<<std::endl;
			break;
		}
	}
}


void Simulation::initiateElementsFromSave(){
    size_t n;
	saveFileToDisplayMesh >> n;
    cout<<"number of elements: "<<n<<std::endl;
    for (size_t i=0; i<n; ++i){
		int shapeType;
		saveFileToDisplayMesh >> shapeType;
		if (shapeType == 1){
			initiatePrismFromSave();
		}
		else{
            std::cerr<<"Error in shape type, corrupt save file! - currShapeType: "<<shapeType<<std::endl;
		}
	}
    cout<<"number of elements: "<<nElements<<std::endl;
}

bool Simulation::readElementDataToContinueFromSave(){    
    /**
     * The order of the element data is : \n
     * [number of elements] \n
     * [ShapeBase#ShapeType] [shape type specific data read by Simulation#readShapeData]
     */
    size_t n;
	saveFileToDisplayMesh >> n;
	if (nElements != n){
        std::cerr<<"The element number from save file and model input are not consistent - cannot continue simulation from save"<<std::endl;
        std::cerr<<"n: "<<n<<" nElements: "<<std::endl;
		return false;
	}
    for (size_t i=0; i<nElements; ++i){
		int shapeType;
		saveFileToDisplayMesh >> shapeType;
		if (Elements[i]->getShapeType() != shapeType){
            std::cerr<<"The element type from save file and model input are not consistent - cannot continue simulation from save"<<std::endl;
			return false;
		}
		if (shapeType == 1){
			bool success = readShapeData(i);
			if (!success){
                std::cerr<<"Error reading shape data, element: "<<i<<std::endl;
				return false;
			}
		}
		else if (shapeType == 4){
			double height = 0.0;
			saveFileToDisplayMesh >> height;
			if (Elements[i]->ReferenceShape->height != height){
                std::cerr<<"The element height from save file and model input are not consistent - cannot continue simulation from save"<<std::endl;
				return false;
			}
			bool success = readShapeData(i);
			if (!success){
                std::cerr<<"Error reading shape data, element: "<<i<<std::endl;
				return false;
			}
		}
		else{
            std::cerr<<"Error in shape type, corrupt save file! - currShapeType: "<<shapeType<<std::endl;
		}

	}
	return true;
}

void Simulation::initiatePrismFromSave(){
    /**
      * Inserts a new prism at order k into elements vector, ShapeBase#NodeIds and reference shape positions
      * ShapeBase#updateShapeFromSave
      */
    std::vector<int> NodeIds(6,0);
    std::unique_ptr<ShapeBase> PrismPnt01 = make_unique<Prism>(NodeIds, Nodes, currElementId);
    PrismPnt01->updateShapeFromSave(saveFileToDisplayMesh);
    Elements.push_back(std::move(PrismPnt01));
    nElements = Elements.size();
	currElementId++;
}

bool Simulation::readShapeData(int i){
    /**
     * The order of the element shape specific information is : \n
     * [ShapeBase#IsAblated] [node ids read by ShapeBase#readNodeIdData] [reference positions read by ShapeBase#readReferencePositionData]
     */
	bool IsAblated;
	saveFileToDisplayMesh >> IsAblated;
	if ( Elements[i]->IsAblated != IsAblated){
        std::cerr<<"The element "<<i<<" ablation from save file and model input are not consistent - cannot continue simulation from save"<<std::endl;
		return false;
	}
	bool success = Elements[i]->readNodeIdData(saveFileToDisplayMesh);
	if (!success){
        std::cerr<<"The element "<<i<<" node ids from save file and model input are not consistent - cannot continue simulation from save"<<std::endl;
		return false;
	}
	success = Elements[i]->readReferencePositionData(saveFileToDisplayMesh);
	if (!success){
        std::cerr<<"The element "<<i<<" reference shape from save file and model input are not consistent - cannot continue simulation from save"<<std::endl;
		return false;
	}
	return true;
}


void Simulation::initiatePrismFromMeshInput(){
    std::vector<int> NodeIds(6,0);
    for (size_t i =0 ;i<6; ++i){
		int savedId;
		inputMeshFile >> savedId;
		NodeIds[i] = savedId;
	}
    std::unique_ptr<ShapeBase> PrismPnt01 = make_unique<Prism>(NodeIds, Nodes, currElementId);
    PrismPnt01->updateReferencePositionMatrixFromMeshInput(inputMeshFile);
    PrismPnt01->checkRotationConsistency3D();
    Elements.push_back(std::move(PrismPnt01));
	nElements = Elements.size();
	currElementId++;
}



void Simulation::reInitiateSystemForces(){
    SystemForces.clear();
    PackingForces.clear();
    for (size_t j=0;j<nNodes;++j){
        SystemForces.push_back(std::array<double,3>{0.0});
        PackingForces.push_back(std::array<double,3>{0.0});
	}
}

void Simulation::updateTensionCompressionFromSave(){
    for(const auto&  itElement : Elements){
        for (size_t j=0; j<6; ++j){
            double S = gsl_matrix_get(itElement->Strain,j,0);
            saveFileToDisplayTenComp.read((char*) &S, sizeof S);
            gsl_matrix_set(itElement->Strain,j,0,S);
        }
	}
}

void Simulation::updateGrowthRedistributionFromSave(){
    for(const auto&  itElement : Elements){
		bool thereIsDistribution = false;
		bool shrinksElement = false;
		saveFileToDisplayGrowthRedistribution.read((char*) &thereIsDistribution, sizeof thereIsDistribution);
        itElement->thereIsGrowthRedistribution = thereIsDistribution;
		saveFileToDisplayGrowthRedistribution.read((char*) &shrinksElement, sizeof shrinksElement);
        itElement->growthRedistributionShrinksElement = shrinksElement;
	}
}

void Simulation::updateNodeBindingFromSave(){
    for(const auto& itNode : Nodes){
        itNode->slaveTo[0] = -1;
        itNode->slaveTo[1] = -1;
        itNode->slaveTo[2] = -1;
        itNode->isMaster[0] = false;
        itNode->isMaster[1] = false;
        itNode->isMaster[2] = false;
	}
    size_t n = 0; //size Of Master-Slave List
	saveFileToDisplayNodeBinding>>n;
    for (size_t i=0; i<n; ++i){
		int dofSlave, dofMaster;
		saveFileToDisplayNodeBinding>>dofSlave;
		saveFileToDisplayNodeBinding>>dofMaster;
		int dim = dofSlave % 3;
		int nodeSlave = (dofSlave - dim)/3;
		int nodeMaster = (dofMaster - dim)/3;
		Nodes[nodeSlave]->slaveTo[dim] = nodeMaster;
		Nodes[nodeMaster]->isMaster[dim] = true;
	}
}

void Simulation::updateGrowthFromSave(){
    for(const auto&  itElement : Elements){
		gsl_matrix* currFg = gsl_matrix_calloc(3,3);
        for (size_t j=0; j<3; ++j){
            for(size_t k=0; k<3; ++k){
                double Fgjk;
                saveFileToDisplayGrowth.read((char*) &Fgjk, sizeof Fgjk);
                gsl_matrix_set(currFg,j,k,Fgjk);
            }
        }
        itElement->setFg(currFg);
        gsl_matrix_free(currFg);
    }
}

void Simulation::updateGrowthRateFromSave(){
    /** For each element, the format is: \n
     * growth rate [r_x] [r_y] [r_z] \n
     * Then the shape information is updated with correct timestep length in ShapeBase#setGrowthRateViaInputTimeMultipliedMagnitude.
     */
    for(const auto&  itElement : Elements){
		double rx=0.0,ry=0.0,rz=0.0;
		saveFileToDisplayGrowthRate.read((char*) &rx, sizeof rx);
		saveFileToDisplayGrowthRate.read((char*) &ry, sizeof ry);
		saveFileToDisplayGrowthRate.read((char*) &rz, sizeof rz);
        itElement->setGrowthRateViaInputTimeMultipliedMagnitude(rx,ry,rz);
	}
}

void Simulation::updatePhysicalPropFromSave(){
	readPhysicalPropToContinueFromSave();
}

void Simulation::readTensionCompressionToContinueFromSave(){
    for(const auto&  itElement : Elements){
        for (size_t j=0; j<6; ++j){
            double S;
            saveFileToDisplayTenComp.read((char*) &S, sizeof S);
            gsl_matrix_set(itElement->Strain,j,0,S);
		}
	}
}

void Simulation::readGrowthRedistributionToContinueFromSave(){
    for(const auto&  itElement : Elements){
		bool thereIsDistribution = false;
		bool shrinksElement = false;
		saveFileToDisplayGrowthRedistribution.read((char*) &thereIsDistribution, sizeof thereIsDistribution);
        itElement->thereIsGrowthRedistribution = thereIsDistribution;
        saveFileToDisplayGrowthRedistribution.read((char*) &shrinksElement, sizeof shrinksElement);
        itElement->growthRedistributionShrinksElement = shrinksElement;
	}
}

void Simulation::readNodeBindingToContinueFromSave(){
	//clear all first:
    for(const auto& itNode : Nodes){
        itNode->slaveTo[0] = -1;
        itNode->slaveTo[1] = -1;
        itNode->slaveTo[2] = -1;
        itNode->isMaster[0] = false;
        itNode->isMaster[1] = false;
        itNode->isMaster[2] = false;
	}
    size_t n = 0; //size Of Master-Slave List
	saveFileToDisplayNodeBinding>>n;
    for (size_t i=0; i<n; ++i){
		int dofSlave, dofMaster;
		saveFileToDisplayNodeBinding>>dofSlave;
		saveFileToDisplayNodeBinding>>dofMaster;
		int dim = dofSlave % 3;
		int nodeSlave = (dofSlave - dim)/3;
		int nodeMaster = (dofMaster - dim)/3;
		Nodes[nodeSlave]->slaveTo[dim] = nodeMaster;
		Nodes[nodeMaster]->isMaster[dim] = true;
	}
}

void Simulation::readGrowthToContinueFromSave(){
	updateGrowthFromSave();
}

void Simulation::readGrowthRateToContinueFromSave(){
	updateGrowthRateFromSave();
}

void Simulation::readPhysicalPropToContinueFromSave(){
    /**
     * For each element, the format is: \n
     * [Young's modulus] [internal viscosity] [zRemodelling] \n
     * The values are assigned to elemetns via functions: ShapeBase#setYoungsModulus,
     * ShapeBase#setViscosity, ShapeBase#setZRemodellingSoFar. \n
     * Then for each node: \n
     * [external viscosity_x] [external viscosity_y] [external viscosity_z]
     */
	double E;
	double internalViscposity;
	double externalViscosity[3];
	double zRemodelling;
    for(const auto&  itElement : Elements){
		saveFileToDisplayPhysicalProp.read((char*) &E, sizeof E);
		saveFileToDisplayPhysicalProp.read((char*) &internalViscposity, sizeof internalViscposity);
		saveFileToDisplayPhysicalProp.read((char*) &zRemodelling, sizeof zRemodelling);
        itElement->setYoungsModulus(E);
        itElement->setViscosity(internalViscposity);
        itElement->setZRemodellingSoFar(zRemodelling);
	}
    for (const auto& itNode : Nodes){
        for (size_t i=0; i<3; ++i){
			saveFileToDisplayPhysicalProp.read((char*) &externalViscosity[i], sizeof externalViscosity[i]);
            itNode->externalViscosity[i] = externalViscosity[i];
		}
	}
}

void Simulation::initiatePrismFromSaveForUpdate(int k){
    /**
      * Inserts a new prism at order k into elements vector, ShapeBase#NodeIds and reference shape positions
      * ShapeBase#updateShapeFromSave
      */
    std::vector<int> NodeIds(6,0);
    std::unique_ptr<ShapeBase> PrismPnt01 = make_unique<Prism>(NodeIds, Nodes, currElementId);
    PrismPnt01->updateShapeFromSave(saveFileToDisplayMesh);
    std::vector<std::unique_ptr<ShapeBase>>::iterator it = Elements.begin();
	it += k;
    Elements.insert(it,std::move(PrismPnt01));
	nElements = Elements.size();
	currElementId++;
}


void Simulation::updateOneStepFromSave(){
	string currline;
	getline(saveFileToDisplayMesh,currline);
	if(saveFileToDisplayMesh.eof()){
		reachedEndOfSaveFile = true;
		return;
	}
    if (implicitPacking){
        resetForces(true);	// reset packing forces
    }
    else{
        resetForces(false);	// do not reset packing forces
    }
    updateNodeNumberFromSave();
	updateNodePositionsFromSave();
	updateElementStatesFromSave();
    for(const auto&  itElement : Elements){
		//This is updating positions from save.
        itElement->updatePositions(Nodes);
	}
	if (TensionCompressionSaved){
		updateTensionCompressionFromSave();
	}
    if (GrowthSaved){
        updateGrowthFromSave();
    }
    if (GrowthRateSaved){
    	updateGrowthRateFromSave();
    }
	if(physicalPropertiesSaved){
		updatePhysicalPropFromSave();
	}
	if (growthRedistributionSaved){
		updateGrowthRedistributionFromSave();
	}
	if (nodeBindingSaved){
		updateNodeBindingFromSave();
	}
	clearNodeMassLists();
	assignNodeMasses();
	assignConnectedElementsAndWeightsToNodes();
	calculateBoundingBox();
	//skipping the footer:
	string currline2;
	getline(saveFileToDisplayMesh,currline2);
	getline(saveFileToDisplayMesh,currline2);
	while (currline.empty() && !saveFileToDisplayMesh.eof()){
		//skipping empty line
		getline(saveFileToDisplayMesh,currline2);
	}
    timestep = timestep + dataSaveInterval;
	currSimTimeSec += dt*dataSaveInterval;
}

void  Simulation::updateNodeNumberFromSave(){
    size_t n;
	saveFileToDisplayMesh>> n;
    size_t currNodeNumber = Nodes.size();
	if (n>currNodeNumber){
        for (size_t i = 0; i<(n-currNodeNumber); ++i){
            std::array<double,3> pos = {0.0, 0.0, 0.0};
            /** This function only initiates the node where necessary, the positions will be read and updated in
             * Simulaition#updateNodePositionsFromSave.
             */
            std::unique_ptr<Node> tmp_nd = make_unique<Node>(i, 3, pos,-1, -1);
            Nodes.push_back(std::move(tmp_nd));
			nNodes = Nodes.size();
		}
	}
	else{
        //shorten the vector until node numbers match
        for (size_t i = 0; i<(currNodeNumber-n); ++i){
            Nodes.pop_back();
		}
	}
	n = Nodes.size();
	if ( n != nNodes){
		//the node number is change, I updated the node list, now I need to fix system forces:
        reInitiateSystemForces();
	}
}

void  Simulation::updateElementStatesFromSave(){
    size_t n;
	saveFileToDisplayMesh >> n;
    size_t currElementNumber = Elements.size();

	while(currElementNumber>n){
        /** If the elements list is longer than necessary, elements will be removed from list via Simulation#removeElementFromEndOfList.
        */
		removeElementFromEndOfList();
		currElementNumber = Elements.size();
	}
    for (size_t i=0; i<currElementNumber; ++i){
        /**
          * Once existing element list size is shorter or equal to the size on save file,
          * the shape properties will be updated from the save files. The format is
          * detailed in Simulartion#readElementDataToContinueFromSave.
          */
		int shapeType;
		saveFileToDisplayMesh >> shapeType;
		int currShapeType = Elements[i]->getShapeType();
		if (shapeType == currShapeType || shapeType == 2 || shapeType ==4){
			if (shapeType ==4){
				double height;
				saveFileToDisplayMesh >> height;
			}
			Elements[i]->updateShapeFromSave(saveFileToDisplayMesh);
		}
		else{
            /**
              * If the element type does not match the one read from save a new element is initiated and inseted here
              * (Simulation#initiatePrismFromSaveForUpdate), followed by
              * an element removed from the end of the list (Simulation#removeElementFromEndOfList) to preserve element
              * number consistency.
              */
			if (shapeType == 1){
				//the new shape is a prism, I will add it now;
				initiatePrismFromSaveForUpdate(i);
			}
			removeElementFromEndOfList();
			currElementNumber = Elements.size();
		}
	}
	while(n>currElementNumber){
        /**
         * If there are more elements saved than I have in my elemetn list new elements will be initiated and updated.
         * (Simulation#initiatePrismFromSaveForUpdate).
         */
		int shapeType;
		saveFileToDisplayMesh >> shapeType;
		if (shapeType == 1){
			int i = Elements.size()-1;
            /**
             * Only the initiation will be carried out in this function, the
             * positions will be updated in Simulation#updatePositions called via Simulation#updateOneStepFromSave.
             */
			initiatePrismFromSaveForUpdate(i);
			currElementNumber = Elements.size();
		}
	}
}

void Simulation::removeElementFromEndOfList(){
    Elements.pop_back();
	nElements = Elements.size();
}

void Simulation::updateNodePositionsFromSave(){
    /** The formatting is detailed in , Simulation#readNodeDataToContinueFromSave,
     * the main difference is that the node numebr is already read once this function is called.
     */
    for (size_t i=0; i<nNodes; ++i){
		saveFileToDisplayMesh >> Nodes[i]->Position[0];
		saveFileToDisplayMesh >> Nodes[i]->Position[1];
		saveFileToDisplayMesh >> Nodes[i]->Position[2];
		saveFileToDisplayMesh >> Nodes[i]->tissuePlacement;
		saveFileToDisplayMesh >> Nodes[i]->tissueType;
		saveFileToDisplayMesh >> Nodes[i]->atCircumference;
	}
}

bool Simulation::initiateMesh(int MeshType, int Row, int Column, float SideLength, float zHeight){
	if ( MeshType == 1 ){
        std::cerr<<"Error: Too many arguments for a single element system"<<std::endl;
		return false;
	}
	if ( MeshType == 2){
        /**
         * This call will initiate the nodes first, with the selected side length and height via
         * Simulation#initiateNodesByRowAndColumn, then the elemetns will be initiated
         * via Simulation#initiateElementsByRowAndColumn accordingly.
         */
		initiateNodesByRowAndColumn(Row,Column,SideLength,zHeight);
		initiateElementsByRowAndColumn(Row,Column);
	}
	else if ( MeshType == 3 ){
        std::cerr<<"Error: Wrong set of arguments for mesh triangulation"<<std::endl;
		return false;
	}
	else if ( MeshType ==4 ){
        std::cerr<<"Error: Too many arguments for reading the mesh from file"<<std::endl;
		return false;
	}
	else {
        std::cerr<<"Error: Mesh Type not recognised"<<std::endl;
		return false;
	}
	return true;
}

bool Simulation::initiateMesh(int MeshType){
	if ( MeshType == 1 ){
        std::cerr<<"Error: Too many arguments for a single element system"<<std::endl;
		return false;
	}
	if ( MeshType == 2){
        std::cerr<<"Error: Too few arguments for mesh by dimensions"<<std::endl;
		return false;
	}
	else if ( MeshType == 3 ){
        std::cerr<<"Error: Wrong set of arguments  for mesh triangulation"<<std::endl;
		return false;
	}
	else if ( MeshType ==4 ){
        /**
         * The mesh information will be read form the file read into Simulation#inputMeshFileName parameter.
         * The nodes and elements will be initiated through functions Simulation#initiateNodesFromMeshInput and
         * Simulation#initiateElementsFromMeshInput. Then if the tissue type weights are recorded, then they will be read in
         * via Simulation#readInTissueWeights, the recorded value is peripodialness weight, ShapeBase#peripodialGrowthWeight.
         * The input maseh file will be closed at teh end of the function.
         */
		const char* name_inputMeshFile = inputMeshFileName.c_str();;
		inputMeshFile.open(name_inputMeshFile, ifstream::in);
		if (!(inputMeshFile.good() && inputMeshFile.is_open())){
            std::cerr<<"Cannot open the input mesh file: "<<name_inputMeshFile<<std::endl;
			return false;
		}
		initiateNodesFromMeshInput();
		initiateElementsFromMeshInput();
        std::cout<<" nNodes: "<<nNodes<<" nEle: "<<nElements<<std::endl;
		bool areTissueWeightsRecorded = checkIfTissueWeightsRecorded();
		if (areTissueWeightsRecorded){
			readInTissueWeights();
		}
		inputMeshFile.close();	
	}
	else {
        std::cerr<<"Error: Mesh Type not recognised"<<std::endl;
		return false;
	}
	return true;
}

bool Simulation::checkIfTissueWeightsRecorded(){
	bool tissueWeightsRecorded;
	saveFileToDisplayMesh >> tissueWeightsRecorded;
	return tissueWeightsRecorded;
}

void Simulation::readInTissueWeights(){
    /**
     * The tissue type weigths for all elemetns are read from mesh input, via the function
     * ShapeBase#setGrowthWeightsViaTissuePlacement. The recorded value is
     * peripodialness weight, ShapeBase#peripodialGrowthWeight.
     */
	double wPeri;
    for (const auto&  itElement : Elements){
		inputMeshFile >> wPeri;
        itElement->setGrowthWeightsViaTissuePlacement(wPeri); //peripodialness weight is recorded
	}
}

void Simulation::clearCircumferenceDataFromSymmetricityLine(){
    /**
      * The input mesh can be specified to be symmetric in x, y or both x&y axes. This means the actual
      * tissue of interest is a larger, but symmetric tissue. Examples are one half of the fly winf, where the
      * tissue is assummed to be symmetric in its long axis, or simulations of a circular tissue pathc, where
      * only one quadrant of the tissue is simulated. The mesh input file most likely will report the borders of
      * the symmetrticity axes to be at circumference. But in real terms, these nodes are not at the circumference,
      * they are not exposed to the outer world, they are part of a continuing tissue. The circumference flags
      * should thenbe cleaned from these nodes prior to simulation intiation. \n
      *
      * First, the tips are detected.
      * The y-symmetry is at the y=0 line and similarly x-symmetry is at the x=0 line. Therefore, any node recorded to be
      * at circumference, but is part of a symmetric tissue and within the relevant distance limit, they should have their
      * Node#atCircumference flag set to false. The nodes at the tip should be an exception, as the nodes at the tips of the tissue,
      * will still be at the circumference, although they are at the limits of the symetricity line. One specific case is the
      * origin in case of x&y symmetry. The origin will be a tissue tip, but will reside at t he centre of the tissue,
      * therefore should not be symmetric.
      *
      */
	if (symmetricY){
		double xTipPos = -1000, xTipNeg = 1000;
        for (const auto& itNode : Nodes){
            if (itNode->Position[0] > xTipPos ){
                xTipPos = itNode->Position[0];
			}
            if (itNode->Position[0] < xTipNeg ){
                xTipNeg = itNode->Position[0];
			}
		}
		double yLimPos = 0.1;
		double yLimNeg = (-1.0)*yLimPos;
        for (const auto& itNode : Nodes){
            double x = itNode->Position[0];
            double y = itNode->Position[1];
			if ( y < yLimPos){
				if ( y  > yLimNeg){
                    itNode->atSymmetricityBorder = true;
                    fixY(itNode.get(),false); //this is for symmetricity, the fixing has to be hard fixing, not with external viscosity under any condition
					bool atTip = false;
					if ( x < xTipPos+yLimPos && x >xTipPos+yLimNeg){
						atTip = true;
					}
					if (!symmetricX){
						//This if clause is checking the x-tip at the negative end.
						//As described above, the node should still be a circumference node,
						//if it is at the line of symmetry, but at the tip as well.
						//BUT, if there is also xSymmetry, this "negative end tip" will
						//be the tip at the x symmetry line. (If I am modelling a quarter of a circle,
						//with x and y symmetry, this will point be the centre of the circle.)
						//Under these conditions, it is not actually at the tip. It should be
						//removed from circumference node list.
						//Therefore, I am carrying the atTip? check for the negative end only if
						//there is no x-symmetry.
						if ( x > xTipNeg+yLimNeg && x < xTipNeg+yLimPos){
							atTip = true;
						}
					}
					if (!atTip){
						//removing node from list:
                        itNode->atCircumference = false;
					}
				}
			}
		}
	}
	if (symmetricX){
		double yTipPos = -1000, yTipNeg = 1000;
        for (const auto& itNode : Nodes){
            if (itNode->Position[0] > yTipPos ){
                yTipPos = itNode->Position[1];
			}
            if (itNode->Position[0] < yTipNeg ){
                yTipNeg = itNode->Position[1];
			}
		}
		double xLimPos = 0.1;
		double xLimNeg = (-1.0)*xLimPos;
        for (const auto& itNode : Nodes){
            double x = itNode->Position[0];
            double y = itNode->Position[1];
			if ( x < xLimPos){
				if ( x  > xLimNeg){
                    itNode->atSymmetricityBorder = true;
                    fixX(itNode.get(),false); //this is for symmetricity, the fixing has to be hard fixing, not with external viscosity under any condition
					//the node is indeed at the border, BUT, I will remove it only if it is not at the tip:
					bool atTip = false;
					if ( y < yTipPos+xLimPos && y >yTipPos+xLimNeg){
						atTip = true;
					}
					if ( y > yTipNeg+xLimNeg && y < yTipNeg+xLimPos){
						atTip = true;
					}
					if (!atTip){
                        itNode->atCircumference = false;
					}
				}
			}
		}
	}
}

void Simulation::getAverageSideLength(double& periAverageSideLength, double& colAverageSideLength){
	double dsumPeri =0.0, dsumCol = 0.0;
	int colCounter =0, periCounter=0;
    double packingDetectionThresholdGridCounter[10][5];
    for (size_t i=0; i<10;++i){
        for (size_t j=0;j<5;++j){
			packingDetectionThresholdGrid[i][j] = 0.0;
			packingDetectionThresholdGridCounter[i][j] = 0.0;
		}
	}
    for (auto &itElement : Elements){
        if (itElement->tissueType==0){ //element belongs to columnar layer
            double currValue = itElement->getApicalSideLengthAverage();
            if(currValue>0){
                /**
                  * There may be alements with fully collapsed areas, where currValue will be zero.
                  * These are not counted in averaging the side length, as the zero length is not relevant for packing
                  * distance calculations.
                  */
                dsumCol += currValue;
                colCounter++;
            }
            if(itElement->tissuePlacement == 0 || itElement->spansWholeTissue){
                currValue = itElement->getBasalSideLengthAverage();
                if(currValue>0){
                    dsumCol += currValue;
                    colCounter++;
                }
            }
            /** The local side length data is recorded in a
             * 10x5 grid, Simulation#packingDetectionThresholdGrid.
             */
            std::array<double,2> reletivePos = itElement->getRelativePosInBoundingBox();
            int relX = floor(reletivePos[0]);
            int relY = floor(reletivePos[1]/2.0);
            if (relX < 0) {relX = 0;}
            if (relX > 9) {relX = 9;}
            if (relY < 0) {relY = 0;}
            if (relY > 4) {relY = 4;}
            packingDetectionThresholdGrid[relX][relY]+=currValue;
            packingDetectionThresholdGridCounter[relX][relY]++;

        }
        else{
            dsumPeri += itElement->getApicalSideLengthAverage();
            periCounter++;
        }
	}
	colAverageSideLength = dsumCol / (double)colCounter;
    for (size_t i=0; i<10;++i){
        for (size_t j=0;j<5;++j){
			if(packingDetectionThresholdGridCounter[i][j]>0){
				packingDetectionThresholdGrid[i][j] /= packingDetectionThresholdGridCounter[i][j];
			}
		}
	}
	if (periCounter>0){
		periAverageSideLength = dsumPeri / (double)periCounter;
	}
}

bool  Simulation::isColumnarLayer3D(){
	bool ColumnarLayer3D = false;
    for (auto &itElement : Elements){
        if (itElement->ShapeDim == 3){
			ColumnarLayer3D  = true;
			break;
		}
	}
	return ColumnarLayer3D;
}

bool Simulation::calculateTissueHeight(){
    /**
     * First check will be if the tissue is made up of 3D elements. If the tissue is 3D, then the height
     * will be calculated accordingly. If the tissue is made up of 2D elements, then the
     * height will simply be the assigned heigt of the elements, ReferenceShapeBase#height. \n
     */
    bool ColumnarLayer3D = isColumnarLayer3D();
	TissueHeight = 0;
	if (ColumnarLayer3D){
        /**
         * The calculation for tissue height will find the first basal node on Simulation#Nodes.
         */
		//Find the first basal node on the Nodes List
		//Find the first element that has this node on the elements list
		//Move apically through the elements until you reach the apical surface - an apical node
		//Find a basal node:
        int currNodeId=0;
		bool foundNode = false;
        for (auto& itNode : Nodes){
            if(itNode->tissueType == 0 && itNode->tissuePlacement == 0){ //Columnar node is basal
				foundNode = true;
                currNodeId = itNode->Id;
				break;
			}
		}
		if (!foundNode){
			return false;
		}
		//Find an element using the basal node, and move on the elements apically, until you reach the apical surface:
		bool foundElement = true;
		TissueHeightDiscretisationLayers = 0;
        double currentH =0;
        /**
         * Then the first element on Simulation#Elements that utilises the found node as a basal node will be identified.
         * Once the element is identified, an apical node of this element will be selected as the current search node,
         * the Simulation#TissueHeight and Simulation#TissueHeightDiscretisationLayers will be incremented up, and
         * the search will continue to find the first element to utilise the new search node as its apicel. The loopo
         * will continue until and apical node is readhed (Node#tissuePlacement == 1).
         */
        while(Nodes[currNodeId]->tissuePlacement != 1 && foundElement){ //while the node is not apical, and I could find the next element
            foundElement = false;
            for (auto &itElement : Elements){
                bool IsBasalOwner = itElement->IsThisNodeMyBasal(currNodeId);
				if (IsBasalOwner){
					foundElement = true;
                    currentH = itElement->getElementHeight();
                    currNodeId = itElement->getCorrecpondingApical(currNodeId); //have the next node
					break;
				}
			}
			TissueHeight += currentH;
			TissueHeightDiscretisationLayers++;
		}
		if (!foundElement){
			return false;
		}
	}
	else{
		TissueHeight = Elements[0]->ReferenceShape->height;
		if (TissueHeight == 0){
            std::cerr<<"The coulmanr layer is 2D, but the tissue height of the elements is not assigned properly, cannot obtain TissueHeight"<<std::endl;
			return false;
		}
	}
	return true;
}

void Simulation::assignInitialZPositions(){
    for (auto &itElement : Elements){
        itElement->setInitialZPosition(boundingBox[0][2], TissueHeight);  //minimum z of the tissue and the tissue height given as input
	}
}

void Simulation::calculateShapeFunctionDerivatives(){
    for (auto &itElement : Elements){
        itElement->calculateElementShapeFunctionDerivatives();
    }
}

void Simulation::fixAllD(Node* currNode, bool fixWithViscosity){
    for (size_t j =0 ; j<currNode->nDim; ++j){
		if(fixWithViscosity){
			currNode->externalViscosity[j] = fixingExternalViscosity[j];
			currNode->externalViscositySetInFixing[j] = true;
		}
		else{
			currNode->FixedPos[j]=true;
		}
	}
}

void Simulation::fixAllD(int i, bool fixWithViscosity){
    for (size_t j =0 ; j<Nodes[i]->nDim; ++j){
		if(fixWithViscosity){
			Nodes[i]->externalViscosity[j] = fixingExternalViscosity[j];
			Nodes[i]->externalViscositySetInFixing[j] = true;
		}
		else{
			Nodes[i]->FixedPos[j]=true;
		}
	}
}

void Simulation::fixX(Node* currNode, bool fixWithViscosity){
	if(currNode->nDim>0){
		if(fixWithViscosity){
			currNode->externalViscosity[0] = fixingExternalViscosity[0];
			currNode->externalViscositySetInFixing[0] = true;
		}
		else{
			currNode->FixedPos[0]=true;
		}
	}
	else{
        std::cerr<<"ERROR: Node : "<<currNode->Id<<" does not have x-dimension"<<std::endl;
	}
}
void Simulation::fixX(int i, bool fixWithViscosity){
	if(Nodes[i]->nDim>0){
		if(fixWithViscosity){
			Nodes[i]->externalViscosity[0] = fixingExternalViscosity[0];
			Nodes[i]->externalViscositySetInFixing[0] = true;
		}
		else{
			Nodes[i]->FixedPos[0]=true;
		}
	}
	else{
        std::cerr<<"ERROR: Node : "<<Nodes[i]->Id<<" does not have x-dimension"<<std::endl;
	}

}

void Simulation::fixY(Node* currNode, bool fixWithViscosity){
	if(currNode->nDim>1){
		if(fixWithViscosity){
			currNode->externalViscosity[1] = fixingExternalViscosity[1];
			currNode->externalViscositySetInFixing[1] = true;
		}
		else{
			currNode->FixedPos[1]=true;
		}
	}
	else{
        std::cerr<<"ERROR: Node : "<<currNode->Id<<" does not have y-dimension"<<std::endl;
	}
}

void Simulation::fixY(int i, bool fixWithViscosity){
	if(Nodes[i]->nDim>1){
		if(fixWithViscosity){
			Nodes[i]->externalViscosity[1] = fixingExternalViscosity[1];
			Nodes[i]->externalViscositySetInFixing[1] = true;
		}
		else{
			Nodes[i]->FixedPos[1]=true;
		}
	}
	else{
        std::cerr<<"ERROR: Node : "<<Nodes[i]->Id<<" does not have y-dimension"<<std::endl;
	}
}

void Simulation::fixZ(Node* currNode, bool fixWithViscosity){
	if(currNode->nDim>2){
		if(fixWithViscosity){
			currNode->externalViscosity[2] = fixingExternalViscosity[2];
			currNode->externalViscositySetInFixing[2] = true;
		}
		else{
			currNode->FixedPos[2]=true;
		}
	}
	else{
        std::cerr<<"ERROR: Node : "<<currNode->Id<<" does not have z-dimension"<<std::endl;
	}
}

void Simulation::fixZ(int i, bool fixWithViscosity){
	if(Nodes[i]->nDim>2){
		if(fixWithViscosity){
			Nodes[i]->externalViscosity[2] = fixingExternalViscosity[2];
			Nodes[i]->externalViscositySetInFixing[2] = true;
		}
		else{
			Nodes[i]->FixedPos[2]=true;
		}
	}
	else{
        std::cerr<<"ERROR: Node : "<<Nodes[i]->Id<<" does not have z-dimension"<<std::endl;
	}
}

void Simulation::zeroForcesOnNode(size_t i){
    /**
     * The forces corresponding to node indexed at input i on Simulation#Nodes will be set to zero.
     */
	double ForceBalance[3];
	ForceBalance[0] = SystemForces[i][0];
	ForceBalance[1] = SystemForces[i][1];
	ForceBalance[2] = SystemForces[i][2];
    for (size_t i=0;i<nNodes;++i){
		SystemForces[i][0]-=ForceBalance[0];
		SystemForces[i][1]-=ForceBalance[1];
		SystemForces[i][2]-=ForceBalance[2];
	}
}

void Simulation::initiateSystemForces(){
    for (size_t j=0;j<nNodes;++j){
		//3 dimensions
        SystemForces.push_back(std::array<double,3>{0.0});
        PackingForces.push_back(std::array<double,3>{0.0});
	}
}


void Simulation::initiateNodesByRowAndColumn(size_t Row, size_t Column, float SideLength, float zHeight){
	dorsalTipIndex = Row;
	ventralTipIndex = 0;
    /**
     * This function will place nodes in space to form the tirangles to genreate a hexagonal mesh of input rows and
     * columns at its widest point. First in a simple geometric calculaion, the height of of the equilateral triangle
     * with the input side length is calculated. Then the nodes will be placed at a h difference in y for each column.
     * The nodes position in x will be shifted by half the side length between columns.
     */
	double sqrt3 = 1.7321;
    float h = sqrt3/2*SideLength; // height of triangle
	vector <double> xPos, yPos;
	vector <int> NodesToFix;
	int toprowcounter = 0;	//number of nodes at the terminating end, the top row. I will need this to add the sideway prisms;
    for (size_t ColCount = 0; ColCount < Column+1; ++ColCount){
		double CurrY = ColCount*h;
		int CurrRowNum = Row + 1 - ColCount;
		double RowOffset = 0.5*SideLength*ColCount;
		for ( int RowCount = 1; RowCount<CurrRowNum+1; ++RowCount){
			double CurrX = RowOffset + RowCount * SideLength;
			xPos.push_back(CurrX);
			yPos.push_back(CurrY);
			if (RowCount ==1 || RowCount == CurrRowNum || ColCount == Column){
				NodesToFix.push_back(xPos.size()-1);
				if(ColCount == Column){
					toprowcounter++;
				}
			}
		}
		if (ColCount>0){
			CurrY = (-1.0)*CurrY;
			for ( int RowCount = 1; RowCount<CurrRowNum+1; ++RowCount){
				double CurrX = RowOffset + RowCount * SideLength;
				xPos.push_back(CurrX);
				yPos.push_back(CurrY);
				if (RowCount ==1 || RowCount == CurrRowNum || ColCount == Column){
					NodesToFix.push_back(xPos.size()-1);
				}
			}
		}
	}
    size_t n =  xPos.size();
    std::array<double,3> pos = {0.0, 0.0, 0.0};
    for (size_t i =0; i< n; ++i){
		pos[0] = xPos[i];
		pos[1] = yPos[i];
		pos[2] = 0.0;
        std::unique_ptr<Node> tmp_nd = make_unique<Node>(i, 3, pos,0,0);
        Nodes.push_back(std::move(tmp_nd));
        nNodes = Nodes.size();
	}
	//Adding the apical level, all will form columnar elements:
    for (size_t i =0; i< n; ++i){
		pos[0] = xPos[i];
		pos[1] = yPos[i];
		pos[2] = zHeight;
        std::unique_ptr<Node> tmp_nd = make_unique<Node>(n+i, 3, pos,1,0);
        Nodes.push_back(std::move(tmp_nd));
        nNodes = Nodes.size();
	}
}

void Simulation::setLinkerCircumference(){
	double maxX = -100, ZofTip = -100;
    for (auto& itNode : Nodes){
        if (itNode->tissueType == 2){
			//The node is linker, is the x position higher than already recorded?
            if (itNode->Position[0]> maxX){
                maxX  = itNode->Position[0];
                ZofTip = itNode->Position[2];
			}
		}
	}
	double thres = 0.2;
    //cout<<"Setting circumference, maxX: "<<maxX<<" Z of tip: "<<ZofTip<<std::endl;
	//Now I have the maxX, and the corresponding z height.
	//Declare all linkers at the basal/or apical surface are the circumferential nodes:
    for (auto & itNode : Nodes){
        if (itNode->tissueType == 2){
			//The node is linker, if it is apical or basal, it can be circumferntial:
            if ( itNode->tissuePlacement == 0 || itNode->tissuePlacement == 1 ){
				//The node is linker, if it is in the range of the height of the tip then it is circumferential
                if ( itNode->Position[2] < ZofTip+thres && itNode->Position[2] > ZofTip-thres ){
                    itNode->atCircumference = true;
				}
			}

		}
	}
}

void Simulation::checkForNodeFixing(){
    /**
     * The Simulation#CircumferentialNodeFix array gives fixing options as booleans for [5options][in 3D] : \n
     * [option = 0] : apical circumference fixing options [fix_x][ix_y][fix_z]. \n
     * [option = 1] : basal  circumference fixing options [fix_x][ix_y][fix_z]. \n
     * [option = 2] : linker apical circumference fixing options [fix_x][ix_y][fix_z]. \n
     * [option = 3] : linker basal circumference fixing options [fix_x][ix_y][fix_z]. \n
     * [option = 4] : ALL circumference fixing options [fix_x][ix_y][fix_z]. \n
     * \n
     * First check point will be for checking if there is circumferential fixing for any part of the tissue.
     */
	//Are there any circumferential node fixing options enabled:
	bool thereIsCircumFix = false;
	if (!thereIsCircumFix){
        for (size_t i=0;i<5; ++i){
            for (size_t j=0;j<3; ++j){
				if (CircumferentialNodeFix[i][j] == true){
					thereIsCircumFix = true;
					break;
				}
			}
			if (thereIsCircumFix){
				break;
			}
		}
	}
    /** If there is any fixing option active in any dimention, then the detailed check will continue.
     */
	if (thereIsCircumFix){
        for (auto& itNode : Nodes){
            if ( itNode->atCircumference){
                /** The node is protected from fixing if it is at peripodial tissue.
                 */
				bool doNotFix = false;
                if (itNode->tissueType == 1){ //the node is peripodial, no not fix this node
					doNotFix = true;
				}
				if (!doNotFix){ 
					//Node is at the circumference, now checking for all possibilities:
					// if i == 0 , I am checking for apical circumference
					// if i == 1 , I am checking for basal  circumference
					// if i == 2 , I am checking for the linker apical circumference
					// if i == 3 , I am checking for the linker basal circumference
					// if i == 4 , I am checking for all    circumference
                    /** Then checking all options against Node#tissuePlacement and Node#tissueType, the node fixing is assigned accordingly,
                     * through functions Node#fixX, Node#fixY, Node#fixZ.
                     */
                    for (size_t i=0;i<5; ++i){
                        if ( (i == 0 && itNode->tissuePlacement == 1 ) ||  //tissuePlacement == 1 is apical
                             (i == 1 && itNode->tissuePlacement == 0 ) ||  //tissuePlacement == 0 is basal
                             (i == 2 && itNode->tissueType == 2 && itNode->tissuePlacement == 1 ) ||  //tissuePlacement == 1 is apical
                             (i == 3 && itNode->tissueType == 2 && itNode->tissuePlacement == 0 ) ||  //tissuePlacement == 0 is basal
							 (i == 4)){										  //tissuePlacement is irrelevant, fixing all
							//The node is at circumference; if
							if (CircumferentialNodeFix[i][0]){
                                fixX(itNode.get(),CircumferentialNodeFixWithHighExternalViscosity[i]);
							}
							if (CircumferentialNodeFix[i][1]){
                                fixY(itNode.get(),CircumferentialNodeFixWithHighExternalViscosity[i]);

							}
							if (CircumferentialNodeFix[i][2]){
                                fixZ(itNode.get(),CircumferentialNodeFixWithHighExternalViscosity[i]);
							}
						}
					}
				}
			}
		}
	}
    /** Once checking the circumference fixing is complete, then fixing whole surfaces are checked,
     * (Simulation#BasalNodeFix and Simulation#ApicalNodeFix) for instance when the tissue is stuck
     * on glass for a pipette aspiration experiment.
     */
	if (BasalNodeFix[0] || BasalNodeFix[1] || BasalNodeFix[2] ){
        for (auto& itNode : Nodes){
            if ( itNode->tissuePlacement == 0){
				if (BasalNodeFix[0]){
                    fixX(itNode.get(), BasalNodeFixWithExternalViscosity);
				}
				if (BasalNodeFix[1]){
                    fixY(itNode.get(), BasalNodeFixWithExternalViscosity);
				}
				if (BasalNodeFix[2]){
                    fixZ(itNode.get(), BasalNodeFixWithExternalViscosity);
				}
			}
		}
	}
	if (ApicalNodeFix[0] || ApicalNodeFix[1]  || ApicalNodeFix[2]){
        for (auto & itNode : Nodes){
            if ( itNode->tissuePlacement == 1){
				if (ApicalNodeFix[0]){
                    fixX(itNode.get(), ApicalNodeFixWithExternalViscosity);
				}
				if (ApicalNodeFix[1]){
                    fixY(itNode.get(), ApicalNodeFixWithExternalViscosity);
				}
				if (ApicalNodeFix[2]){
                    fixZ(itNode.get(), ApicalNodeFixWithExternalViscosity);
				}
			}
		}
	}
    /** The tissue can be fixed basally through tissue domain specification, specifically for th notum,
     * which is set by tissue relative position and the parameters
     * Siumulation#notumFixingRange and Simulation#NotumNodeFix.
     */
	if (NotumNodeFix[0] || NotumNodeFix[1]  || NotumNodeFix[2]){
        for (auto& itNode : Nodes){
            if ( itNode->tissuePlacement == 0){
				double boundingBoxXMin = boundingBox[0][0];
				double boundingBoxLength = boundingBoxSize[0];
                double relativeXPos = (itNode->Position[0] - boundingBoxXMin)/boundingBoxLength;
				if ( relativeXPos >= notumFixingRange[0] && relativeXPos <= notumFixingRange[1]){
					if (NotumNodeFix[0]){
                        fixX(itNode.get(), NotumNodeFixWithExternalViscosity);
					}
					if (NotumNodeFix[1]){
                        fixY(itNode.get(), NotumNodeFixWithExternalViscosity);
					}
					if (NotumNodeFix[2]){
                        fixZ(itNode.get(), NotumNodeFixWithExternalViscosity);
					}
				}
			}
		}
	}
}

void Simulation::induceClones(){
    /**
      * Growth of a subset of elemetns can be manipulated as mutant clones. These are set with
      * positional information from model input file, and the mutation growth influence type, which can be
      * an absolute growth value input, or scaling the existing growth rate.
      */
    for (size_t i=0;i<numberOfClones; ++i){
		double inMicronsX = cloneInformationX[i]*boundingBoxSize[0] + boundingBox[0][0];
		double inMicronsY = cloneInformationY[i]*boundingBoxSize[1] + boundingBox[0][1];
		double inMicronsRadius = cloneInformationR[i];
        //cout<<" clone: "<<i<<" position: "<<inMicronsX<<" "<<inMicronsY<<" radius: "<<inMicronsRadius<<std::endl;
		double r2 = inMicronsRadius*inMicronsRadius;
		double growthRateORFold = cloneInformationGrowth[i];
        for(const auto& itElement : Elements){
            if (!itElement->isECMMimicing){
                std::array<double,3> c = itElement->getCentre();
				double dx = inMicronsX - c[0];
				double dy = inMicronsY - c[1];
				double d2 = dx*dx + dy*dy;
				if (d2 < r2){
                    //cout<<"mutating element: "<<(*itElement)->Id<<std::endl;
					if (cloneInformationUsingAbsolueGrowth[i]){
                        itElement->mutateElement(0,growthRateORFold); //the mutation is absolute, using an absolute value
					}
					else{
                        itElement->mutateElement(growthRateORFold,0); //the mutation is not absolute, using relative values
					}
				}
			}
		}
	}
}

void Simulation::initiateElementsByRowAndColumn(size_t Row, size_t Column){
	int xinit1 = 0;
	int xinit2 = xinit1+Row+1;
	int xinit3 = 0;
	int xinit4 = xinit1+2*(Row+1)-1;
    size_t n = Nodes.size() /2.0;
    //initialising the tissue elements:
    for (size_t ColCount = 0; ColCount < Column; ++ColCount){
        size_t CurrRowNum = Row + 1 - ColCount;
        for (size_t RowCount = 0; RowCount<CurrRowNum-1; ++RowCount ){
            std::vector<int> NodeIds(6,0);
			NodeIds[0] = xinit1+RowCount;
			NodeIds[1] = xinit1+RowCount+1;
			NodeIds[2] = xinit2+RowCount;
			NodeIds[3] = NodeIds[0] + n;
			NodeIds[4] = NodeIds[1] + n;
			NodeIds[5] = NodeIds[2] + n;
            std::unique_ptr<ShapeBase> PrismPnt01 = make_unique<Prism>(NodeIds, Nodes, currElementId);
            Elements.push_back(std::move(PrismPnt01));
			nElements = Elements.size();
			currElementId++;

			NodeIds[0] = xinit3+RowCount;
			NodeIds[1] = xinit4+RowCount;
			NodeIds[2] = xinit3+RowCount+1;
			NodeIds[3] = NodeIds[0] + n;
			NodeIds[4] = NodeIds[1] + n;
			NodeIds[5] = NodeIds[2] + n;
            std::unique_ptr<ShapeBase> PrismPnt02 = make_unique<Prism>(NodeIds, Nodes, currElementId);
            Elements.push_back(std::move(PrismPnt02));
			nElements = Elements.size();
			currElementId++;
		}
        for (size_t RowCount = 0; RowCount<CurrRowNum-2; ++RowCount ){
            std::vector<int> NodeIds(6,0);
			NodeIds[0] = xinit2+RowCount;
			NodeIds[1] = xinit1+RowCount+1;
			NodeIds[2] = xinit2+RowCount+1;
			NodeIds[3] = NodeIds[0] + n;
			NodeIds[4] = NodeIds[1] + n;
			NodeIds[5] = NodeIds[2] + n;
            std::unique_ptr<ShapeBase> PrismPnt01 = make_unique<Prism>(NodeIds, Nodes, currElementId);
            Elements.push_back(std::move(PrismPnt01));
            nElements = Elements.size();
			currElementId++;

			NodeIds[0] = xinit4+RowCount;
			NodeIds[1] = xinit4+RowCount+1;
			NodeIds[2] = xinit3+RowCount+1;
			NodeIds[3] = NodeIds[0] + n;
			NodeIds[4] = NodeIds[1] + n;
			NodeIds[5] = NodeIds[2] + n;
            std::unique_ptr<ShapeBase> PrismPnt02 = make_unique<Prism>(NodeIds, Nodes, currElementId);
            Elements.push_back(std::move(PrismPnt02));
			nElements = Elements.size();
			currElementId++;

		}
		xinit1 = xinit2;
		xinit2 = xinit4 + CurrRowNum-1;
		xinit3 = xinit4;
		xinit4 = xinit2 + CurrRowNum-2;
	}
}

void Simulation::calculateSystemCentre(){
    for (size_t i = 0; i< nNodes; ++i){
        for (size_t j =0; j<Nodes[i]->nDim; ++j){
			SystemCentre[j] += Nodes[i]->Position[j];
		}
	}
	SystemCentre[0]= SystemCentre[0]/nNodes;
	SystemCentre[1]= SystemCentre[1]/nNodes;
	SystemCentre[2]= SystemCentre[2]/nNodes;
}

void Simulation::addSoftPeriphery(std::vector<double>& fractions){
    /**
     * The elemetns within the range of soft periphery are identified first with their
     * distance from the nearest circumference node. This calculation will be noisy as the
     * side length of elements and node positions sill influence the cut-off, but this is
     * an acceptable noise given the noise in biological systems. \n
     * The array of booleans states:
     * [applied ot apical side] [applied to basal side] [applied to columnar tissue] [applied to peripodial tissue]
     */
	double t2 = softDepth*softDepth;
	double midlineFraction = 0.0;
	double currSoftnessFraction = softnessFraction;
	double tissueTypeMultiplier = 1;
	if (softPeripheryBooleans[0] && softPeripheryBooleans[1]){
		midlineFraction = softnessFraction;
	}
	else if (softPeripheryBooleans[0] || softPeripheryBooleans[1]){
		midlineFraction = 1.0 - (1.0-softnessFraction)*0.5;
	}

    for(const auto& itNode : Nodes){
        if ( itNode->atCircumference ){
			//this node is at circumference, I will calculate the distance of all elements to this node
			//if an element is close enough, update the fraction matrix set above:
            for(const auto& itElement : Elements){
				bool applyToThisElement = true;
                if ( itElement->tissuePlacement == 1 && !softPeripheryBooleans[0]){
					//element is apical, the softness is NOT applied to apical
					applyToThisElement = false;
				}
                if (itElement->tissuePlacement == 0 && !softPeripheryBooleans[1]){
					//element is basal, the softness is NOT applied to apical
					applyToThisElement = false;
				}
                if (itElement->tissueType == 0 && !softPeripheryBooleans[2]){
					//element is columnar, the softness is NOT applied to columnar
					applyToThisElement = false;
				}
                if (itElement->tissueType == 1 && !softPeripheryBooleans[3]){
					//element is peripodial, the softness is NOT applied to peripodial
					applyToThisElement = false;
				}
				if ( applyToThisElement ){
					//scaling the softness fraction if necessary:
					currSoftnessFraction  = softnessFraction;
                    if (itElement->tissuePlacement == 2 ){ //element is on the midline
						currSoftnessFraction= midlineFraction;
					}
					//scaling the tissue type multiplier if necessary (important for linker elements:
                    if (itElement->tissueType == 2 ){ // element is on the linker zone
						if (softPeripheryBooleans[2] && softPeripheryBooleans[3]){
							//softness is applied to both peripodial and columnar zones, the linkers should not be scaling anything
							tissueTypeMultiplier = 1.0;
						}else if (softPeripheryBooleans[2] ){
							//softness only applied to columnar layer:
                            tissueTypeMultiplier = itElement->getColumnarness();
						}
						else if (softPeripheryBooleans[3] ){
							//softness only applied to peripodial layer:
                            tissueTypeMultiplier = itElement->getPeripodialness();
						}
					}
					else{ //element is not linker
						tissueTypeMultiplier = 1;
					}
					currSoftnessFraction  = 1.0 - (1.0-currSoftnessFraction)*tissueTypeMultiplier;
                    std::array<double,3> c = itElement->getCentre();
                    double dx = c[0]-itNode->Position[0];
                    double dy = c[1]-itNode->Position[1];
					//double dz = c[2]-(*itNode1)->Position[2];
					double d = dx*dx + dy*dy;
					if (d<t2){
						d = pow(d,0.5);
						double f = currSoftnessFraction + (1.0 - currSoftnessFraction)*d/softDepth;
                        if ((softnessFraction< 1.0 && f < fractions[itElement->Id]) ||(softnessFraction> 1.0 && f > fractions[itElement->Id])){
                            fractions[itElement->Id] = f;
						}
					}
				}
			}
		}
	}
}

void Simulation::assignPhysicalParameters(){
    /** This funcition starts by generating the array of stiffness fractions to allow for perturbation implementation at the same time with the
     * setting up of all elastic properties.
     */
    std::vector<double> fractions;
    for (size_t i=0; i<nElements; ++i){
        fractions.push_back(1.0);
    }
	if(softPeriphery){
        /** If there is a soft periphery implemented in the system, this is scaled here by modigying the fractions array via function
         * Simulation#addSoftPeriphery.
         */
		addSoftPeriphery(fractions);
	}
    for(const auto& itElement : Elements){
        /** If there is noise on the physical proerties, this is randomly assigned first.
         */
		double r = (rand() % 200) / 100.0;	//random number between 0.00 and 2.00
		r = r - 1.0; 						//random number between -1.00 and 1.00
		float noise1 = r*noiseOnPysProp[0];	//percent noise on current element
		r = (rand() % 200) / 100.0;
		r = r - 1.0;
		float noise2 = r*noiseOnPysProp[1];
        /** The input Young's moduli are scaled for the stiffness scaling of the current element.
         * The noise is added on at the same step, and the Poission's ratio is set to zero for ECM mimicking elements.
         * Then the elastic properties of the elemetn are set and necessary tensors updated via ShapeBase#setElasticProperties and
         * ShapeBase#setViscosity functions.
         */
        if (itElement->tissueType == 0){ //Element is on the columnar layer
            double currEApical 	= fractions[itElement->Id] * EApical*(1 + noise1/100.0);
            double currEBasal	= fractions[itElement->Id] * EBasal*(1 + noise1/100.0);
            double currEMid		= fractions[itElement->Id] * EMid*(1 + noise1/100.0);
            double currEECM		= fractions[itElement->Id] * EColumnarECM*(1 + noise1/100.0);
            double currPoisson = poisson*(1 + noise2/100);
            if(itElement->isECMMimicing){
				//the Poisson ratio is zero so that the ECM layer will not thin!
				currPoisson = 0;
			}
            itElement->setElasticProperties(currEApical,currEBasal,currEMid,currEECM,currPoisson);
            itElement->setViscosity(discProperApicalViscosity,discProperBasalViscosity,discProperMidlineViscosity);
		}
        else if (itElement->tissueType == 1){ //Element is on the peripodial membrane
            double currEApical = fractions[itElement->Id] * PeripodialElasticity*(1 + noise1/100.0);
			double currEBasal = currEApical;
			double currEMid = currEApical;
            double currEECM = fractions[itElement->Id] * EPeripodialECM*(1 + noise1/100.0);
			double currPoisson = poisson*(1 + noise2/100);
            if(itElement->isECMMimicing){
				//the Poisson ratio is zero so that the ECM layer will not thin!
				currPoisson = 0;
			}
            itElement->setElasticProperties(currEApical,currEBasal,currEMid,currEECM, currPoisson);
            itElement->setViscosity(peripodialApicalViscosity,peripodialBasalViscosity,peripodialMidlineViscosity);
		}
        else if (itElement->tissueType == 2 ){ //Element is on the linker Zone,
			//The elastic properties of the linker zone ECM are always based on the
			//peripodialness factor.
            double currPeripodialECME = fractions[itElement->Id] * EPeripodialECM*(1 + noise1/100.0);
            double currColumnarECME = fractions[itElement->Id] * EColumnarECM*(1 + noise1/100.0);
            double periWeight 		= itElement->getPeripodialness();
            double colWeight = itElement->getColumnarness();
			double currEECM		= colWeight * currColumnarECME + periWeight * currPeripodialECME;
			if (BaseLinkerZoneParametersOnPeripodialness){
				//I will weight the values:
                double currPeripodialE 	= fractions[itElement->Id] * PeripodialElasticity * (1 + noise1/100.0);
                double currEApical 		= fractions[itElement->Id] * EApical * (1 + noise1/100.0);
                double currEBasal 		= fractions[itElement->Id] * EBasal * (1 + noise1/100.0);
                double currEMid 		= fractions[itElement->Id] * EMid * (1 + noise1/100.0);
				double currPoisson      = poisson*(1 + noise2/100);
				currEApical = colWeight * currEApical + periWeight * currPeripodialE;
				currEBasal  = colWeight * currEBasal  + periWeight * currPeripodialE;
				currEMid    = colWeight * currEMid    + periWeight * currPeripodialE;
                itElement->setElasticProperties(currEApical,currEBasal,currEMid,currEECM, currPoisson);
				double currViscApical  = colWeight * discProperApicalViscosity  + periWeight * peripodialApicalViscosity;
				double currViscBasal   = colWeight * discProperBasalViscosity   + periWeight * peripodialBasalViscosity;
                itElement->setViscosity(currViscApical,currViscBasal); //There is no midline in the linker zone.
			}
			else{
				//I have the inputs provided:
                double currEApical 	= fractions[itElement->Id] * LinkerZoneApicalElasticity*(1 + noise1/100.0);
                double currEBasal	= fractions[itElement->Id] * LinkerZoneBasalYoungsModulus*(1 + noise1/100.0);
                double currEMid		= fractions[itElement->Id] * 0.5 * (LinkerZoneApicalElasticity+LinkerZoneBasalYoungsModulus)*(1 + noise1/100.0);
				double currPoisson  = poisson*(1 + noise2/100);


                if(itElement->isECMMimicing){
					//the Poisson ratio is zero so that the ECM layer will not thin!
					currPoisson = 0;
				}
                itElement->setElasticProperties(currEApical,currEBasal,currEMid,currEECM, currPoisson);
                itElement->setViscosity(linkerZoneApicalViscosity,linkerZoneBasalViscosity,linkerZoneMidlineViscosity);
			}
		}
	}
    /** Once all elemental properties are set, the nodal viscosities are asssigned through their owenr elements.
     */
    for(const auto& itNode : Nodes){
		double r = (rand() % 200) / 100.0;
		r = r - 1.0;
		float noise3 = r*noiseOnPysProp[1];
		noise3 = (1 + noise3/100.0);
        if (itNode->tissueType == 2 ){ //linker, the viscosity can be averaged, or based on indivudual set of parameters:
			if (BaseLinkerZoneParametersOnPeripodialness){
				//take average of the two:
				double currExternalViscApical = 0.5*(externalViscosityDPApical + externalViscosityPMApical);
				double currExternalViscBasal  = 0.5*(externalViscosityDPBasal  + externalViscosityPMBasal);
                itNode->setExternalViscosity(currExternalViscApical,currExternalViscBasal, extendExternalViscosityToInnerTissue);
			}
			else{
                itNode->setExternalViscosity(externalViscosityLZApical*noise3, externalViscosityLZBasal*noise3, extendExternalViscosityToInnerTissue);
			}
        }
        else if (itNode->tissueType == 0){ //disc proper
            itNode->setExternalViscosity(externalViscosityDPApical*noise3, externalViscosityDPBasal*noise3, extendExternalViscosityToInnerTissue);
        }
        else if (itNode->tissueType == 1){
            itNode->setExternalViscosity(externalViscosityPMApical*noise3, externalViscosityPMBasal*noise3, extendExternalViscosityToInnerTissue);
		}
	}
}

void Simulation::checkForZeroExternalViscosity(){
    /**
     * Among all nodes, if there is at least one node with external viscosity in a given dimension, then the
     * viscosity is not zero in that dimension.
     */
    for (const auto& itNode : Nodes){
        for (size_t i=0; i<3; ++i){
            if (itNode->externalViscosity[i] > 0){
				zeroExternalViscosity[i] = false;
			}
		}
	}
    /**
     * If after checking all nodes, there is at least one dimension of x, y & z that has zero viscosity, at
     * least one degree of freedom in that dimension should be fixed to have a unique displacement
     * solution to the system in the given dimension. In other words, if no node feels external resistanc ein x, then
     * the forces can fix the relative movement of all nodes in x, but nothing prevents the tissue rigid body displacement
     * in x, therefore resulting in infinite number of solutions.
     */
    if (zeroExternalViscosity[0] || zeroExternalViscosity[1] || zeroExternalViscosity[2]){
    	//At least one of the dimensions have zero viscosity
    	if(nNodes < 3) {
            std::cerr<<"Cannot run zero viscosity simulation with less than 3 nodes, run simulation at your own risk!"<<std::endl;
    	}        
        if (!symmetricY){
        	bool sufficientXFix = false;
        	bool sufficientYFix = false;
        	bool sufficientZFix = false;
            for (size_t a=0; a<3; ++a){
        		//if zeroExternalViscosity[0] is false, it means there is already  sufficient x fix
        		//OR, if I spotted sufficientXFix in previous nodes, it persists,
        		//OR, there is circumferential fix on x dimension
        		sufficientXFix = sufficientXFix || !zeroExternalViscosity[0] || CircumferentialNodeFix[a][0];
        		sufficientYFix = sufficientYFix || !zeroExternalViscosity[1] || CircumferentialNodeFix[a][1];
        		sufficientZFix = sufficientZFix || !zeroExternalViscosity[2] || CircumferentialNodeFix[a][2];
        	}
			//if there are no z fixe on the circumference, check the whole surfaces:
        	if (!sufficientXFix){
        		sufficientXFix = sufficientXFix || !zeroExternalViscosity[0] || BasalNodeFix[0] ;
				sufficientXFix = sufficientXFix || !zeroExternalViscosity[0] || ApicalNodeFix[0] ;
			}
        	if (!sufficientYFix){
        		sufficientYFix = sufficientYFix || !zeroExternalViscosity[1] || BasalNodeFix[1] ;
        		sufficientYFix = sufficientYFix || !zeroExternalViscosity[1] || ApicalNodeFix[1] ;
        	}
        	if (!sufficientZFix){
				sufficientZFix = sufficientZFix || !zeroExternalViscosity[2] || BasalNodeFix[2];
				sufficientZFix = sufficientZFix || !zeroExternalViscosity[2] || ApicalNodeFix[2];
			}

        	if (!sufficientXFix){
        		Nodes[ventralTipIndex]->FixedPos[0] = true;
        	}
        	if (!sufficientYFix){
        		Nodes[ventralTipIndex]->FixedPos[1] = true;
        	}
        	if (!sufficientZFix){
        		Nodes[ventralTipIndex]->FixedPos[2] = true;
        	}
        	if (!sufficientYFix){
        		Nodes[dorsalTipIndex]->FixedPos[1] = true;
        	}
			if (!sufficientZFix){
				Nodes[dorsalTipIndex]->FixedPos[2] = true;
			}
			if (!sufficientZFix){
				//if there is symmetricity, then the mid-line nodes will be fixed in y, and I do not need to fix the third node.
				// in fact, fixing the position of the third node in z will cause problems.
				if (dorsalTipIndex != 1 && 	ventralTipIndex!= 1 ){
					Nodes[1]->FixedPos[2] = true;
				}
				else if (dorsalTipIndex != 2 && ventralTipIndex!= 2 ){
					Nodes[2]->FixedPos[2] = true;
				}
				else if (dorsalTipIndex != 3 && ventralTipIndex!= 3 ){
					Nodes[3]->FixedPos[2] = true;
				}
			}
        }
        else{
        	//There is symmetricity, then there is y fix, if there are any node fixes, on the surfaces in z, then I do not need to do fix z:
        	bool sufficientZFix = false;
            for (size_t a=0; a<3; ++a){
        		sufficientZFix = sufficientZFix || !zeroExternalViscosity[2] || CircumferentialNodeFix[a][2];
			}
        	//if there are no z fixe on the circumferenece, check the whole surfaces:
			if (!sufficientZFix){
				sufficientZFix = sufficientZFix || !zeroExternalViscosity[2] || BasalNodeFix[2];
				sufficientZFix = sufficientZFix || !zeroExternalViscosity[2] || ApicalNodeFix[2];
			}
        	bool sufficientXFix = false;
            for (size_t a=0; a<3; ++a){
				sufficientXFix = sufficientXFix || !zeroExternalViscosity[0] || CircumferentialNodeFix[a][0] ;
			}
			//if there are no x fixes on the circumferenece, check the whole surfaces:
			if (!sufficientXFix){
				sufficientXFix = sufficientXFix || !zeroExternalViscosity[0] || BasalNodeFix[0] ;
				sufficientXFix = sufficientXFix || !zeroExternalViscosity[0] || ApicalNodeFix[0] ;
			}
			//there is no circumference or surface fixing of suffieint nature, then I should fix the nodes:
			if (!sufficientXFix){
              	Nodes[ventralTipIndex]->FixedPos[0] = true;
			}
			Nodes[ventralTipIndex]->FixedPos[1] = true;
			if (!sufficientZFix){
				Nodes[ventralTipIndex]->FixedPos[2] = true;
				Nodes[dorsalTipIndex]->FixedPos[2] = true;
			}
        }
    }
}

void Simulation::fixNode0InPosition(double x, double y, double z){
	double dx = Nodes[0]->Position[0]-x;
	double dy = Nodes[0]->Position[1]-y;
	double dz = Nodes[0]->Position[2]-z;
    for (const auto& itNode : Nodes){
        itNode->Position[0] -=dx;
        itNode->Position[1] -=dy;
        itNode->Position[2] -=dz;
	}
    for(const auto& itElement : Elements){
        itElement->updatePositions(Nodes);
	}
}

void Simulation::manualPerturbationToInitialSetup(bool deform, bool rotate){
    /**
     * This function allows the developer to add manual perturbations to initial setup.
     * The simple positonal perturbations are self-explanatory, the function is used
     * to test the response of any system to impulse perturbations. It also allows for
     * rigid body rotation adn translations.
     */
    //cout<<" inside manual perturbation, timestep "<<timestep<<std::endl;
	if(timestep==0){
        double scaleX = 2.0;
        double scaleY = 1.0;
        double scaleZ = 1.0;

        double PI = 3.14159265359;
        double tetX = 0 *PI/180.0;
        double tetY = 45 *PI/180.0;
        double tetZ = 0 *PI/180.0;
    	int axisToRotateOn = 1; //0: rotate around x axis, 1: around y-axis, and 2: around z-axis.
        if(deform){
            for (const auto& itNode : Nodes){
                itNode->Position[0] *=scaleX;
                itNode->Position[1] *=scaleY;
                itNode->Position[2] *=scaleZ;
            }
        }
        if (rotate){
            std::cerr<<"Rotating system"<<std::endl;
            double R[3][3] =  {{1,0,0},{0,1,0},{0,0,1}};
            if (axisToRotateOn == 0){
            	double c = cos(tetX);
            	double s = sin(tetX);
            	//Rx = {{1,0,0},{0,c,-s},{0,s,c}};
            	R[1][1] = c;
				R[1][2] = -s;
				R[2][1] = s;
				R[2][2] = c;
            }
            else if(axisToRotateOn == 1){
            	double c = cos(tetY);
            	double s = sin(tetY);
            	//Ry = {{c,0,s},{0,1,0},{-s,0,c}};
            	R[0][0] = c;
				R[0][2] = s;
				R[2][0] = -s;
				R[2][2] = c;
            }
            else if(axisToRotateOn == 2){
            	double c = cos(tetZ);
            	double s = sin(tetZ);
            	//Rz = {{c,-s,0},{s,c,0},{0,0,1}};
            	R[0][0] = c;
				R[0][1] = -s;
				R[1][0] = s;
				R[1][1] = c;
            }
            for (const auto&  itNode : Nodes){
                double x = itNode->Position[0]*R[0][0] + itNode->Position[1]*R[0][1] + itNode->Position[2]*R[0][2];
                double y = itNode->Position[0]*R[1][0] + itNode->Position[1]*R[1][1] + itNode->Position[2]*R[1][2];
                double z = itNode->Position[0]*R[2][0] + itNode->Position[1]*R[2][1] + itNode->Position[2]*R[2][2];
                itNode->Position[0]=x;
                itNode->Position[1]=y;
                itNode->Position[2]=z;
            }
        }
        for(const auto& itElement : Elements){
            itElement->updatePositions(Nodes);
        }
    }
}

void Simulation::updateGrowthRotationMatrices(){
    /**
     * Growth rotation matrices are updated for each element via ShapeBase#CalculateGrowthRotationByF.
     */
    for(const auto & itElement : Elements){
        itElement->CalculateGrowthRotationByF();
    }
}


std::array<double,3> Simulation::calculateNewNodePosForPeripodialNodeAddition(int nodeId0, int nodeId1, int nodeId2, double sideThickness){
    std::array<double,3> vec1{0.0};
	if (symmetricY && nodeId1 == -100){
		vec1[0] = -1.0;
	}
	else if (symmetricY && nodeId2 == -100){
		vec1[0] =  1.0;
	}
	else{
        std::array<double,3> vec2{0.0};
        for (size_t j=0; j<Nodes[nodeId0]->nDim; j++){
			vec1[j] = (Nodes[nodeId0]->Position[j] - Nodes[nodeId1]->Position[j]);
			vec2[j] = (Nodes[nodeId0]->Position[j] - Nodes[nodeId2]->Position[j]);
		}
        (void) Elements[0]->normaliseVector3D(vec1);
        (void) Elements[0]->normaliseVector3D(vec2);
		vec1[0] += vec2[0];
		vec1[1] += vec2[1];
		vec1[2] += vec2[2];
		double vec1Mag2 = vec1[0]*vec1[0] + vec1[1]*vec1[1] + vec1[2]*vec1[2];
		if (vec1Mag2 < 1E-5){
			//the two nodes are linear, the resulting vector is of zero length.
			//I will rotate one of the vectors 90 degrees and it will be pointing in the correct orientation, the actual direction can be fixed in the below loop:
			//this vector is already normalised, I am skipping the normalisation step
			vec1[0] = -1.0*vec2[1];
			vec1[1] = vec2[0];
		}
		else{
            (void) Elements[0]->normaliseVector3D(vec1);
		}
		//now I have the vector to point out from the base node 0. BUT, this will point outwards only if the tissue curvature is convex at all points
		//I need to check if it actually is pointing out, as the experimentally driven tissues can be concave at points.
		//the cross produc of vec[2] and the vector to the cell centre should have the opposite sign with the corss product of my orientation vector and vector 2.
        std::array<double,3> vecCentre{0.0};
		vecCentre[0] =  SystemCentre[0] - Nodes[nodeId0]->Position[0];
		vecCentre[1] =  SystemCentre[1] - Nodes[nodeId0]->Position[1];
		vecCentre[2] =  SystemCentre[2] - Nodes[nodeId0]->Position[2];
        (void) Elements[0]->normaliseVector3D(vecCentre);
        std::array<double,3> cross1 = Elements[0]->crossProduct3D(vec2,vecCentre);
        (void) Elements[0]->normaliseVector3D(cross1);
        std::array<double,3> cross2 = Elements[0]->crossProduct3D(vec2,vec1);
        (void) Elements[0]->normaliseVector3D(cross2);
		double dotp = Elements[0]->dotProduct3D(cross1,cross2);
		if (dotp >0 ){
			//the vectors are pointing to the same direction! Need to rotate vec1 180 degrees:
			vec1[0] *= -1.0;
			vec1[1] *= -1.0;
			vec1[2] *= -1.0;
		}
	}
    std::array<double,3> pos{0.0};
	pos[0] = Nodes[nodeId0]->Position[0] + vec1[0]*sideThickness;
	pos[1] = Nodes[nodeId0]->Position[1] + vec1[1]*sideThickness;
	pos[2] = Nodes[nodeId0]->Position[2];
    return pos;
}


std::array<double,3> Simulation::calculateNewNodePosForPeripodialNodeAddition(int nodeId0, int nodeId1, double sideThickness){
    std::array<double,3> vec0{0.0};
    std::array<double,3> midpoint{0.0};
    for (size_t j=0; j<Nodes[nodeId0]->nDim; j++){
		vec0[j] = (Nodes[nodeId1]->Position[j] - Nodes[nodeId0]->Position[j])/2.0;
		midpoint[j] = Nodes[nodeId0]->Position[j] + vec0[j];
	}
    (void) Elements[0]->normaliseVector3D(vec0);
	//The list is sorted counter-cock-wise, to point out, I will rotate normalised vector v0 -90 degrees on z axis:
	// (x,y,z) -> (y,-x,z);
	// then I will add this vector to the calculated mid point to gt the new node's position.
    std::array<double,3> pos {0.0};
    pos[0] = midpoint[0] + vec0[1]*sideThickness;
	pos[1] = midpoint[1] - vec0[0]*sideThickness;
	pos[2] = midpoint[2];
    return pos;
}

bool Simulation::checkIfThereIsPeripodialMembrane(){
	bool Success = true;
    for(const auto& itNode : Nodes){
        if (itNode->tissueType == 1){
			thereIsPeripodialMembrane = true;
			//I have not added the peripodial membrane as yet,
			//if there is one, it came from the input mesh
			//Then I want to set up the circumferencial input properly.
			setLinkerCircumference();
			break;
		}
	}
	return Success;
}

void Simulation::checkForExperimentalSetupsBeforeIteration(){
	if (stretcherAttached){
        //moveClampedNodesForStretcher();
	}
}

void Simulation::checkForExperimentalSetupsWithinIteration(){
}

void Simulation::checkForExperimentalSetupsAfterIteration(){

}

void Simulation::calculateStiffnessChangeRatesForActin(int idOfCurrentStiffnessPerturbation){
    /**
      * If the actin layer has perturbations, the elements are checked if they are influenced via
      * ShapeBase#isActinStiffnessChangeAppliedToElement and
      * rates are set in this function via ShapeBase#calculateStiffnessPerturbationRate.
      */
    startedStiffnessPerturbation[idOfCurrentStiffnessPerturbation] = true;
    //not using range based loops here to ensure openMP comaptibility
    #pragma omp parallel for
    for(std::vector<std::unique_ptr<ShapeBase>>::iterator itElement = Elements.begin(); itElement <Elements.end(); ++itElement){
        bool applyToThisElement =  (*itElement)->isActinStiffnessChangeAppliedToElement(ThereIsWholeTissueStiffnessPerturbation[idOfCurrentStiffnessPerturbation], ThereIsApicalStiffnessPerturbation[idOfCurrentStiffnessPerturbation], ThereIsBasalStiffnessPerturbation[idOfCurrentStiffnessPerturbation], ThereIsBasolateralWithApicalRelaxationStiffnessPerturbation[idOfCurrentStiffnessPerturbation],ThereIsBasolateralStiffnessPerturbation[idOfCurrentStiffnessPerturbation], stiffnessPerturbationEllipseBandIds[idOfCurrentStiffnessPerturbation], numberOfStiffnessPerturbationAppliesEllipseBands[idOfCurrentStiffnessPerturbation]);
        if (applyToThisElement){
            (*itElement)->calculateStiffnessPerturbationRate(ThereIsBasolateralWithApicalRelaxationStiffnessPerturbation[idOfCurrentStiffnessPerturbation], stiffnessPerturbationBeginTimeInSec[idOfCurrentStiffnessPerturbation],stiffnessPerturbationEndTimeInSec[idOfCurrentStiffnessPerturbation], stiffnessChangedToFractionOfOriginal[idOfCurrentStiffnessPerturbation]);
        }
    }
}

void Simulation::updateStiffnessChangeForActin(int idOfCurrentStiffnessPerturbation){
    /**
      * If the actin layer has perturbations, the elements are checked if they are influenced via
      * ShapeBase#isActinStiffnessChangeAppliedToElement, the ShapeBase#stiffnessMultiplier is altered via
      * ShapeBase#updateStiffnessMultiplier and elastic property tensors are updated
      * via ShapeBase#updateElasticProperties.
      */
    //not using range based loops here to ensure openMP comaptibility
    #pragma omp parallel for
    for(std::vector<std::unique_ptr<ShapeBase>>::iterator itElement = Elements.begin(); itElement <Elements.end(); ++itElement){
        bool applyToThisElement = (*itElement)->isActinStiffnessChangeAppliedToElement(ThereIsWholeTissueStiffnessPerturbation[idOfCurrentStiffnessPerturbation], ThereIsApicalStiffnessPerturbation[idOfCurrentStiffnessPerturbation], ThereIsBasalStiffnessPerturbation[idOfCurrentStiffnessPerturbation], ThereIsBasolateralWithApicalRelaxationStiffnessPerturbation[idOfCurrentStiffnessPerturbation],  ThereIsBasolateralStiffnessPerturbation[idOfCurrentStiffnessPerturbation], stiffnessPerturbationEllipseBandIds[idOfCurrentStiffnessPerturbation], numberOfStiffnessPerturbationAppliesEllipseBands[idOfCurrentStiffnessPerturbation]);
        if (applyToThisElement){
            (*itElement)->updateStiffnessMultiplier(dt);
            (*itElement)->updateElasticProperties();
        }
    }
}

void Simulation::checkStiffnessPerturbation(){
    /**
     * For all input stiffness perturbaitons, the perturbation activity time is checked from the Simulation#currSimTimeSec
     * being between the selected Simulation#stiffnessPerturbationBeginTimeInSec and Simulation#stiffnessPerturbationEndTimeInSec.
     * If this is the first time the perturbation is being called (checked via Simulation#startedStiffnessPerturbation flag) then
     * the rate is calculated via Simulation#calculateStiffnessChangeRatesForActin. Then stiffnesses are uodated via
     * Simulation#updateStiffnessChangeForActin.
     */
	int n = stiffnessPerturbationBeginTimeInSec.size();
	for (int idOfCurrentStiffnessPerturbation=0; idOfCurrentStiffnessPerturbation<n; ++idOfCurrentStiffnessPerturbation){
		if (currSimTimeSec >=stiffnessPerturbationBeginTimeInSec[idOfCurrentStiffnessPerturbation] && currSimTimeSec <stiffnessPerturbationEndTimeInSec[idOfCurrentStiffnessPerturbation]){
			if (startedStiffnessPerturbation[idOfCurrentStiffnessPerturbation] == false){
				calculateStiffnessChangeRatesForActin(idOfCurrentStiffnessPerturbation);
			}
			updateStiffnessChangeForActin(idOfCurrentStiffnessPerturbation);
		}
	}
}

void Simulation::updateChangeForExplicitECM(int idOfCurrentECMPerturbation){
    /**
      * If the ECM has perturbationson elasticity, the elements are checked if they are influenced via
      * ShapeBase#isECMChangeAppliedToElement, the ShapeBase#stiffnessMultiplier is altered via
      * ShapeBase#updateStiffnessMultiplier and elastic property tensors are updated
      * via ShapeBase#updateElasticProperties.
      */
    //not using range based loops here to ensure openMP comaptibility
    #pragma omp parallel for
    for(std::vector<std::unique_ptr<ShapeBase>>::iterator itElement = Elements.begin(); itElement <Elements.end(); ++itElement){
        bool applyToThisElement = (*itElement)->isECMChangeAppliedToElement(changeApicalECM[idOfCurrentECMPerturbation], changeBasalECM[idOfCurrentECMPerturbation], ECMChangeEllipseBandIds[idOfCurrentECMPerturbation], numberOfECMChangeEllipseBands[idOfCurrentECMPerturbation]);
        if (applyToThisElement){
             (*itElement)->updateStiffnessMultiplier(dt);
             (*itElement)->updateElasticProperties();
        }
    }
}

void Simulation::updateChangeForViscosityBasedECMDefinition(int idOfCurrentECMPerturbation){
    /**
      * If the ECM has perturbations on viscosity, the nodes are checked if they are influenced via
      * Node#tissuePlacement and Simulation#changeBasalECM - Simulation#changeApicalECM boolean arrays,
      * as well as the marker ellipse Ids. The visocisties of each node are updated via
      * Node#ECMViscosityChangePerHour and the time step Simulation#dt (in sec). Once new viscosity is
      * calculated, it is capped to be positive and at most equal to the maximum external viscosity.
      */
   //not using range based loops here to ensure openMP comaptibility
   #pragma omp parallel for
   for (std::vector<std::unique_ptr<Node>>::iterator itNode = Nodes.begin(); itNode<Nodes.end(); ++itNode){
        if(((*itNode)->tissuePlacement == 0 && changeBasalECM[idOfCurrentECMPerturbation] ) || ((*itNode)->tissuePlacement == 1 && changeApicalECM[idOfCurrentECMPerturbation] )){
            if((*itNode)->insideEllipseBand){
                for (int ECMReductionRangeCounter =0; ECMReductionRangeCounter<numberOfECMChangeEllipseBands[idOfCurrentECMPerturbation]; ++ECMReductionRangeCounter){
                    if ((*itNode)->coveringEllipseBandId == ECMChangeEllipseBandIds[idOfCurrentECMPerturbation][ECMReductionRangeCounter]){
                        for (size_t i =0; i<(*itNode)->nDim; ++i){
                            double viscosityChange = (*itNode)->ECMViscosityChangePerHour[i]/3600*dt;
                            double newViscosity = (*itNode)->externalViscosity[i] - viscosityChange;
                            //avoiding setting negative viscosity!
                            if (newViscosity>(*itNode)->maximumExternalViscosity[i]){
                                (*itNode)->externalViscosity[i] =(*itNode)->maximumExternalViscosity[i];
                            }
                            else if(newViscosity<(*itNode)->minimumExternalViscosity[i]){
                                (*itNode)->externalViscosity[i] =(*itNode)->minimumExternalViscosity[i];
                            }
                            else{
                                (*itNode)->externalViscosity[i] = newViscosity;
                            }
                        }
                    }
                }
            }
        }
	}
}

void Simulation::calculateChangeRatesForECM(int idOfCurrentECMPerturbation){
    /**
      * If the ECM has perturbations, the elements are checked if they are influenced via
      * ShapeBase#isECMChangeAppliedToElement and
      * rates are set in this function via ShapeBase#calculateStiffnessPerturbationRate.
      */
    changedECM[idOfCurrentECMPerturbation] = true; //this will not be used for emergent ecm perturbations
	//this is the first time step I am changing the ECM stiffness.
	//I need to calculate rates first.
	//If there is explicit ECM, I will calculate the young modulus change via elements.
	//If there is no explicit ECM, the viscosity reflects the ECM stiffness, and I will change the viscosity on a nodal basis.
	if( thereIsExplicitECM){
        //not using range based loops here to ensure openMP comaptibility
        #pragma omp parallel for
        for(std::vector<std::unique_ptr<ShapeBase>>::iterator itElement = Elements.begin(); itElement <Elements.end(); ++itElement){
            bool applyToThisElement = (*itElement)->isECMChangeAppliedToElement(changeApicalECM[idOfCurrentECMPerturbation], changeBasalECM[idOfCurrentECMPerturbation], ECMChangeEllipseBandIds[idOfCurrentECMPerturbation], numberOfECMChangeEllipseBands[idOfCurrentECMPerturbation]);
			if (applyToThisElement){
				//the first input is used for checking basolateral stiffenning combined with apical relaxation
				//the ECM does not have such options. Will give the boolean as false and continue.
                (*itElement)->calculateStiffnessPerturbationRate(false, ECMChangeBeginTimeInSec[idOfCurrentECMPerturbation],ECMChangeEndTimeInSec[idOfCurrentECMPerturbation], ECMStiffnessChangeFraction[idOfCurrentECMPerturbation]);
			}
		}
		//now I need to check for the viscosity based calculation:
        //not using range based loops here to ensure openMP comaptibility
        /**
         * For each node, the applicability of ECM change is checked via Node#isECMChangeAppliedToNode, and the rate
         * Node#ECMViscosityChangePerHour per dimension is calculated via the change fraction Simulation#ECMViscosityChangeFraction, the
         * total time of applied change obtined from Simulation#ECMChangeBeginTimeInSec and Simulation#ECMChangeEndTimeInSec.
         */
        #pragma omp parallel for
        for (std::vector<std::unique_ptr<Node>>::iterator itNode = Nodes.begin(); itNode<Nodes.end(); ++itNode){
            bool applyToThisNode = (*itNode)->isECMChangeAppliedToNode(changeApicalECM[idOfCurrentECMPerturbation], changeBasalECM[idOfCurrentECMPerturbation], ECMChangeEllipseBandIds[idOfCurrentECMPerturbation], numberOfECMChangeEllipseBands[idOfCurrentECMPerturbation]);
			if(applyToThisNode){
				double timeDifferenceInHours = (ECMChangeEndTimeInSec[idOfCurrentECMPerturbation] - ECMChangeBeginTimeInSec[idOfCurrentECMPerturbation])/3600;
				for (int i=0; i<3; ++i){
                    (*itNode)->ECMViscosityChangePerHour[i] = (*itNode)->initialExternalViscosity[i]*(1-ECMViscosityChangeFraction[idOfCurrentECMPerturbation])/timeDifferenceInHours;
                    if ((*itNode)->ECMViscosityChangePerHour[i]<0){
                        (*itNode)->maximumExternalViscosity[i] = (*itNode)->initialExternalViscosity[i]*ECMViscosityChangeFraction[idOfCurrentECMPerturbation];
					}
					else{
                        (*itNode)->minimumExternalViscosity[i] = (*itNode)->initialExternalViscosity[i]*ECMViscosityChangeFraction[idOfCurrentECMPerturbation];
					}
				}
			}
		}
	}
	else{
        /**
         * The viscosity based ECM perturbation can be applied for ECM definitions beyond explicit ECM, and this is purely on nodal viscosity basis.
         */
		double timeDifferenceInHours = (ECMChangeEndTimeInSec[idOfCurrentECMPerturbation] - ECMChangeBeginTimeInSec[idOfCurrentECMPerturbation])/3600;
        //not using range based loops here to ensure openMP comaptibility
        #pragma omp parallel for
        for (std::vector<std::unique_ptr<Node>>::iterator itNode = Nodes.begin(); itNode<Nodes.end(); ++itNode){
            bool applyToThisNode = (*itNode)->isECMChangeAppliedToNode(changeApicalECM[idOfCurrentECMPerturbation], changeBasalECM[idOfCurrentECMPerturbation], ECMChangeEllipseBandIds[idOfCurrentECMPerturbation], numberOfECMChangeEllipseBands[idOfCurrentECMPerturbation]);
			if(applyToThisNode){
				for (int i=0; i<3; ++i){
                    (*itNode)->ECMViscosityChangePerHour[i] = (*itNode)->initialExternalViscosity[i]*(1-ECMViscosityChangeFraction[idOfCurrentECMPerturbation])/timeDifferenceInHours;
                    if ((*itNode)->ECMViscosityChangePerHour[i]<0){
                        (*itNode)->maximumExternalViscosity[i] = (*itNode)->initialExternalViscosity[i]*ECMViscosityChangeFraction[idOfCurrentECMPerturbation];
					}
					else{
                        (*itNode)->minimumExternalViscosity[i] = (*itNode)->initialExternalViscosity[i]*ECMViscosityChangeFraction[idOfCurrentECMPerturbation];
					}
				}
			}
		}
	}
}

void Simulation::updateECMRenewalHalflifeMultiplier(int idOfCurrentECMPerturbation){
    /**
     * If there is perturbation on ECM renewal halflife, the type of perturbation emergence is checked,
     * if the perturbation initiation is emergend based on topology (active at fold grove basal sides), then
     * the rate is calculated from the total application time.
     */
    if(thereIsExplicitECM){
		if (ECMChangeTypeIsEmergent[idOfCurrentECMPerturbation]){
			double totalTimeChange = ECMChangeEndTimeInSec[idOfCurrentECMPerturbation] - ECMChangeBeginTimeInSec[idOfCurrentECMPerturbation];
			double incrementPerSec = (ECMRenewalHalfLifeTargetFraction[idOfCurrentECMPerturbation] - 1.0) / totalTimeChange;
            //not using range based loops here to ensure openMP comaptibility
            #pragma omp parallel for
            for(std::vector<std::unique_ptr<ShapeBase>>::iterator itElement = Elements.begin(); itElement <Elements.end(); ++itElement){
                if ((*itElement)->isECMMimicing){
                    bool applyToThisElement = (*itElement)->isECMChangeAppliedToElement(changeApicalECM[idOfCurrentECMPerturbation], changeBasalECM[idOfCurrentECMPerturbation], ECMChangeEllipseBandIds[idOfCurrentECMPerturbation], numberOfECMChangeEllipseBands[idOfCurrentECMPerturbation]);
					if (applyToThisElement){
                        (*itElement)->plasticDeformationHalfLifeMultiplier +=  incrementPerSec*dt;
                        if (incrementPerSec<0 && (*itElement)->plasticDeformationHalfLifeMultiplier < ECMRenewalHalfLifeTargetFraction[idOfCurrentECMPerturbation]){
                            (*itElement)->plasticDeformationHalfLifeMultiplier = ECMRenewalHalfLifeTargetFraction[idOfCurrentECMPerturbation];
						}
                        if (incrementPerSec>0 && (*itElement)->plasticDeformationHalfLifeMultiplier > ECMRenewalHalfLifeTargetFraction[idOfCurrentECMPerturbation]){
                            (*itElement)->plasticDeformationHalfLifeMultiplier = ECMRenewalHalfLifeTargetFraction[idOfCurrentECMPerturbation];
						}
					}
				}
			}
		}
		else{
			if (currSimTimeSec>ECMChangeBeginTimeInSec[idOfCurrentECMPerturbation]){
				//perturbation on ECM renewal half life started
				double currECMRenewaHalfLifeMultiplier = 1.0;
				if (currSimTimeSec>=ECMChangeEndTimeInSec[idOfCurrentECMPerturbation]){
					currECMRenewaHalfLifeMultiplier = ECMRenewalHalfLifeTargetFraction[idOfCurrentECMPerturbation];
				}
				else{
					double totalTimeChange = ECMChangeEndTimeInSec[idOfCurrentECMPerturbation] - ECMChangeBeginTimeInSec[idOfCurrentECMPerturbation];
					double currTimeChange = currSimTimeSec - ECMChangeBeginTimeInSec[idOfCurrentECMPerturbation];
					currECMRenewaHalfLifeMultiplier = 1 + (ECMRenewalHalfLifeTargetFraction[idOfCurrentECMPerturbation] - 1.0) * currTimeChange / totalTimeChange;
				}
                //not using range based loops here to ensure openMP comaptibility
                #pragma omp parallel for
                for(std::vector<std::unique_ptr<ShapeBase>>::iterator itElement = Elements.begin(); itElement <Elements.end(); ++itElement){
                   if ((*itElement)->isECMMimicing){
                        bool applyToThisElement = (*itElement)->isECMChangeAppliedToElement(changeApicalECM[idOfCurrentECMPerturbation], changeBasalECM[idOfCurrentECMPerturbation], ECMChangeEllipseBandIds[idOfCurrentECMPerturbation], numberOfECMChangeEllipseBands[idOfCurrentECMPerturbation]);
						if (applyToThisElement){
                            (*itElement)->plasticDeformationHalfLifeMultiplier =  currECMRenewaHalfLifeMultiplier;
						}
					}
				}
			}
		}
	}
}

void Simulation::checkECMChange(){
    /**
     * If there is explicit ECM, then first the compartment based perturbations are checked.
     * Then each perturbation based on tissue placement, such as apical, basal or emergnt with markers, are checked.
     * The stiffness, viscosity and half life updates are carried out via ShapeBase#updateStiffnessMultiplier,
     * ShapeBase#updateChangeForViscosityBasedECMDefinition, ShapeBase#updateChangeForExplicitECM and ShapeBase#updateECMRenewalHalflifeMultiplier.
     */
    if( thereIsExplicitECM){
        //not using range based loops here to ensure openMP comaptibility
        #pragma omp parallel for
        for(std::vector<std::unique_ptr<ShapeBase>>::iterator itElement = Elements.begin(); itElement <Elements.end(); ++itElement){
			bool updateStiffness = false;
            if( (*itElement)->isECMMimicing && (*itElement)->tissuePlacement == 0 && (*itElement)->tissueType ==0 ){//columnar basal ecmmimicking element
				//check notum:
                if (notumECMChangeFraction != 1.0 && (*itElement)->compartmentType == 2){//notum:
					if (currSimTimeSec>= notumECMChangeInitTime && currSimTimeSec< notumECMChangeEndTime){
                        double currentFraction = 1.0 + (notumECMChangeFraction-1)*(*itElement)->compartmentIdentityFraction;
                        (*itElement)->calculateStiffnessPerturbationRate(false, notumECMChangeInitTime,notumECMChangeEndTime, currentFraction);
						updateStiffness=true;
					}
				}
				//check hinge:
                if (hingeECMChangeFraction != 1.0 && (*itElement)->compartmentType == 1){//hinge:
					if (currSimTimeSec>= hingeECMChangeInitTime && currSimTimeSec< hingeECMChangeEndTime){
                        double currentFraction = 1.0 + (hingeECMChangeFraction-1)*(*itElement)->compartmentIdentityFraction;
                        (*itElement)->calculateStiffnessPerturbationRate(false, hingeECMChangeInitTime,hingeECMChangeEndTime, currentFraction);
						updateStiffness=true;
					}
				}
				//check pouch:
                if (pouchECMChangeFraction != 1.0 && (*itElement)->compartmentType == 0){//pouch:
					if (currSimTimeSec>= pouchECMChangeInitTime && currSimTimeSec< pouchECMChangeEndTime){
                        double currentFraction = 1.0 + (pouchECMChangeFraction-1)*(*itElement)->compartmentIdentityFraction;
                        (*itElement)->calculateStiffnessPerturbationRate(false, pouchECMChangeInitTime,pouchECMChangeEndTime, currentFraction);
						updateStiffness=true;
					}
				}
			}
			if (updateStiffness){
                (*itElement)->updateStiffnessMultiplier(dt);
                (*itElement)->updateElasticProperties();
			}
		}
	}
	//Now go through ellipses, ellipses overwrite the compartment based definitions
	int n = ECMChangeBeginTimeInSec.size();
	for (int idOfCurrentECMPerturbation=0; idOfCurrentECMPerturbation<n; ++idOfCurrentECMPerturbation){
		if (ECMChangeTypeIsEmergent[idOfCurrentECMPerturbation]){
			if (currSimTimeSec >=ECMChangeBeginTimeInSec[idOfCurrentECMPerturbation]){
				calculateChangeRatesForECM(idOfCurrentECMPerturbation);
				if( thereIsExplicitECM){
					updateECMRenewalHalflifeMultiplier(idOfCurrentECMPerturbation);
					updateChangeForExplicitECM(idOfCurrentECMPerturbation);
				}
				updateChangeForViscosityBasedECMDefinition(idOfCurrentECMPerturbation);
			}
		}
		else{
			if (currSimTimeSec >=ECMChangeBeginTimeInSec[idOfCurrentECMPerturbation] && currSimTimeSec <ECMChangeEndTimeInSec[idOfCurrentECMPerturbation]){
				if (changedECM[idOfCurrentECMPerturbation] == false){
					calculateChangeRatesForECM(idOfCurrentECMPerturbation);
				}
				if( thereIsExplicitECM){
					updateECMRenewalHalflifeMultiplier(idOfCurrentECMPerturbation);
					updateChangeForExplicitECM(idOfCurrentECMPerturbation);
				}
				updateChangeForViscosityBasedECMDefinition(idOfCurrentECMPerturbation);
			}
		}
	}
}

void Simulation::checkEllipseAllocationWithCurvingNodes(){
    /**
     * The emergent curved regions of the tissue are marked via ellipse ban ids.
     * These bands can be defined by the user, as such marker ellipse band Ids 100, 101 and 102 are reserved
     * for emergent band initiation. \n
     * - If a surface is on the apical fold initiation surface, then it acquires Id 100. \n
     * - If a surface is on the basal fold initiation, then is acquires Id 101. \n
     * - If there are emergent perturbations, and initiation of perturbation is conditional (
     * tissue shape change perturbation starts only after the ECM is reduced below a certain fraction of its original value), then
     * the elements initially tagged 100 are updated to tag 102 once htis threshold (Simulation#shapeChangeECMLimit) is reached. The
     * function Simulation#checkForEllipseIdUpdateWithECMDegradation
     * will carry out this operation.
     */
	int selectedEllipseBandId =100;
    for (const auto& itNode : Nodes){
        if (itNode->onFoldInitiation){
            if(itNode->coveringEllipseBandId != 100 && itNode->coveringEllipseBandId != 101){
                if (itNode->tissuePlacement == 1){ //apical collapse: ECM relaxation, cell shortening, volume redistribution to shrink top
					selectedEllipseBandId = 100;
				}
				else{ //basal collapse, volume redistribution to shrink bottom
					selectedEllipseBandId = 101;
				}
                size_t n = itNode->connectedElementIds.size();
                for (size_t i=0; i<n; ++i){
                    int currElementId = itNode->connectedElementIds[i];
					if (Elements[currElementId]->coveringEllipseBandId == 100 || Elements[currElementId]->coveringEllipseBandId == 101 || Elements[currElementId]->isECMMimimcingAtCircumference){
						continue;
					}
					bool changeEllipseId = Elements[currElementId]->hasEnoughNodesOnCurve(Nodes);
					if (changeEllipseId){
						Elements[currElementId]->coveringEllipseBandId = selectedEllipseBandId;
						Elements[currElementId]->assignEllipseBandIdToWholeTissueColumn(TissueHeightDiscretisationLayers,Nodes,Elements);
					}
				}
			}
		}
	}
}

void Simulation::checkForLeftOutElementsInEllipseAssignment(){
    /**
     * When elements are assigned emergent marker identities, their nodes will follow.
     * As a result, in densly folding regions, an element can have all its nodes assigned to a specific marker
     * by its neighbours, without itself being assigned to it. These elements themselves will be tagged via their nodes.
     */
    for (const auto& itNode : Nodes){
        if (itNode->checkOwnersforEllipseAsignment){
            int currentEllipseId = itNode->coveringEllipseBandId;
            int n = itNode->connectedElementIds.size();
			for (int i=0; i<n ; ++i){
                int currElementId = itNode->connectedElementIds[i];
				if (Elements[currElementId]->coveringEllipseBandId ==currentEllipseId ){
					continue;
				}
                const std::vector<int>& nodeIds = Elements[currElementId]->getNodeIds();
				bool hasNodeOutsideEllipse = false;
                for (auto currNodeId : nodeIds){
					if (Nodes[currNodeId]->coveringEllipseBandId != currentEllipseId){
						hasNodeOutsideEllipse = true;
						break;
					}
				}
				if (!hasNodeOutsideEllipse){
					//the element does not have any nodes outside the ellipse
					//should assign to ellipse
					Elements[currElementId]->coveringEllipseBandId = currentEllipseId;
					Elements[currElementId]->assignEllipseBandIdToWholeTissueColumn(TissueHeightDiscretisationLayers,Nodes,Elements);
				}
			}
		}
	}
}


void Simulation::updateEllipseWithCollapse(){
    /**
     * If there is an apical collapse of nodes, then this elemetn must be on a curving surface, thus the emergent marker ids are updated.
     * The element is not checked if it is already assigned to a fold via emergent mariking (ShapeBase#coveringEllipseBandId is 100 or 101).
     * Lateral ECM elements are not modified. If the element is not already assigned to an emergent marker group (on fold initiation), then
     * its collapse status is checked via ShapeBase#checkForCollapsedNodes.
     */
    for (const auto& itEle : Elements){
        if (itEle->coveringEllipseBandId == 100 || itEle->coveringEllipseBandId == 101 || itEle->isECMMimimcingAtCircumference){
			continue;
		}
        if(itEle->tissuePlacement == 0 || itEle->tissuePlacement == 1 || (itEle->tissuePlacement == 2 && itEle->spansWholeTissue ) ){
            itEle->checkForCollapsedNodes(TissueHeightDiscretisationLayers, Nodes, Elements);
		}
	}
}

void Simulation::checkForEllipseIdUpdateWithECMDegradation(){
    /**
     * If there are emergent marker allocations, and initiation of a perturbation is conditional (
     * tissue shape change perturbation starts only after the ECM is reduced below a certain fraction of its original value), then
     * the elements initially tagged 100 are updated to tag 102 once this threshold (Simulation#shapeChangeECMLimit) is reached.
     */
    for (const auto& itEle : Elements){
        if (itEle->isECMMimicing && itEle->coveringEllipseBandId == 100){
            if (itEle->stiffnessMultiplier<shapeChangeECMLimit){
                itEle->coveringEllipseBandId = 102;
                itEle->assignEllipseBandIdToWholeTissueColumn(TissueHeightDiscretisationLayers,Nodes,Elements);
			}
		}
	}
}

void Simulation::updateOnFoldNodesFromCollapse(){
    /**
     * When an element collapses its nodes, it is a sign that it is residing on a fold. These elements
     * assign their nodes to be on folds, by flagging the boolean Node#onFold.
     */
	if (checkedForCollapsedNodesOnFoldingOnce){
		return;
	}
	//update detection threshold  grid:
	double periAverageSideLength = 0,colAverageSideLength = 0;
	getAverageSideLength(periAverageSideLength, colAverageSideLength);
	if (thereIsPeripodialMembrane){
		colAverageSideLength = (periAverageSideLength+colAverageSideLength)/2.0;
	}
    for (size_t i=0;i<10;++i){
        for(size_t j=0;j<5;++j){
			packingDetectionThresholdGrid[i][j] = 0.4 * packingDetectionThresholdGrid[i][j];
		}
	}
	//check collapsed nodes, run this code only on initiation!
    for (auto& itNode : Nodes){
        if ( itNode->collapsedWith.size()>0) {
			bool checkForFold = true;
            int collapsedId = itNode->collapsedWith[0];
            if (collapsedId == itNode->Id){
                if (itNode->collapsedWith.size()>1){
                    collapsedId = itNode->collapsedWith[1];
				}
				else{
					checkForFold=false;
				}
			}
			if (checkForFold){
                assignFoldRegionAndReleasePeripodial(itNode.get(), Nodes[collapsedId].get());
			}
		}
	}
	checkedForCollapsedNodesOnFoldingOnce = true;
}

void Simulation::checkForEmergentEllipseFormation(){
    /**
     * The emergent marker ellipse Ids, reserved 100, 101 and 102, are assigned to
     * curved surfaces of hte tissue. This information can be obtained via
     * collapse of an elemental surface, or adhesion of two nodes. Once the
     * identity is assigned, it can be updated to mark ECM relaxation threshold.
     * The respnsible functions are Simulation#updateEllipseWithCollapse,
     * Simulation#updateOnFoldNodesFromCollapse, Simulation#checkEllipseAllocationWithCurvingNodes,
     *  Simulation#checkForLeftOutElementsInEllipseAssignment and Simulation#checkForEllipseIdUpdateWithECMDegradation.
     */
	updateEllipseWithCollapse();
	updateOnFoldNodesFromCollapse();

	checkEllipseAllocationWithCurvingNodes();
	checkForLeftOutElementsInEllipseAssignment();
	checkForEllipseIdUpdateWithECMDegradation();
}

void Simulation::artificialRelax(){
    /**
     * This function relaxes all accumulated forces on elements of the tissue by transfering the current elastic deformation
     * gradient onto growth gradient. The relaxation can optionally be extended to ECM, or be applied to only
     * cellular material of the tissue.
     */
    for (const auto& itEle : Elements){
		bool relaxElement = false;
        if (itEle->isECMMimicing){
			if( relaxECMInArtificialRelaxation ){
				relaxElement = true;
			}
		}
        else if(itEle->tissueType == 0 ){
			relaxElement = true;
		}
		if (relaxElement){
            itEle->relaxElasticForces();
		}
	}
}


bool Simulation::runOneStep(){
    /**
     * The core function in the Simulation class is the function to update the whole system from one time step to the next.
     */
    bool Success = true;
    for (size_t i=0; i< nElements; ++i){
		if (Elements[i]->tissuePlacement == 1 && Elements[i]->tissueType == 0){
			Elements[i]->calculateApicalArea();
		}
	}
    manualPerturbationToInitialSetup(false,false); //bool deform, bool rotate, optional perturbations at first time step, see the function documentation,
    /**
     * The process starts by resetting all forces of the current time step via Simulation#resetForces. Then emergent marker initiation is checked,
     * (Simulation#checkForEmergentEllipseFormation) as this could initiate physical property changes and should be carried out before the relavent updates. \n
     */
    resetForces(true); // reset the packing forces together with all the rest of the forces here
    checkForEmergentEllipseFormation();
    /**
     * Then the relative position updates of the tissue are carreid out.  First the tissue long axis is aligned on
     * the x axis (Simulation#alignTissueDVToXPositive), then the bounsing box is calculated (ShapeBase#calculateBoundingBox),
     * and the relative positions of apical (top layer) elements are obtained (ShapeBase#calculateRelativePosInBoundingBox). The
     * relative positions of the midline and bottom elements are updated to be the same as the apical layer, as growth should
     * not differ in the z-axis, the relative positions should read the same grid point from growth maps.
     */

    alignTissueDVToXPositive();
    calculateBoundingBox();
    calculateDVDistance();
    for(const auto& itElement : Elements){
        itElement->calculateRelativePosInBoundingBox(boundingBox[0][0],boundingBox[0][1],boundingBoxSize[0],boundingBoxSize[1]);
    }
    updateRelativePositionsToApicalPositioning();
    /**
     * Any experimental setup, that requires an update prior to entring the NR iterations for the
     * time step are carried out by the function Simulation#checkForExperimentalSetupsBeforeIteration. This is folllowed by
     * any existing perturbation on physical properties, through functions Simulation#checkStiffnessPerturbation and
     * Simulation#checkECMChange. \n
     */
    checkForExperimentalSetupsBeforeIteration();
    if (ThereIsStiffnessPerturbation) {
    	checkStiffnessPerturbation();
    }
    if (thereIsECMChange) {
    	checkECMChange();
    }
    /**
     * Once the physical proerties are updated, the growth of each element is calculated. For growth functions based on reading
     * gorwht maps (the main attribute utilised in morphogenesis simulations), the relative position in usecan be updated less
     * frequently then each step, this update is checked in Simulation#checkForPinningPositionsUpdate. Then the rigid body rotations
     * around z axis are extracted to be eliminated through ShapeBase#updateGrowthRotationMatrices. For all growth functions and
     * induced mutant clones the growth is calculated in Simulation#calculateGrowth. Simularly, shape changes are updated via
     * Simulation#calculateShapeChange. \n
     */
    if(nGrowthFunctions>0 || nShapeChangeFunctions >0 || numberOfClones >0 ){
    	checkForPinningPositionsUpdate();
        updateGrowthRotationMatrices();
        if(nGrowthFunctions>0 || numberOfClones > 0){
        	calculateGrowth();
        }
        if(nShapeChangeFunctions>0){
        	calculateShapeChange();
        }
    }
    /**
     * Based on user preferences, the volumes of the elements can be conserved individually (default), or the volume can be
     * conserved through the column, but volume exchange is allowed between elements of the same column, this is stated
     * by the boolean Simulation#conservingColumnVolumes and operation is handled by the function Simulation#conserveColumnVolume. \n
     */
    if (conservingColumnVolumes){
    	conserveColumnVolume();
    }
    /**
     * Following potential volume exchange, the remodelling (plastic deformation) is updated. All explicit ECM elements (ShapeBase#isECMMimicing = true)
     * do have remodelling by default. Remodelling of cellular elements (all the rest) can be defiend by the user preferences
     * and is defined in Simulation#thereIsPlasticDeformation. The remodelling function is Simulation#updatePlasticDeformation.
     */
    if(thereIsPlasticDeformation || thereIsExplicitECM){
    	updatePlasticDeformation();
    }
    /**
     * The growth, shape change, volume exchange and remodelling all induce theri effects eventually through growth of the elements. Once
     * all these options are covered, the elemetns are grown in ShapeBase#growShapeByFg, which updated the growth deformation gradient of the element. \n
     */
    for(const auto &itElement: Elements){
        itElement->growShapeByFg();
    }
    /**
    * As the shape size is changed, nodal masses (volumes) (Simulation#updateNodeMasses) and
    * exposed surfaces  (Simulation#updateNodeViscositySurfaces) are updated, then the
    * weight of each elements contribution on the nodes are updated (Simulation#updateElementToConnectedNodes). Then
    * the collapse nodes positions are updated, if the collapse is being carried out in stages. \n
    */
    updateNodeMasses();
    updateNodeViscositySurfaces();
    updateElementToConnectedNodes(Nodes);
    updatePositionsOfNodesCollapsingInStages();
    /**
     * Then packing to all posssible surfaces are detected, while the actual packing forces are calculated in the NR iterations as
     * they are position dependent. This initial check creates a list of potential packing nodes, to check only these through
     * the NR itaration, rather than checking all nodes at all times. Packing is checked against self contact (Simulation#detectPacingNodes)
     * then to enclosing surfaces (Simulation#detectPacingToEnclosingSurfacesNodes) and pipette (Simulation#detectPackingToPipette)
     * if they are implemented. \n
     */
    detectPacingNodes();
    if (encloseTissueBetweenSurfaces){
    	detectPacingToEnclosingSurfacesNodes();
    }
    if (PipetteSuction){
        detectPackingToPipette();
    }
    if (addingRandomForces){
    	calculateRandomForces();
    }
    /**
     * Then collapse of elemental surfaces are checked (Simulation#checkEdgeLenghtsForBindingPotentiallyUnstableElementss),
     * and the NR solver is notified if there are changes made.
     */
    bool thereIsBinding = false;
    if (thereNodeCollapsing){
    	thereIsBinding = checkEdgeLenghtsForBindingPotentiallyUnstableElements();
    }
    if (thereIsBinding){
		NRSolver->boundNodesWithSlaveMasterDefinition = true;
	}
    /**
     * After the collapse check, the nodal adhesion is updated (Simulation#adhereNodes), makinguse of the detected packing, as any node couple close enough
     * to adhere must be close enough to activate self-contact. \n
     */
    if (thereIsAdhesion){
    	thereIsBinding = adhereNodes();
    }
    if (thereIsBinding){
    	NRSolver->boundNodesWithSlaveMasterDefinition = true;
    }
    /**
     * The actual positional updates are delegated to teh newton-Raphson solver object via function Simulation#updateStepNR. \n
     */
    updateStepNR();
    /**
     * Once the positions are updated, the bounding box is calculated again, the elements are checked agains flipping (Simulation#checkFlip).
     */
    calculateBoundingBox();
    //cout<<" step: "<<timestep<<" Pressure: "<<SuctionPressure[2]<<" Pa, maximum z-height: "<<boundingBox[1][2]<<" L/a: "<<(boundingBox[1][2]-50)/(2.0*pipetteRadius)<<std::endl;
    Success = checkFlip();
    /**
     * If there is aftificail relaxation defined by the user, this is done at the end of the time step via Simulation#artificialRelax.
     */
    if (thereIsArtificaialRelaxation && artificialRelaxationTime == currSimTimeSec){
		artificialRelax();
	}
    /**
     * Then the time step is updated, and the data is saved if necessary in Simulation#processDisplayDataAndSave.
     */
    timestep++;
    currSimTimeSec += dt;
    if (Success){
    	processDisplayDataAndSave();
    }
    return Success;
}

void Simulation::assignIfElementsAreInsideEllipseBands(){
    for(const auto& itElement : Elements){
        if ((itElement->tissueType == 0 || itElement->tissueType == 2)){
            //The element is either a columnar or a linker element.
			//Peripodial elements are not counted in the ellipses, as the perturbations of
			// interest are not applied to them			
            itElement->checkIfInsideEllipseBands(nMarkerEllipseRanges, markerEllipseBandXCentres,markerEllipseBandR1Ranges, markerEllipseBandR2Ranges, Nodes);
		}
	}
}

void Simulation::updateRelativePositionsToApicalPositioning(){
    /**
     * All the elements of a single column are stored in ShapeBase#elementsIdsOnSameColumn vector for each element.
     * Going through the apical (top layer) elements of the tissue, the relative position of the apical element
     * is assigned to all the elements on the same column.
     */
    for( const auto& itElement : Elements){
        if ( itElement->tissuePlacement == 1 && itElement->tissueType == 0){//apical element of columnar layer
            std::array<double,2> ReletivePos = itElement->getRelativePosInBoundingBox();
            for (size_t i=0; i < TissueHeightDiscretisationLayers-2;++i){
                int idOfElementBelowThisApicalElement = itElement->elementsIdsOnSameColumn[i];
				Elements[idOfElementBelowThisApicalElement]->setRelativePosInBoundingBox(ReletivePos[0],ReletivePos[1]);
			}
		}
	}
}

void Simulation::checkForPinningPositionsUpdate(){
	if (GridGrowthsPinnedOnInitialMesh){
		for (int i=0; i<nGrowthPinning; ++i){
			if (currSimTimeSec >= growthPinUpdateTime[i] &&  !growthPinUpdateBools[i]){
				updatePinningPositions();
				growthPinUpdateBools[i] = true;
			}
		}
	}
}

void Simulation::updatePinningPositions(){
    for( const auto& itElement : Elements){
        itElement->setInitialRelativePosInBoundingBox();
	}
}

bool Simulation::checkFlip(){
    for(const auto& itElement : Elements){
        if (itElement->isFlipped){
			//there is a flipped element:
            outputFile<<"There is A flipped element: "<<itElement->Id<<std::endl;
            std::cerr<<"There is A flipped element: "<<itElement->Id<<std::endl;
			return false;
		}
	}
	return true;
}
void Simulation::wrapUpAtTheEndOfSimulation(){
    alignTissueDVToXPositive();
}


void Simulation::clearProjectedAreas(){
    for (const auto& currNode : Nodes){
        currNode->zProjectedArea = 0.0;
    }
}

void Simulation::correctzProjectedAreaForMidNodes(){
   for (const auto& currNode : Nodes){
        if (currNode->tissuePlacement == 2 ){ // the node is on midlayer
            /**
             * For the nodes in mid-line nodes the area is added from apical and basal surfaces of elemetns on both sides.
             * This is corrected in this function. !! A more efficient approach would be to calculate these areas only once!!
             */
            currNode->zProjectedArea /= 2.0;
        }
    }
}

void Simulation::calculateZProjectedAreas(){
    clearProjectedAreas();
    for(const auto& itElement : Elements){
        itElement->calculateZProjectedAreas();
        itElement->assignZProjectedAreas(Nodes);
    }
    correctzProjectedAreaForMidNodes();
}

void Simulation::updatePlasticDeformation(){
    //not using range based loops here to ensure openMP comaptibility
    #pragma omp parallel for
    for(std::vector<std::unique_ptr<ShapeBase>>::iterator itElement = Elements.begin(); itElement <Elements.end(); ++itElement){
        if (thereIsExplicitECM &&  (*itElement)->isECMMimicing){
            /**
             * The ezplicitely defined ECM elements are always subject to non-volume econserving plastic deformation.
             * The remodelling is calculated by ShaopeBase#calculatePlasticDeformation3D.
             */
            (*itElement)->calculatePlasticDeformation3D(false,dt,ECMRenawalHalfLife, 0.1, 10.0);
        }
        else if (thereIsPlasticDeformation){
            /**
             * If there is user prefered remodelling on elements other than the ECM then the ShapeBase#tissueType is
             * chacked agains the parameters Simulation#plasticDeformationAppliedToColumnar and
             * Simulation#plasticDeformationAppliedToPeripodial. The parameters of the ShapeBase#calculatePlasticDeformation3D
             * are again defined by the user, such as the Simulation#volumeConservedInPlasticDeformation, Simulation#plasticDeformationHalfLife
             * Simulation#zRemodellingLowerThreshold and Simulation#zRemodellingUpperThreshold. See function documentation for further details.
             */
             if(( ((*itElement)->tissueType == 0 || (*itElement)->tissueType == 2) && plasticDeformationAppliedToColumnar) || ( ((*itElement)->tissueType == 1 || (*itElement)->tissueType == 2) && plasticDeformationAppliedToPeripodial))
            {
                (*itElement)->calculatePlasticDeformation3D(volumeConservedInPlasticDeformation,dt,plasticDeformationHalfLife, zRemodellingLowerThreshold, zRemodellingUpperThreshold);
            }
        }
		else{
            /**
             * If there is no remodelling, then the elements plastic deformation increment is set to identity in ShapeBase#setPlasticDeformationIncrement.
             */
            (*itElement)->setPlasticDeformationIncrement(1.0,1.0,1.0);
		}
	}
}

void Simulation::conserveColumnVolume(){
    //not using range based loops here to ensure openMP comaptibility
    #pragma omp parallel for
    for(std::vector<std::unique_ptr<ShapeBase>>::iterator itElement = Elements.begin(); itElement <Elements.end(); ++itElement){
        std::vector <int> elementIdsForRedistribution;
        if ((*itElement)->tissueType == 0 && (*itElement)->tissuePlacement == 1){
            for (size_t i=0; i < TissueHeightDiscretisationLayers;++i){
                int idOfElementOnSameColumn = (*itElement)->elementsIdsOnSameColumn[i];
    			if (Elements[idOfElementOnSameColumn]->tissueType != 0){
    				//must be columnar
					continue;
				}
    			if (Elements[idOfElementOnSameColumn]->isECMMimicing){
    				//exclude ECM
    				continue;
    			}
    			if (Elements[idOfElementOnSameColumn]->isActinMimicing){
    				//exclude top actin layer if there is one
					continue;
				}
    			elementIdsForRedistribution.push_back(idOfElementOnSameColumn);
    		}
            /**
             * The column-vise volume conservation requires volume redistribution amongst cell layer of the tissue,
             * and a diffusion based volume exchange where the most compressed elements give volume to least compressed.
             */
    		//now I have a list of elements that I would like to redistribute the volumes of:
    		//get sum of ideal volumes:
    		const int n = elementIdsForRedistribution.size();
    		double ratioOfCurrentVolumeToIdeal[n];
    		double sumIdealVolumes = 0;
    		double sumCurrentVolumes = 0;
    		for (int i=0; i < n;++i){
    			double idealVolume = Elements[elementIdsForRedistribution[i]]->GrownVolume;
    			sumIdealVolumes += idealVolume;
    			//get current volume from average detF:
    			double currentVolume = Elements[elementIdsForRedistribution[i]]->getCurrentVolume();
    			if (currentVolume < 10E-10){
    				currentVolume = idealVolume;
    			}
    			ratioOfCurrentVolumeToIdeal[i]= currentVolume/idealVolume;
    			sumCurrentVolumes += currentVolume;
    			if (elementIdsForRedistribution[i] == 7009 || elementIdsForRedistribution[i] == 5207 || elementIdsForRedistribution[i] == 4622 || elementIdsForRedistribution[i] == 2820){
                    std::cout<<" redistributing element: "<<elementIdsForRedistribution[i]<<" idealVolume: "<<idealVolume<<" currentVolume: "<<currentVolume<<" ratio: "<<ratioOfCurrentVolumeToIdeal[i]<<std::endl;
    			}
    		}
    		double redistributionFactors[n];
    		for (int i=0; i < n;++i){
    			//calculate redistribution factors:
    			redistributionFactors[i] = sumIdealVolumes/sumCurrentVolumes *ratioOfCurrentVolumeToIdeal[i];
    			if (elementIdsForRedistribution[i] == 7009 || elementIdsForRedistribution[i] == 5207 || elementIdsForRedistribution[i] == 4622 || elementIdsForRedistribution[i] == 2820){
                    std::cout<<" redistributing element: "<<elementIdsForRedistribution[i]<<" redist factor: "<<redistributionFactors[i]<<std::endl;
    			}
    			//calculateGrowth:
    			double uniformGrowthIncrement = pow(redistributionFactors[i],1.0/3.0);
    			gsl_matrix* columnarFgIncrement = gsl_matrix_calloc(3,3);
    			gsl_matrix* peripodialFgIncrement = gsl_matrix_calloc(3,3);
    			gsl_matrix_set_identity(columnarFgIncrement);
    			gsl_matrix_set_identity(peripodialFgIncrement);
    			gsl_matrix_set(columnarFgIncrement,0,0,uniformGrowthIncrement);
    			gsl_matrix_set(columnarFgIncrement,1,1,uniformGrowthIncrement);
    			gsl_matrix_set(columnarFgIncrement,2,2,uniformGrowthIncrement);
    			Elements[elementIdsForRedistribution[i]]->updateGrowthIncrement(columnarFgIncrement,peripodialFgIncrement);
    			gsl_matrix_free(columnarFgIncrement);
    			gsl_matrix_free(peripodialFgIncrement);
    		}
    	}
    }
}

void Simulation::updateStepNR(){
    int iteratorK = 0;
    /**
     * The iteration will be carried out for 20 steps, if upon 20 trials the simulation have not converged, an error will be
     * generated. The iteration starts by clearing all force, displacement and Jacobian matrices of the NR solver in
     * NewtonRaphsonSolver#setMatricesToZeroAtTheBeginningOfIteration. Then the vector containing the positions of teh nodes at the
     * end of the previous time step "n", \f$ \boldsymbol{u_n} \f$ is calculated in NewtonRaphsonSolver#constructUnMatrix. The positions
     * of the current iteration "k", \f$ \boldsymbol{u_k} \f$ are initiated equal to \f$ \boldsymbol{u_n} \f$ via
     * NewtonRaphsonSolver#initialteUkMatrix. \n
     */
    int maxIteration =20;
    bool converged = false;

    NRSolver->setMatricesToZeroAtTheBeginningOfIteration();
    NRSolver->constructUnMatrix(Nodes);
    NRSolver->initialteUkMatrix();
    while (!converged){
        /**
         * While the system has not converged, all the forces are rest at the beginning of each iteration with Simulation#resetForces.
         * Then the matrices to be cumulated from scratch in each iteration are reset in NewtonRaphsonSolver#setMatricesToZeroInsideIteration().
         * The displacement matrix per time step Simulation#dt is calculated for use in external viscous forces in
         * NewtonRaphsonSolver#calculateDisplacementMatrix. Then all the internal elemental and nodal forces are calculated in NewtonRaphsonSolver#calculateForcesAndJacobianMatrixNR.
         * The forces are moved from the elements and nodes to the system vectors in NewtonRaphsonSolver#writeForcesTogeAndgvInternal, the elastic and viscous
         * terms of the elemental Jacobians are mapped and added onto the system Jacobian in  NewtonRaphsonSolver#writeImplicitElementalKToJacobian.
         * The external forces are calculated in NewtonRaphsonSolver#calculateExternalViscousForcesForNR and their derivatives are added to
         * the system Jacobian in NewtonRaphsonSolver#addImplicitKViscousExternalToJacobian. \n
         */
        std::cout<<"iteration: "<<iteratorK<<std::endl;
        if (implicitPacking){
            resetForces(true);	// reset packing forces
        }
        else{
            resetForces(false);	// do not reset packing forces
        }
        NRSolver->setMatricesToZeroInsideIteration();
        NRSolver->calculateDisplacementMatrix(dt);
        NRSolver->calculateForcesAndJacobianMatrixNR(Nodes, Elements, dt);
	    //Writing elastic Forces and elastic Ke:
		NRSolver->writeForcesTogeAndgvInternal(Nodes, Elements, SystemForces);
	    NRSolver->writeImplicitElementalKToJacobian(Elements);
        NRSolver->calculateExternalViscousForcesForNR(Nodes);
	    NRSolver->addImplicitKViscousExternalToJacobian(Nodes,dt);
        /**
         * If the packing forces calculated implicitely (default) they are updated via Simulation#calculatePackingForcesImplicit3D and
         * the corresponding Jacobian update is carried out with Simulation#calculatePackingJacobian3D. The packing for enclosing surfaces and the pipette
         * are also carried out here (Simulation#calculatePackingForcesToEnclosingSurfacesImplicit3D,
         * Simulation#calculatePackingToEnclosingSurfacesJacobian3D, Simulation#calculatePackingToPipetteForcesImplicit3D
         * and Simulation#calculatePackingToPipetteJacobian3D).
         */
        if (implicitPacking){
            calculatePackingForcesImplicit3D();
            calculatePackingJacobian3D(NRSolver->K);
	        if (encloseTissueBetweenSurfaces){
	        	calculatePackingForcesToEnclosingSurfacesImplicit3D();
	        	calculatePackingToEnclosingSurfacesJacobian3D(NRSolver->K);
	        }
            if (PipetteSuction){
                calculatePackingToPipetteForcesImplicit3D();
                calculatePackingToPipetteJacobian3D(NRSolver->K);
            }
	    }
        /**
         * All the internal forces and the external viscous resistance forces are collated in NewtonRaphsonSolver#gSum
         * via NewtonRaphsonSolver#calculateSumOfInternalForces. The external forceas from the pipette suction if set up, all packing,
         * and random forces if assigned are collated in NewtonRaphsonSolver#gExt, and these are added to system forces in
         * NewtonRaphsonSolver#addExernalForces(). \n
         */
        NRSolver->calculateSumOfInternalForces();
        if (PipetteSuction && timestep >= PipetteInitialStep){
			calculateZProjectedAreas();
            addPipetteForces(NRSolver->gExt);
		}
        //packing can come from both encapsulation and tissue-tissue packing. I add the forces irrespective of adhesion.
		addPackingForces(NRSolver->gExt);

		if (addingRandomForces){
			addRandomForces(NRSolver->gExt);
		}
        NRSolver->addExernalForces();
        checkForExperimentalSetupsWithinIteration();
        /**
         * Once all forces and their derivatives are collated in system forces and Jacobian, then degrees of freedom fixing is
         * reflected on these matrices in  NewtonRaphsonSolver#calcutateFixedK and NewtonRaphsonSolver#calculateBoundKWithSlavesMasterDoF.
         */
    	NRSolver->calcutateFixedK(Nodes);
        NRSolver->calculateBoundKWithSlavesMasterDoF();
        /** Then NewtonRaphsonSolver#solveForDeltaU function arranges the matrices and solves for the incremental displacements via
         * PARDISO sparse matrix solver. The convergence is checked via the norm of incremental displacements in
         * NewtonRaphsonSolver#checkConvergenceViaDeltaU, the position of the nodes in the iteration, \f $ \boldsymbol{u_k} \f $ are updated
         * in NewtonRaphsonSolver#updateUkInIteration. The nodal and elemental positions are updated (Simulation#updateElementPositionsinNR,
         * Simulation#updateNodePositionsNR) with the new node positions of the iteration.
         */
        NRSolver->solveForDeltaU();
        converged = NRSolver->checkConvergenceViaDeltaU();
        NRSolver->updateUkInIteration();
        updateElementPositionsinNR(NRSolver->uk);
        updateNodePositionsNR(NRSolver->uk);
        iteratorK ++;
        if (!converged && iteratorK > maxIteration){
            std::cerr<<"Error: did not converge!!!"<<std::endl;
            converged = true;
        }
    }
    checkForExperimentalSetupsAfterIteration();
    //Now the calculation is converged, I update the node positions with the latest positions uk:
    updateNodePositionsNR(NRSolver->uk);
     //Element positions are already up to date.
    std::cout<<"finished run one step"<<std::endl;
    if (PipetteSuction){
    	//find the max z:
    	double zMax = -10000;
    	int idMax = -10;
        for (const auto& itNode : Nodes){
            if (itNode->Position[2] > zMax){
                zMax = itNode->Position[2];
                idMax = itNode->Id;
    		}
    	}
        std::cout<<"Pipette suction: "<<SuctionPressure[2]<<" max suction: "<<zMax<<" from node "<<idMax<<std::endl;
    }
}

void Simulation::calculateRandomForces(){
	randomForces.clear();
	randomForces=RandomGenerator::Obj().getNormRV( randomForceMean,randomForceVar, 3*nNodes );
	//making the sum of forces zero:
	double sumRandomForceX = 0.0;
	double sumRandomForceY = 0.0;
	double sumRandomForceZ = 0.0;

	vector<double>::iterator itDouble;
	for (itDouble =randomForces.begin(); itDouble < randomForces.end(); itDouble = itDouble+3){
		sumRandomForceX+= (*itDouble);
		sumRandomForceY+= (*(itDouble+1));
		sumRandomForceZ+= (*(itDouble+2));
	}
	sumRandomForceX /= nNodes;
	sumRandomForceY /= nNodes;
	sumRandomForceZ /= nNodes;

	for (itDouble =randomForces.begin(); itDouble < randomForces.end(); itDouble = itDouble+3){
		(*itDouble) -= sumRandomForceX;
		(*(itDouble+1)) -= sumRandomForceY;
		(*(itDouble+2)) -= sumRandomForceZ;
	}
}

void Simulation::addRandomForces(gsl_matrix* gExt){
        for (size_t j=0; j<3*nNodes; ++j){
			double F = randomForces[j];
			F += gsl_matrix_get(gExt,j,0);
			gsl_matrix_set(gExt,j,0,F);
		}
}

void Simulation::updateNodePositionsNR(gsl_matrix* uk){
    /**
     * The nodal positions vector contains each position of each node in a single vector, and indexing is carried out accordingly.
     */
    int dim = 3;
    for (size_t i = 0; i<nNodes; ++i){
    	for (int j=0; j<dim; ++j){
            Nodes[i]->Position[j]=gsl_matrix_get(uk,dim*i+j,0);
        }
    }
}

void Simulation::updateElementPositionsinNR(gsl_matrix* uk){
    int dim = 3;
    for(const auto& itElement : Elements){
        const vector<int>& nodeIds = itElement->getNodeIds();
        size_t nNodes= itElement->getNodeNumber();
        for (size_t j=0; j<nNodes; ++j){
            double x = gsl_matrix_get(uk,dim*nodeIds[j],0);
            double y = gsl_matrix_get(uk,dim*nodeIds[j]+1,0);
            double z = gsl_matrix_get(uk,dim*nodeIds[j]+2,0);
            itElement->Positions[j][0] = x;
            itElement->Positions[j][1] = y;
            itElement->Positions[j][2] = z;
        }
    }
}

void Simulation::calculatePackingForcesToEnclosingSurfacesImplicit3D(){
    /**
     * See the documentation for Simulation#calculatePackingForcesImplicit3D. The methodology is the same,
     * with rigid wall positions in two z coordinates instead of a second node position.
     */
    size_t nPositive=nodesPackingToPositiveSurface.size();
    size_t nNegative=nodesPackingToNegativeSurface.size();
	double multiplier = packingMultiplier;
	double sigmoidSaturation = sigmoidSaturationForPacking;
    for(size_t i = 0 ; i<nNegative; ++i){
		int id0 = nodesPackingToNegativeSurface[i];
		double dz = Nodes[id0]->Position[2] - zEnclosementBoundaries[0];
		if (initialWeightPackingToNegativeSurface[i]>0){
			dz *= -1.0;
		}
		double mass = Nodes[id0]->mass;
		double Fz = multiplier * mass / (1 + exp(sigmoidSaturation / packingToEnclosingSurfacesThreshold * (-1.0*dz)));
		Fz *= initialWeightPackingToNegativeSurface[i];
		PackingForces[id0][2] += Fz;
	}
    for(size_t i = 0 ; i<nPositive; ++i){
		int id0 = nodesPackingToPositiveSurface[i];
		double dz = Nodes[id0]->Position[2] - zEnclosementBoundaries[1];
		if (initialWeightPackingToPositiveSurface[i]>0){
			dz *= -1.0;
		}
		double mass = Nodes[id0]->mass;
		double Fz = multiplier * mass / (1 + exp(sigmoidSaturation / packingToEnclosingSurfacesThreshold * (-1.0*dz)));
		Fz *= initialWeightPackingToPositiveSurface[i];
		PackingForces[id0][2] += Fz;
	}
}

void Simulation::calculatePackingToEnclosingSurfacesJacobian3D(gsl_matrix* K){
    /**
     * See the documentation for Simulation#calculatePackingJacobian3D. The methodology is the same,
     * with rigid wall positions in two z coordinates instead of a second node position.
     */
    size_t nPositive=nodesPackingToPositiveSurface.size();
    size_t nNegative=nodesPackingToNegativeSurface.size();
	double sigmoidSaturation = sigmoidSaturationForPacking;
	double multiplier = packingMultiplier;
    //not using range based loops here to ensure openMP comaptibility
    #pragma omp parallel for
    for(size_t i = 0 ; i<nPositive; ++i){
		int id0 =  nodesPackingToPositiveSurface[i];
		//sigmoid test:
		double dz = Nodes[id0]->Position[2] - zEnclosementBoundaries[1];
		double mass = Nodes[id0]->mass;
		if (initialWeightPackingToPositiveSurface[i]>0){
			dz *= -1.0;
		}
		double sigmoidz =  1 / (1 + exp(sigmoidSaturation/ packingToEnclosingSurfacesThreshold * (-1.0 * dz) ));
		double dFzdz0 = sigmoidz * (1 - sigmoidz) * multiplier * mass * initialWeightPackingToPositiveSurface[i] * (sigmoidSaturation/packingToEnclosingSurfacesThreshold);
		if (initialWeightPackingToPositiveSurface[i]>0){
			dFzdz0 *= -1.0;
		}
		//z values:
        addValueToMatrix(K,3*id0+2,3*id0+2,-dFzdz0);
	}
    //not using range based loops here to ensure openMP comaptibility
    #pragma omp parallel for
    for(size_t i = 0 ; i<nNegative; ++i){
		int id0 =  nodesPackingToNegativeSurface[i];
		//sigmoid test:
		double dz = Nodes[id0]->Position[2] - zEnclosementBoundaries[0];
		double mass = Nodes[id0]->mass;
		if (initialWeightPackingToNegativeSurface[i]>0){
			dz *= -1.0;
		}
		double sigmoidz =  1 / (1 + exp(sigmoidSaturation/ packingToEnclosingSurfacesThreshold * (-1.0 * dz) ));
		double dFzdz0 = sigmoidz * (1 - sigmoidz) * multiplier * mass * initialWeightPackingToNegativeSurface[i] * (sigmoidSaturation/packingToEnclosingSurfacesThreshold);
		if (initialWeightPackingToNegativeSurface[i]>0){
			dFzdz0 *= -1.0;
		}
		//z values:
        addValueToMatrix(K,3*id0+2,3*id0+2,-dFzdz0);
	}
}

bool Simulation::areNodesOnNeighbouingElements(int masterNoeId, int slaveNodeId){
	bool neigElements = false;
    size_t nOwnersMaster = Nodes[masterNoeId]->connectedElementIds.size();
    size_t nOwnersSlave = Nodes[slaveNodeId]->connectedElementIds.size();
    for (size_t i =0; i<nOwnersMaster; ++i){
        int idOwnerMaster = Nodes[masterNoeId]->connectedElementIds[i];
        const vector<int>&  nodeIdsMasterOwner = Elements[idOwnerMaster]->getNodeIds();
        for (size_t j =0; j<nOwnersSlave; ++j){
			int idOwnerSlave = Nodes[slaveNodeId]->connectedElementIds[j];
			int counter = 0;
            const vector<int> nodeIdsSlaveOwner = Elements[idOwnerSlave]->getNodeIds();
            for (auto currMasterNodeId : nodeIdsMasterOwner ){
                for (auto currSlaveNodeId : nodeIdsSlaveOwner){
                    if (currMasterNodeId == currSlaveNodeId){
						counter++;
					}
					if (counter>2){
						neigElements = true;
						return neigElements;
					}
				}
			}
		}
	}
	return neigElements;
}

bool Simulation::isAdhesionAllowed(int masterNodeId, int slaveNodeId){
    /**
     * Two nodes are being attempted to adhere. Certain criteria must be met before operation
     * can proceed. \n
     * If any of the nodes is already in the process of being moved due to an adhesion collapse
     * a new adhesion cannot be formed at this step.
     */
	if(Nodes[slaveNodeId]->positionUpdateOngoing || Nodes[masterNodeId]->positionUpdateOngoing){
		return false;
	}
    /**
     * The nodes must be of the same tissue type, peripodial elements cannot bind to columner elements.
     */
	if(Nodes[slaveNodeId]->tissueType != Nodes[masterNodeId]->tissueType){
        return false;
	}
    /**
     * The nodes at circumference can only adhere with each other, there are biological resoans with the structure of the matrix for this.
     */
    if(  (Nodes[slaveNodeId]->atCircumference && !Nodes[masterNodeId]->atCircumference)
		||(!Nodes[slaveNodeId]->atCircumference && Nodes[masterNodeId]->atCircumference) )
	{
		return false;
	}
    /**
     * If the nodes are collapsed on adhesion, then the adhesion can only be one to one, and an alrady adhered node
     * cannot adhere a another node. The default adhesion id is "-1" indicating no adhesion.
     */
    if(collapseNodesOnAdhesion){
		if (Nodes[slaveNodeId]->adheredTo > -1 || Nodes[masterNodeId]->adheredTo > -1) {
            //std::cout<<"one is already adhered"<<std::endl;
			return false;
		}
	}
    /**
     * If this couple have already been adhered, no need to continue the operation. This is done by checking the
     * slave node id in the list of collapsed nodes for master id (Node#collapsedWith)
     */
    if (binary_search(Nodes[masterNodeId]->collapsedWith.begin(), Nodes[masterNodeId]->collapsedWith.end(),slaveNodeId)){
		//the couple is already collapsed:
        return false;
	}
    /**
     * If this node couple belong to neightbouring elements, do not collapse.
     */
    bool nodeIsOnNeigElement = areNodesOnNeighbouingElements(masterNodeId,slaveNodeId);
	if (nodeIsOnNeigElement){
		return false;
	}
	bool collapedNodeIsNeig = false;
	collapedNodeIsNeig = Nodes[masterNodeId]->isNeigWithMyCollapsedNodes(slaveNodeId,Nodes);
	if (collapedNodeIsNeig){
        //std::cout<<"node are collapsed with neighbours (master)"<<std::endl;
		return false;
	}
	collapedNodeIsNeig = Nodes[slaveNodeId]->isNeigWithMyCollapsedNodes(masterNodeId,Nodes);
	if (collapedNodeIsNeig){
        //std::cout<<"node are collapsed with neighbours (slave)"<<std::endl;
		return false;
	}
	return true;
}

double Simulation::distanceSqBetweenNodes(int id0, int id1){
	double dx = Nodes[id0] ->Position[0] - Nodes[id1] ->Position[0];
	double dy = Nodes[id0] ->Position[1] - Nodes[id1] ->Position[1];
	double dz = Nodes[id0] ->Position[2] - Nodes[id1] ->Position[2];
	double dSq = dx*dx + dy*dy + dz*dz;
	return dSq;
}

bool Simulation::checkForElementFlippingUponNodeCollapse(std::vector<int> &newCollapseList, std::array<double,3> avrPos){
    /**
     * If the collapse of two nodes will move the nodes (there is collapse) the owner elemetns are first checked to see if
     * this movement will cause element flipping. The adhesion is avoided if it will.
     */
	int nCollapsedList = newCollapseList.size();
	bool elementWillCollapse = false;
	for (int nodeIterator=0; nodeIterator<nCollapsedList;++nodeIterator){
		int currNodeId = newCollapseList[nodeIterator];
		int nElements=  Nodes[currNodeId]->connectedElementIds.size();
		for (int i=0; i<nElements;++i){
            //std::cout<<"checking element : "<< Nodes[currNodeId]->connectedElementIds[i]<<" from node : "<< currNodeId<<std::endl;
			elementWillCollapse = Elements[Nodes[currNodeId]->connectedElementIds[i]]->isElementFlippedInPotentialNewShape(currNodeId, avrPos[0], avrPos[1], avrPos[2]);
			if (elementWillCollapse){
                std::cout<<" element "<<Nodes[currNodeId]->connectedElementIds[i]<<" will flip if adhered"<<std::endl;
				return false;
			}
		}
	}
	return true;
}


void Simulation::updatePositionsOfNodesCollapsingInStages(){
    for(const auto& itNode : Nodes){
        if (itNode->positionUpdateOngoing){
            std::array<double,3> avrPos = {0.0, 0.0, 0.0};
            std::array<bool,3> fix = {false, false, false};

            int n = itNode->collapsedWith.size();
            for (auto collapsedNodeId : itNode->collapsedWith){
                for (size_t j=0; j<3; ++j){
                    avrPos[j] += Nodes[collapsedNodeId]->Position[j]/n;
				}
			}
            for (auto collapsedNodeId : itNode->collapsedWith){
                for (size_t j=0; j<3; ++j){
                    if(Nodes[collapsedNodeId]->FixedPos[j]){
                        avrPos[j] = Nodes[collapsedNodeId]->Position[j];
						fix[j] = true;
					}
				}
			}
            itNode->updatePositionTowardsPoint(avrPos,fix);
		}

	}
}

bool Simulation::adhereNodes(){
    for(const auto& itNode : Nodes){
         itNode->clearDuplicatesFromCollapseList();
	}
    size_t n = pacingNodeCouples0.size();
	int dim = 3;
	bool thereIsBinding = false;
	vector<bool> adhesionAllowedList;
    /**
     * The node couples idedntified as potentially packing to each other Simulation#pacingNodeCouples0 and Simulation#pacingNodeCouples1
     * are checked for adhesion. First if the adhesion between the pair is feasible is checked via Simulation#isAdhesionAllowed and recorded in a
     * vector (adhesionAllowedList).
     */
    for(size_t nodeCoupleIterator = 0 ; nodeCoupleIterator<n; ++nodeCoupleIterator){
		bool adhesionAllowed = isAdhesionAllowed(pacingNodeCouples0[nodeCoupleIterator], pacingNodeCouples1[nodeCoupleIterator]);
		adhesionAllowedList.push_back(adhesionAllowed);
	}
    for(size_t nodeCoupleIterator = 0 ; nodeCoupleIterator<n; ++nodeCoupleIterator){
		if (!adhesionAllowedList[nodeCoupleIterator]){
			continue;
		}
		int masterNodeId = pacingNodeCouples0[nodeCoupleIterator];
		int slaveNodeId = pacingNodeCouples1[nodeCoupleIterator];
		if (collapseNodesOnAdhesion){
            /**
             * If the adhesion collapses nodes, then the adhesion is one-to-one, and should be carried out between the closest couples.
             */
            for (size_t nodeCoupleIteratorFromHereOn = nodeCoupleIterator+1; nodeCoupleIteratorFromHereOn<n; ++nodeCoupleIteratorFromHereOn){
				if (masterNodeId == pacingNodeCouples0[nodeCoupleIteratorFromHereOn] ||
					masterNodeId == pacingNodeCouples1[nodeCoupleIteratorFromHereOn] ||
					slaveNodeId  == pacingNodeCouples0[nodeCoupleIteratorFromHereOn] ||
					slaveNodeId  == pacingNodeCouples1[nodeCoupleIteratorFromHereOn]
					){
					//one of the node couples occurs somewhere else, is this occurance viable?
					if (adhesionAllowedList[nodeCoupleIteratorFromHereOn]){
                        /**
                         * If the adhesion is allowed for this pair, then all the coming pairs are checked to see if there is a pair involving these two
                         * nodes that is closer then this one.
                         */
						//yes this adhesion is possible, compare distances:
						double dSqOriginal = distanceSqBetweenNodes(masterNodeId, slaveNodeId);
						double dSqNext = distanceSqBetweenNodes(pacingNodeCouples0[nodeCoupleIteratorFromHereOn],pacingNodeCouples1[nodeCoupleIteratorFromHereOn]);
						if (dSqNext < dSqOriginal){
							//the next couple I will reach is closer, adhere them (or another one that may be further down the list and even closer)
                            /**
                             * If a consequent pair I will reach is closer, they should be adhered (or they will also be surpassed by yet another
                             * closer pair). Declare this pari cannot adhere, and continue.
                             */
                            adhesionAllowedList[nodeCoupleIterator] = false;
							break;
						}
						else{
                            /**
                             * Alternatively, if a consequent pair is further away, then they should be flagged as not feasible to avoid
                             * re-calculation. The current pair will be adhered and there will be no need to check the following one. \n
                             */
                            adhesionAllowedList[nodeCoupleIteratorFromHereOn] = false;
                        }
					}
				}
			}
			if (!adhesionAllowedList[nodeCoupleIterator]){
				continue;
			}
		}
        /**
         * If the current pair is still adhering, then the average position in the middle will be calculated and the nodes will be collapsed.
         */
        if (collapseNodesOnAdhesion){
			vector<int> newCollapseList;
            std::array <double,3> avrPos = {0.0, 0.0, 0.0};
            std::array <bool,3> fix = {false, false, false};
			Nodes[slaveNodeId]->getNewCollapseListAndAveragePos(newCollapseList, avrPos, fix,Nodes, masterNodeId);
			Nodes[slaveNodeId]->adheredTo = masterNodeId;
			Nodes[masterNodeId]->adheredTo = slaveNodeId;
			if (collapseNodesOnAdhesion){
                    Nodes[slaveNodeId]->collapseOnNodeInStages(newCollapseList, avrPos, fix, Nodes);
			}
		}
		else{
			//Adhering the nodes
			Nodes[slaveNodeId]->adheredTo = masterNodeId;
			Nodes[masterNodeId]->adheredTo = slaveNodeId;
		}
        /**
         * The pair will be flagged as adhered in Simulation#pacingNodeCouplesHaveAdhered, and packing forces will not be calculated for them.
         */
		pacingNodeCouplesHaveAdhered[nodeCoupleIterator] = true;
        /**
         * If any of the pair have fixed dimensions, this will be reflected on the other.
         */
        for (size_t i =0 ; i<3; ++i){
			//if the dimension is fixed in space, fix the other and move on.
			if (Nodes[masterNodeId]->FixedPos[i]){
				Nodes[slaveNodeId]->FixedPos[i]=true;
			}
			else if (Nodes[slaveNodeId]->FixedPos[i]){
				Nodes[masterNodeId]->FixedPos[i]=true;
			}
			if (!Nodes[slaveNodeId]->FixedPos[i]){
				int dofmaster = masterNodeId*dim+i;
                int dofslave  = slaveNodeId*dim+i;
                /**
                 * If the slave is alrady slave to another node (Node#slaveTo is not set to -1, the default non-slave indice), then the master of this processed couple is set to
                 * be the slave of the existing master of the slave, and all three degree of freedoms are bound to each other.
                 */
				if (Nodes[slaveNodeId]->slaveTo[i] > -1){
                    /**
                     * One possiblity is that the potential slave is already a slave to the potential master, then nothing is needed to be done.
                     */
                    if(Nodes[slaveNodeId]->slaveTo[i] == masterNodeId || Nodes[slaveNodeId]->slaveTo[i] == Nodes[masterNodeId]->slaveTo[i] ){
						pacingNodeCouplesHaveAdhered[nodeCoupleIterator] = true;
						continue;
					}
					dofslave = Nodes[slaveNodeId]->slaveTo[i]*dim+i;
					slaveNodeId = Nodes[slaveNodeId]->slaveTo[i];
				}
                /**
                 * If the potential master node is already a slave to another node, then the potential master is moved to
                 * the original master of the potential master DoF in function NewtonRaphsonSolver#checkMasterUpdate.
                 */
                NRSolver->checkMasterUpdate(dofmaster,masterNodeId);
                /**
                 * In this last check point, if the original couple was flipped such that the potential master vas the slave of the potential slave,
                 * then now both master and slave nodes are equal to the initial potential slave id ( as it has been the master prior to this step).
                 * Then again, nothing needs to be implemented, slave master coupling is already in place.
                 */
                if (dofmaster != dofslave){
					bool continueAddition =  NRSolver->checkIfCombinationExists(dofslave,dofmaster);
					if (continueAddition){
                        /**
                         * If the couple is not implemented, then the check for the status of the slave is carried out, if the slave is already master of
                         * other nodes, the coupling is checked and corrected in function NewtonRaphsonSolver#checkIfSlaveIsAlreadyMasterOfOthers. If the
                         * function updates the node couple, then the potential slave was a master, and now is a slave to the potential master.
                         * All slaves of the potential master are moved on to the potential master. \n
                         */
						bool madeChange = NRSolver->checkIfSlaveIsAlreadyMasterOfOthers(dofslave,dofmaster);
						if (madeChange){
                            for (size_t nodeIt = 0 ; nodeIt<nNodes; ++nodeIt){
								if(Nodes[nodeIt]->slaveTo[i]==slaveNodeId){
									Nodes[nodeIt]->slaveTo[i]=masterNodeId;
								}
							}
						}
						vector <int> fixDOF;
						fixDOF.push_back(dofslave);
						fixDOF.push_back(dofmaster);
						NRSolver->slaveMasterList.push_back(fixDOF);
						Nodes[slaveNodeId]->slaveTo[i] = masterNodeId;
						Nodes[masterNodeId]->isMaster[i] = true;
                        //std::cout<<" adhereing nodes: "<<masterNodeId<<" "<<slaveNodeId<<" in dof "<<i<<std::endl;
						thereIsBinding = true;
						pacingNodeCouplesHaveAdhered[nodeCoupleIterator] = true;

					}
				}
			}
		}
	}
	return thereIsBinding;
}

void Simulation::calculatePackingForcesImplicit3D(){
    /**
     * This function calculates the packing forces between nodes, inside the Newton-Raphson iteration steps.
     * The list of nodes that can potentially pack are detected in function , and recorded in Simulation#pacingNodeCouples0 and Simulation#pacingNodeCouples1.
     * The applied packing force is a function of the distance between nodes, and is calculated
     * with an inverse logic function:
     *
      \f[ f\left( d \right) = \frac{L}{1+ e^{-k\left(d-d_{0}\right)}}
        \f]
     *
     * Here, the amplitude \f$ L \f$ is defined such that the force will scale with the average mass of the two nodes.
     * The steepness of the curve is the force profile will approach to zero as the distance between nodes approaches
     * to packing threshold. It is defined the sigmoid saturation term and the packing threshold as detected by the
     * current average side lengths of mesh elements as below. The sigmoid saturation is set to 5, as this is the approximate
     * saturation distance of the standard logistic function.
     *
     *  \f[ -k = \frac{2\:sigmoid\:saturation}{packing\:threshold}
        \f]
     *
     * The distance is shifted with distance \f$ d_{0} \f$ to move the mid point of the function to approximately 60 per cent of
     * the packing threshold distance.
     *
     * \image html packingForce.png
     *
     * Then the forces on each node i and j become:
     * \f[
     *   \mathbf{F_{i}}\left( d \right) = f\left( d \right)\: \mathbf{e_{i}}\:\:,\:\: F_{j}\left( d \right) = f\left( d \right)\: \mathbf{e_{j}}=-f\left( d \right)\: \mathbf{e_{i}}=-\mathbf{F_{i}}
     *   \f]
     *
     *   where the distance \f$ d \f$ is \f$ ||\mathbf{x_{i}}-\mathbf{x_{j}} || \f$, and the normal is \f$  \mathbf{e_{i}} = \left( \mathbf{x_{i}}-\mathbf{x_{j}} \right) / ||\mathbf{x_{i}}-\mathbf{x_{j}} || \f$.
     *
     * Procedure:
     */

    double zeroThreshold = 1E-3;
    /**
     * - Go through all the node couples that are selected to be packing in function
     */
    size_t n = pacingNodeCouples0.size();
    for(size_t node_counter = 0 ; node_counter<n; ++node_counter){
        if (pacingNodeCouplesHaveAdhered[node_counter]){
            continue;
        }
        int id0 = pacingNodeCouples0[node_counter];
        int id1 = pacingNodeCouples1[node_counter];
        double multiplier = packingMultiplier;
        double sigmoidSaturation = sigmoidSaturationForPacking;
        double distanceShiftScale  = 0.6*packingThreshold;
        /**
         * - For each node pair id0 and id1, calculate the distance  \f$ d \f$ and the unit normal between them
         */
        double dx = Nodes[id0]->Position[0] - Nodes[id1]->Position[0];
        double dy = Nodes[id0]->Position[1] - Nodes[id1]->Position[1];
        double dz = Nodes[id0]->Position[2] - Nodes[id1]->Position[2];
        double d = pow((dx*dx + dy*dy + dz*dz),0.5);
        double averageMass = 0.5 *( Nodes[id0]->mass + Nodes[id1]->mass );
        double normal[3]= {1.0,0.0,0.0};
        /**
         * - If the distance between the node pair is zero, then set the distance to the zero threshold value set in the function, to avoid division by zero.
         */
        if (d > zeroThreshold){
            normal[0] = dx/d;
            normal[1] = dy/d;
            normal[2] = dz/d;
        }
        else{
            d = zeroThreshold;
        }
        /**
         * - Shift the distance by \f$ d_{0} \f$, such that the sigmoid is 0.5 at this distance.
         * Currently \f$ d_{0} \f$ is set to 60 percent of Simulation#packingThreshold.
         * A shift of minimum of 50 per cent is necessary, as it will saturate the
         * sigmoid at distance zero. A higher shift will make the saturation earlier,
         * keeping a distance between nodes.
         */
        double shifted_d = d - distanceShiftScale; //shift the distance by the selected percentage,
        /**
         * - calculate the sigmoid value
         */
        double sigmoid = 1 / (1 + exp(sigmoidSaturation/ packingThreshold * 2.0 * shifted_d));
        /**
         * - Scale with average mass to obtain force per node. There is an additional multiplier that can be used to
         * scale the force under specific conditions if desired.
         * - Assign the forces in x, y and z  directions with the normal
         */
        double Fmagnitude = multiplier * averageMass * sigmoid;
        double Fx = Fmagnitude*normal[0];
        double Fy = Fmagnitude*normal[1];
        double Fz = Fmagnitude*normal[2];


        bool printForces = false;
        if (printForces){
            std::cout<<" id0-id1: "<<id0<<"-"<<id1<<" F : "<<Fx<<" "<<Fy<<" "<<Fz<<" d: "<<d<<" packingThreshold: "<<packingThreshold<<" sigmoid: "<<sigmoid;
            std::cout<<" node0pos: "<<Nodes[id0]->Position[0]<<" "<<Nodes[id0]->Position[1]<<" "<<Nodes[id0]->Position[2];
            std::cout<<" node1pos: "<<Nodes[id1]->Position[0]<<" "<<Nodes[id1]->Position[1]<<" "<<Nodes[id1]->Position[2]<<std::endl;
        }
        /**
         * - With the direction of the normal calculation, the algorithm will calculate \f$  \mathbf{F_{0}} \f$.
         * 	\f$  \mathbf{F_{1}} \f$ is in the opposite direction. The forces are recorded on the Simulation#PackingForces vector.
         */
        PackingForces[id0][0] += Fx;
        PackingForces[id0][1] += Fy;
        PackingForces[id0][2] += Fz;
        PackingForces[id1][0] -= Fx;
        PackingForces[id1][1] -= Fy;
        PackingForces[id1][2] -= Fz;
        if (std::isnan(Fx)){
              std::cout<<" packing force Fx is nan for nodes "<<pacingNodeCouples0[node_counter]<<" - "<<pacingNodeCouples1[node_counter]<<std::endl;
        }
        if (std::isnan(Fy)){
              std::cout<<" packing force Fy is nan for nodes "<<pacingNodeCouples0[node_counter]<<" - "<<pacingNodeCouples1[node_counter]<<std::endl;
        }
        if (std::isnan(Fz)){
              std::cout<<" packing force Fz is nan for nodes "<<pacingNodeCouples0[node_counter]<<" - "<<pacingNodeCouples1[node_counter]<<std::endl;
        }
    }
}

void Simulation::calculatePackingJacobian3D(gsl_matrix* K){
    /**
     * This function calculates the derivatives of packing forces between nodes with respect to nodal positions, and fills in hte Jacobian
     * inside the Newton-Raphson iteration steps. The list of nodes that can potentially
     * pack are detected in function , and recorded in and recorded in Simulation#pacingNodeCouples0 and Simulation#pacingNodeCouples1..
     * The applied packing force is calculated in function Simulation#calculatePackingForcesImplicit3D.
     *
     * The derivatives are calculated as
     * \f[
     *   \frac{d\mathbf{F_{i}}}{d\mathbf{x_{i}}} =
     *   		\frac{f\left( d \right)}{d}\left( \mathbf{I} -  \mathbf{e_{i}} \mathbf{e_{i}^{T}}\right)
     *   		+\frac{df\left( d \right)}{dd} \mathbf{e_{i}} \mathbf{e_{i}^{T}}
     *   \f]
     * where distance \f$ d \f$ is \f$ ||\mathbf{x_{i}}-\mathbf{x_{j}} || \f$,
     * the normal is \f$  \mathbf{e_{i}} = \left( \mathbf{x_{i}}-\mathbf{x_{j}} \right) / ||\mathbf{x_{i}}-\mathbf{x_{j}} || \f$,
     * and  \f$  \mathbf{I} \f$ is the identity matrix. The function \f$ f\left( d \right) \f$, linking the distance between the nodes is an inverse sigmoid function,
     * detailed in Simulation#calculatePackingForcesImplicit3D. Its derivative is then:
     * \f[
     * \frac{df\left( d \right)}{dd} =
     * mass\:\frac{-2\:sigmoid\:saturation}{packing\:threshold} f\left( d \right) \left( 1 - f\left( d \right) \right)
     * \f]
     * Procedure:
     */
    double zeroThreshold = 1E-3;
    size_t n = pacingNodeCouples0.size();
    /**
     * - Go through all the node couples that are selected to be packing in function
     */
    for(size_t node_counter = 0 ; node_counter<n; ++node_counter){
        if (pacingNodeCouplesHaveAdhered[node_counter]){
            continue;
        }
        /**
         * - For each node pair id0 and id1, calculate the distance  \f$ d \f$ and the unit normal between them
         */
        int id0 = pacingNodeCouples0[node_counter];
        int id1 = pacingNodeCouples1[node_counter];
        double multiplier = packingMultiplier;
        double sigmoidSaturation = sigmoidSaturationForPacking;
        double distanceShiftScale  = 0.6*packingThreshold;

        double dx = Nodes[id0]->Position[0] - Nodes[id1]->Position[0];
        double dy = Nodes[id0]->Position[1] - Nodes[id1]->Position[1];
        double dz = Nodes[id0]->Position[2] - Nodes[id1]->Position[2];
        double d = pow((dx*dx + dy*dy + dz*dz),0.5);
        double averageMass = 0.5 *( Nodes[id0]->mass + Nodes[id1]->mass );
        /**
         * - If the distance between the node pair is zero, then set the distance to the zero threshold value set in the function, to avoid division by zero.
         */
        double normal[3]= {1.0,0.0,0.0};
        if (d > zeroThreshold){
            normal[0] = dx/d;
            normal[1] = dy/d;
            normal[2] = dz/d;
        }
        else{
            d = zeroThreshold;
        }
        /**
         * - Shift the distance by \f$ d_{0} \f$, such that the sigmoid is 0.5 at this distance.
         * Currently \f$ d_{0} \f$ is set to 60 percent of Simulation#packingThreshold. And calculate the sigmoid function value.
         */
        double shifted_d = d - distanceShiftScale; //shift the distance by the selected percentage,
        double sigmoid = 1 / (1 + exp(sigmoidSaturation/ packingThreshold * 2.0 * shifted_d));
        //F0 = sigmoid * mass
        //F1 = -sigmoid *mass
        //dFi / dXi = - dFi / dxj
        //dFj / dXj = - dFj / dxi
        /**
         * - Shift the distance by \f$ d_{0} \f$, such that the sigmoid is 0.5 at this distance.
         * Currently \f$ d_{0} \f$ is set to 60 percent of Simulation#packingThreshold. And calculate the sigmoid function value.
         * Calculate the derivative, leave the multiplication by the amplitude (mass) to the last stage.
         */
        double dSdXi = (-sigmoidSaturation/ packingThreshold * 2.0) * sigmoid * (1-sigmoid);
        gsl_matrix* normalMat = gsl_matrix_calloc(3,1);
        gsl_matrix_set(normalMat, 0,0,normal[0]);
        gsl_matrix_set(normalMat, 1,0,normal[1]);
        gsl_matrix_set(normalMat, 2,0,normal[2]);
        gsl_matrix* norm_normT = gsl_matrix_calloc(3,3);
        gsl_blas_dgemm (CblasNoTrans, CblasTrans,1.0, normalMat, normalMat, 0.0, norm_normT);
        gsl_matrix* dFidXi = gsl_matrix_calloc(3,3);
        gsl_matrix_set_identity(dFidXi);

        gsl_matrix_sub(dFidXi,norm_normT);  //(I - norm*norm^T)
        gsl_matrix_scale(dFidXi,sigmoid/d);	//sigmoid/d * (I - norm*norm^T)
        gsl_matrix_scale(norm_normT,dSdXi); // norm*norm^T * dSdXi
        gsl_matrix_add(dFidXi,norm_normT);  // sigmoid/d * (I - norm*norm^T) + norm*norm^T * dSdXi
        gsl_matrix_scale(dFidXi,multiplier * averageMass); //scale by mass and multiplier when needed

        //dFidxi = -dFidxj
        //dFj/dxj = dFidxi = -dFj/dxi

        bool printForces = false;
        if (printForces){
            std::cout<<" id0-id1: "<<id0<<std::endl;
            Elements[0]->displayMatrix(dFidXi,"dFidXi");
        }
        //bool printForces = true;
        /**
         * - Add the resulting 3 by 3 derivative matrices are added into the Jacobian in the form:
         * \f[
         *  - \frac{d\mathbf{F_{0}}}{d\mathbf{x_{0}}} \:\: \rightarrow
         *
         * K
         *  \begin{bmatrix}
         *  	      & \vdots          & \vdots   & \vdots            &        \\
                \dots & x_{3id0,3id0}   & \dots    & x_{3id0,3id0+2}   & \dots  \\
                \dots & \vdots          & \dots    & \vdots            & \dots  \\
                \dots & x_{3id0+2,3id0} & \dots    & x_{3id0+2,3id0+2} & \dots  \\
                      & \vdots          & \vdots   & \vdots             &
            \end{bmatrix},
            \f]
            \f[
            \frac{d\mathbf{F_{0}}}{d\mathbf{x_{0}}} \:\: \rightarrow
                \begin{bmatrix}
                x_{3id0,3id1}   & \dots\\
                \dots           & x_{3id0+2,3id1+2}
                \end{bmatrix},

            -\frac{d\mathbf{F_{0}}}{d\mathbf{x_{0}}} \:\: \rightarrow
                \begin{bmatrix}
                x_{3id1,3id1}   & \dots\\
                \dots           & x_{3id1+2,3id1+2}
                \end{bmatrix},

            \frac{d\mathbf{F_{0}}}{d\mathbf{x_{0}}} \:\: \rightarrow
                \begin{bmatrix}
                x_{3id1,3id0}   & \dots\\
                \dots           & x_{3id1+2,3id0+2}
                \end{bmatrix}.

            \f]
         */
        for (size_t i=0; i<3; ++i){
            for (size_t j=0; j<3; ++j){
                double derivativeij = gsl_matrix_get(dFidXi,i,j);
                addValueToMatrix(K,3*id0+i,3*id0+j,-derivativeij);

                //dFidxi = -dFidxj
                addValueToMatrix(K,3*id0+i,3*id1+j,derivativeij);

                //dFj/dxj = dFidxi
                addValueToMatrix(K,3*id1+i,3*id1+j,-derivativeij);

                //dFj/dxi = -dFidxi
                addValueToMatrix(K,3*id1+i,3*id0+j,derivativeij);
                            }
        }
        gsl_matrix_free(normalMat);
        gsl_matrix_free(dFidXi);
        gsl_matrix_free(norm_normT);
    }
}

void Simulation::addValueToMatrix(gsl_matrix* K, int i, int j , double value){
	value += gsl_matrix_get(K,i,j);
	gsl_matrix_set(K,i,j,value);
}

void Simulation::processDisplayDataAndSave(){
	if (saveData && ( (int) (currSimTimeSec/dt) )% dataSaveInterval == 0){ //timestep % dataSaveInterval == 0){
		saveStep();
    }
}

void Simulation::updateNodeMasses(){
    for(const auto& itNode : Nodes){
        itNode->mass = 0;
    }
    for( const auto&  itElement : Elements){
        itElement->assignVolumesToNodes(Nodes);
    }
}

void Simulation::updateNodeViscositySurfaces(){
    //not using range based loops here to ensure openMP comaptibility
    #pragma omp parallel for
    for(std::vector<std::unique_ptr<ShapeBase>>::iterator itElement = Elements.begin(); itElement <Elements.end(); ++itElement){
        (*itElement)->calculateViscositySurfaces();
	}
    //calculated the areas, now assigning them
    for(const auto& itNode : Nodes){
        itNode->viscositySurface = 0;
	}
    for( const auto& itElement : Elements){
        itElement->assignViscositySurfaceAreaToNodes(Nodes);
    }
}

void 	Simulation:: updateElementToConnectedNodes(const std::vector <std::unique_ptr<Node>>& Nodes){
    for(auto& itNode : Nodes){
        if (itNode->mass > 0){//an ablated node will have this as zero
            size_t n = itNode->connectedElementIds.size();
            for (size_t i=0; i<n; ++i){
                itNode->connectedElementWeights[i] = Elements[itNode->connectedElementIds[i]]->VolumePerNode/itNode->mass;
			}
	    }
    }
}

void Simulation::fillInNodeNeighbourhood(){
    for (const auto& itEle : Elements){
        itEle->fillNodeNeighbourhood(Nodes);
	}
}

void Simulation::setBasalNeighboursForApicalElements(){
    /** The basal neighbours are set with ShapeBase#setBasalNeigElementId function.
     */
    //not using range based loops here to ensure openMP comaptibility
    #pragma omp parallel for
    for(std::vector<std::unique_ptr<ShapeBase>>::iterator itElement = Elements.begin(); itElement <Elements.end(); ++itElement){
        (*itElement)->setBasalNeigElementId(Elements);
	}
}

void Simulation::fillInElementColumnLists(){
    //not using range based loops here to ensure openMP comaptibility
    #pragma omp parallel for
    for(std::vector<std::unique_ptr<ShapeBase>>::iterator itElement = Elements.begin(); itElement <Elements.end(); ++itElement){
        if ( (*itElement)->tissueType ==0 ){//Check only the columnar layer elements.
            if ((*itElement)->tissuePlacement == 0){
				//start from the basal element and move up
				//then you can copy the element numbers to the other elements on the list
                (*itElement)->constructElementStackList(TissueHeightDiscretisationLayers, Elements);
			}
		}
	}
}

void Simulation::detectPacingToEnclosingSurfacesNodes(){
    /**
     * The current disctance of packing surfaces are stored in Simulation#zEnclosementBoundaries. The z distance of
     * nodes to these surfaces are checked against the threshold for packing detection Simulation#packingDetectionToEnclosingSurfacesThreshold.
     * The actual forces of packing are not calculate4d here, this is for detection of potentially packing nodes.
     */
	packingToEnclosingSurfacesThreshold = 3;  //pack to the boundary at 3 microns distance
	packingDetectionToEnclosingSurfacesThreshold = 1.2 * packingToEnclosingSurfacesThreshold;
	double t2 = packingDetectionToEnclosingSurfacesThreshold*packingDetectionToEnclosingSurfacesThreshold;	//threshold square for rapid calculation
	nodesPackingToPositiveSurface.clear();
	nodesPackingToNegativeSurface.clear();
	initialWeightPackingToPositiveSurface.clear();
	initialWeightPackingToNegativeSurface.clear();
	//calculate the current boundaries:
	if (currSimTimeSec < initialTimeToEncloseTissueBetweenSurfacesSec){
		return;
	}
	else if (currSimTimeSec>=finalTimeToEncloseTissueBetweenSurfacesSec){
		zEnclosementBoundaries[0] = finalZEnclosementBoundaries[0];
		zEnclosementBoundaries[1] = finalZEnclosementBoundaries[1];
	}
	else{
		double totalTime = finalTimeToEncloseTissueBetweenSurfacesSec - initialTimeToEncloseTissueBetweenSurfacesSec;
		double currTimeDiff = currSimTimeSec - initialTimeToEncloseTissueBetweenSurfacesSec;
		double zNegDifference = finalZEnclosementBoundaries[0] - initialZEnclosementBoundaries[0];
		double zPosDifference = finalZEnclosementBoundaries[1] - initialZEnclosementBoundaries[1];
		zEnclosementBoundaries[0] = initialZEnclosementBoundaries[0] + currTimeDiff*zNegDifference / totalTime;
		zEnclosementBoundaries[1] = initialZEnclosementBoundaries[1] + currTimeDiff*zPosDifference / totalTime;

	}
    std::cout<<"curr time in sec: "<<currSimTimeSec<<" z enclosement boundaries: "<<zEnclosementBoundaries[0]<<" "<<zEnclosementBoundaries[1]<<" initial boundaries: "<<initialZEnclosementBoundaries[0]<<" "<<initialZEnclosementBoundaries[1]<<" final boundaries: "<<finalZEnclosementBoundaries[0]<<" "<<finalZEnclosementBoundaries[1]<<std::endl;

	//go over nodes to detect packing:
    for (const auto& itNode : Nodes){
        if (itNode->mass >0){ //node is not ablated
			bool checkAgainstPositiveSurface = false;
			bool checkAgainstNegativeSurface = false;
            if (itNode->atCircumference){
				//Node is at the ciscumference, with sufficient rotation, it can pack anywhere
				//chack against both surfaces.
				checkAgainstPositiveSurface = true;
				checkAgainstNegativeSurface = true;
			}
			else {
				//basal node of the columnar layer is checked against the negative surface:
                if ( itNode->tissuePlacement == 0){
					//node is basal
					//if it is columnar layer, check against negative surface:
                    if( itNode->tissueType ==0 ){
						checkAgainstNegativeSurface = true;
					}
                    else if(itNode->tissueType ==1 ){
						//if it is peripodial layer, check against positive surface:
						checkAgainstPositiveSurface = true;
					}
				}
				//if tissue placement is apical, and the tissue tyoe is columnar,
				//check against the positive surface:
                if( itNode->tissuePlacement == 1 && itNode->tissueType ==0 ){	//Node is apical, check against positive border
					checkAgainstPositiveSurface = true;
				}
			}
			if( checkAgainstNegativeSurface ){
                std::array<double,3> pos = itNode->getCurrentPosition();
				double dz = pos[2] - zEnclosementBoundaries[0];
				double d2 = dz*dz;
				if (d2<t2){
					//node is close enough to pack to surface:
                    nodesPackingToNegativeSurface.push_back(itNode->Id);
					double d = pow (d2,0.5);
					initialWeightPackingToNegativeSurface.push_back(dz/d);
				}
			}
			if( checkAgainstPositiveSurface){
                std::array<double,3> pos = itNode->getCurrentPosition();
                double dz = pos[2] - zEnclosementBoundaries[1];
				double d2 = dz*dz;
				if (d2<t2){
					//node is close enough to pack to surface:
                    nodesPackingToPositiveSurface.push_back(itNode->Id);
					double d = pow (d2,0.5);
					initialWeightPackingToPositiveSurface.push_back(dz/d);
				}
			}
		}
	}
}


void 	Simulation::assignFoldRegionAndReleasePeripodial(Node* NodeMaster, Node* NodeSlave ){
    size_t masterXGridIndex = floor( (NodeMaster->Position[0] - boundingBox[0][0])/boundingBoxSize[0] );
    size_t masterYGridIndex = floor( (NodeMaster->Position[1] - boundingBox[0][1])/boundingBoxSize[1] );
    /**
     * If two nodes are adhered they are on an indenting surface and therefore on a fold initiation surface.
     * The Node#onFoldInitiation booleans are set to true. Then all the nodes falling in between these two nodes
     * are also declared to be on fold initiation.
     */
    NodeMaster->onFoldInitiation = true;
	NodeSlave->onFoldInitiation = true;
    //std::cout<<"assigning to fold initiation, nodes "<<(*itNode)->Id <<" and "<<(*itNodeSlave)->Id<<std::endl;
	//check for other nodes in the vicinity:
	//if I am on the apical side, I will check the borders from basal sides
	//If I am on the basal side, I will chack from apical sides:
    int masterCorrespondingX = NodeMaster->Position[0];
    int slaveCorrespoindingX = NodeSlave->Position[0];
    int masterCorrespondingY = NodeMaster->Position[1];
    int slaveCorrespoindingY = NodeSlave->Position[1];

    for (auto itInt : NodeMaster->immediateNeigs){
        if(Nodes[itInt]->tissuePlacement != NodeMaster->tissuePlacement){
            int commonOwnerId = NodeMaster->getCommonOwnerId(Nodes[itInt].get());
			if (commonOwnerId > -1){
				//there is common owner
                bool nodesConnected = Elements[commonOwnerId]->areNodesDirectlyConnected(NodeMaster->Id, Nodes[itInt]->Id);
				if (nodesConnected){
                    masterCorrespondingX = Nodes[itInt]->Position[0];
                    masterCorrespondingY = Nodes[itInt]->Position[1];
					break;
				}
			}
		}
	}
    for (auto itInt : NodeSlave->immediateNeigs){
        if(Nodes[itInt]->tissuePlacement != NodeSlave->tissuePlacement){
            int commonOwnerId = NodeSlave->getCommonOwnerId(Nodes[itInt].get());
			if (commonOwnerId > -1){
                bool nodesConnected = Elements[commonOwnerId]->areNodesDirectlyConnected(NodeSlave->Id, Nodes[itInt]->Id);
				if (nodesConnected){
                    slaveCorrespoindingX = Nodes[itInt]->Position[0];
                    slaveCorrespoindingY = Nodes[itInt]->Position[1];
					break;
				}
			}
		}
	}
    double xmin = min(masterCorrespondingX, slaveCorrespoindingX);
    double xmax = max(masterCorrespondingX, slaveCorrespoindingX);
    double ymin = min(masterCorrespondingY, slaveCorrespoindingY);
    double ymax = max(masterCorrespondingY, slaveCorrespoindingY);
    double zmax = max(NodeMaster->Position[2], NodeSlave->Position[2]);
    double zmin = min(NodeMaster->Position[2], NodeSlave->Position[2]);
    ymin -= 4.0*packingDetectionThresholdGrid[masterXGridIndex][masterYGridIndex];
    ymax += 4.0*packingDetectionThresholdGrid[masterXGridIndex][masterYGridIndex];
    for (auto& itNodeInBetween :  Nodes){
        if(itNodeInBetween->tissuePlacement == NodeMaster->tissuePlacement ){
        if(!itNodeInBetween->onFoldInitiation ){
            //for apical curves, chack for zmax, for basal , check for z min
            if( (itNodeInBetween->tissuePlacement == 1 && itNodeInBetween->Position[2] <= zmax ) ||(itNodeInBetween->tissuePlacement == 0 && itNodeInBetween->Position[2] >= zmin )   ){
                if(itNodeInBetween->Position[0] >= xmin && itNodeInBetween->Position[0] <= xmax){
                    if(itNodeInBetween->Position[1] >= ymin && itNodeInBetween->Position[1] <= ymax){
                        itNodeInBetween->onFoldInitiation = true;
                    }
                }
            }
        }}
    }
}

void Simulation::detectPacingNodes(){
    /**
     * The current distance of node-node packing is dependent on the local average side length and
     * is sored in a 10x5 grid Simulation#packingDetectionThresholdGrid.
     */
    double periAverageSideLength = 0,colAverageSideLength = 0;
	getAverageSideLength(periAverageSideLength, colAverageSideLength);

	if (thereIsPeripodialMembrane){
		colAverageSideLength = (periAverageSideLength+colAverageSideLength)/2.0;
	}
	double packingDetectionThresholdGridSq[10][5];
    for (size_t i=0;i<10;++i){
        for(size_t j=0;j<5;++j){
            packingDetectionThresholdGrid[i][j] = 0.4 * packingDetectionThresholdGrid[i][j];
			packingDetectionThresholdGridSq[i][j] = packingDetectionThresholdGrid[i][j]*packingDetectionThresholdGrid[i][j];
		}
	}
	packingThreshold = 0.4*colAverageSideLength;  //0.6 was old value
	packingDetectionThreshold = 1.2 * packingThreshold; //0.9 * packingThreshold;
	packingMultiplier = 5000; //1000 in original, increasing 5 fold to test run54 on 16 april 2018
	sigmoidSaturationForPacking = 5;
	double t2 = packingDetectionThreshold*packingDetectionThreshold;	//threshold square for rapid calculation
	pacingNodeCouples0.clear();
	pacingNodeCouples1.clear();
	pacingNodeCouplesHaveAdhered.clear();
     //not using range based loops here to ensure openMP comaptibility
    #pragma omp parallel for
    for(std::vector<std::unique_ptr<ShapeBase>>::iterator itElement = Elements.begin(); itElement <Elements.end(); ++itElement){
        if ((*itElement)->tissueType == 0 && ((*itElement)->tissuePlacement == 1 || (*itElement)->spansWholeTissue)){
			//columnar element at the apical surface or spans whole tissue
            (*itElement)->calculateApicalNormalCurrentShape();
		}
	}
    for (auto itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
        bool NodeHasPacking = (*itNode)->checkIfNodeHasPacking();
        if (NodeHasPacking){
            std::array<double,3> pos = (*itNode)->getCurrentPosition();
            int masterXGridIndex = floor( (pos[0] - boundingBox[0][0])/boundingBoxSize[0] );
            int masterYGridIndex = floor( (pos[1] - boundingBox[0][1])/boundingBoxSize[1] );
            t2 = packingDetectionThresholdGridSq[masterXGridIndex][masterYGridIndex];
            for (auto itNodeSlave=itNode+1; itNodeSlave<Nodes.end(); ++itNodeSlave){
                bool SlaveNodeHasPacking = (*itNodeSlave)->checkIfNodeHasPacking();
                if (SlaveNodeHasPacking){
                    if ((*itNode)->tissuePlacement == (*itNodeSlave)->tissuePlacement){
                        //nodes can pack, are they connected?
                        bool neigbours = (*itNode)->isMyNeig((*itNodeSlave)->Id);
                        if (neigbours){
                            continue;
                        }
                        //check if the master node is in the collapsed list?
                        //the vector is sorted, as I sort it to remove duplicates
                        if (binary_search((*itNode)->collapsedWith.begin(), (*itNode)->collapsedWith.end(),(*itNodeSlave)->Id)){
                            //the couple is already collapsed:
                            continue;
                        }
                        //check if the collapsed nodes of master are neigs of slave
                        neigbours = (*itNode)->isNeigWithMyCollapsedNodes((*itNodeSlave)->Id,Nodes);
                        if (neigbours){
                            continue;
                        }
                        neigbours = (*itNodeSlave)->isNeigWithMyCollapsedNodes((*itNode)->Id,Nodes);
                        if (neigbours){
                            continue;
                        }
                        //the nodes can potentially pack, are they close enough?
                        std::array<double,3> posSlave = (*itNodeSlave)->getCurrentPosition();
                        double dx = pos[0] - posSlave[0];
                        double dy = pos[1] - posSlave[1];
                        double dz = pos[2] - posSlave[2];
                        double d2 = dx*dx + dy*dy + dz*dz;
                        if (d2<t2){
                            //close enough for packing , add to list:
                            pacingNodeCouples0.push_back((*itNode)->Id);
                            pacingNodeCouples1.push_back((*itNodeSlave)->Id);
                            pacingNodeCouplesHaveAdhered.push_back(false);
                        }
                        bool checkForCurvedRegions = false;

                        double curveIdentificationThresholdSq = packingDetectionThresholdGrid[masterXGridIndex][masterYGridIndex];
                        curveIdentificationThresholdSq *= curveIdentificationThresholdSq;
                        curveIdentificationThresholdSq *=9;//0.625;
                        if (d2<curveIdentificationThresholdSq){
                            if (thereIsEmergentEllipseMarking && (!(*itNode)->onFoldInitiation || !(*itNodeSlave)->onFoldInitiation) ){
                                //only check if the nodes are not covered by emergent ellipse bands
                                checkForCurvedRegions=true;
                            }
                        }
                        if (checkForCurvedRegions){
                            //std::cout<<" check for curved regions active"<<std::endl;
                            int elementIdMaster = (*itNode)->connectedElementIds[0];
                            int elementIdSlave = (*itNodeSlave)->connectedElementIds[0];

                            double dotP = Elements[elementIdSlave]->dotProduct3D(Elements[elementIdMaster]->apicalNormalCurrentShape,Elements[elementIdSlave]->apicalNormalCurrentShape);
                            //std::cout<<"dotp: "<<dotP<<std::endl;
                            if (dotP <0){
                                //normals point opposite directions
                                //do they face each other?
                                bool releaseNodes = false;
                                std::array<double,3> connectingVec = {0.0};
                                for (int dimIter=0; dimIter<3; ++dimIter){
                                    connectingVec[dimIter] = (*itNode)->Position[dimIter] - (*itNodeSlave)->Position[dimIter];
                                }
                                double dotPmaster = Elements[elementIdMaster]->dotProduct3D(connectingVec,Elements[elementIdMaster]->apicalNormalCurrentShape);
                                if (dotPmaster >0 ){
                                    releaseNodes = true;
                                }
                                else{
                                    double dotPslave = Elements[elementIdSlave]->dotProduct3D(connectingVec,Elements[elementIdSlave]->apicalNormalCurrentShape);
                                    if (dotPslave>0){
                                        releaseNodes = true;
                                    }
                                }
                                if (!releaseNodes){
                                    continue;
                                }
                                assignFoldRegionAndReleasePeripodial((*itNode).get(),(*itNodeSlave).get());
                            }
                        }
                    }
                }
            }
        }
    }
}

void Simulation::addPackingForces(gsl_matrix* gExt){
    /**
     * The packing forces calculated in Simulation#calculatePackingForcesImplicit3D are added to the
     * system force vector. The vector containd all dimensions of all nodes in one column, therefore indexing is carried out
     * accordingly.
     */
    for (size_t j=0; j<nNodes; ++j){
        double Fx = PackingForces[j][0];
		double Fy = PackingForces[j][1];
		double Fz = PackingForces[j][2];
		int indice = j*3;
		Fx += gsl_matrix_get(gExt,indice,0);
		gsl_matrix_set(gExt,indice,0,Fx);
		Fy += gsl_matrix_get(gExt,indice+1,0);
		gsl_matrix_set(gExt,indice+1,0,Fy);
		Fz += gsl_matrix_get(gExt,indice+2,0);
		gsl_matrix_set(gExt,indice+2,0,Fz);
	}
}

void Simulation::updateElementPositions(){
    for(const auto& itElement : Elements){
        itElement->updatePositions(Nodes);
	}
}


void Simulation::updateElementPositionsSingle(size_t i ){
	Elements[i]->updatePositions(Nodes);
}

void Simulation::assignTips(){
	double xTips[2] ={ 10000, -10000};
	double yTips[2] ={ 10000, -10000};
    for (size_t i =0; i<nNodes; ++i){
		if (Nodes[i]->tissuePlacement == 0 && Nodes[i]->tissueType == 0 ) { //taking the basal nodes of the columnar layer only
			if (Nodes[i]->Position[0] < xTips[0] ){
				dorsalTipIndex = i;
				xTips[0] = Nodes[i]->Position[0];
			}
			if (Nodes[i]->Position[1] < yTips[0] ){
				anteriorTipIndex = i;
				yTips[0] = Nodes[i]->Position[1];
			}
			if (Nodes[i]->Position[0] > xTips[1] ){
				ventralTipIndex = i;
				xTips[1] = Nodes[i]->Position[0];
			}
			if (Nodes[i]->Position[1] > yTips[1] ){
				posteriorTipIndex = i;
				yTips[1] = Nodes[i]->Position[1];
			}
		}
	}
	//I have identified the tips assuming the tissue is represented in full.
	//If there is symmetry in the input, then I will need to correct this.
	//In the case of x symmetry, the minimum x will be correct in DV, but not the maximum x.
	//In the case of y symmetry, the maximum y will be correct in AP. but not the minimum.
	if (symmetricX){
		double delta = 0.02;
		//the ventralTipIndex is correct, it is at the minimum of x axis. The dorsalTipIndex is at x=0. but there are a series of nodes at x=0,
		//and the selected one should be at y=0 too.
        for (size_t i =0; i<nNodes; ++i){
				if (Nodes[i]->tissuePlacement == 0 && Nodes[i]->tissueType == 0 ) { //taking the basal nodes of the columnar layer only
					if (Nodes[i]->Position[0]< delta &&  Nodes[i]->Position[0]> -delta){
						//this node is at x = 0, is it at y = 0 too?
						if (Nodes[i]->Position[1]< delta &&  Nodes[i]->Position[1]> -delta){
							dorsalTipIndex = i;
							break;
						}
					}
				}
		}
	}
	if (symmetricY){
		double delta = 0.02;
		//the posteriorTipIndex is correct, it is at the maximum of y axis. The anteriorTipIndex is at y=0. but there are a series of nodes at y=0,
		//and the selected one should be at x=0 too.
        for (size_t i =0; i<nNodes; ++i){
				if (Nodes[i]->tissuePlacement == 0 && Nodes[i]->tissueType == 0 ) { //taking the basal nodes of the columnar layer only
					if (Nodes[i]->Position[1]< delta &&  Nodes[i]->Position[1]> -delta){
						//this node is at x = 0, is it at y = 0 too?
						if (Nodes[i]->Position[0]< delta &&  Nodes[i]->Position[0]> -delta){
							anteriorTipIndex = i;
							break;
						}
					}
				}
		}
	}
    std::cout<<"Tip node indexes - dorsalTipIndex: "<<dorsalTipIndex<<" ventralTipIndex: "<<ventralTipIndex<<std::endl;
    std::cout<<"Tip node indexes - anteriorTipIndex: "<<anteriorTipIndex<<" posteriorTipIndex: "<<posteriorTipIndex<<std::endl;
}

void Simulation::alignTissueDVToXPositive(){
    /**
     * The bounding box calculateion is relies on the fact that the tissue is lying parallel to the x axis in its initial long axis.
     * Bring the long axis defined by Simulation#ventralTipIndex and Simulation#dorsalTipIndex Node#Ids parallel to the x axis with a rigid
     * body rotation.
     */
	if (!symmetricX){
        std::array<double,3> u = {0.0, 0.0, 0.0};
		//For simulations with no external viscosity, the position of Dorsal tip is fixed, Ventral tip is fixed in y and z, another node is fixed in z
        for (size_t i=0;i<3;++i){
			u[i] = Nodes[ventralTipIndex]->Position[i] - Nodes[dorsalTipIndex]->Position[i];
		}
        (void) Elements[0]->normaliseVector3D(u);
        std::array<double,3> v = {1.0, 0.0, 0.0};
        double c, s;
		Elements[0]->calculateRotationAngleSinCos(u,v,c,s);
        std::array<double,3> rotAx  = {0.0, 0.0, 1.0};
        std::array<double,9> rotMat = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        Elements[0]->calculateRotationAxis(u,v,rotAx,c);	//calculating the rotation axis that is perpendicular to both u and v
		Elements[0]->constructRotationMatrix(c,s,rotAx,rotMat);
        for(size_t i=1;i<nNodes;++i){
            for(size_t j = 0; j< Nodes[i]->nDim; ++j){
				u[j] = Nodes[i]->Position[j] - Nodes[0]->Position[j];
			}
			Elements[0]->rotateVectorByRotationMatrix(u,rotMat);
            for(size_t j = 0; j< Nodes[i]->nDim; ++j){
				Nodes[i]->Position[j] = u[j] + Nodes[0]->Position[j];
			}
		}
        for(auto &currElement : Elements){
            currElement->updatePositions(Nodes);
		}
	}
}


void Simulation::calculateDVDistance(){
	double d[3];
	if (symmetricX){
		for (int i=0;i<3;++i){
			d[i] = Nodes[ventralTipIndex]->Position[i];
			d[i] *=d[i];
		}
	}
	else{
		for (int i=0;i<3;++i){
			d[i] = Nodes[dorsalTipIndex]->Position[i] - Nodes[ventralTipIndex]->Position[i];
			d[i] *=d[i];
		}
	}
	double dmag = d[0]+d[1]+d[2];
	dmag = pow(dmag,0.5);
	outputFile<<"time: "<<currSimTimeSec<<" DV distance is: "<<dmag<<" ";
    std::cout<<"time: "<<currSimTimeSec<<" DV distance is: "<<dmag<<" ";
	if (symmetricY){
		for (int i=0;i<3;++i){
			d[i] = Nodes[posteriorTipIndex]->Position[i];
			d[i] *=d[i];
		}
	}
	else{
		for (int i=0;i<3;++i){
			d[i] = Nodes[anteriorTipIndex]->Position[i] - Nodes[posteriorTipIndex]->Position[i];
			d[i] *=d[i];
		}
	}
	dmag = d[0]+d[1]+d[2];
	dmag = pow(dmag,0.5);
    outputFile<<" AP distance is: "<<dmag<<std::endl;
    std::cout<<" AP distance is: "<<dmag<<std::endl;
}

void Simulation::resetForces(bool resetPacking){
    for (auto& currNodeForces : SystemForces){
		//3 dimensions
        currNodeForces.fill(0.0);
	}
    if (resetPacking){
        for (auto& currNodePackingForces : PackingForces){
            currNodePackingForces.fill(0.0);
        }
    }
}


void Simulation::calculateApicalSize(){
	double sizeLim[2][2] = {{0.0,0.0},{0.0,0.0}};
	bool found[2][2] = {{false,false},{false,false}};
    for (const auto& itNode : Nodes){
        if (itNode->tissueType == 0 && itNode->tissuePlacement == 1){ //checking only columnar apical layer nodes
			for (int i=0; i<2; ++i){
                if ( itNode->Position[i] < sizeLim[0][i] ){
                    sizeLim[0][i] = itNode->Position[i];
					found[0][i] = true;
				}
                else if(itNode->Position[i]>sizeLim[1][i]){
                    sizeLim[1][i] = itNode->Position[i];
					found[1][i] = true;
				}
			}
		}
	}
	if (!found[0][0] && !found[0][1] && !found[1][0] && !found[1][1]){
        std::cerr<<" error in apical bounding  box calculation! Found? :"<<found[0][0]<<" "<<found[0][1]<<" "<<found[1][0]<<" "<<found[1][1]<<std::endl;
	}
	double DV = sizeLim[1][0] - sizeLim[0][0];
	double AP = sizeLim[1][1] - sizeLim[0][1];
    outputFile<<"At time: "<<currSimTimeSec<<" apical bounding box size: "<<DV<<" "<<AP<<std::endl;
}

void Simulation::calculateBoundingBox(){
	boundingBox[0][0] =  100000.0;	//lower left x
	boundingBox[0][1] =  100000.0;	//lower left y
	boundingBox[0][2] =  100000.0;	//lower z
	boundingBox[1][0] = -100000.0;	//upper right x
	boundingBox[1][1] = -100000.0;	//upper right y
	boundingBox[1][2] = -100000.0;	//upper z
	bool found[2][3] = {{false,false,false},{false,false,false}};
    for (const auto& itNode : Nodes){
		//Do not count node if it is part of an explicit ECM:
        if(!itNode->allOwnersECMMimicing){
            for (size_t i=0; i<itNode->nDim; ++i){
                if ( itNode->Position[i] < boundingBox[0][i] ){
                    boundingBox[0][i] = itNode->Position[i];
					found[0][i] = true;
				}
                else if(itNode->Position[i]>boundingBox[1][i]){
                    boundingBox[1][i] = itNode->Position[i];
					found[1][i] = true;
				}
			}
		}
	}
	if (symmetricY){
		boundingBox[0][1] = (-1.0)*boundingBox[1][1]; //if there is Y symmetricity, then the bounding box is extended to double the size in y
	}
	for (int i=0; i<3; ++i){
		boundingBoxSize[i] = boundingBox[1][i] - boundingBox[0][i];
	}
	if (!found[0][0] && !found[0][1] && !found[0][2] && !found[1][0] && !found[1][1] && !found[1][2]){
        std::cerr<<" error in bounding box calculation! Found? :"<<found[0][0]<<" "<<found[0][1]<<" "<<found[0][2]<<" "<<found[1][0]<<" "<<found[1][1]<<" "<<found[1][2]<<std::endl;
	}
}

void Simulation::saveStep(){
    outputFile<<"Saving step: "<< timestep<<" this is :"<<currSimTimeSec<<" sec ("<<currSimTimeSec/3600<<" hr )"<<std::endl;
	writeSaveFileStepHeader();
	writeNodes();
	writeElements();
	writeSaveFileStepFooter();
	writeTensionCompression();
    writeGrowth();
	writePhysicalProp();
	writeGrowthRedistribution();
	writeNodeBinding();
}

void Simulation::writeSaveFileStepHeader(){
    saveFileMesh<<"=============== TIME: ";
	saveFileMesh.precision(6);
	saveFileMesh.width(10);
	saveFileMesh<<currSimTimeSec;
    saveFileMesh<<"==================================================="<<std::endl;
}

void Simulation::writeSaveFileStepFooter(){
	saveFileMesh<<"=============== END OF TIME: ";
	saveFileMesh.precision(6);
	saveFileMesh.width(10);
	saveFileMesh<<currSimTimeSec;
    saveFileMesh<<"============================================"<<std::endl;
}

void Simulation::writeNodes(){
    saveFileMesh<<nNodes<<std::endl;
    for (size_t i = 0; i<nNodes; ++i){
		if(Nodes[i]->nDim==2){
			saveFileMesh.precision(10);saveFileMesh.width(20);
			saveFileMesh<<Nodes[i]->Position[0];
			saveFileMesh.precision(10);saveFileMesh.width(20);
			saveFileMesh<<Nodes[i]->Position[1];
			saveFileMesh.precision(10);saveFileMesh.width(20);
			saveFileMesh<<0.0;

		}
		else{
			int ndim = Nodes[i]->nDim;
			for (int j=0;j<ndim;++j){
				saveFileMesh.precision(10);saveFileMesh.width(20);
				saveFileMesh<<Nodes[i]->Position[j];
			}

		}
		saveFileMesh.width(4);
		saveFileMesh<<Nodes[i]->tissuePlacement;
		saveFileMesh.width(4);
		saveFileMesh<<Nodes[i]->tissueType;
		saveFileMesh.width(4);
		saveFileMesh<<Nodes[i]->atCircumference;
        saveFileMesh<<std::endl;
	}
}

void Simulation::writeElements(){
    saveFileMesh<<nElements<<std::endl;
    for (size_t i = 0; i<nElements; ++i){
		int shapetype = Elements[i]->getShapeType();
		saveFileMesh.width(4);
		saveFileMesh<<shapetype;
		if(shapetype == 4 ){
			double height = Elements[i]->getElementHeight();
			saveFileMesh.precision(5);saveFileMesh.width(12);
			saveFileMesh<<height;
		}
		saveFileMesh.width(4);
		saveFileMesh<<Elements[i]->IsAblated;
        size_t nodeNumber = Elements[i]->getNodeNumber();
        const vector<int>&  NodeIds = Elements[i]->getNodeIds();
        for (auto currNodeId : NodeIds ){
            saveFileMesh.width(9);saveFileMesh<<currNodeId;
		}
        size_t dim  = Elements[i]->getDim();
        const vector<std::array<double,3>> refPos = Elements[i]->getReferencePos();
        for (size_t j = 0; j<nodeNumber; ++j ){
            for (size_t k = 0; k<dim; ++k ){
				saveFileMesh.precision(4);saveFileMesh.width(15);
				saveFileMesh<<refPos[j][k];
			}
		}
        saveFileMesh<<std::endl;
	}
}

void Simulation::writeTensionCompression(){
    for(const auto& itElement : Elements){
        for (size_t j=0; j<6; ++j){
            double S = gsl_matrix_get(itElement->Strain,j,0);
            saveFileTensionCompression.write((char*) &S, sizeof S);
        }
    }
	saveFileTensionCompression.flush();
}



void Simulation::writeGrowthRedistribution(){
    for(const auto& itElement : Elements){
        bool thereIsDistribution = itElement->thereIsGrowthRedistribution;
        bool shrinksElement = itElement->growthRedistributionShrinksElement;
		saveFileGrowthRedistribution.write((char*) &thereIsDistribution, sizeof thereIsDistribution);
		saveFileGrowthRedistribution.write((char*) &shrinksElement, sizeof shrinksElement);
	}
	saveFileGrowthRedistribution.flush();
}

void Simulation::writeNodeBinding(){
	vector<int> slaveDOFs, masterDOFs;
    int counter = 0;
    for (const auto& itNode : Nodes){
        for (size_t i=0; i<3; ++i){
            if (itNode->slaveTo[i]>0){
                int slavedof = 3*itNode->Id + i;
                int masterdof = 3*itNode->slaveTo[i] + i;
				slaveDOFs.push_back(slavedof);
				masterDOFs.push_back(masterdof);
				counter++;
			}
		}
	}
    saveFileNodeBinding<<counter<<" "<<std::endl;
    for(int i=0;i<counter;++i){
		int slave  =slaveDOFs[i];
		int master = masterDOFs[i];
		saveFileNodeBinding<<slave<<" ";
        saveFileNodeBinding<<master<<" "<<std::endl;
        //std::cout<<"slaveDOFs :"<<slaveDOFs[i]<<" masterDOFs "<<masterDOFs[i]<<std::endl;
	}
}

void Simulation::writeGrowth(){
    for(const auto& itElement : Elements){
        gsl_matrix* currFg = itElement->getFg();
        const std::array<double,3>& growthRate = itElement->getGrowthRate();
        for (int j=0; j<3; ++j){
            for (int k=0; k<3; ++k){
                double Fgjk = gsl_matrix_get(currFg,j,k);
                saveFileGrowth.write((char*) &Fgjk, sizeof Fgjk);
            }
            saveFileGrowthRate.write((char*) &growthRate[j], sizeof growthRate[j]);
        }
        gsl_matrix_free(currFg);
    }
    saveFileGrowth.flush();
    saveFileGrowthRate.flush();
}

void Simulation::writePhysicalProp(){
    for(const auto& itElement : Elements){
        double E = itElement->getYoungModulus();
        double interalVisc	= itElement->getInternalViscosity();
        double zRemodellingSoFar = itElement->getZRemodellingSoFar();
		saveFilePhysicalProp.write((char*) &E, sizeof E);
		saveFilePhysicalProp.write((char*) &interalVisc, sizeof interalVisc);
		saveFilePhysicalProp.write((char*) &zRemodellingSoFar, sizeof zRemodellingSoFar);
	}
    for (const auto& itNode : Nodes){
        saveFilePhysicalProp.write((char*) &itNode->externalViscosity[0], sizeof itNode->externalViscosity[0]);
        saveFilePhysicalProp.write((char*) &itNode->externalViscosity[1], sizeof itNode->externalViscosity[1]);
        saveFilePhysicalProp.write((char*) &itNode->externalViscosity[2], sizeof itNode->externalViscosity[2]);
	}
	saveFilePhysicalProp.flush();
}

void Simulation::calculateGrowth(){
	cleanUpGrowthRates();
	for (int i=0; i<nGrowthFunctions; ++i){
		if (GrowthFunctions[i]->Type == 1){
            calculateGrowthUniform(GrowthFunctions[i].get());
		}
		else if(GrowthFunctions[i]->Type == 2){
            calculateGrowthRing(GrowthFunctions[i].get());
		}
		else if(GrowthFunctions[i]->Type == 3){
            calculateGrowthGridBased(GrowthFunctions[i].get());
		}
	}
    for(auto& itElement : Elements){
        if (itElement->isMutated){
            itElement->updateGrowthByMutation(dt);
		}
	}
}

void Simulation::calculateShapeChange(){
	cleanUpShapeChangeRates();
	for (int i=0; i<nShapeChangeFunctions; ++i){
		if (ShapeChangeFunctions[i]->Type == 1){
            calculateShapeChangeUniform(ShapeChangeFunctions[i].get());
		}
		if (ShapeChangeFunctions[i]->Type == 2){
            calculateShapeChangeMarkerEllipseBased(ShapeChangeFunctions[i].get());
		}
	}
}


void Simulation::cleanUpGrowthRates(){
    for(const auto& currElement : Elements){
        currElement->setGrowthRate(dt,0.0,0.0,0.0);
        currElement->updateGrowthIncrementFromRate();
	}
}

void Simulation::cleanUpShapeChangeRates(){
    for(const auto& currElement : Elements){
        currElement->setShapeChangeRate(0.0,0.0,0.0,0.0,0.0,0.0);
        currElement->setShapeChangeInrementToIdentity();
	}
}

void Simulation::calculateShapeChangeMarkerEllipseBased (GrowthFunctionBase* currSCF){
	if(currSimTimeSec >= currSCF->initTime && currSimTimeSec < currSCF->endTime ){
			gsl_matrix* columnarShapeChangeIncrement = gsl_matrix_calloc(3,3);
            std::array<double,3> growthRates = currSCF->getShapeChangeRateRate();
            for(const auto& itElement : Elements){
                bool appliedToElement = itElement->isShapeChangeAppliedToElement(currSCF->appliedEllipseBandIds, currSCF->applyToBasalECM, currSCF->applyToLateralECM, currSCF->applyTissueApical, currSCF->applyTissueBasal, currSCF->applyTissueMidLine );
				if (appliedToElement){
					gsl_matrix_set_identity(columnarShapeChangeIncrement);
                    itElement->calculateShapeChangeIncrementFromRates(dt, growthRates[0],growthRates[1],growthRates[2],columnarShapeChangeIncrement);
                    itElement->updateShapeChangeIncrement(columnarShapeChangeIncrement);
				}
	    	}
			gsl_matrix_free(columnarShapeChangeIncrement);
		}
}

void Simulation::calculateShapeChangeUniform (GrowthFunctionBase* currSCF){
	if(currSimTimeSec >= currSCF->initTime && currSimTimeSec < currSCF->endTime ){
        std::array<double,3> maxValues = currSCF->getGrowthRate();
        for(const auto& itElement : Elements){
			if ( currSCF->applyToColumnarLayer){
                if (itElement->tissueType == 0){ //columnar layer, grow directly
                    itElement->updateShapeChangeRate(maxValues[0],maxValues[1],maxValues[2],0,0,0);
				}
                else if (itElement->tissueType == 2){ //Linker zone, need to weight the growth
                    double weight = itElement->getColumnarness();
                    itElement->updateShapeChangeRate(weight*maxValues[0],weight*maxValues[1],weight*maxValues[2],0,0,0);
				}
			}
			if ( currSCF->applyToPeripodialMembrane){
                if (itElement->tissueType == 1){ //peripodial membrane, grow directly
                    itElement->updateShapeChangeRate(maxValues[0],maxValues[1],maxValues[2],0,0,0);
				}
                else if (itElement->tissueType == 2){ //Linker zone, need to weight the growth
                    double weight = itElement->getPeripodialness();
                    itElement->updateShapeChangeRate(weight*maxValues[0],weight*maxValues[1],weight*maxValues[2],0,0,0);
				}
			}
		}
	}
}

void Simulation::calculateGrowthUniform(GrowthFunctionBase* currGF){
	if(currSimTimeSec >= currGF->initTime && currSimTimeSec < currGF->endTime ){
		gsl_matrix* columnarFgIncrement = gsl_matrix_calloc(3,3);
		gsl_matrix* peripodialFgIncrement = gsl_matrix_calloc(3,3);
        std::array<double,3> growthRates = currGF->getGrowthRate();
        for(const auto& itElement : Elements){
			gsl_matrix_set_identity(columnarFgIncrement);
			gsl_matrix_set_identity(peripodialFgIncrement);
            if (currGF->applyToBasalECM || currGF->applyToLateralECM){
                if (currGF->applyToBasalECM){
                    if (itElement->isECMMimicing && itElement->tissuePlacement == 0){
                        itElement->calculateFgFromRates(dt, growthRates[0],growthRates[1],growthRates[2], currGF->getShearAngleRotationMatrix(), columnarFgIncrement, 0, currGF->zMin, currGF->zMax);
                    }
                }
                if (currGF->applyToLateralECM){
                    if (itElement->isECMMimimcingAtCircumference && !itElement->tissuePlacement == 0){ //do not grow the basal element twice
                        itElement->calculateFgFromRates(dt, growthRates[0],growthRates[1],growthRates[2], currGF->getShearAngleRotationMatrix(), columnarFgIncrement, 0, currGF->zMin, currGF->zMax);
                    }
                }
            }
            if (!thereIsExplicitECM || !itElement->isECMMimicing ){
                //There is either no explicit ECM definition, or the element is not ECM mimicing.
                //If there is explicit ECM, the basal elements should not grow, all others should proceed as usual
                //If there is no explicit ecm, then all should proceed as usual.
                if (currGF->applyToColumnarLayer){
                    itElement->calculateFgFromRates(dt, growthRates[0],growthRates[1],growthRates[2], currGF->getShearAngleRotationMatrix(), columnarFgIncrement, 0, currGF->zMin, currGF->zMax);
                }
                if (currGF->applyToPeripodialMembrane){
                    itElement->calculateFgFromRates(dt, growthRates[0],growthRates[1],growthRates[2], currGF->getShearAngleRotationMatrix(), peripodialFgIncrement, 1, currGF->zMin, currGF->zMax);
                }
            }
            itElement->updateGrowthIncrement(columnarFgIncrement,peripodialFgIncrement);
    	}
		gsl_matrix_free(columnarFgIncrement);
		gsl_matrix_free(peripodialFgIncrement);
	}
}

void Simulation::calculateGrowthRing(GrowthFunctionBase* currGF){
	if(currSimTimeSec >= currGF->initTime && currSimTimeSec < currGF->endTime ){
		float centre[2];
		currGF->getCentre(centre[0], centre[1]);
		float innerRadius = currGF->getInnerRadius();
		float outerRadius = currGF->getOuterRadius();
        std::array<double,3> maxValues = currGF->getGrowthRate();
		float innerRadius2 = innerRadius*innerRadius;
		float outerRadius2 = outerRadius*outerRadius;
		gsl_matrix* columnarFgIncrement = gsl_matrix_calloc(3,3);
		gsl_matrix* peripodialFgIncrement = gsl_matrix_calloc(3,3);
        for(const auto& itElement : Elements){
            std::array<double,3> Elementcentre = itElement->getCentre();
			//the distance is calculated in the x-y projection
			double d[2] = {centre[0] - Elementcentre[0], centre[1] - Elementcentre[1]};
			double dmag2 = d[0]*d[0] + d[1]*d[1];
			if (dmag2 > innerRadius2 && dmag2 < outerRadius2){
				//the element is within the growth zone.
				float distance = pow(dmag2,0.5);
				//calculating the growth rate: as a fraction increase within this time point
				double sf = (1.0 - (distance - innerRadius) / (outerRadius - innerRadius) );
				double growthscale[3] = {maxValues[0]*sf,maxValues[1]*sf,maxValues[2]*sf};
				gsl_matrix_set_identity(columnarFgIncrement);
				gsl_matrix_set_identity(peripodialFgIncrement);
                if (currGF->applyToBasalECM || currGF->applyToLateralECM){
                    if (currGF->applyToBasalECM){
                        if (itElement->isECMMimicing && itElement->tissuePlacement == 0){
                            itElement->calculateFgFromRates(dt, growthscale[0],growthscale[1],growthscale[2], currGF->getShearAngleRotationMatrix(), columnarFgIncrement, 0, currGF->zMin, currGF->zMax);
                        }
                    }
                    if (currGF->applyToLateralECM){
                        if (itElement->isECMMimimcingAtCircumference && !itElement->tissuePlacement == 0){ //do not grow the basal element twice
                            itElement->calculateFgFromRates(dt, growthscale[0],growthscale[1],growthscale[2], currGF->getShearAngleRotationMatrix(), columnarFgIncrement, 0, currGF->zMin, currGF->zMax);
                        }
                    }
                }
                if (!thereIsExplicitECM || !itElement->isECMMimicing ){
                    //There is either no explicit ECM definition, or the element is not ECM mimicing.
                    //If there is explicit ECM, the basal elements should not grow, all others should proceed as usual
                    //If there is no explicit ecm, then all should proceed as usual.
                    if (currGF->applyToColumnarLayer){
                        itElement->calculateFgFromRates(dt, growthscale[0],growthscale[1],growthscale[2], currGF->getShearAngleRotationMatrix(), columnarFgIncrement, 0, currGF->zMin, currGF->zMax);
                    }
                    if (currGF->applyToPeripodialMembrane){
                        itElement->calculateFgFromRates(dt, growthscale[0],growthscale[1],growthscale[2], currGF->getShearAngleRotationMatrix(), peripodialFgIncrement, 1, currGF->zMin, currGF->zMax);
                    }
                }
                itElement->updateGrowthIncrement(columnarFgIncrement,peripodialFgIncrement);
			}
		}
		gsl_matrix_free(columnarFgIncrement);
		gsl_matrix_free(peripodialFgIncrement);
	}
}

void Simulation::setUpECMMimicingElements(){
   /**
    * First all nodes are assigned as if they are owne explicitly by ECM mimicking elements (Node#allOwnersECMMimicing = true), then
    * the flags is cleared through the elemtns.
    */
    for(const auto& itNode : Nodes){
        itNode->allOwnersECMMimicing = true;
	}

    for(const auto& itElement : Elements){
        if ( itElement->tissuePlacement == 0 ){
            if(itElement->tissueType == 0){
                /**
                 * The basal elements (ShapeBase#tissuePlacement = 0) are ECM, their identity is set accordingly.
                 */
                itElement->setECMMimicing(true);
			}
		}
        else if ( itElement->tissuePlacement == 1 ){
            if(itElement->tissueType == 1){
                /**
                 * On the peripodial tissue, conversly apical elements  (ShapeBase#tissuePlacement = 1) are ECM.
                 */
                itElement->setECMMimicing(true);
			}
		}
        else if ( itElement->tissuePlacement == 2 && itElement->spansWholeTissue ){
            /**
             * If the elements are spanning the whole tissue, (ShapeBase#tissuePlacement = 2), then they are ECM. This
             * scenario simulates a layer of ECM only, without any cell material.
             */
            itElement->setECMMimicing(true);
		}
        if (!itElement->isECMMimicing){
            /**
             * At the end of element ECM status update, if the element is NOT part of the ECM, then its nodes are not
             * owned exclusively by ECM elements, therefore their flag Node#allOwnersECMMimicing is set false.
             */
			//the element is not ecm mimicking, its nodes will not be ecm mimicking:
            const vector<int>& nodeIds = itElement->getNodeIds();
            for (auto currentNodeId : nodeIds){
                Nodes[currentNodeId]->allOwnersECMMimicing = false;
			}
		}
	}
    /**
     * If there is explicit ECM, then all basal elements will be ECM, then I need another marker to identify cell material at he most
     * basal side, attached to ECM. This is done in Simulation#assigneElementsAtTheBorderOfECM.
     */
	assigneElementsAtTheBorderOfECM();
}


void Simulation::assigneElementsAtTheBorderOfECM(){
    /**
     * If there is explicit ECM, then all basal elements will be ECM, then I need another marker to identify cell material at he most
     * basal side, attached to ECM. If an element is ECM, all its nodes are recorded. Then all the elements owning this node, but are not
     * ECM themselves become new basal elements, and their status is recorded in ShapeBase#atBasalBorderOfECM.
     */
	vector<int> nodeListToCheck;
    for(const auto& itElement : Elements){
        if ( itElement->isECMMimicing  && !itElement->isECMMimimcingAtCircumference){ //I do not want to track the elements at the side borders
            int n = itElement->getNodeNumber();
			for (int i=0; i<n; ++i){
				bool alreadyRecorded = false;
				int nVector = nodeListToCheck.size();
				for (int j=0 ;j< nVector; ++j){
                    if (nodeListToCheck[j] ==itElement->NodeIds[i] ){
						alreadyRecorded = true;
						break;
					}
				}
				if (!alreadyRecorded){
                    nodeListToCheck.push_back(itElement->NodeIds[i]);
				}
			}
		}
	}
	int nVector = nodeListToCheck.size();
    for(const auto& itElement : Elements){
        if (!itElement->isECMMimicing && !itElement->atBasalBorderOfECM){
			for (int j=0 ;j< nVector; ++j){
                if (itElement->DoesPointBelogToMe(nodeListToCheck[j])){
                    itElement->atBasalBorderOfECM = true;
					break;
				}
			}
		}
	}
}

void Simulation::assigneElementsAtTheBorderOfActin(){
    /**
     * If there is explicit actin, then all apical elements will be actin rich, I want another marker to identify non-actin material
     * at the most apical (top)  side of the tissue.
     * If an element is actin rich, all its nodes are recorded. Then all the elements owning this node, but are not
     * actin-rich themselves become new basal elements, and their status is recorded in ShapeBase#atApicalBorderOfActin.
     */
    vector<int> nodeListToCheck;
    for(const auto& itElement : Elements){
        if ( itElement->isActinMimicing  ){
            int n = itElement->getNodeNumber();
			for (int i=0; i<n; ++i){
				bool alreadyRecorded = false;
				int nVector = nodeListToCheck.size();
				for (int j=0 ;j< nVector; ++j){
                    if (nodeListToCheck[j] == itElement->NodeIds[i] ){
						alreadyRecorded = true;
						break;
					}
				}
				if (!alreadyRecorded){
                    nodeListToCheck.push_back(itElement->NodeIds[i]);
				}
			}
		}
	}
	int nVector = nodeListToCheck.size();
    for(const auto& itElement : Elements){
        if (!itElement->isActinMimicing && !itElement->isECMMimicing && !itElement->atBasalBorderOfECM){
			for (int j=0 ;j< nVector; ++j){
                if (itElement->DoesPointBelogToMe(nodeListToCheck[j])){
                    itElement->atApicalBorderOfActin = true;
					break;
				}
			}
		}
	}
}

void Simulation::setUpActinMimicingElements(){
    /**
     * All elements on the apical side (ShapeBase#tissuePlacement = 1) are actin-rich elements.
     */
    for(const auto& itElement : Elements){
        if ( itElement->tissuePlacement == 1 ){
			//apical elements are mimicing actin:
            if (!itElement->isECMMimimcingAtCircumference){
                itElement->setActinMimicing(true);
			}
		}
        if ( itElement->tissuePlacement == 2 && itElement->spansWholeTissue ){
            /**
             * If the elements are spanning the whole tissue, (ShapeBase#tissuePlacement = 2), then they are actin-rich. This
             * scenario simulates a layer of actin mesh only, without any softer material.
             */
            if (!itElement->isECMMimimcingAtCircumference){
                itElement->setActinMimicing(true);
			}
		}
	}
    /**
     * If there is explicit actin, then all apical elements will be actin rich, I want another marker to identify non-actin material
     * at the most apical (top)  side of the tissue. This is done via Simulation#assigneElementsAtTheBorderOfActin.
     */
	assigneElementsAtTheBorderOfActin();
}

void Simulation::calculateGrowthGridBased(GrowthFunctionBase* currGF){
    int nGridX = currGF->getGridX();
	int nGridY = currGF->getGridY();
	if(currSimTimeSec >= currGF->initTime && currSimTimeSec < currGF->endTime ){
		gsl_matrix* columnarFgIncrement = gsl_matrix_calloc(3,3);
		gsl_matrix* peripodialFgIncrement = gsl_matrix_calloc(3,3);
        for(const auto& itElement : Elements){
			gsl_matrix_set_identity(columnarFgIncrement);
			gsl_matrix_set_identity(peripodialFgIncrement);

			int IndexX = 0.0, IndexY = 0.0;
			double FracX = 1.0,  FracY = 1.0;
			if (GridGrowthsPinnedOnInitialMesh){
                itElement->getInitialRelativePositionInTissueInGridIndex(nGridX, nGridY, IndexX, IndexY, FracX, FracY);
			}
			else{
                itElement->getRelativePositionInTissueInGridIndex(nGridX, nGridY, IndexX, IndexY, FracX, FracY);
			}
            if (currGF->applyToBasalECM || currGF->applyToLateralECM){
                if (currGF->applyToBasalECM){
                    if (itElement->isECMMimicing && itElement->tissuePlacement == 0){
                        itElement->calculateFgFromGridCorners(gridGrowthsInterpolationType, dt, currGF, columnarFgIncrement, 0, IndexX,  IndexY, FracX, FracY); 	//sourceTissue is 0 for columnar Layer
                    }
                }
                if (currGF->applyToLateralECM){
                    if (itElement->isECMMimimcingAtCircumference && !itElement->tissuePlacement == 0){ //do not grow the basal element twice
                        itElement->calculateFgFromGridCorners(gridGrowthsInterpolationType, dt, currGF, columnarFgIncrement, 0, IndexX,  IndexY, FracX, FracY); 	//sourceTissue is 0 for columnar Layer
                    }
                }
            }
            if (!thereIsExplicitECM || !itElement->isECMMimicing ){
                //There is either no explicit ECM definition, or the element is not ECM mimicing.
                //If there is explicit ECM, the basal elements should not grow, all others should proceed as usual
                //If there is no explicit ecm, then all should proceed as usual.
                if (currGF->applyToColumnarLayer){
                    itElement->calculateFgFromGridCorners(gridGrowthsInterpolationType, dt, currGF, columnarFgIncrement, 0, IndexX,  IndexY, FracX, FracY); 	//sourceTissue is 0 for columnar Layer

                }
                if (currGF->applyToPeripodialMembrane){
                    itElement->calculateFgFromGridCorners(gridGrowthsInterpolationType, dt, currGF, peripodialFgIncrement, 1, IndexX,  IndexY, FracX, FracY); 	//sourceTissue is 1 for peripodial membrane
                }
            }
            itElement->updateGrowthIncrement(columnarFgIncrement,peripodialFgIncrement);
		}
		gsl_matrix_free(columnarFgIncrement);
		gsl_matrix_free(peripodialFgIncrement);
	}
}

void Simulation::setStretch(){
    /**
     * The clap orientation is set in Simulation#DVClamp boolean, where if true, the tissue is stretched in x axis, and stretched in
     * y axis otherwise.
     */
    for (int i=0;i<3;i++){
		leftClampForces[i] = 0.0;
		rightClampForces[i] = 0.0;
	}
	vector <int> clampedNodeIds;
	double distance = 0;

    /**
     * The initial distance is calculated to set the stretched distance at maximum stretch Simulation#StretchStrain.
     */

	if (DVClamp){
		distance = fabs(Nodes[ventralTipIndex]->Position[0] - Nodes[dorsalTipIndex]->Position[0]);
		distanceIndex = 0; //if the clamp is on DV axis, then the direction of interest is x, index is 0;
	}
	else{
		distance = fabs(Nodes[anteriorTipIndex]->Position[1] - Nodes[posteriorTipIndex]->Position[1]);
		distanceIndex = 1; //if the clamp is on AP axis, then the direction of interest is y, index is 1.
	}
    /**
     * For each node falling under a clamp, their Ids are recorded to move with clamp, their movement is fixed otherwise.
     */
    for (const auto& currNode : Nodes){
        if (currNode->Position[distanceIndex]> StretchMax || currNode->Position[distanceIndex] < StretchMin){
            currNode->FixedPos[0]=1;
            currNode->FixedPos[1]=1;
            currNode->FixedPos[2]=1;
            clampedNodeIds.push_back(currNode->Id);
		}
	}
	setUpClampBorders(clampedNodeIds);
	//the distance that is to be moved:
	distance *= StretchStrain;
    /**
     * Then with the time steps of stretch identified by the user input, the stretch distance at each time step to move the clamped nodes is
     * obtained.
     */
	double StretchTimeSteps = (StretchEndTime - StretchInitialTime)/dt;
	StretchDistanceStep = 0.5* (distance / StretchTimeSteps);
}

void Simulation::setUpClampBorders(std::vector<int>& clampedNodeIds){
	int n = clampedNodeIds.size();
    for(const auto& itElement : Elements){
        const vector<int> nodeIds = itElement->getNodeIds();
		bool hasClampedNode = false;
		bool hasNonClampedNode = false;
		bool leftHandSide = false;
		vector <int> clampedBorderNodes;
        for (auto currNodeId : nodeIds){
			bool lastNodeWasClamped = false;
			for (int k=0; k<n; k++){
                if (currNodeId == clampedNodeIds[k]){
					hasClampedNode = true;
					lastNodeWasClamped = true;
					//check if node is recorded:
					bool alreadyRecorded = false;
					int nLeftClamp = leftClampBorder.size();
					for (int m=0; m<nLeftClamp; m++){
                        if (leftClampBorder[m] == currNodeId){
							alreadyRecorded = true;
							break;
						}
					}
					if (!alreadyRecorded){
						int nRightClamp = rightClampBorder.size();
						for (int m=0; m<nRightClamp; m++){
                            if (rightClampBorder[m] == currNodeId){
								alreadyRecorded = true;
								break;
							}
						}
					}
					if (!alreadyRecorded){
                        clampedBorderNodes.push_back(currNodeId);
					}
                    if (Nodes[currNodeId]->Position[distanceIndex] < StretchMin){
						leftHandSide = true;
					}
					break;
				}
			}
			if (lastNodeWasClamped == false){
				hasNonClampedNode = true;
			}
		}
		if (hasClampedNode && hasNonClampedNode){
            std::cout<<" Element "<<itElement->Id<<" is at the border"<<std::endl;
			int nClamp = clampedBorderNodes.size();
			if(leftHandSide){
                std::cout<<"Element is on left hand side"<<std::endl;
				for (int k=0; k<nClamp; k++){
					leftClampBorder.push_back(clampedBorderNodes[k]);
				}
			}
			else {
                std::cout<<"Element is on right hand side"<<std::endl;
				for (int k=0; k<nClamp; k++){
					rightClampBorder.push_back(clampedBorderNodes[k]);
				}
			}
		}
	}int nLeftClamp = leftClampBorder.size();
	for (int k=0; k<nLeftClamp; k++){
        std::cout<<"left clamp border nodes: "<<leftClampBorder[k]<<std::endl;
	}
	int nRightClamp = rightClampBorder.size();
	for (int k=0; k<nRightClamp; k++){
        std::cout<<"right clamp border nodes: "<<rightClampBorder[k]<<std::endl;
	}
}

void Simulation::setupPipetteExperiment(){
	pipetteInnerRadiusSq = pipetteInnerRadius*pipetteInnerRadius;
    effectLimitsInZ[0] = pipetteCentre[2] - pipetteDepth;
    effectLimitsInZ[1] = pipetteCentre[2] + pipetteDepth;
    /**
     * When setting the effectiveness boundaries of the suction pipette, the boundary inside the tube needs to be extended,
     * such that the nodes will not stop being suced in after they go inside a certain length into the tube.
     * This means increasing the upper boundary for apical suction, and reducing the lower boundary for basal suction.
     */
    if(ApicalSuction){
		effectLimitsInZ[1] += 1000;
	}
	else{
		effectLimitsInZ[0] -= 1000;
	}
    //Now I am sticking the other side of the tissue to a surface
    /**
     * If the other side of the tissue is attached to a surface, the nodesa re fixed accordingly.
     */
	if (TissueStuckOnGlassDuringPipetteAspiration){
		if (ApicalSuction){
            for (const auto& currNode : Nodes){
				//fix basal nodes of columnar layer:
                if(currNode->tissuePlacement == 0 && currNode->tissueType == 0) {
                    fixAllD(currNode.get(), false); //this is fixing with adhesives, should be a hard fix at all times
				}
			}
		}
		else{
             for (const auto& currNode : Nodes){
				//fix apical nodes and all peripodial membrane nodes:
                if(currNode->tissuePlacement == 1 || currNode->tissueType == 1) {
                    fixAllD(currNode.get(), false); //this is fixing with adhesives, should be a hard fix at all times
				}
			}
		}
	}
}

void Simulation::addPipetteForces(gsl_matrix* gExt){
	for (int i=0; i<nPipetteSuctionSteps; ++i){
		if (currSimTimeSec >= pipetteSuctionTimes[i]){
			SuctionPressure[2] = pipetteSuctionPressures[i];
		}
		else{
			break;
		}
	}
	int dim = 3;
    for (size_t i=0; i<nNodes; ++i){
        /**
         * Each node is checked agains the pipette boundaries. If the node is apical and the suction
         * is apical OR the suction is basal and the node is basal, the position is checked.
         */
		if( ( ApicalSuction && Nodes[i]->tissuePlacement == 1) || ( !ApicalSuction && Nodes[i]->tissuePlacement == 0)){
			if (Nodes[i]->Position[2]> effectLimitsInZ[0] &&  Nodes[i]->Position[2]< effectLimitsInZ[1]){
				double dx = pipetteCentre[0] - Nodes[i]->Position[0];
				double dy = pipetteCentre[1] - Nodes[i]->Position[1];
				double d2 = dx*dx+dy*dy ;
				if (d2 < pipetteInnerRadiusSq){
                    /**
                     * The suction pressure of the pipette is applied to the Node#zProjectedArea and added as an
                     * external force to the nodal force vector.
                     */
					double multiplier = 1.0;
					bool scalePressure = false;
					if(scalePressure){ multiplier = (1 - d2/pipetteInnerRadiusSq);}
					double LocalPressure[3] = { multiplier*SuctionPressure[0], multiplier*SuctionPressure[1],multiplier*SuctionPressure[2]};
					//node is within suciton range
					if(!Nodes[i]->FixedPos[0]){
						double value = gsl_matrix_get(gExt,i*dim,0);
						value +=LocalPressure[0]*Nodes[i]->zProjectedArea;
						gsl_matrix_set(gExt,i*dim,0,value);
					}
					if(!Nodes[i]->FixedPos[1]){
						double value = gsl_matrix_get(gExt,i*dim+1,0);
						value +=LocalPressure[1]*Nodes[i]->zProjectedArea;
						gsl_matrix_set(gExt,i*dim+1,0,value);
					}
					if(!Nodes[i]->FixedPos[2]){
						double value = gsl_matrix_get(gExt,i*dim+2,0);
						value +=LocalPressure[2]*Nodes[i]->zProjectedArea;
						gsl_matrix_set(gExt,i*dim+2,0,value);
					}
				}
			}
		}
	}
}

void Simulation::detectPackingToPipette(){
    /**
     * The current tip of the pipette is defined as the range starting Simulation#pipetteInnerRadius going for Simulation#pipetteThickness, and
     * the centre of the pipette tip is recorded in Simulation#pipetteCentre. Any node falling within this threshold on the corect surface
     * can potentially pack to the pipette.
     *
     */
    TransientZFixListForPipette.clear();
    std::array<double,2> pipLim {0.95*pipetteInnerRadius, 1.2*(pipetteInnerRadius+pipetteThickness)};
    std::array<double,2> pipLimSq {pipLim[0] *pipLim[0], pipLim[1]* pipLim[1]};
    double zPackThreshold = 2.0;
    std::array<double,2> zLimits{-1.0*zPackThreshold,zPackThreshold};
    for (size_t currId=0; currId<nNodes;currId++){
        if( ( ApicalSuction && Nodes[currId]->tissuePlacement == 1) || ( !ApicalSuction && Nodes[currId]->tissuePlacement == 0)){
             std::array<double,3> distanceToPipetCentre{Nodes[currId]->Position[0]-pipetteCentre[0], Nodes[currId]->Position[1]-pipetteCentre[1], Nodes[currId]->Position[2]-pipetteCentre[2]};
            if (distanceToPipetCentre[2] < zLimits[1] && distanceToPipetCentre[2]> zLimits[0]){
                double xydist2 = distanceToPipetCentre[0]*distanceToPipetCentre[0]+distanceToPipetCentre[1]*distanceToPipetCentre[1];
                if (xydist2>pipLimSq[0] && xydist2<pipLimSq[1]){
                    TransientZFixListForPipette.push_back(currId);
                }
            }
        }
    }
}

void Simulation::calculatePackingToPipetteForcesImplicit3D(){
    double zeroThreshold = 1E-3;
    double zPackThreshold = 2.0;
    /**
     * See the documentation for Simulation#calculatePackingForcesImplicit3D. The methodology is the same,
     * with pipeete tip centre used instead of the second node.
     */
    size_t n = TransientZFixListForPipette.size();
    for(size_t node_counter = 0 ; node_counter<n; ++node_counter){
        int id0 = TransientZFixListForPipette[node_counter];
        double multiplier = packingMultiplier;
        double sigmoidSaturation = sigmoidSaturationForPacking;
        double distanceShiftScale  = 0.6*zPackThreshold;
        double d = fabs(Nodes[id0]->Position[2] - pipetteCentre[2]);
        double averageMass = Nodes[id0]->mass;
        double normalZ= -1.0;//normal from pipette to node
        //if node is apical, normal points to -1
        if (Nodes[node_counter]->tissuePlacement == 1){
            normalZ = 1.0;
        }
        if (d < zeroThreshold){
            d = zeroThreshold;
        }
        double shifted_d = d - distanceShiftScale; //shift the distance by the selected percentage,
        double sigmoid = 1 / (1 + exp(sigmoidSaturation/ zPackThreshold * 2.0 * shifted_d));
        double Fmagnitude = multiplier * averageMass * sigmoid;
        double Fz = Fmagnitude*normalZ;
        PackingForces[id0][2] += Fz;
        if (std::isnan(Fz)){
              std::cerr<<" packing force from pipette is nan for node "<<TransientZFixListForPipette[node_counter]<<std::endl;
        }
    }
}

void Simulation::calculatePackingToPipetteJacobian3D(gsl_matrix* K){
    double zeroThreshold = 1E-3;
    double zPackThreshold = 2.0;
    /**
     * See the documentation for Simulation#calculatePackingJacobian3D. The methodology is the same,
     * with pipette tip centre used instead of the second node.
     */
    size_t n = TransientZFixListForPipette.size();
     for(size_t node_counter = 0 ; node_counter<n; ++node_counter){
         int id0 = TransientZFixListForPipette[node_counter];
        double multiplier = packingMultiplier;
        double sigmoidSaturation = sigmoidSaturationForPacking;
        double distanceShiftScale  = 0.6*zPackThreshold;
        double d = fabs(Nodes[id0]->Position[2] - pipetteCentre[2]);
        double averageMass = Nodes[id0]->mass;
        double normalZ= -1.0;//normal from pipette to node
        //if node is apical, normal points to -1
        if (Nodes[node_counter]->tissuePlacement == 1){
            normalZ = 1.0;
        }
        if (d < zeroThreshold){
            d = zeroThreshold;
        }
        double shifted_d = d - distanceShiftScale; //shift the distance by the selected percentage,
        double sigmoid = 1 / (1 + exp(sigmoidSaturation/ zPackThreshold * 2.0 * shifted_d));
        double dSdXi = (-sigmoidSaturation/ zPackThreshold * 2.0) * sigmoid * (1-sigmoid);
        gsl_matrix* normalMat = gsl_matrix_calloc(3,1);
        gsl_matrix_set(normalMat, 0,0,0.0);
        gsl_matrix_set(normalMat, 1,0,0.0);
        gsl_matrix_set(normalMat, 2,0,normalZ);
        gsl_matrix* norm_normT = gsl_matrix_calloc(3,3);
        gsl_blas_dgemm (CblasNoTrans, CblasTrans,1.0, normalMat, normalMat, 0.0, norm_normT);
        gsl_matrix* dFidXi = gsl_matrix_calloc(3,3);
        gsl_matrix_set_identity(dFidXi);

        gsl_matrix_sub(dFidXi,norm_normT);  //(I - norm*norm^T)
        gsl_matrix_scale(dFidXi,sigmoid/d);	//sigmoid/d * (I - norm*norm^T)
        gsl_matrix_scale(norm_normT,dSdXi); // norm*norm^T * dSdXi
        gsl_matrix_add(dFidXi,norm_normT); // sigmoid/d * (I - norm*norm^T) + norm*norm^T * dSdXi
        gsl_matrix_scale(dFidXi,multiplier * averageMass); //scale by mass and multiplier when needed
        for (size_t i=0; i<3; ++i){
            for (size_t j=0; j<3; ++j){
                double derivativeij = gsl_matrix_get(dFidXi,i,j);
                addValueToMatrix(K,3*id0+i,3*id0+j,-derivativeij);
            }
        }
        gsl_matrix_free(normalMat);
        gsl_matrix_free(dFidXi);
        gsl_matrix_free(norm_normT);
    }
}

void Simulation::updateElementVolumesAndTissuePlacements(){
    for (auto &itElement : Elements){
        std::cout<<"updating element: "<<itElement->Id<<std::endl;
        itElement->updateElementVolumesAndTissuePlacementsForSave(Nodes);
	}
}

void Simulation::updateElasticPropertiesForAllNodes(){
    /**
     * Elastic property update function ShapeBase#updateElasticProperties is called on all elements.
     */
    for (auto &itElement : Elements){
        itElement->updateElasticProperties();
	}
}

void Simulation::clearNodeMassLists(){
    for (const auto& currNode : Nodes){
        currNode->connectedElementIds.size();
        currNode->connectedElementIds.clear();
        currNode->connectedElementWeights.clear();
        currNode->mass=0.0;
	}
}

void Simulation::setupYsymmetricity(){
    /**
     * If a node is at the symmetriciy zone, its degree of freedom is fixed on the symmetricity axis. There will be elements of exactly equal strengt
     * actin on the node and cancelling the forces in the symmetric axis.
     */
	double yLimPos = 0.1;
	double yLimNeg = (-1.0)*yLimPos;
    for (size_t i=0; i<nNodes; ++i){
		if (Nodes[i]->Position[1]< yLimPos){
			if (Nodes[i]->Position[1] > yLimNeg){
				Nodes[i]->atSymmetricityBorder = true;
                fixY(Nodes[i].get(),false); //this is for symmetricity, the fixing has to be hard fixing, not with external viscosity under any condition
			}
		}
	}		
}

void Simulation::setupXsymmetricity(){
    /**
     * If a node is at the symmetriciy zone, its degree of freedom is fixed on the symmetricity axis. There will be elements of exactly equal strengt
     * actin on the node and cancelling the forces in the symmetric axis.
     */
    double xLimPos = 0.1;
	double xLimNeg = (-1.0)*xLimPos;
    for (size_t i=0; i<nNodes; ++i){
		if (Nodes[i]->Position[0]< xLimPos){
			if (Nodes[i]->Position[0] > xLimNeg){
				Nodes[i]->atSymmetricityBorder = true;
                fixX(Nodes[i].get(),false); //this is for symmetricity, the fixing has to be hard fixing, not with external viscosity under any condition
			}
		}
	}
}

void Simulation::clearScaleingDueToApikobasalRedistribution(){
    //not using range based loops here to ensure openMP comaptibility
    #pragma omp parallel for
    for(std::vector<std::unique_ptr<ShapeBase>>::iterator itElement = Elements.begin(); itElement <Elements.end(); ++itElement){
        (*itElement)->thereIsGrowthRedistribution = false;
        (*itElement)->growthRedistributionScale = 0.0;
    }
}
