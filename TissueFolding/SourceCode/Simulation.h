#ifndef Simulation_H
#define Simulation_H

#include <stdio.h>
#include <iostream>
#include <ostream>
#include <fstream>
#include <vector>

#include "ModelInputObject.h"
#include "ShapeBase.h"
#include "Node.h"
#include "GrowthFunctionBase.h"
#include "GrowthFunctionTypes.h"
#include "NewtonRaphsonSolver.h"

//#include <omp.h>

class ModelInputObject;
using namespace std;

class Simulation : public std::enable_shared_from_this<Simulation> {
private:
    int currElementId;                                       ///< The next element ID to be used in a new element
    std::shared_ptr<ModelInputObject> ModInp;                ///< The ModelInputObject pointer, used to read the user input file.
    std::ofstream saveFileMesh;                              ///< Output file to save the mesh coordinate data.
    std::ofstream saveFileTensionCompression;                ///< Output file to save the strains.
    std::ofstream saveFileGrowth;                            ///< Output file to save the growth.
    std::ofstream saveFileGrowthRate;                        ///< Output file to save the growth rates.
    std::ofstream saveFilePhysicalProp;                      ///< Output file to save the physical properties of elements.
    std::ofstream saveFileSpecificType;                      ///< Output file to save the specific node types, such as the ECM, actin, marker ellipses.
    std::ofstream saveFileGrowthRedistribution;              ///< Output file to save the volume redistribution.
    std::ofstream saveFileNodeBinding;                       ///< Output file to save node binding.
    std::ofstream saveFileSimulationSummary;                 ///< Output file to save the simulation inputs summary.
    std::ifstream saveFileToDisplayMesh;                     ///< Input file to save the mesh coordinate data.
    std::ifstream inputMeshFile;                             ///< Input file to save the strains.
    std::ifstream saveFileToDisplayTenComp;                  ///< Input file to save the strains.
    std::ifstream saveFileToDisplayGrowth;                   ///< Input file to save the growth.
    std::ifstream saveFileToDisplayGrowthRate;              ///< Input file to save the growth rates.
    std::ifstream saveFileToDisplayPhysicalProp;            ///< Input file to save the physical properties of elements.
    std::ifstream saveFileToDisplaySpecificNodeTypes;       ///< Input file to save the specific node types, such as the ECM, actin, marker ellipses.
    std::ifstream saveFileToDisplaySimSum;                  ///< Input file to save the simulation inputs summary.
    std::ifstream saveFileToDisplayGrowthRedistribution;    ///< Input file to save the volume redistribution.
    std::ifstream saveFileToDisplayNodeBinding;             ///< Input file to save node binding.
    bool ContinueFromSave;                                  ///< Model input boolean stating if the simualitons is continuing from stretch.

    bool TensionCompressionSaved;                           ///< Boolean stating if the strains are saved.
    bool GrowthSaved;                                       ///< Boolean stating if the growth is saved.
    bool GrowthRateSaved;                                   ///< Boolean stating if the growth rates are saved.
    bool ForcesSaved;                                       ///< Boolean stating if the forces are saved.
    bool physicalPropertiesSaved;                           ///< Boolean stating if the pshyscal properties of the elements are saved.
    bool PackingSaved;                                      ///< Boolean stating if the packing is  saved.
    bool growthRedistributionSaved;                         ///< Boolean stating if the volume redistribution is saved.
    bool nodeBindingSaved;                                  ///< Boolean stating if the node binding is saved.
    bool specificElementTypesRecorded;                      ///< Boolean stating if the specific node types, such as the ECM, actin, marker ellipses are saved.
    bool collapseAndAdhesionSaved;                          ///< Boolean stating if the strains are saved.
    int	 nCircumferencialNodes;                             ///< Boolean stating if the collapse and ahdesion are saved.
    int  dorsalTipIndex;                                    ///< The Node#Id for the dorsal tip (-x) of the tissue on xy plane.
    int  ventralTipIndex;                                   ///< The Node#Id for the ventral tip (+x) of the tissue on xy plane.
    int  anteriorTipIndex;                                  ///< The Node#Id for the dorsal tip (+y) of the tissue on xy plane.
    int  posteriorTipIndex;                                 ///< The Node#Id for the dorsal tip (-y) of the tissue on xy plane.
    double StretchDistanceStep;                             ///< The step of the stretcher movement in microns per time step.
    double boundingBoxSize[3];                              ///< The size of the bounding  box of the tissue in 3D, [x][y][z] in microns.
    int growthRotationUpdateFrequency;                      ///< The time step frequency to update the planar rotations of the elements.
    bool checkedForCollapsedNodesOnFoldingOnce;             ///< The flag to state the folding regions node binding is checked once, necessary to safety check at the beginning of simulations continuing from save.
    bool boundLateralElements;                              ///< Te flog to state if the lateral elements are bound in their degrees of freesdom to avoid boundary buckling.

    bool readModeOfSim(int& i, int argc, char **argv);              ///< User input reading function, reading mode if the simulations, can be DisplaySave, SimulationOnTheGo, ContinueFromSave. input tag is -mode
    bool readParameters(int& i, int argc, char **argv);             ///< User input reading function, reading modelinput file, the input tag is -i, should give path to modelinput file.
    bool readOutputDirectory(int& i, int argc, char **argv);        ///< User input reading function, reading output directory, the input tag is -od, shoud give path to output file
    bool readSaveDirectoryToDisplay(int& i, int argc, char **argv); ///< User input reading function, reading the directory to obtain the save files to continue the simulation from, input tag is -od
    bool openFilesToDisplay();                                      ///< This function opens the files to display a saved simulation.
    bool readSystemSummaryFromSave();                               ///< This function reads the simulation parameters from the summary file, to display or continue simulation from saved simulation.
    bool readSpecificNodeTypesFromSave();                           ///< This function reads the specific node types, such as the ECM, actin, marker ellipses from file.
    void initiateNodesFromSave();                                   ///< This function reads coordinates and initiates nodes from save file.
    bool readNodeDataToContinueFromSave();                          ///< This function reads nodal coordinates through the simulation to continue simulation from saved simulation.
    void initiateNodesFromMeshInput();                              ///< This function reads coordinates and initiates nodes from a mesh file.
    void initiateElementsFromSave();                                ///< This function reads shape properties, nodal data, reference shape, and initiates elements from save file.
    bool readElementDataToContinueFromSave();                       ///< This function reads shape properties, nodal data, reference shape through the simulation to continue simulation from saved simulation.
    void initiateElementsFromMeshInput();                           ///< This function reads shape properties, nodal data, reference shape, and initiates elements from mesh file.
    void initiatePrismFromSave();                                   ///< This function initiates a new prism and stores a unique pointer to it in Simulation#Elements from save data.
    bool readShapeData(int i);                                      ///< This function reads the element data.
    void initiatePrismFromMeshInput();                              ///< This function initiates a new prism and stores a unique pointer to it in Simulation#Elements from mesh input data.
    void initiatePrismFromSaveForUpdate(int k);                     ///< This function initiates a new prism and stores a unique pointer to it in Simulation#Elements from continuing simulation data.
    void removeElementFromEndOfList();                              ///< This function remves one element form the end of the vector of elements, Simulation#Elements
    void updateNodeNumberFromSave();                                ///< This function updates the node number of the sytem according to data read from save file.
    void updateNodePositionsFromSave();                             ///< This function updates the node positions of the sytem according to data read from save file.
    void updateElementStatesFromSave();                             ///< This function updates the element stated of the sytem according to data read from save file.
    void updateTensionCompressionFromSave();                        ///< This function updates system strains from save file to display.
    void updateGrowthFromSave();                                    ///< This function updates system growths from save file to display.
    void updateGrowthRateFromSave();                                ///< This function updates system growth rates from save file to display.
    void updatePhysicalPropFromSave();                              ///< This function updates system physical properties from save file to display.
    void updateNodeBindingFromSave();                               ///< This function updates node binding information from save file to display.
    void updateGrowthRedistributionFromSave();                      ///< This function updates volume redistribution from save file to display.
    void readTensionCompressionToContinueFromSave();                ///< This function updates system strains from save file to continue from save.
    void readGrowthToContinueFromSave();                            ///< This function updates system growths from save file to continue from save.
    void readGrowthRateToContinueFromSave();                        ///< This function updates system gorwth rates from save file to continue from save.
    void readPhysicalPropToContinueFromSave();                      ///< This function updates system physical properties from save file to continue from save.
    void readGrowthRedistributionToContinueFromSave();              ///< This function updates system volume redistribution from save file to continue from save.
    void readNodeBindingToContinueFromSave();                       ///< This function updates system node binding data from save file to continue from save.
    void readCollapseAndAdhesionToContinueFromSave();               ///< This function updates adhered and collapsed nodes from save file to continue from save.
    bool readFinalSimulationStep();                                 ///< This function reads a saved simulation all at once up to the end point.
    void reInitiateSystemForces();                                  ///< This function re-initiates arrays to store system forces in the scenario where system size may change (remeshing)
    bool checkInputConsistency();                                   ///< This function checks the consistency of model inputs.
    void setDefaultParameters();                                    ///< This function sets the default parameters to the simulation parameters.
    bool openFiles();                                               ///< This function opens the save files.
    void initiateSystemForces();                                    ///< This function initiates system force arrays.
    bool initiateMesh(int MeshType, int Row, int Column, float SideLength, float zHeight);      ///< This function initiates a mesh from a selected number of rows and columns of elements.
    void initiateNodesByRowAndColumn(size_t Row, size_t Column,  float SideLength, float zHeight); ///< This function initiates nodes of the sytem for a mesh definition by number of raws and columns of elements
    void initiateElementsByRowAndColumn(size_t Row, size_t Column); ///< This initiates elements for a mesh definition by number of row and colmn of elements.
    bool initiateMesh(int MeshType);                                ///< This function initiates the simulation mesh depending on the input mesh type.
    void readInTissueWeights();                                     ///< This function reads in the tissue type weight of elements, peripodialness is recorded.
    bool checkIfTissueWeightsRecorded();                            ///< This function checks if the tissue type weights are recorded.
    void clearCircumferenceDataFromSymmetricityLine();              ///< This function clears the circumference identity from nodes and elements that are at a symmetricity boundary.
    void getAverageSideLength(double& periAverageSideLength, double& colAverageSideLength); ///< This function calculates the average side length of elements, necessary for packing/binding/collapse thresholds.
    int getElementArrayIndexFromId(int currId);                     ///< This function will return the ShapeBase#Id of the element from its index on the vector Simulation#Elements
    int getNodeArrayIndexFromId(int currId);                        ///< This function will return the Node#Id of the node from its index on the vector Simulation#Nodes
    bool isColumnarLayer3D();                                       ///< This function will check ifthe input mesh tissue is 3D.
    bool checkIfThereIsPeripodialMembrane();                        ///< This function will check ifthe input mesh contains a peripodial membrane.
    void setLinkerCircumference();                                  ///< This function will set the circumference mid line band of linker elements.
    bool calculateTissueHeight();                                   ///< This funciton will calculate the tissue height of the input mesh.
    void assignInitialZPositions();                                 ///< This function sets the initial relative z height of the elements on the tissue height.
    std::array<double,3> calculateNewNodePosForPeripodialNodeAddition(int nodeId0, int nodeId1, double sideThickness);  ///< This function calculates the position of a new node to construct a new element at tissue border, for peripodial membrane addition.
    std::array<double,3> calculateNewNodePosForPeripodialNodeAddition(int nodeId0, int nodeId1, int nodeId2, double sideThickness); ///< This function calculates the position of a new node to construct a new element at tissue border, for peripodial membrane addition.

    void checkForNodeFixing();                                      ///< This function checks for fixed nodes with the model input parametrs.
    void checkForNodeBinding();                                     ///< This function checks for node binding as boundary condition with the model input parametrs.
    bool bindCircumferenceXY();                                     ///< This function binds the x&y positional degrees of freedom of the circumference nodes to the base column on each column of nodes. This avoids bending of the boundary.
    bool areNodesToCollapseOnLateralECM(int slaveNodeId, int masterNodeId); ///< This funciton check if the nodes assigned to collapse are on the lateral ECM, the nodeas are bound but not collapsed at the side ECM.
    bool checkEdgeLenghtsForBindingPotentiallyUnstableElements();   ///< This function checks the eedge lengths of elemetns to collapse too small edges, avoiding element flipping.
    bool areNodesOnNeighbouingElements(int masterNoeId, int slaveNodeId);   ///< This function checks if the two nodes with the input IDs are on neightbouring elements.
    void manualAdhesion(int masterNodeId,int slaveNodeId);          ///< This function forces adhesion between the two input nodes.
    double distanceSqBetweenNodes(int id0, int id1);                ///< Helper funciton calculates the square of the distance between the two nodes with input node IDs.
    bool isAdhesionAllowed(int masterNodeId, int slaveNodeId);      ///< This function checks if the two nodes can adhere.
    bool checkForElementFlippingUponNodeCollapse(std::vector<int> &newCollapseList, std::array<double,3>  avrPos); ///< This function checks all elements on the potential collapse list against flipping.
    bool adhereNodes();                                             ///< This function adheres the nodes on the proximity lists if there is no rule against the adhesion.
    void updatePositionsOfNodesCollapsingInStages();                ///< This function brings the adhred node couples closer to each other in stages, to avoid too large movements and causing instabilities ata  single step.
    void induceClones();                                            ///< This function induces the growth mutant clones as specified in the model input file
    void assignPhysicalParameters();                                ///< This function assigns the physical properties of the elemetns.
    void checkForZeroExternalViscosity();                           ///< Checking if the system has zero external viscosity in any degree of freedon (x,y,z), and fixing node degrees of freedom as necessary to converge to one unique solution.
    void calculateShapeFunctionDerivatives();                       ///< Calculates the shape function derivatives.
    void assignNodeMasses();                                        ///< This function assigns masses to nodes.
    void updateNodeMasses();                                        ///< This function updates masses of the nodes.
    void assignElementalSurfaceAreaIndices();                       ///< This function assigns the indices of externally exposed surfaces in elemetns.
    void updateNodeViscositySurfaces();                             ///< This function calculates the surfaces feeling viscosity ojn teh elements and assign the surfaces to constructing nodes.
    void updateElementToConnectedNodes(const std::vector <std::unique_ptr<Node>>& Nodes);   ///< This function updates the connectivity of the elemetns.
    void assignConnectedElementsAndWeightsToNodes();                ///< This function assigns the wighted contribution of elemental properties on the nodes.
    void fixAllD(int i, bool fixWithViscosity);                     ///< Fixes 3D degrees of freedom of node, can be rigid fixing or fixing with increased viscosity.
    void fixAllD(Node* currNode, bool fixWithViscosity);            ///< Fixes 3D degrees of freedom of node, can be rigid fixing or fixing with increased viscosity.
    void fixX(int i, bool fixWithViscosity);                        ///< Fixes x degree of freedom of node, can be rigid fixing or fixing with increased viscosity.
    void fixX(Node* currNode, bool fixWithViscosity);               ///< Fixes x degree of freedom of node, can be rigid fixing or fixing with increased viscosity.
    void fixY(int i, bool fixWithViscosity);                        ///< Fixes y degree of freedom of node, can be rigid fixing or fixing with increased viscosity.
    void fixY(Node* currNode, bool fixWithViscosity);               ///< Fixes y degree of freedom of node, can be rigid fixing or fixing with increased viscosity.
    void fixZ(int i, bool fixWithViscosity);                        ///< Fixes z degree of freedom of node, can be rigid fixing or fixing with increased viscosity.
    void fixZ(Node* currNode, bool fixWithViscosity);               ///< Fixes z degree of freedom of node, can be rigid fixing or fixing with increased viscosity.
    void zeroForcesOnNode(size_t i);                                   ///< Sets the forces on the node to zero.
    void processDisplayDataAndSave();                               ///< This function brings all displayable parameters up to date, for correct display purposes.
    void saveStep();                                                ///< Save one step of simulation
    void writeSimulationSummary();                                  ///< Write the summary of the simulation into Save_Summary in the output directory specified by the user input
    void writeSpecificNodeTypes();                                  ///< Write the specific node types of the simulation to the output directory specified by the user input
    void writeMeshFileSummary();                                    ///< Write the summary of the simulation into Save_Summary in the output directory specified by the user input
    void writeGrowthRatesSummary();                                 ///< Write the summary of the mesh file to the simulation summary file, Save_Summary in the output directory specified by the user input
    void writeECMSummary();                                         ///< Write the summary of the ECM type into SaveFile_Summary in the output directory specified by the user input
    void writeECMProperties();                                      ///< Write the summary of the ECM physical properties into Save_Summary in the output directory specified by the user input
    void writeActinSummary();                                       ///< Write the summary of the actin type into SaveFile_Summary in the output directory specified by the user input
    void writeExperimentalSummary();                                ///< Write the summary of the experiments (stretcher, pipette aspiration) into SaveFile_Summary in the output directory specified by the user input
    void writePipetteSumary();                                      ///< Write the summary of the pipette aspiration experiment into Save_Summary in the output directory specified by the user input
    void writeSaveFileStepHeader();                                 ///< Write the time step header, marking beginning of time step, at each saved frame.
    void writeNodes();                                              ///< Write the nodal information to save file at each time step.
    void writeElements();                                           ///< Write the elemental information to save file at each time step.
    void writeSaveFileStepFooter();                                 ///< Write the time step footer, marking end of time step, at each saved frame.
    void writeTensionCompression();                                 ///< Write the strains of the simulation to the output directory specified by the user input
    void writeGrowth();                                             ///< Write the growths of the simulation to the output directory specified by the user input.
    void writePhysicalProp();                                       ///< Write the physical properties of the simulation to the output directory specified by the user input
    void writeNodeBinding();                                        ///< Write the node binding information of the simulation to the output directory specified by the user input
    void writeGrowthRedistribution();                               ///< Write the volume redistribution of the simulation to the output directory specified by the user input
    void calculateGrowth();                                         ///< Calculates the growths of each element from input growth functions as specified in the model input file.
    void calculateShapeChange();                                    ///< Calculates the shape change of each element from input growth functions as specified in the model input file.
    void cleanUpGrowthRates();                                      ///< Set growths increments of each element to zero at the beginning of each time step.
    void assignIfElementsAreInsideEllipseBands();                   ///< Assigns if elemetns are inside marker ellipses.
    void checkForPinningPositionsUpdate();                          ///< Checks if the relative positions of the elemetns for growth rate readouts are updated.
    void updateRelativePositionsToApicalPositioning();              ///< Update the relative positions of all elemetns to their apical counterpart.
    void updatePinningPositions();                                  ///< Update the relative positions of teh elements, and update the pinning reference frame.
    void updateGrowthRotationMatrices();                            ///< Update the growth rotation matrices of all elements, to eliminte the rigid body rotation around z axis in growth orientation calculation.
    void cleanUpShapeChangeRates();                                 ///< Set shape change increments of each element to zero at the beginning of each time step.
    void calculateGrowthUniform(GrowthFunctionBase* currGF);        ///< Calculates growth increments applied on the elements from a uniform growth function.
    void calculateGrowthRing(GrowthFunctionBase* currGF);           ///< Calculates growth increments applied on the elements from a ring shaped growth function.
    void calculateGrowthGridBased(GrowthFunctionBase* currGF);      ///< Calculates growth increments applied on the elements from a input growth map grid.
    void calculateShapeChangeUniform (GrowthFunctionBase* currSCF); ///< Calculates shape change increments applied on the elements from a uniform growth function.
    void calculateShapeChangeMarkerEllipseBased (GrowthFunctionBase* currSCF); ///< Calculates shape change increments applied on the elements from a shape function based on spacial marking of certain elemetns. Marking can be emergent at fold initiation (curvature based) or input in modelinpiut file.
    void setStretch();                                              ///< Set up the stretcher experiment.
    void setUpClampBorders(std::vector<int>& clampedNodeIds);       ///< Set up the clamping borders of the stretcher device.
    void setupPipetteExperiment();                                  ///< Set up the pipette aspiration experiment.
    void addPipetteForces(gsl_matrix *gExt);                        ///< Add suction forces of the pipette.
    void setUpAFM();                                                ///< Set up the Atomic Force Microscopy experiment.
    void moveAFMBead();                                             ///< Push the bead of the Atomic Force Microscopy setup towrds the tissue incrementally.
    void fillInNodeNeighbourhood();                                 ///< Fill in the neighbourhood information on the nodes.
    void setBasalNeighboursForApicalElements();                     ///< Set basal neighbours of all apical elements.
    void conserveColumnVolume();                                    ///< Conserve the volume of the column of elements in the mesh, alowing transfer of material from more compressed elemetns to their more relaxed neightbours.
    void fillInElementColumnLists();                                ///< Prepare arrays of element IDs marking element residing on the same colum on the tissue.
    void updateElementVolumesAndTissuePlacements();                 ///< Update volumes and tissue placements of the elements from save.
    void updateElasticPropertiesForAllNodes();                      ///< When a simulation is read from save, the modified physical properties will be read.  The elasticity calculation parameters (lambda, mu, and related matrices) should be updated after reading the files. This function is called upon finishing reading a save file, to ensure the update.

    void clearNodeMassLists();                                      ///< Set all assigned masses of the nodes to zero.
    void manualPerturbationToInitialSetup(bool deform, bool rotate); ///< Specify manual perturbation at the beginning of simulation. such as stretch, fixing or tilting of the tissue.x
    void addSoftPeriphery(std::vector<double>& fractions);          ///< Add  soft periphery to the tissue with the geometric prperties specified in model input file.
    void setupYsymmetricity();                                      ///< Set up symettric mesh around y-axis.
    void setupXsymmetricity();                                      ///< Set up symettric mesh around x-axis.
    void setUpECMMimicingElements();                                ///< Set up ECM mimicking elements.
    void setUpActinMimicingElements();                              ///< Set up actin elements.
    void assigneElementsAtTheBorderOfECM();                         ///< Mark elements at the border of the ECM (bottom & sides - basal & lateral).
    void assigneElementsAtTheBorderOfActin();                       ///< Mark elements at the border of the actin (top - apical).
    void clearScaleingDueToApikobasalRedistribution();              ///< Set volume redistribution scaling to zero.

public:
    std::ofstream outputFile;                                   ///< The output file to direct eroors and wornings, progress summary such as last converged time step.
    bool displayIsOn;                                           ///< The simulation is displayed.
    bool DisplaySave;                                           ///< The package is displaying a saved simulation from.
    bool reachedEndOfSaveFile;                                  ///< The simulation has reached the end of the read of save file
    float dt;                                                   ///< The time step increment in seconds, such that Simulation#currSimTimeSec = Simulation#dt * Simulation#timestep.
    int timestep;                                               ///< The time step count, such that  Simulation#currSimTimeSec = Simulation#dt * Simulation#timestep.
    double currSimTimeSec;                                      ///< The current time of the simulation in time seconds, obtained as Simulation#currSimTimeSec = Simulation#dt * Simulation#timestep.
    double SimLength;                                           ///< The desiered length of the current simulation in seconds
    std::string saveDirectory;                                  ///< The full path to the directory where the results of the simulation are saved.
    std::string saveScreenshotsDirectory;                       ///< The full path to the directory where the screenshots of the simulation are saved.
    std::string saveDirectoryToDisplayString;                   ///< The full path to the directory from which the saved simulation is being read.
    std::string inputMeshFileName;                              ///< The path to the input mesh file.
    std::string outputFileString;                               ///< The full path to the output file (Simulation#outputFile) directory.
    bool saveImages;                                            ///< Boolean stating if the screenshot images should be saved.
    bool saveData;                                              ///< Boolean stating if the simulation data should be saved.
    //std::string name_saveFile;                                ///< The name of the save file
    int imageSaveInterval;                                      ///< The interval to save screenshot images, in time step counts.
    int dataSaveInterval;                                       ///< The interval to save simalation data, in time step counts.
    double EApical;                                             ///< Young's modulus of the apical section of tissue
    double EBasal;                                              ///< Young's modulus of the basal section of tissue
    double EMid;                                                ///< Young's modulus of the mid-line of the tissue
    double EColumnarECM;                                        ///< Young's modulus of the columnar extracellular matrix (ECM)
    double EPeripodialECM;                                      ///< Young's modulus of the peripodial extracellular matrix (ECM)
    double poisson;                                             ///< Poisson ratio of the tissue.
    double discProperApicalViscosity;                           ///< The apical internal viscosity of the columnar layer.
    double discProperBasalViscosity;                            ///< The basal internal viscosity of the columnar layer.
    double discProperMidlineViscosity;                          ///< The mid-line internal viscosity of the columnar layer.
    std::array<int,4> noiseOnPysProp;                           ///< The desired noise percentage on the physical properties [Simulation#EApical] [Simulation#EBasal][Simulation#EMidLine][Simulation#poisson]
    bool zeroExternalViscosity[3];                              ///< The boolean stating if there is zero external viscosity on any of the 3 dimensions
    bool extendExternalViscosityToInnerTissue;                  ///< The boolean stating if the external viscosity will be extended to the inner layers.
    double externalViscosityDPApical;                           ///< The apical external viscosity of the columnar layer.
    double externalViscosityDPBasal;                            ///< The basal external viscosity of the columnar layer.
    double externalViscosityPMApical;                           ///< The apical external viscosity of the peripodial layer.
    double externalViscosityPMBasal;                            ///< The basal external viscosity of the peripodial layer.
    double externalViscosityLZApical;                           ///< The apical external viscosity of the linker zone between columnar and peripodial layer.
    double externalViscosityLZBasal;                            ///< The basal external viscosity of the linker zone between columnar and peripodial layer.
    int MeshType;                                               ///< The mesh type of the tissue; 2: a hexagonal mesh is generated with Simulation#Row and Simulation#Column number of elemetns; 4: an input file is specified.
    int Row;                                                    ///< The number of rows of elemetns for initiating a hexagonal mesh for Simulation#MeshType = 2.
    int Column;                                                 ///< The number of columns of elemetns for initiating a hexagonal mesh for Simulation#MeshType = 2.
    float SideLength;                                           ///< The side length of an element for initiating a hexagonal mesh for Simulation#MeshType = 2.
    float zHeight;                                              ///< The z height of an element for initiating a hexagonal mesh for Simulation#MeshType = 2.
    bool ApicalNodeFixWithExternalViscosity;                    ///< The boolean stating if the apical nodes should be fixed with external viscosity;
    bool BasalNodeFixWithExternalViscosity;                     ///< The boolean stating if the basal nodes should be fixed with external viscosity;
    bool CircumferentialNodeFixWithHighExternalViscosity[5];    ///< The boolean array stating if the circumferential nodes should be fixed with external viscosity ([apical][basal][linker apcial][linker basal] [all])
    bool NotumNodeFixWithExternalViscosity;                     ///< The boolean stating if the notum nodes should be fixed with external viscosity;
    double fixingExternalViscosity[3];                          ///< The array of viscosity values to fix nodes by external viscosity [x][y][z]
    bool ApicalNodeFix[3];                                      ///< The boolean array stating if the apical nodes are fixed, rigid or via viscosity [x][y][z] (accompanying boolean array to specify the fixing method Simulation#ApicalNodeFixWithExternalViscosity)
    bool BasalNodeFix[3];                                       ///< The boolean array stating if the basal nodes are fixed, rigid or via viscosity [x][y][z] (accompanying boolean array to specify the fixing method Simulation#BasalNodeFixWithExternalViscosity)
    bool NotumNodeFix[3];                                       ///< The boolean array stating if the notum nodes are fixed, rigid or via viscosity [x][y][z] (accompanying boolean array to specify the fixing method Simulation#CircumferentialNodeFixWithHighExternalViscosity)
    double notumFixingRange[2];                                 ///< The boundaries of the notum zone in relative x [minX][maxX].
    /**
     * @brief The boolean array stating if the circumference should be fixed, rigid or via viscosity: \n
     * row 0: apical circumferece x,y,z ; row 1: basal circumference x,y,z; row 2: linker apical circumference x,y,z, row 3: linker basal circumference x,y,z, row 4: all circumference x,y,z
     */
    bool CircumferentialNodeFix[5][3];
    double PeripodialElasticity;                                ///< Young's modulus of the peropodial tissue.
    double peripodialApicalViscosity;                           ///< Apical viscosity of the peripodial tissue.
    double peripodialBasalViscosity;                            ///< Basal viscosity of the peripodial tissue.
    double peripodialMidlineViscosity;                          ///< Mid-line viscosity of the peripodial tissue.
    //double currViscMidline;                                     ///<
	//linker zone parameters
    bool BaseLinkerZoneParametersOnPeripodialness;              ///< The boolean stating if the physical properties should be calculated based on the tisseu type weights for linker zones (ShapeBase#columnarGrowthWeight)
    double LinkerZoneApicalElasticity;                          ///< Young's modulus of the apical section of linker tissue between columnar and peripodial
    double LinkerZoneBasalYoungsModulus;                        ///< Young's modulus of the basal section of linker tissue between columnar and peripodial
    double linkerZoneApicalViscosity;                           ///< Internal viscosity of the apical section of linker tissue between columnar and peripodial
    double linkerZoneBasalViscosity;                            ///< Internal viscosity of the basal section of linker tissue between columnar and peripodial
    double linkerZoneMidlineViscosity;                          ///< Internal viscosity of the mid-line section of linker tissue between columnar and peripodial

    int nGrowthFunctions;                                       ///< The number of growth functions active in the simulation.
    bool GridGrowthsPinnedOnInitialMesh;                        ///< The boolean stating if the grid based gorwth is pinned on the mesh for prolonged periods, the timings stored in Simulation#growthPinUpdateTime
    int nGrowthPinning;                                         ///< The number of times growth pinning positions will be updated, the relative positions of elements in tissue.
    std::vector<int> growthPinUpdateTime;                       ///< The vector storing the growth rate pinning update times, as Simulation#timestep counters.
    std::vector<bool> growthPinUpdateBools;                     ///< The vector storint the booleans to check the pinning update has been carried out.
    int gridGrowthsInterpolationType;                           ///< The type of interpolation done on growth map grid, 0 = no interpolation, step function, 1 = linear interpolation (default = 1).
    std::vector<std::unique_ptr<GrowthFunctionBase>> GrowthFunctions;   ///< The vector containing the unique pointers to the growth functions active in the simulation
    int nShapeChangeFunctions;                                  ///< The number of shape change functions active in the simulation.
    std::vector<std::unique_ptr<GrowthFunctionBase>> ShapeChangeFunctions; ///< The vector containing the unique pointers to the shape change functions active in the simulation
    double shapeChangeECMLimit;                                 ///< The threshold of ECM density upon which emergent shape change will be activated (reduction of ECM strength inducing shape change).

    bool thereIsPlasticDeformation;                             ///< The boolean stating that there is plastic deformation (remodelling) in the tissue.
    bool plasticDeformationAppliedToPeripodial;                 ///< The boolean stating that plastic deformation (remodelling) is applied to peripodial tissue.
    bool plasticDeformationAppliedToColumnar;                   ///< The boolean stating that plastic deformation (remodelling) is applied to columnar tissue.
    bool volumeConservedInPlasticDeformation;                   ///< The boolean stating if volume is conserved during plastic deformation (remodelling).
    double plasticDeformationHalfLife;                          ///< The half life or plastic deformation (remodelling).
    double zRemodellingLowerThreshold;                          ///< The lower threshold of z remodelling, below which an element cannot be shrunk in z, fraction.
    double zRemodellingUpperThreshold;                          ///< The lower threshold of z remodelling, below which an element cannot be extended in z, fraction.

    std::vector <std::unique_ptr<Node>> Nodes;                  ///< The vector storing the unique pointers ot the nodes of the simulation.
    std::vector <std::unique_ptr<ShapeBase>> Elements;          ///< The vector storing the unique pointers ot the elements of the simulation.
    size_t nElements;                                           ///< The number of elements of the simulation.
    size_t nNodes;                                              ///< The number of nodes of the simulation.
    std::vector<std::array<double,3>> SystemForces;             ///< The vector storing the system forces applied on each node [Simulation#nNodes][3D]
    std::vector<std::array<double,3>> PackingForces;            ///< The vector storing the packing forces applied on each node [Simulation#nNodes][3D]
    bool addingRandomForces;                                    ///< The boolean stating if random forces are added as noise.
    std::vector <double> randomForces;                          ///< The vector storing the random forces applied on each node, stored in vector form [Simulation#nNodes * 3D].
    double randomForceMean;                                     ///< The mean of the random forces
    double randomForceVar;                                      ///< The variance of the random forces
    double SystemCentre[3];                                     ///< The geometric centre of the tissue mesh.
    bool needPeripodialforInputConsistency;                     ///< The boolean stating if the simulation inputs require peripodial tissue.
    bool thereIsPeripodialMembrane;                             ///< The boolean stating if there is a peripodial membrane in the simulation
    bool symmetricY;                                            ///< The boolean stating if there is y-axis symettricity in tissue
    bool symmetricX;                                            ///< The boolean stating if there is x-axis symettricity in tissue
    bool conservingColumnVolumes;                               ///< The boolean stating if the simulation is conserving the volume of each elemental column rather than each element.
    bool thereIsArtificaialRelaxation;                          ///< The boolean stating if there is artifical relaxation of forces in tissue at a selected time point stored in Simulation#artificialRelaxationTime
    bool relaxECMInArtificialRelaxation;                        ///< The boolean stating if the ECM should be relaxed with the tissue in artifical relaxation of forces at a selected time point stored in Simulation#artificialRelaxationTime
    double artificialRelaxationTime;                            ///< The time of artifical relaxation of forces, in seconds.
    bool stretcherAttached;                                     ///< The boolean stating if there is a stretcher experimental setup attached to tissue.
    std::vector <int> leftClampBorder;                          ///< The Node#Id vector for the nodes at the border of stretcher experiment left side clamp.
    std::vector <int> rightClampBorder;                         ///< The Node#Id vector for the nodes at the border of stretcher experiment right side clamp.
    double leftClampForces[3];                                  ///< The total force vector at the stretcher experiment left side clamp.
    double rightClampForces[3];                                 ///< The total force vector at the stretcher experiment right side clamp.
    bool DVClamp;                                               ///< Boolean stating if the tissue should be clapmed in the DV axis (long axis).
    int distanceIndex;                                          ///< The index of dimension for stretcher clamp position, 0 (x) id attached on long axis, 1 (y) otherwise.
    int StretchInitialTime;                                     ///< The activation time of the stretcher, timestep counter unit.
    int StretchEndTime;                                         ///< The end time of gradual stretching stretcher, timestep counter unit.
    double StretchMin;                                          ///< The minimum relative position in the selected direction (Simulation#distanceIndex), below which the nodes are clamped and move with the stretcher.
    double StretchMax;                                          ///< The maximum relative position in the selected direction (Simulation#distanceIndex), below which the nodes are clamped and move with the stretcher.
    double StretchStrain;                                       ///< The stretch applied in the stretcher device.
    bool PipetteSuction;                                        ///< The boolean if there is a pipette aspiration experiment in the simulation.
    bool ApicalSuction;                                         ///< The boolean if the pipette aspiration experiment is attached to apical surface, if not, it is attached to basal surface.
    bool TissueStuckOnGlassDuringPipetteAspiration;             ///< The boolean stating if the tissue is stuck on glass on the non-aspirated surface
    std::vector <int> TransientZFixListForPipette;              ///< The vector of Node#Id marking the nodes packling to the tip of the pipette.
    int PipetteInitialStep;                                     ///< The initial time step of pipette aspiration
    int nPipetteSuctionSteps;                                   ///< The number of steps the pressure of the pipette aspiration experimetnis increased, times and pressures saved in Simulation#pipetteSuctionTimes and Simulation#pipetteSuctionPressures, respectively.
    std::vector<double> pipetteSuctionTimes;                    ///< The time points where the pipette aspiration pressure is increased, the pressures are stored in Simulation#pipetteSuctionPressures.
    std::vector<double> pipetteSuctionPressures;                ///< The pressures steps of the pipette aspiration experiment, the times are stored in Simulation#pipetteSuctionTimes.
    double pipetteCentre[3];                                    ///< The centre of the pipette tip [x][y][z]
    double pipetteDepth;                                        ///< The z-depth of the suction pressure effectiveness from the pipette tip.
    double pipetteInnerRadius;                                  ///< The inner radius of the pipette.
    double pipetteThickness;                                    ///< The thickness of the pipette.
    double pipetteInnerRadiusSq;                                ///< The pipette inner radius (Simulation#pipetteInnerRadius) squared.
    double effectLimitsInZ[2];                                  ///< The effective suction limits of the pipette in z, from tip to infinity in the opposing direction.
    double SuctionPressure[3];                                  ///< Current suction pressure of the pipette suction experiment.
    double TissueHeight;                                        ///< The height of the tissue
    size_t TissueHeightDiscretisationLayers;                    ///< The number of elemetns the tissue height is discretisized to.
    double boundingBox[2][3];                                   ///< The bounding box of the tissue [lower left corner [x][y][z]] [upper right corner[x][y][z]]
    std::vector <int> pacingNodeCouples0;                       ///< The vector storing the Node#Id list of packing node couples, coupled with Simulation#pacingNodeCouples1.
    std::vector <int> pacingNodeCouples1;                       ///< The vector storing the Node#Id list of packing node couples, coupled with Simulation#pacingNodeCouples0.
    std::vector <bool> pacingNodeCouplesHaveAdhered;            ///< The boolean vector stating if the potentially packing node couple have adheres (couple stored in Simulation#///< The vector storing the Node#Id list of packing node couples, coupled with Simulation#pacingNodeCouples0 and Simulation#///< The vector storing the Node#Id list of packing node couples, coupled with Simulation#pacingNodeCouples1).
    //std::vector <double> initialWeightPointx;                   ///< initial weights of packing nodes in x.
    //std::vector <double> initialWeightPointy;                   ///< initial weights of packing nodes in y.
    //std::vector <double> initialWeightPointz;                   ///< initial weights of packing nodes in z.

    int	nMarkerEllipseRanges;                                   ///< Number of marker ellipses used in perturbation activation
    std::vector<double> markerEllipseBandXCentres;              ///< The vecotr storing the x centres of the marker ellipse bands used in perturbation activation
    std::vector<double> markerEllipseBandR1Ranges;              ///< The vecotr storing the inner radia of the marker ellipse bands used in perturbation activation
    std::vector<double> markerEllipseBandR2Ranges;              ///< The vecotr storing the outer radia of the marker ellipse bands used in perturbation activation
    bool 	thereIsECMChange;                                   ///< There is perturnaton on ECM
    std::vector <int> numberOfECMChangeEllipseBands;            ///< Number of marker ellipses used in perturbation of the ECM
    std::vector< vector<int> > ECMChangeEllipseBandIds;         ///< The band Ids of marker ellipses used in perturbation of the ECM
    std::vector <double> ECMChangeBeginTimeInSec;               ///< The time point where ECM change is activated.
    std::vector <double> ECMChangeEndTimeInSec;                 ///< The time point where ECM change is ended.
    std::vector <double>	ECMStiffnessChangeFraction;         ///< The final ECM stiffness as a fraction of the initial value at the end of perturbation.
    std::vector <double> ECMRenewalHalfLifeTargetFraction;      ///< The final ECM renewal halflife as a fraction of the initial value at the end of perturbation.
    std::vector <double> ECMViscosityChangeFraction;            ///< The final ECM external viscosity as a fraction of the initial value at the end of perturbation.
    std::vector <bool> 	changeApicalECM;                        ///< Boolean stating if the apical ECM should be affected by perturbation
    std::vector <bool> 	changeBasalECM;                         ///< Boolean stating if the basal ECM should be affected by perturbation
    std::vector <bool> 	ECMChangeTypeIsEmergent;                ///< Boolean stating if the ECM perturbation is emergent with fold initiation, rather tahn input ellipse markers.
    std::vector <bool> changedECM;                              ///< The boolaen vector stating if the ECM perturbation has been set up.
    double notumECMChangeInitTime;                              ///< The time point where ECM change to the notum tissue domain is activated.
    double notumECMChangeEndTime;                               ///< The time point where ECM change to the notum tissue domain is ended.
    double notumECMChangeFraction;                              ///< The final ECM stiffness as a fraction of the initial value at the end of perturbation to the notum tissue domain is ended.
    double hingeECMChangeInitTime;                              ///< The time point where ECM change to the hinge tissue domain is activated.
    double hingeECMChangeEndTime;                               ///< The time point where ECM change to the hinge tissue domain is ended.
    double hingeECMChangeFraction;                              ///< The final ECM stiffness as a fraction of the initial value at the end of perturbation to the hinge tissue domain is ended.
    double pouchECMChangeInitTime;                              ///< The time point where ECM change to the pouch tissue domain is activated.
    double pouchECMChangeEndTime;                               ///< The time point where ECM change to the pouch tissue domain is activated.
    double pouchECMChangeFraction;                              ///< The final ECM stiffness as a fraction of the initial value at the end of perturbation to the pouch tissue domain is ended.
	

    bool ThereIsStiffnessPerturbation;                          ///< The bolean stating if there is stiffness perturbation in the system
    std::vector <bool> ThereIsApicalStiffnessPerturbation;      ///< The vector storing if the shape change function is applied to apical surface
    std::vector <bool> ThereIsBasalStiffnessPerturbation;       ///< The vector storing if the shape change function is applied to basal surface
    std::vector <bool> ThereIsWholeTissueStiffnessPerturbation; ///< The vector storing if the shape change function is applied to whole tissue
    std::vector <bool> ThereIsBasolateralStiffnessPerturbation; ///< The vector storing if the shape change function is applied to basal and lateral tissue
    std::vector <bool> ThereIsBasolateralWithApicalRelaxationStiffnessPerturbation;  ///< The vector storing if the shape change function is applied to basal and lateral tissue with inverse compansation on the apical layer
    std::vector <double> stiffnessChangedToFractionOfOriginal;  ///The vector storing the maximum stiffness change as a fraction of the original stiffness.
    std::vector <double> stiffnessPerturbationBeginTimeInSec;   ///The vector storing the stiffness perturbation initiation time points in sec.
    std::vector <double> stiffnessPerturbationEndTimeInSec;     ///The vector storing the stiffness perturbation end time points in sec.
    std::vector <int> numberOfStiffnessPerturbationAppliesEllipseBands; ///< The number of marker ellipse band ids the perturbation is applied to.
    std::vector< std::vector<int> > stiffnessPerturbationEllipseBandIds; ///< The IDs of the marker ellipse band ids the perturbation is applied to.
    std::vector <bool> startedStiffnessPerturbation;

    double packingDetectionThreshold;                           ///< The threshold to detect potentially packing nodes, scaled to average side length.
    double packingDetectionThresholdGrid[10][5];                ///< The threshold to detect potentially packing nodes as a grid, scaled to local average side length.
    double packingThreshold;                                    ///< The threshold to detect packing nodes, scaled to average side length.
    double packingMultiplier;                                   ///< A multiplier stored to scale packing functions where necessary, current value is unity.
    double sigmoidSaturationForPacking;                         ///< The sigmoid saturation for packing, current value is 5, see packing function for details.

    bool 	softPeriphery;                  ///< Boolen stating if there is soft periphery.
    double 	softDepth;                      ///< The range of the soft periphery on xy plane.
    double 	softnessFraction;               ///< The fraction of stiffness perturbation at the soft periphery.
    bool 	softPeripheryBooleans[4];       ///< Boolen array stating if soft peripheri is applicable to [applyToApical]  [applyToBasal]  [applyToColumnar]  [applyToPeripodial]
    bool 	thereIsAdhesion;                ///< Boolean stating if there is adhesion in the system.
    bool    collapseNodesOnAdhesion;        ///< Boolean stating if the nodes should be collapsed in their degrees of freedom upon adheison.
    bool    savedStepwise;                  ///< Boolean stating if collapse is stepwise rather than immediate, to avoid instabilities.
    bool	thereIsEmergentEllipseMarking;  ///< Boolean stating if the emerging folding regions will be marked for further perturbation.
    bool    thereNodeCollapsing;            ///< Boolean stating if the nodes of elements at hogh risk of flipping will be collapsed in their degrees of freedom.

    bool implicitPacking;                   ///< Boolean stating if the hard-wall packing will be calculated implicitly, inside NR iterations.

    bool thereIsExplicitECM;    ///< Boolean stating if there is explcite ECM definition in the system.
    bool addLateralECMManually; ///< Boolean stating if the side ECM will be added during simulation time rather than being a mexh input.
    double lateralECMThickness; ///< The thickness of the side ECM to be added.
    bool thereIsExplicitActin; ///< Boolean stating if there is explcite actin definition in the system.
    double ECMRenawalHalfLife; ///< The half life for ECM renewal inside plastic deformation

    size_t numberOfClones;      ///< Number of mutant clones
    std::vector<double> cloneInformationX;  ///< The vector storing x positions of clone centres.
    std::vector<double> cloneInformationY;  ///< The vector storing y positions of clone centres.
    std::vector<double> cloneInformationR;  ///< The vector storing z positions of clone centres.
    std::vector<bool> cloneInformationUsingAbsolueGrowth; ///< The vector of booleans defining the growth perturbation type of clones, true if absolute growth is provided, false if scaled.
    std::vector<double> cloneInformationGrowth; ///< The vector storing the absoule growth rates OR growth fold changes of clones.

    bool encloseTissueBetweenSurfaces; ///< Boolean stating if there is external enclosing surfaces.
    double initialZEnclosementBoundaries[2]; ///< Initial z position of the enclosing surfaces [top][bottom]
    double finalZEnclosementBoundaries[2];   ///< Final z position of the enclosing surfaces [top][bottom]
    double initialTimeToEncloseTissueBetweenSurfacesSec; ///< Initial time point to start enclosement.
    double finalTimeToEncloseTissueBetweenSurfacesSec;   ///< Final time point to apply enclosement.
    double zEnclosementBoundaries[2];                   ///< Current z position of the enclosing surfaces [top][bottom]
    double packingToEnclosingSurfacesThreshold;         ///< The distance threshold for detection of nodes potentially packing to the surface.
    double packingDetectionToEnclosingSurfacesThreshold; ///< The distance threshold for packing to nodes to surfaces.
    std::vector <int> nodesPackingToPositiveSurface;     ///< Node ids packing to enclosing surface on the top.
    std::vector <int> nodesPackingToNegativeSurface;     ///< Node ids packing to enclosing surface at the bottom.
    std::vector <double> initialWeightPackingToPositiveSurface; ///< Weights of packing to state directionality, for nodes packing to surface on top.
    std::vector <double> initialWeightPackingToNegativeSurface; ///< Weights of packing to state directionality, for nodes packing to surface at the bottom.

    bool thereIsCircumferenceXYBinding;     ///< The boolean stating if circumference is bound in xy plane, all nodes of the same column have their x&y degrees of freedom fixed on the bottom node. This avoids boundary buckling.

    std::unique_ptr<NewtonRaphsonSolver> NRSolver;  ///< The pointer to the newton raphson solver

    Simulation();                                       ///< Constructor
    ~Simulation();                                      ///< Destructor
    void assignTips();                                  ///< This function assigns the nodes marking the tips of the tissue in xy plane.
    bool readExecutableInputs(int argc, char **argv);   ///< This function reads the input executables from user input.
    bool initiateSystem();                              ///< This function initiates the system with model inputs.
    bool initiateSavedSystem();                         ///< This function initiates a saved system from save files.
    void calculateSystemCentre();                       ///< This function calculates the tissue geometric centre.
    void resetForces(bool resetPacking);                ///< This function resets all forces to zero at the beginning of each iteration. Depending on the packing calculation type (implicit/explicit) the boolean states if teh packing forces should be reset or not.
    void calculateApicalSize();                         ///< This function calculates the xy-plane bounding box of the apical side of the tissue
    void calculateBoundingBox();                        ///< This function calculates the 3D bounding box of the tissue
    void calculateZProjectedAreas();                    ///< This function calculates the z-projected areas of all elemetns and nodes.
    void correctzProjectedAreaForMidNodes();            ///< This function corrects the z-projected areas for mid-l;ayer nodes, as their areas are counted twice, by both elements.
    void clearProjectedAreas();                         ///< This function clears the z-projected areas of all elemetns and nodes.
    void checkForExperimentalSetupsBeforeIteration();   ///< This macro function checks for all experimental setup updates at the beginning of a time step, prior to iterations.
    void checkForExperimentalSetupsWithinIteration();   ///< This macro function checks for all experimental setup updates at each N-R iteration.
    void checkForExperimentalSetupsAfterIteration();    ///< This macro function checks for all experimental setup updates at the end of iterations, post convergence.

    void checkForEmergentEllipseFormation();            ///< Check for emergent marker elipse bands with fold initiation in the tissue.
    void checkForEllipseIdUpdateWithECMDegradation();   ///< Check for emergence of ellipse bands as a result of ECM loss.
    void checkEllipseAllocationWithCurvingNodes();      ///< Check the ellipse band allocation of nodes residing on fold initiation.
    void updateEllipseWithCollapse();                   ///< Update marker ellipses with collapsing nodes.
    void checkForLeftOutElementsInEllipseAssignment();    ///< Check if elements at border of two emergent folds are missed in marking, and correct as needed.

    void checkECMChange();                                                              ///< Check if there is ECM perturbation.
    void updateChangeForExplicitECM(int idOfCurrentECMPerturbation);                    ///< Update physical properties of the ECM with the ECM perturbaiton function of input ID.
    void updateECMRenewalHalflifeMultiplier(int idOfCurrentECMPerturbation);            ///< Update ECM renewal half-life with the ECM perturbaiton function of input ID.
    void updateChangeForViscosityBasedECMDefinition(int idOfCurrentECMPerturbation);    ///< Update ECM viscosity with the ECM perturbaiton function of input ID.
    void calculateChangeRatesForECM(int idOfCurrentECMPerturbation);                    ///< Calculate the physical property change rates for ECM perturbation with the ECM perturbaiton function of input ID.

    void checkStiffnessPerturbation();                                                  ///< Check for stiffness perturbation on the tissue.
    void updateStiffnessChangeForActin(int idOfCurrentStiffnessPerturbation);           ///< Update stiffness of actin  with the stiffness perturbaiton function of input ID.
    void calculateStiffnessChangeRatesForActin(int idOfCurrentStiffnessPerturbation);   ///< Calculate the stiffness change rates for actin with the stiffness perturbaiton function of input ID.


    void updateOnFoldNodesFromCollapse();                   ///< Mark nodes as on fold, due to collapsed elements.
    void assignFoldRegionAndReleasePeripodial(Node* NodeMAster, Node* NodeSlave); ///< Nodes are assigned to be on folding regions if they lie in between two adhered nodes.
    void artificialRelax();                                 ///< Artificially relax all deformation on all elements at selected time point.
    bool runOneStep();                                      ///< Run the simulation for one time step.
    void updateOneStepFromSave();                           ///< Update the simulation for time step from saved files.
    void updatePlasticDeformation();                        ///< Update the plastic deformation (remodelling) of all nodes.
    void updateStepNR();                                    ///< Update the positions with solving for the displacements with the N-R iterations.
    void updateElementPositionsinNR(gsl_matrix* uk);        ///< Update elemental positions during one iteration of the NR numerical solving for the displacements.
    void updateNodePositionsNR(gsl_matrix* uk);             ///< Update nodal positions during one iteration of the NR numerical solving for the displacements.
    void calculateRandomForces();                           ///< Calculate the random forces on elements.
    void addRandomForces(gsl_matrix* gExt);                 ///< Add random forces to system forces.

    void detectPacingNodes();                                           ///< Detect nodes potentially packing to each other
    void calculatePackingForcesImplicit3D();                            ///< Calculate packing forces on nodes due to packing to each other.
    void calculatePackingJacobian3D(gsl_matrix* K);                     ///< Update the system jaconbian with the derivatives of forces due to packing to each other, with respect to positions.
    void detectPackingToPipette();                                      ///< Detect nodes potentially packing to the pipette wall in pipette aspiration experiment.
    void calculatePackingToPipetteForcesImplicit3D();                   ///< Calculate packing forces on nodes due to packing to pipette wall in pipette aspiration experiment.
    void calculatePackingToPipetteJacobian3D(gsl_matrix* K);            ///< Update the system jaconbian with the derivatives of forces due to packing to the pipette wall in pipette aspiration experiment, with respect to positions.
    void detectPacingToEnclosingSurfacesNodes();                        ///< Detect nodes potentially packing to the enclosing surfaces in z.
    void calculatePackingForcesToEnclosingSurfacesImplicit3D();         ///< Calculate packing forces on nodes due to packing to enclosing surfaces.
    void calculatePackingToEnclosingSurfacesJacobian3D(gsl_matrix* K);  ///< Update the system jaconbian with the derivatives of forces due to packing to enclosing surfaces, with respect to positions.
    void addPackingForces(gsl_matrix* gExt);                            ///< Add al packing forces to external foces vector.


    void addValueToMatrix(gsl_matrix* K, int i, int j, double value);   ///< Helper function, adds value to matrix K, in indices (i,j).
    void updateElementPositions();                                      ///< Update element positions, from the nodal positions.
    void updateMasterSlaveNodesInBinding();                             ///< Update node degrees of freedom binding information stored in N-R solver to the saved data.
    void updateElementPositionsSingle(size_t i );                       ///< Update the position of element at index i, in the Simulation#Elements vector.
    void alignTissueDVToXPositive();                                    ///< Align tissue DV axis to x axis, the DV axis is defiend by Simulation#ventralTipIndex and Simulation#dorsalTipIndex, as assigned in Simulation#assignTips.
    bool checkFlip();                                                   ///< Check if any of the elements have flipped, and report error accordingly.
    void wrapUpAtTheEndOfSimulation();                                  ///< Correct tissue alignement once more before simulation is finalised to have the correct final mesh.
    void calculateDVDistance();                                         ///< Calclate the long (DV) axis tip to tip distance (not contour length), the DV axis is defiend by Simulation#ventralTipIndex and Simulation#dorsalTipIndex, as assigned in Simulation#assignTips.
    void fixNode0InPosition(double x, double y, double z);              ///< Fix the node 0 in 3D space, in case of zero external viscosity on all axes, checked in Simulation#checkForZeroExternalViscosity.

};

#endif
