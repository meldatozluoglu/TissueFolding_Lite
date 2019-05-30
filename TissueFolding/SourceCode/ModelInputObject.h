/*
 * ModelInputObject.h
 *
 *  Created on: 5 Aug 2014
 *      Author: melda
 */

#ifndef MODELINPUTOBJECT_H_
#define MODELINPUTOBJECT_H_

#include <stdio.h>
#include <iostream>
#include <memory>

class Simulation;

/*! ModelInputObject class */
class ModelInputObject {
private:
    void printErrorMessage(std::string currentInput, std::string sourceFuction, std::string expectedInput); ///< This function writes the error message
    bool checkFileStatus(std::ifstream &file, std::string fileName);        ///< This function checks the health of given file.
    bool readPysicalProperties(std::ifstream &file);                        ///< This function reads the physical parameters of the tissue from file.
    bool readSaveOptions(std::ifstream &file);                              ///< This function reads the save options of the simulation from file.
    bool readMeshParameters(std::ifstream& file);                           ///< This function reads the mesh structure related parameters of the tissue from file.
    bool readPeripodialMembraneParameters(std::ifstream& file);             ///< This function reads the peripodial membrane related parameters of the tissue from file.
    bool readLinkerZoneParameters(std::ifstream& file);                     ///< This function reads the linker zone related parameters of the tissue from file.
    bool readExternalViscosityParameters(std::ifstream& file);              ///< This function reads the external viscosity setup for the whole tissue.
    bool readNodeFixingParameters(std::ifstream& file);                     ///< This function reads the inputs relating to fixing the nodes of the tissue from file.
    bool readNodeBindingParameters(std::ifstream& file);
    bool readTimeParameters(std::ifstream &file);                           ///< This function reads the time related parameters of the simulation from file.
    bool readMeshType2(std::ifstream& file);                                ///< This function reads the mesh structure details for a mesh input as columns and rows of prisms.
    bool readMeshType4(std::ifstream& file);                            	///< This function reads the mesh structure details for a mesh given as a pre-assembled mesh in a separate input file.
    bool readGrowthOptions(std::ifstream& file);                    		///< This function reads the growth functions and related parameters of the simulation from file.
    bool readGrowthType1(std::ifstream& file);                              ///< This function reads the uniform growth parameters from file (UniformGrowthFunction).
    bool readGrowthType2(std::ifstream& file);                              ///< This function reads the ring growth parameters from file (RingGrowthFunction).
    bool readGrowthType3(std::ifstream& file);                          	///< This function reads the grid based growth parameters from file (GridBasedGrowthFunction). It will utilise a separate input file storing the growth rates.
    bool readShapeChangeOptions(std::ifstream& file);                       ///< This function reads the active shape change functions and related parameters of the simulation from file.
    bool readPlasticDeformationOptions(std::ifstream& file);            	///< This function reads the parametrs for plastic deformation, as a response to strains and stresses in the tissue.
    bool readShapeChangeType1(std::ifstream& file);                         ///< This function reads the shape change  parameters from file (UniformShapeChange).
    bool readShapeChangeType2(std::ifstream& file);                         ///< This function reads the shape change  parameters from file (markerEllipseBased).
    bool readMarkerEllipseBandOptions(std::ifstream& file);
    bool readShapeChangeType3(std::ifstream& file);                         ///< This function reads the shape change  parameters from file (GridBasedShapeChange).
    bool readStretcherSetup(std::ifstream& file);                           ///< This function reads the stretcher experimental setup parameters of the simulation from file.
    bool readPipetteSetup(std::ifstream& file);                             ///< This function reads the pipette aspiration experimental setup parameters of the simulation from file.
    bool readECMPerturbation(std::ifstream& file);
    bool readStiffnessPerturbation(std::ifstream& file);
    bool readExplicitECMOptions(std::ifstream& file);                       ///< This function reads parameters relating to the definitin of an explicit ECM layer.
    bool readAdhesionOptions(std::ifstream& file);
    bool readNodeCollapseOptions(std::ifstream& file);
    bool readExplicitActinOptions(std::ifstream& file);
    bool readColumnViseVolumeConservationOptions(std::ifstream& file);
    bool readMutationOptions(std::ifstream& file);
    bool readartificialRelaxationOptions(std::ifstream& file);
    bool readEnclosementOptions(std::ifstream& file);

public:
    std::shared_ptr<Simulation> Sim;        ///< The pointer to the simulation object, for which the parameters are being read from the modelinput file.
    std::string parameterFileName;          ///< The name (including path) of the file containing the model input parameters.
    ModelInputObject();                     ///< The constructor of the ModelInputObject.
	~ModelInputObject();
    bool readParameters();                  ///< This is the main funciton reading the parameters from file
};

#endif /* MODELINPUTOBJECT_H_ */
