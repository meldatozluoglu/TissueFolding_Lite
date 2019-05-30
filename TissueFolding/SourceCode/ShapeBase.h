#ifndef ShapeBase_H
#define ShapeBase_H


#include <iostream>
#include <ostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include "Node.h"
#include "ReferenceShapeBase.h"
#include "GrowthFunctionBase.h"
#include "GrowthFunctionTypes.h"

class ShapeBase{
private:
    void ParentErrorMessage(std::string functionName);				///<Error message displayed when a virtual function is called through the ShapeBase(parent), while it should have been called through a child (eg. Prism). For functions taking in string inputs.
    bool ParentErrorMessage(std::string functionName, bool returnValue); 	///<Error message displayed when a virtual function is called through the ShapeBase(parent), while it should have been called through a child (eg. Prism). For functions taking in string and bool inputs, returning bool.
    double ParentErrorMessage(std::string functionName, double returnValue);	///<Error message displayed when a virtual function is called through the ShapeBase(parent), while it should have been called through a child (eg. Prism). For functions taking in string and double inputs, returning double.
    int ParentErrorMessage(std::string functionName, int returnValue);		///<Error message displayed when a virtual function is called through the ShapeBase(parent), while it should have been called through a child (eg. Prism). For functions taking in string and int inputs, returning int.
protected:
    int                         ShapeType;                      ///< The integer defining the type of the shape, Prisms shape type = 1;
    size_t                      nNodes;                         ///< The number of nodes of the element, it is based on ShapeBase#ShapeType
    size_t                      nDim;                           ///< The number of dimensions for the positions of each of the nodes of the element
    std::array<int,3>           IdentifierColour;				///< The unique identifier colour of the element, this is used for "picking" in the visual interface.
    std::array<double,3>        GrowthRate;                     ///< Growth rate recording for display purposes only. The recorded growth rate in x, y, and z  coordinates, does not record shear deformation induced in growth. Recorded in exponential form through time step, converted to rate per hour for display within the visual interface
    gsl_matrix*                 growthIncrement;				///< The matrix (3,3) representing the incremental growth in current time step. Reset to identity at the beginning of each time step, updated in growth functions, and utilised to update Fg.
    gsl_matrix*                 plasticDeformationIncrement;	///< The matrix (3,3) representing the incremental plastic deformation (treated as growth) in current time step. Set in plastic deformation calculation at each step, and utilised to update Fg.
    gsl_matrix*                 shapeChangeIncrement;           ///< The matrix (3,3) representing the incremental shape change in current time step. Reset to identity at the beginning of each time step, updated in shape change functions, and utilised to update Fg.
    double                      zRemodellingSoFar;				///< The z remodelling that have been applied to elemetn up to the current time step. This parameter is used to limit extreme thinning or elongation of elements.
    double                      columnarGrowthWeight;			///< The fraction defining how close to the columnar layer the element is. 1.0 for columnar layer, 0.0 for peripodial membrane elements, and scaled according to position in the elements surrounding the lumen.
    double                      peripodialGrowthWeight;         ///< The fraction defining how close to the peripodial membrane the element is. 0.0 for columnar layer, 1.0 for peripodial membrane elements, and scaled according to position in the elements surrounding the lumen.
    std::array<double,6>        ShapeChangeRate;				///< Shape change rate of the elements, only orthagonal shape changes are allowed (x, y, z). Shape changes will be scaled to conserve volume, thus three values will not be independent.
    bool                        rotatedGrowth;					///< The boolean stating if the element has rotated from the growth axis, hence the calculated growth requires further rotation to follow tissue axes.
    std::array<double,3>        relativePosInBoundingBox;		///< The relative position on x-y plane, within the bounding box of the tissue(x,y).
    std::array<double,3>        initialRelativePosInBoundingBox;///< The relative position on x-y plane, within the bounding box of the tissue(x,y) at the beginning of simulation. This is used when growth rates are pinned to the initial structure of the tissue.
    double                      initialRelativePositionInZ;     ///< The relative position on z-height of tissue, taken not in z direction but in tissue layers, 0 being on the apical surface and 1 being on the basal surface.
    size_t                      numberOfGaussPoints;            ///< The number of Gauss points used in numerical deforamtion calculation.
    std::array<std::array<double,3>,6>	gaussPoints;            ///< The array contianing all the Gauss points for element. Set up is for 6, can work for any number as loops are kept indexed up to numberOfGaussPoints.
    std::array<double,6>        gaussWeights;                   ///< The array for storing the weights of each Gauss point for element.
    std::vector<gsl_matrix*> 	ShapeFuncDerivatives;           ///< The array of matrices for shape function derivatives. The array stores a ShapeBase#nDim by ShapeBase#nNodes matrix for each gauss point (there are 3 Gauss points for prisms).
    std::vector<gsl_matrix*> 	ShapeFuncDerStacks;                    		///< The array of matrices of shape function derivatives in stacked format for ease of matrix operations. The array stores a (ShapeBase#nDim * ShapeBase#nDim) by (ShapeBase#nDim * ShapeBase#nNodes) matrix for each gauss point (there are 3 Gauss points for prisms).
    std::vector<gsl_matrix*> 	InvdXdes;                              		///< The array stores inverse of the matrix for derivatives of world coordinates with respect to barycentric coordinates (dX / de). The array stores an ShapeBase#nDim by ShapeBase#nDim  matrix for each gauss point (there are 3 Gauss points for prisms).
    std::array<double,6>     	detdXdes;                              		///< The array stores the determinants of the matrices for derivatives of world coordinates with respect to barycentric coordinates (dX / de). The array stores a double value for each gauss point (there are 3 Gauss points for prisms).
    std::vector<gsl_matrix*> 	Bmatrices;                             		///< The array stores the B matrix for the calculation of stiffness matrix, see for ShapeBase#calculateBTforNodalForces calculation. The array stores an ShapeBase#nNodes by (ShapeBase#nDim*ShapeBase#nNodes)  matrix for each Gauss point (there are 3 Gauss points for prisms).
    std::vector<gsl_matrix*> 	FeMatrices;					///< The array stores the elastic part of the deformation matrix. The array stores an ShapeBase#nDim by ShapeBase#nDim  matrix for each Gauss point (there are 6 Gauss points for prisms).
    std::vector<gsl_matrix*> 	invJShapeFuncDerStack;				///< The array stores the shape function derivatives multiplied by the inverse Jacobian stack, for each Gauss point. See ShapeBase#calculateBTforNodalForces for calculation.
    std::vector<gsl_matrix*> 	invJShapeFuncDerStackwithFe;			///< See ShapeBase#calculateInvJShFuncDerSWithFe for calculation.
    std::vector<gsl_matrix*> 	elasticStress;					///< The array of matrices for elastic stress of the element. The array stores a 6 by 6 matrix for each Gauss point (there are 6 Gauss points for prisms).
    std::vector<gsl_matrix*> 	viscousStress;					///< The array of matrices for internal viscous stress of the element. The array stores a 6 by 6 matrix for each Gauss point (there are 6 Gauss points for prisms).
    gsl_matrix*              	TriPointF;                             		///< The deformation matrix of the element resulting from iteration over all Gauss points. The dimensions of the matrix is ShapeBase#nDim by ShapeBase#nDim.
    gsl_matrix*              	ElementalElasticSystemForces;          		///< The matrix stores the elemental elastic forces. The dimensions of the matrix is ShapeBase#nNodes by ShapeBase#nDim.
    gsl_matrix*          	ElementalInternalViscousSystemForces;  		///< The matrix stores the elemental internal viscous forces. The dimensions of the matrix is ShapeBase#nNodes by ShapeBase#nDim.
    std::array<double,6> 	detFs;                              		///< The array stores the determinant of the deformation matrix for each Gauss point.
    double 			ZProjectedBasalArea;				///< The z-projected area of the basal surface of the element.
    double 			ZProjectedApicalArea;			    	///< The z-projected area of the apical surface of the element.
    double 			BasalArea;					///< The area of the basal surface of the element.
    double 			ApicalArea;					///< The area of the apical surface of the element.
    double			exposedLateralAreaApicalSide;			///< The area of the element on a linker position, and has lateral sides exposed to outside of the tissue, on the apical side, therefore should feel external viscosity.
    double			exposedLateralAreaBasalSide;			///< The area of the element on a linker position, and has lateral sides exposed to outside of the tissue, on the basal side, therefore should feel external viscosity.

    bool 			elementHasExposedApicalSurface;			///< The boolean stating if the element has any apical surface exposed to the environment
    bool 			elementHasExposedBasalSurface;			///< The boolean stating if the element has any basal surface exposed to the environment
    int 			exposedApicalSurfaceNodeIds[3];			///< The int array of size 3, listing the node IDs of element that form the exposed apical surface. The IDs are the node IDs on the element (0-5 for prism), not the actual Node#Id.
    int 			exposedBasalSurfaceNodeIds[3];			///< The int array of size 3, listing the node IDs of element that form the exposed basal surface. The IDs are the node IDs on the element (0-5 for prism), not the actual Node#Id.
    int 			exposedLateralAreaApicalSideNodeIds[4];		///< The int array of size 4, listing the node IDs of element that form the lateral surface exposed apically. The IDs are the node IDs on the element (0-5 for prism), not the actual Node::Id.
    int 			exposedLateralAreaBasalSideNodeIds[4];		///< The int array of size 4, listing the node IDs of element that form the lateral surface exposed basally. The IDs are the node IDs on the element (0-5 for prism), not the actual Node::Id.
    size_t 			nLateralSurfaceAreaNodeNumber;			///< Number of nodes that form the lateral surfaces for the element.
    size_t 			nSurfaceAreaNodeNumber;				///< Number of nodes that form the apical/basal surfaces for the element.

    double 			stiffnessPerturbationRateInSec;         	///< The rate at which the stiffness of the element will be perturbed, used with the model inputs from "Stiffness_Perturbation:" header in model input file
    double 			minimumValueOfStiffnessMultiplier;		///< The lower bound of stiffness modification multiplier, exists to prevent elements reaching unintended zero or negative stiffness values.
    double 			maximumValueOfStiffnessMultiplier;		///< The upper bound of stiffness modification multiplier, exists to prevent elements reaching unrealistic hard stiffness values.
    double 			mutationGrowthRatePerSec;			///< The growth rate set by a mutant clone covering this element.
    double 			mutationGrowthFold;				///< The rate of fold change in growth rate set by a mutant clone covering this element.
    void 			setShapeType(std::string TypeName);		///< The function sets the type of the shape.
    void 			readNodeIds(const std::vector<int>& inpNodeIds);///< The function sets the Node#Id array that constructs the shape.
    void 			setPositionMatrix(const std::vector<std::unique_ptr<Node>>& Nodes);	///< The function sets the ShapeBase#Positions matrix to define the locations of each constructing node.
    void 			setTissuePlacement(const std::vector<std::unique_ptr<Node>>& Nodes);	///< The function sets the placement of the element within the tissue
    void 			setTissueType(const std::vector<std::unique_ptr<Node>>& Nodes);		///< The function sets the tissue type of the element
    void 			setReferencePositionMatrix();						///< The function sets the RefereneceShapeBase#Positions matrix to define the reference positions of the element.
    void 			setIdentificationColour();						///< The function sets the unique ShapeBase#IdentifierColour colour for the element, which is used in element picking from the user interface.
    void 			rotateReferenceElementByRotationMatrix(std::array<double,9> rotMat);	///< The function rotates the reference of the element (ShapeBase#ReferenceShape) by input rotation matrix, provided as a double pointer of 9 doubles.
    bool 			InvertMatrix(boost::numeric::ublas::matrix<double>& input, boost::numeric::ublas::matrix<double>& inverse); ///< The function takes the first input matrix, and writes the inverse on the second input. False is returned if the matrix is not inverted. Input format is ublas matrices (slow).
    bool 			InvertMatrix(gsl_matrix* input, gsl_matrix* inverse); 			///< The function takes the first input matrix, and writes the inverse on the second input. False is returned if the matrix is not inverted. Input format is gsl matrices (fast).

    void 			updateNodeIdsFromSave(std::ifstream& file);				///< The function reads the ShapeBase#NodeIds of the current shape from save file provided as input.
    void 			updateReferencePositionMatrixFromSave(std::ifstream& file); 		///< The function reads and updates the ShapeBase#ReferenceShape positions (ReferenceShapeBase#Positions) of the current shape from save file provided as input.
    virtual void 		calculateReferenceVolume(){ParentErrorMessage("calculateReferenceVolume");}  ///<Virtual function of the ShapeBase class to calculate volume of the ShapeBase#ReferenceShape

    bool                calculateGrowthStrainsRotMat(double* v);				///< The function calculates the rotation matrix to apply on growth strains to align growth with the current x axis of the tissue.
    void                calculateForces3D(const std::vector<std::unique_ptr<Node>>& Nodes,  gsl_matrix* displacementPerDt); ///< The function calculates the viscous and elastic forces generated by the element.
    gsl_matrix* 		calculateEForNodalForcesKirshoff(gsl_matrix* C);			///< This function calculates the green strains for a Kirshoff material model.
    gsl_matrix* 		calculateCauchyGreenDeformationTensor(gsl_matrix* Fe);			///< This function calculates the Caucy-Green deformation tensor, from the elastic part of the deformation gradient
    gsl_matrix* 		calculateSForNodalForcesKirshoff(gsl_matrix* E);			///< This function calculates the Secons order Piola-Kirshoff stress tensor for Kirshoff material model.
    gsl_matrix* 		calculateSForNodalForcesNeoHookean(gsl_matrix* invC, double lnJ);	///< This function calculates the Secons order Piola-Kirshoff stress tensor for Neo-Hookean material model.v
    void 			updateLagrangianElasticityTensorNeoHookean(gsl_matrix* invC,double lnJ, int pointNo);	///< This function calcualtes the Lagrangian elasticity tensor for Neo-Hookean material model.
    gsl_matrix* 		calculateCompactStressForNodalForces(double detFe,gsl_matrix* Fe, gsl_matrix* S, gsl_matrix *Stress); ///< This function calculates elemental stress in Voigt notation.
    gsl_matrix* 		calculateInverseJacobianStackForNodalForces(gsl_matrix* Jacobian);	///< This function calculates the stack matrix of inverse Jacobians, used to calculate the nodal forces.
    gsl_matrix* 		calculateBTforNodalForces(gsl_matrix* InvJacobianStack, gsl_matrix* ShapeFuncDerStack, gsl_matrix *B, gsl_matrix* invJShFuncDerS); ///< This function calculates the B matrix, to calculate the nodal force.
    void			calculateInvJShFuncDerSWithFe(gsl_matrix * currFe, gsl_matrix * InvDXde, gsl_matrix* ShapeFuncDerStack, gsl_matrix *invJShFuncDerSWithF); ///< This function calculates the collated matrix from inverse jaconians and shape function derivatives.
    gsl_matrix* 		calculateVelocityGradientTensor(gsl_matrix* B, gsl_matrix* displacementPerDt);  ///< This function calculates the velocity gradient tensor.
    gsl_matrix* 		constructElementalDisplacementMatrix(gsl_matrix* displacement); 		///< This function will assemble elemental node displacement matrix from the input displacement matrix for the whole system.
    gsl_matrix* 		calculateRateOfDeformationTensor(gsl_matrix* l);				///< This function will calculate rate of deformation tensor from velocity gradient tensor
    void 			calculateViscousStress(gsl_matrix* d, gsl_matrix* viscousStress);		///< This function will calculate internal viscous stress of the element from rate of deformation matrix.
    void 			calculateViscousForces(gsl_matrix*  gv, gsl_matrix*  BTdetFdetdXde, gsl_matrix* viscousStress); ///< This function  will calculate the elemental viscous forces from viscous stress.

    void    			consturctBaTBb(gsl_matrix* B, gsl_matrix* BaT, gsl_matrix* Bb, int a, int b);	///< This function constructs nodal subrange of B matrix for node couple a & b.
    void   			calculateElasticKIntegral1(gsl_matrix* currElementalK,int pointNo);		///< This function calcultes the first part of the integral for the stiffness matirx, the elastic part of the system Jacobian.
    void			calculateElasticKIntegral2(gsl_matrix* currElementalK,int pointNo);             ///< This function calcultes the second part of the integral for the stiffness matirx, the elastic part of the system Jacobian.
    void			calculateViscousKIntegral1(gsl_matrix* currElementalK, gsl_matrix* paranthesisTermForKv1, int pointNo);	///< This function calcultes the first part of the integral for the internal viscous part of the system Jacobian.
    void			calculateViscousKIntegral2(gsl_matrix* currElementalK,int pointNo);             ///< This function calcultes the second part of the integral for the viscous part of the system Jacobian.
    void			calculateVelocityGradient( gsl_matrix* velocityGradient, gsl_matrix* displacementPerDt, int pointNo);	///< This function calculates the velocity gradient.
    void			calculateOuterProduct(gsl_matrix* a, gsl_matrix* b, gsl_matrix* outerProduct);	///< Calculates the outer product, maths helper function.
    gsl_matrix* 		calculateSymmetricisedTensorProduct(gsl_matrix* a, gsl_matrix* b);          ///< Calculates the symetricised tensor product, maths helper function.

    bool 			disassembleRotationMatrixForZ(gsl_matrix* rotMat);			///< This function extracts the z rotation from a rotation matrix.
    bool 			calculate3DRotMatFromF(gsl_matrix* rotMat);					///< This function dissects the deformation gradient of the element into the rigid body rotation and deformation.

    gsl_matrix* 		D;                                                      ///< elasticity tensor for Kirshoff material
    gsl_matrix* 		CoeffMat;                                               ///< The coefficient matrix relating the shape function derivative stack to the Voigt notation of elemental stress nad strain.
    //need to construct double array of size: D81[nGaussPoints][3][3][3][3];
    std::vector<std::array<std::array<std::array<std::array<double,3>,3>,3>,3>> D81;				///<Lagrangian elasticity tensor, vector for the number of Gauss points in simulaiton.
    double 			E;                                                          ///< Young's modulus of the element.
    double			v;                                                          ///< Poisson's ratio of the element.
    double 			internalViscosity;                                          ///< Current internal viscosity of the element.
    double 			originalInternalViscosity;                                  ///< The internal viscosity of the element at the beginning of the simulation, prior to physical property perturbations.
    double 			lambda;                                                     ///< Lame's second parameter, driven from Young's modulus and Poisson's ratio of the element
    double			mu;                                                         ///< Sheer modulus of the element.
    gsl_matrix* 	InvFg;                                                      ///< Inverse of growth matrix
    gsl_matrix* 	Fsc;                                                        ///< Shape change matrix
    gsl_matrix* 	InvFsc;                                                     ///< Inverse of shape change matrix
    gsl_matrix* 	TriPointKe;                                                 ///< Current elastic part of the Jacobian (stiffness matrix) of the system, averaged over all Gauss Points.
    gsl_matrix* 	TriPointKv;                                                 ///< Current viscous part of the Jacobian of the system, averaged over all Gauss Points.
public:
    double 			stiffnessMultiplier;                                        ///< Current stiffness multiplier of the element, initially 1.0, modulated by stiffness perturbations						///< The double for the multiplier that will define Young's modulus stress stiffening.
    gsl_matrix*		remodellingPlaneRotationMatrix;                             ///< The rotation matrix converting the xyz coordinate system to the plane of remodelling for the lateral elements.
    gsl_matrix* 	Fg;                                                         ///< Growth deformation gradient
    int 			Id;                                                         ///< The unique ID of the element, without remodelling, equal to its indes on the Simulation#Elements vector.
    int				ShapeDim;                                                   ///< The dimension of the shape in workd coordiantes (2D vs 3D).
    std::vector<int> NodeIds;                                                   ///< The vector storing the unique IDs (Node#Id) of nodes constructing this element. Their order is consistent for a given shape type.
    virtual ~ShapeBase(){                                                       /// The ShapeBase destructor. This destructor should not be called uner healthy conditions.
        //while deleting a ShapeBase* that happens to point a child, this destructor will be called after the child destructor
    }
    std::vector<std::array<double,3>> 	Positions;								///< The vector storing the positions of the nodes constructing the element.
    std::unique_ptr<ReferenceShapeBase> ReferenceShape;                         ///< The pointer to the reference shape object that defines th reference shape of this element.
    gsl_matrix* 	Strain;                                                     ///< The gsl_matrix pointer, storing the address of the current strains on the element.

    bool 			isFlipped;                                                  ///< Boolean stating if the element is flipped. The simulation will be stopped if there are flipped elements.
    bool 			IsChangingShape;
    //bool			willBeRefined;
    int 			tissuePlacement;                                            ///< 1 -> apical, 0 -> basal, 2->middle, 3 -> lateral
    int 			tissueType;                                                 ///< The tissue type is 000 for columnar layer, 1 for peripodial membrane, and 2 for linker zone
    bool			spansWholeTissue;                                           ///< Boolean staing is the element spans the whole tissue. This is used to identify mid-layer tagged tissues (tissuePlacement = 2), that should still have apical abd basal responses
    int 			compartmentType;                                            ///< integer identifying the compartment of the tissue in DV axis, 0 pouch, 1 hinge, 2 notum
    double 			compartmentIdentityFraction;                                ///< The weight defining the constibution of each compartment to the physical identity of this element.
    bool			isECMMimicing;                                              ///< Boolean stating if the element is an ECM element.
    bool			isECMMimimcingAtCircumference;                              ///< Boolean stating if the element is an ECM element at the circumference of the tissue.
    bool			atBasalBorderOfECM;                                         ///< Boolean stating if the element is at the basal border of the cellular layer, linking to ECM elements.
    bool			isActinMimicing;                                            ///< Boolean stating if the element is forming the actin dense layer on the apical surface.
    bool			atApicalBorderOfActin;                                      ///< Boolean stating if the element is at the apical border of the soft cellular layer, linking to actin dense layer.
    bool			IsAblated;                                                  ///< Boolean stating if the element is laser ablated, therefore dead.
    bool			atSymetricityBoundary;                                      ///< Boolean stating if the element is at the symmetricity boundary.
    //bool			IsClippedInDisplay;	
    //bool 			IsXSymmetricClippedInDisplay;
    //bool			IsYSymmetricClippedInDisplay;
    double 			CurrShapeChangeToAdd[3];                                    ///< The current shape change to be applied to the element, in form of 3D rates.
    double  		GrownVolume;                                                ///< Current volume of the element after growth.
    double  		VolumePerNode;                                              ///< Volume per node of the element.
    bool 			capElement;                                                 ///< Boolean stating if the element is capped at its remodelling due to restrictions in z remodelling (avoiding too thin ot too thick layers due to numerical error).
    std::vector<int> elementsIdsOnSameColumn;                                   ///< The vector storing the unique element IDs of each element that is on the same columnar region of the tissue, i.e. the elements share apical/basal surfaces.
    int 			basalNeigElementId;                                         ///< This is recorded only for apical nodes of the columnar layer. If not recorded, id is -1.
    bool 			insideEllipseBand;                                          ///< Boolean stating if the element is marked by any identifier bands for physical perturbation.
    int 			coveringEllipseBandId;                                      ///< The unique ID of the covering perturbation band.

    double 			emergentShapeLongAxis[2];                                   ///< The long axis of the emergent shape. This is necessary for analysis of emergent growth orientations.
    double 			emergentShapeShortAxis[2];                                  ///< The short axis of the emergent shape. This is necessary for analysis of emergent growth orientations.

    double			plasticDeformationHalfLifeMultiplier;                       ///< The multiplier to modify the remodelling half-life upon physical property perturbation.
    bool 			isMutated;                                                  ///< Boolean stating if the element is mutated.

    bool 			thereIsGrowthRedistribution;                                ///< Boolean stating if there is redistribution of growth among mesh elements.
    bool 			growthRedistributionShrinksElement;                         ///< Boolean stating if the growth distribution is taking material out of thes element to redistribute ot others.
    double 			growthRedistributionScale;                                  ///< The extent of the redistribution of volume.
    bool			RotatedElement;                                             ///< The boolean stating if the elemetn has rigid body rotation.
    gsl_matrix* 	GrowthStrainsRotMat;                                        ///< The rotation matrix needed to correct for the rigid body rotations of the element.

    std::array<double,3> 	apicalNormalCurrentShape;							///< The apical normal of the current shape.
    int                         getId();                                        ///< The function returns the unique ID of the element.
    std::string                 getName();                                      ///< The function returns the name of the element.
    int                         getShapeType();									///< The function returns the shape type of the element.
    size_t                      getNodeNumber();								///< The function returns the number of nodes o the element.
    const std::vector<int>& 	getNodeIds();									///< The function returns the vector of node IDs.
    int                         getNodeId(int i);								///< The function returns the input i^{th} node's ID.
    size_t                      getDim();                                       ///< The function returns the dimensions of the node, ShapeBase#Dim.
    std::array<int,3>           getIdentifierColour();							///< The unique [r,g,b] identifier colour of the element, utilised in picking in the user interface.
    std::array<double,3>        getCentre();									///< This function returns the centre of the element in world spave.
    double                      getPeripodialness();							///< This function returns the relative influence of the peripodial physical characteristics to this element.
    double                      getColumnarness();								///< This function returns the relative influence of the columnar  physical characteristics ot this element.
    void                        getRelativePositionInTissueInGridIndex(int nGridX, int nGridY, int& IndexX, int& IndexY, double& FracX, double& FracY); 	///< Provides the relative position within the bounding box of the tissue, and calculates which point on the growth maps should be read.
    void                        getInitialRelativePositionInTissueInGridIndex(int nGridX, int nGridY, int& IndexX, int& IndexY, double& FracX, double& FracY); 	///< The relative positions of the iitial configuration of the element within the bounding box of the tissue.
    double                      getStiffnessMultiplier();							///< This function returns the current stiffness multiplier as a result of perturbations to physical properties.
    double                      getCurrentVolume();								///< This function returns the current volume of the element.
    double                      getElementalElasticForce(int nodeIndex, int dimIndex);				///< This function returns the calculated elastic force for the node nodeIndex, in coordinate dimIndex.
    void                        setElementalElasticForce(int nodeIndex, int dimIndex, double value);		///< This functions sets the calculated elastic force for the node nodeIndex, dimension dimIndex, to the input value.

    gsl_matrix*             getCurrentFe();									///< This function returns the current elastic part of the deformation gradient.
    double                  getApicalArea();								///< This function returns the current apical area of the element.
    void                    relaxElasticForces();								///< This function relaxes all teh accumulated elastic forces in the system.
    bool                    isGrowthRateApplicable(int sourceTissue, double& weight, double zmin, double zmax);	///< The function checks if the element if affected by the current growth functions.
    void                    updateGrowthWillBeScaledDueToApikobasalRedistribution(bool thisFunctionShrinksApical, std::vector<int>& ellipseBandIdsForGrowthRedistribution); ///< This function decide if the growth will be redistirbuted in the height of the tissue.
    void                    scaleGrowthForZRedistribution( double& x, double& y, double& z);			///< This function will modify the incremental growth deformation gradient of the element to reflect the volume redistribution in the height of the tissue.
    void                    calculateFgFromRates(double dt, double x, double y, double z, gsl_matrix* rotMat, gsl_matrix* increment, int sourceTissue, double zMin, double zMax);	///< This fucntion will calculate the incremental growth deformation gradient change for the current time step, from input growth rates
    void                    calculateFgFromGridCorners(int gridGrowthsInterpolationType, double dt, GrowthFunctionBase* currGF, gsl_matrix* increment, int sourceTissue, int IndexX, int IndexY, double FracX, double dFracY); ///< This fucntion will calculate the incremental growth deformation gradient change for the current time step by reading it from the grid, and interpolating on 4 corners.
    gsl_matrix*             getGrowthIncrement();								///< This function will return the current growth deformation gradient increment
    void                    updateGrowthIncrement(gsl_matrix* columnar, gsl_matrix* peripodial);		///< This function will update the elemental growth deformation gradient from the current growth deformation gradient increment.
    void                    updateGrowthByMutation(double dt);						///< This function will update the growth growth deformaton gradient increment of the element due to a mutation.
    void                    scaleGrowthIncrement(double multiuplier);					///< This function will scale the growth growth deformaton gradient increment by input double.
    void                    calculateShapeChangeIncrementFromRates(double dt, double rx, double ry, double rz, gsl_matrix* increment); ///< This function will calculate the current shape change deformation gradient increment due to elemental active shape change, from input rates.
    void                    updateShapeChangeIncrement(gsl_matrix* columnarShapeChangeIncrement);		///< This function updates the growth increment of the element with the current shape change increment
    void                    calculateRelativePosInBoundingBox(double boundingBoxXMin, double boundingBoxYMin, double boundingBoxLength, double boundingBoxWidth); ///< This function will calculate the relative positions of the element in the xy-plane bounding box of the tissue.
    void                    mutateElement(double growthFold, double growthRatePerHour);			///< This function will set the element as a mutant, the mutant growth rates will be set.
    void                    updateReferencePositionMatrixFromInput(double** input);				///< This function will update the reference position matrix. Not used under healthy, continuous simulations, to preserve continuity of the mesh.
    void                    displayRelativePosInBoundingBox();						///< Helper function to display the relative position of the element in the bounding box of the tissue.
    std::array<double,2>	getRelativePosInBoundingBox();							///< This function will return the relative position of the element in the xy bounding box of the tissue
    void 			setRelativePosInBoundingBox(double x, double y);				///< This functoin will set the relative position of the element in the bounding box of the tissue to the input coordinates.
    void			setInitialRelativePosInBoundingBox();						///< This function sets the initial relative position in bounding box of the tissue to current reference position.
    void 			setInitialZPosition(double zMax, double TissueHeight);				///< This fucntion sets the initial relative z position of the tissue in tissue height.
    std::array<double,2>	getInitialRelativePosInBoundingBox();						///< This function will return the initial relative position of the element in the xy bounding box of the tissue
    void 			convertRelativePosToGridIndex(std::array<double,2> relpos, int& indexX, int &indexY, double &fracX, double &fracY, int nGridX, int nGridY); ///< This function will convert the relative position of the tissue in xy plane bounding box to growth map grid indices.
    void 			getStrain(int type, float &StrainMag);						///< This function will return the selected strain component of the element.
    void 			getNodeBasedPysProp(int type, int NodeNo, const std::vector<std::unique_ptr<Node>>& Nodes, float& PysPropMag); ///< This function will return the selected physical properties of the element on a nodal basis.
    void 			getPysProp(int type, float &PysPropMag, double dt);				///< This function will return the selected physical properties of the element.
    double			getInternalViscosity();								///< This function will return the internal viscosity of the element.
    double			getOriginalInternalViscosity();							///<  This function will return the internal viscosity of the element prior to any perturbations.
    //void   			updateInternalViscosityTest();
    double 			getYoungModulus();								///< This function will return the Young's modulus of the element
    double 			getPoissonRatio();								///< This function will return the Poissons's ratio of the element
    const std::array<double,3>& getGrowthRate();								///< This function will return the current growth rate of the element.
    const std::array<double,6>& getShapeChangeRate();								///< This function will return the current shape change rate of the element.
    const std::vector<std::array<double,3>>& getReferencePos();							///< This function will return the reference positions of the element.
    void    			getPos(gsl_matrix* Pos);							///< This function will write the position of the element into input matrix.
    gsl_matrix* 		getFg();									///< This function will return the growth component of the deformation gradient.
    gsl_matrix* 		getInvFg();									///< This function calculates the inverse of the growth deformation gradient matrix.
    gsl_matrix* 		getFsc();									///< This function will return the shape change component of the deformation gradient.
    gsl_matrix* 		getInvFsc();									///< This function calculates the inverse of the shape change deformation gradient matrix. 
    gsl_matrix* 		getFe();									///< This function will return the elastic component of the deformation gradient.
    double 			getZRemodellingSoFar();								///< This function will return the z remodelling applied to the element so far, to cap the z remodelling.
    void 			setZRemodellingSoFar(double zRemodellingSoFar); 				///< This function will set the z remodelling applied to the element so far to input value. This is needed during saved input reading.
    void 			displayName();									///< Helper function, display the name of the element
    void			displayNodeIds();								///< Helper function, display the Ids of the nodes of the element
    void 			displayPositions();								///< Helper function, display the nodal positions of the element
    void 			displayReferencePositions();							///< Helper function, display the nodal positions of the reference element
    void 			displayIdentifierColour();							///< Helper function, display the unique identifier colour, for picking in user interface.
    void  		  	setFg(gsl_matrix* currFg);							///< This function sets the current grwoth deformation gradient matrix equal to input matrix 
    void 			setGrowthWeightsViaTissuePlacement (double periWeight);				///< This function sets the weight for growth rate scaling depending on tissue type
    void 			setYoungsModulus(double E);							///< This function sets the Young's modulus of the shape to inout double.
    virtual void 		setElasticProperties(double /*EApical*/,double /*EBasal*/, double /*EMid*/, double /*EECM*/, double /*v*/){ParentErrorMessage("setElasticProperties");} ///< This is the parent virtual function for setting up the elasticity properties of the shape depending on its tissue type placement.
    virtual void 		checkEdgeLenghtsForBinding(std::vector<int>& /*masterIds*/, std::vector<int>& /*slaveIds*/){ParentErrorMessage("checkEdgeLenghtsForBinding");} ///< The virt
    void 			setViscosity(double viscosityApical,double viscosityBasal, double viscosityMid);///< This function sets the viscosity of the element depending on its placement in the tissue.
    void 			setViscosity(double viscosityApical,double viscosityBasal);			///< This function sets the viscosity of the element depending on its placement in the tissue.
    void 			setViscosity(double viscosity);							///< This function sets the viscosity of the element depending on its placement in the tissue.
    double			calculateEmergentShapeOrientation();						///< This fucntion calculates the orientation of the emergent shape of an element in the xy plane of the tissue.
    bool 			isActinStiffnessChangeAppliedToElement(bool ThereIsWholeTissueStiffnessPerturbation, bool ThereIsApicalStiffnessPerturbation, bool ThereIsBasalStiffnessPerturbation, bool ThereIsBasolateralWithApicalRelaxationStiffnessPerturbation, bool ThereIsBasolateralStiffnessPerturbation, std::vector <int> &stiffnessPerturbationEllipseBandIds, int numberOfStiffnessPerturbationAppliesEllipseBands ); ///< This function decides if the actin stiffness perturbation is applied to this element.
    bool 			isECMChangeAppliedToElement(bool changeApicalECM, bool changeBasalECM, std::vector<int> &ECMChangeEllipseBandIds, int numberOfECMChangeEllipseBands);				///< This function decides if the ECM perturbation is applied to this element
    bool 			isShapeChangeAppliedToElement(std::vector<int> &ellipseBandIds, bool applyBasalECM, bool applyToLateralECM, bool applyApically, bool applyBasally, bool applyMidLayer );	///< This function decides if the shape change perturbation is applied to this element
    void 			calculateStiffnessPerturbationRate(bool ThereIsBasolateralWithApicalRelaxationStiffnessPerturbation, double stiffnessPerturbationBeginTimeInSec, double stiffnessPerturbationEndTimeInSec, double stiffnessChangedToFractionOfOriginal); ///< This function will calciulate the stiffness perturbation rate.
    void 			updateStiffnessMultiplier(double dt); ///< The function will update the actin multiplier as a result of stiffness perturbations.
    virtual std::array<double,3>  calculateBasalNormal(){ParentErrorMessage("calculateBasalNormal"); std::array<double,3> dummy = {0.0, 0.0, 0.0}; return dummy;} 	///< The virtual function of the parent for basal normal calculation. The value is dependent on node topology of the element and defined for eac individual child class.
    virtual void 		calculateApicalNormalCurrentShape(){ParentErrorMessage("calculateApicalNormal");}							///< The virtual function of the parent for apical normal calculation. The value is dependent on node topology of the element and defined for eac individual child class.
    void 			calculateCurrentGrowthIncrement(gsl_matrix* resultingGrowthIncrement, double dt, double growthx, double growthy, double growthz, gsl_matrix* ShearAngleRotationMatrix); ///< The function calculates the current growth increment from the input of growth rate and the orientation rotation matrix.
    void 			updateShapeChangeRate(double x, double y, double z, double xy, double yz, double xz);
    virtual void 		calculateElementShapeFunctionDerivatives(){ParentErrorMessage("calculateElementShapeFunctionDerivatives");}	///< The virtual function of the parent for calculation of shape function derivatives. This is topology dependent and implemented in each child shape type.
    virtual void 		calculateCurrNodalForces(gsl_matrix */*gslcurrge*/, gsl_matrix */*gslcurrgv*/, gsl_matrix */*gslcurrF*/, gsl_matrix* /*displacementPerDt*/, int /*pointNo*/){ParentErrorMessage("calculateCurrNodalForces");} ///< The virtual function to calculate nodal force. This is topology dependent and implemented in each child shape type.
    virtual void 		calculateCurrTriPointFForRotation(gsl_matrix */*currF*/,int /*pointNo*/){ParentErrorMessage("calculateCurrTriPointFForRotation");}
    virtual void 		calculateApicalArea(){ParentErrorMessage("calculateApicalArea");}                           ///< The virtual function to calculate apical area of element. This is topology dependent and implemented in each child shape type.
    virtual void 		calculateBasalArea(){ParentErrorMessage("calculateBasalArea");}                             ///< The virtual function to calculate basal area of element. This is topology dependent and implemented in each child shape type.
    double 			calculateCurrentGrownAndEmergentVolumes();                                                      ///< This is the function to calculate hte current ideal volume of the element and its current apparent volume.
    virtual void 		updateElasticProperties(){ParentErrorMessage("updateElasticProperties");}                   ///< This functions updates elastic propertiesand their dependent tensors upon alteration of a physical property.
    void 			writeInternalForcesTogeAndgv(gsl_matrix* ge, gsl_matrix* gvInternal, std::vector<std::array<double,3>>& SystemForces, const std::vector<std::unique_ptr<Node>>& Nodes); ///< This function writes the elemental elastic and viscous forces to the system scale force vector.
    void 			calculateForces(const std::vector<std::unique_ptr<Node>>& Nodes, gsl_matrix* displacementPerDt);
    void 			updatePositions(const std::vector<std::unique_ptr<Node>>& Nodes);                               ///< This function updates the position array of the element from the updated nodal posiitons.
    void 			setGrowthRate(double dt, double rx, double ry, double rz);                                      ///< This function sets the growth of the element from the tome step and the input rates
    void 			setGrowthRateViaInputTimeMultipliedMagnitude(double x, double y, double z);                     ///< This function sets the growht rate to pre-calculated rates given as input.
    void 			updateGrowthIncrementFromRate();                                                                ///< This function fills th egrowht increment matirx from the current growth rate matrix
    double  		calculateVolumeForInputShapeStructure(std::vector<std::array<double,3>> shapePositions, size_t nTriangularFaces, std::vector<std::array<int,3>>  triangularFaces, std::array<double,3> midPoint ); ///< This function calculates the volume of the shape, it is generalised such that the shape is defined as an array of triengles forming a convex hull.
    void 			calculatePrincipalStrains3D(double& e1, double &e2,  double &e3, gsl_matrix* eigenVec);         ///< This function calculates the principal components of the strains through eigen values and eigen vectors in 3D.
    void 			calculatePrincipalStrains2D(double& e1, double &e2,  double &e3, gsl_matrix* eigenVec);         ///< This function calculates the principal components of the strains through eigen values and eigen vectors in 2D.
    //void 			calculatePrincipalStrainAxesOnXYPlane(double& e1, double &e2, double& tet);			///< This function princial 
   // bool			checkIfXYPlaneStrainAboveThreshold(double thres);
   // bool 			calculateIfInsideActiveStripe(double initialPoint,double endPoint, double stripeSize1, double stripeSize2);
    void 			setShapeChangeRate(double x, double y, double z, double xy, double yz, double xz);                  ///< This function sets the shape change rate to pre-calculated rates given as input.
    void  			setShapeChangeInrementToIdentity();                                                                 ///< This function sets the shape change deformation gradient increment to identity
    void 			updateElementVolumesAndTissuePlacementsForSave(const std::vector<std::unique_ptr<Node>>& Nodes);    ///< This function calculates the reference volume, tissue placement and tissue type from nodal information
    bool 			readNodeIdData(std::ifstream& file);                                                                ///< This function reads in the node Ids for the element from save file
    bool			readReferencePositionData(std::ifstream& file);                                                     ///< This function reads the reference element positions from save file

    bool 			areanyOfMyNodesAtCircumference(const std::vector<std::unique_ptr<Node> > &Nodes);                   ///< This function check if the element owns any node at tissue circumference

    virtual void 		 checkHealth(){ParentErrorMessage("checkHealth");}                                              ///< The virtual function in parent to check if element is flipped, implemented for each child as it is node topology dependent.
    void    			writeKelasticToMainKatrix(gsl_matrix* K);                                                       ///< This function writes the elemental elastic component of the Jacobian to system Jacobian.
    void    			writeKviscousToMainKatrix(gsl_matrix* K);                                                       ///< This function writes the elemental viscous component of the Jacobian to system Jacobian.
    void    			calculateImplicitKElastic();                                                                    ///< This function calculates the elemental elastic component of the Jacobian for implicit NR itaration.
    void 			calculateImplicitKViscous(gsl_matrix* displacementPerDt, double dt);                                ///< This function calculates the elemental viscous component of the Jacobian for implicit NR itaration.
    void			calculateForceFromStress(int nodeId, gsl_matrix* Externalstress, gsl_matrix* ExternalNodalForces);  ///< This function calculates the elemental nodal forces from an input external stress matrix


    void 			updateShapeFromSave(std::ifstream& file);                                           ///< This function updates the element geometric poroerties from save file.
    void 			displayMatrix(boost::numeric::ublas::matrix<double>& mat, std::string matname);     ///< Helper function, displays the input blas (double) matrix with the input name.
    void 			displayMatrix(boost::numeric::ublas::matrix<int>& mat, std::string matname);        ///< Helper function, displays the input blas (int) matrix with the input name.
    void 			displayMatrix(boost::numeric::ublas::vector<double>& vec, std::string matname);     ///< Helper function, displays the input blas (double) vector with the input name.
    void 			displayMatrix(gsl_matrix* mat, std::string matname);                                ///< Helper function, displays the input gsl matrix with the input name.
    void 			displayMatrix(gsl_vector* mat, std::string matname);                                ///< Helper function, displays the input gsl vector with the input name.
    void 			createMatrixCopy(gsl_matrix *dest, gsl_matrix* src);                                ///< Helper function, creates a copy of the gsl matrix on new memory locaiton.
    double  		calculateMagnitudeVector3D(std::array<double,3> v);                                 ///< Helper algebraic function, calculates norm of the vector defined in the array<double,3>;
    void			normaliseVector3D(gsl_vector* v);                                                   ///< Helper algebraic function, normalises the input gsl vector (the input vector is modified)
    double			normaliseVector3D(std::array<double,3>& v);                                         ///< Helper algebraic function, normalises the input array<double,3> (the input vector is modified)
    double			getNormVector3D(gsl_vector* v);                                                     ///< Helper algebraic function, calculates norm of the vector defined in the gsl vector, the vector is not modified.
    //double 			determinant3by3Matrix(double* rotMat);                                              ///< Helper algebraic function, calculates determinant of 3by3 matrix stored in the double array pointed by the input pointer
    double 			determinant3by3Matrix(boost::numeric::ublas::matrix<double>& Mat);                  ///< Helper algebraic function, calculates determinant of 3by3 boost matrix
    double 			determinant3by3Matrix(gsl_matrix* Mat);                                             ///< Helper algebraic function, calculates determinant of 3by3 gsl matrix
    double 			determinant2by2Matrix(boost::numeric::ublas::matrix<double>& Mat);                  ///< Helper algebraic function, calculates determinant of 2by2 boost matrix
    void			calculateRotationAngleSinCos(std::array<double,3>& u, std::array<double,3>& v, double& c, double& s);       /// Helper algebraic function calculates the sine and cosine of the rotation angle needed to align vector u onto v.
    void    		calculateRotationAxis(const std::array<double,3>& u, const std::array<double,3> &v, std::array<double,3>& rotAx, double c); /// Helper algebraic function calculates the rotation axis needed to align vector u onto v.
    void			constructRotationMatrix(double c, double s, double* rotAx, double* rotMat);         ///< Helper algebraic functio calculates the rotation matrix from input sine, cosine of the rotation angle and the rotation axis. Writes the matrix into the input rotMat
    void			constructRotationMatrix(double c, double s, std::array<double,3>& rotAx, std::array<double,9>& rotMat); ///< Helper algebraic functio calculates the rotation matrix from input sine, cosine of the rotation angle and the rotation axis. Writes the matrix into the input rotMat
    void			rotateVectorByRotationMatrix(double* u,double* rotMat);                             ///< Helper algebraic function rotates the input vector v by rotation matrix rotMat
    void			rotateVectorByRotationMatrix(double* u,gsl_matrix* rotMat);                         ///< Helper algebraic function rotates the input vector v by rotation matrix rotMat
    void			rotateVectorByRotationMatrix(std::array<double,3>& u, std::array<double,9> rotMat); ///< Helper algebraic function rotates the input vector v by rotation matrix rotMat
    void    		CalculateGrowthRotationByF();                                                       ///< This function calculates the rigid body rotation of the element around the z axis of the tissue from hte fecormation gradient.
    void 			calculateTriPointFForRatation();                                                    ///< This function calculates the current deformaiton gradient as averaged at all Gauss points, for rigid body rotation extraction
    void 			setPlasticDeformationIncrement(double xx, double yy, double zz);                    ///< This function sets diagonal of the plastic deformation gradient increment from input values
    void 			growShapeByFg();                                                                    ///< This function updates the current growth deformaiton gradient with the growt/shape change/plastic deformation increments and their respective rotations.
    void 			changeShapeByFsc(double dt);                                                        ///< This function calculates the shape change increment from shape change rates
    void			checkIfInsideEllipseBands(int nMarkerEllipseRanges, std::vector<double> markerEllipseBandXCentres, std::vector<double> markerEllipseBandR1Ranges, std::vector<double> markerEllipseBandR2Ranges, const std::vector<std::unique_ptr<Node>>& Nodes); ///< This function checks if the element is inside any marker bands for perturbatins.
    bool			checkZCappingInRemodelling(bool volumeConserved, double zRemodellingLowerThreshold, double zRemodellingUpperThreshold, gsl_matrix* increment, gsl_matrix* eigenVec);    ///< This function checks if the remodelling of the element in z axis have reached the specified cap.
    void			assignSoftHinge(double lowHingeLimit, double highHingeLimit,double softnessLevel);  ///< This function modulates the stiffness of the hinge domain of the tissue with the input level. The domain is defined in relative x position boundaries.

    void			calculatePlasticDeformation3D(bool volumeConserved, double dt, double plasticDeformationHalfLife, double zRemodellingLowerThreshold, double zRemodellingUpperThreshold); ///< This function calculates the plastic deformation (remodelling) from the current elastic deformation gradient.
    //void 			addMigrationIncrementToGrowthIncrement(gsl_matrix* migrationIncrement);
    void 			displayDebuggingMatrices();                                                                     ///< This function displays a selected set of matricex for debugging purposes.
    virtual double 		getApicalSideLengthAverage(){return ParentErrorMessage("getApicalSideLengthAverage",0.0);}  ///< The virtual function of the parent to calculate average apical side length, dependent on nodal topology, defined in each child.
    virtual double 		getBasalSideLengthAverage(){return ParentErrorMessage("getBasalSideLengthAverage",0.0);}    ///< The virtual function of the parent to calculate average basal side length, dependent on nodal topology, defined in each child.
    virtual int 		getCorrecpondingApical(int /*currNodeId*/){return ParentErrorMessage("getCorrecpondingApical", -100);}  ///< The virtual function of the parent to obtain the corresponding apical node of a basal node, dependent on nodal topology, defined in each child.
    virtual bool 		IsThisNodeMyBasal(int /*currNodeId*/){return ParentErrorMessage("IsThisNodeMyBasal", false);}   ///< The virtual function of the parent to check if the input node ID is a basal node of the element, dependent on nodal topology, defined in each child.
    virtual bool 		IsThisNodeMyApical(int /*currNodeId*/){return ParentErrorMessage("IsThisNodeMyApical", false);} ///< The virtual function of the parent to check if the input node ID is an apical node of the element, dependent on nodal topology, defined in each child.
    virtual double 		getElementHeight(){return ParentErrorMessage("getElementHeight", 0.0);}                         ///< //< The virtual function of the parent to calculate z height of the element, dependent on nodal topology, defined in each child.
    virtual void		constructElementStackList(const int /*discretisationLayers*/, const std::vector<std::unique_ptr<ShapeBase>>& /*elementsList*/){ParentErrorMessage("constructElementStackList");}
    virtual void 		checkRotationConsistency3D(){ParentErrorMessage("checkRotationConsistency3D");}                                         ///< The virtual function of the parent to check if the two input nodes of the element are directly connected on one of the elemental surfaces, dependent on nodal topology, defined in each child.
    virtual bool 		areNodesDirectlyConnected(int /*node0*/, int /*node1*/){return ParentErrorMessage("areNodesDirectlyConnected",false);} ///< The virtual function of the parent to check if the rotation of the nodes of the element are consistent, dependent on nodal topology, defined in each child.
    bool 			DoesPointBelogToMe(int IdNode);                             ///< This function checks if the input node belogs to the element.
    void 			assignVolumesToNodes(const std::vector<std::unique_ptr<Node>>& Nodes);  ///< This function distributes element's total volume among its owner nodes.
    void 			calculateZProjectedAreas();                                         ///< This function calculated the z-projected (to world xy plane) apical and basal areas of the element.
    void 			assignZProjectedAreas(const std::vector<std::unique_ptr<Node>>& Nodes); ///< This function distributes element's z-projected areas among its owner nodes.
    void 			assignElementToConnectedNodes(const std::vector<std::unique_ptr<Node>>& Nodes); ///< This function assigns the element to the nodes it owns, necessary to construc the owner and connectivity list of nodes.
    void 			setECMMimicing(bool IsECMMimicing);                         ///< This funciton sets the element as an ECM mimicking element (distinct domain in terms of physical characteristics).
    void 			setActinMimicing(bool isActinMimicing);                     ///< This funciton sets the element as an actin mimicking element (distinct domain in terms of physical characteristics).
    virtual void 		assignExposedSurfaceAreaIndices(const std::vector<std::unique_ptr<Node>>& /*Nodes*/){ParentErrorMessage("assignExposedSurfaceAreaIndices");} ///< The virtual function on parent assigns the nodes of the surfaces that are exposed to external world, dependent on topology, defined for each child.
    void 			calculateViscositySurfaces();                               ///< This function calls for the assignment of exposed surfaces if the element has viscosity
    void 			assignViscositySurfaceAreaToNodes(const std::vector<std::unique_ptr<Node>>& Nodes); ///< This function distributes elemenbt's exposed surfaces to nodes.

    void 			calculateEmergentRotationAngles();                          ///< This function calculates the emergent rotation of the element in xy plane for display purposes
    void 			updateReferencePositionMatrixFromMeshInput(std::ifstream& file);    ///< This function updates the reference position of the element from save file
    void			fillNodeNeighbourhood(const std::vector<std::unique_ptr<Node>>& Nodes); ///< This function fills in the node neightbourhood, needed for constrction of the connectivity of nodes
    double 			dotProduct3D(std::array<double,3>& u, std::array<double,3>& v); ///< Helper algebraic function, calculates dot product of two arrays <double,3>
    void			crossProduct3D(gsl_vector* u, gsl_vector* v, gsl_vector* cross); ///< Helper algebraic function, calculates cross product of two gsl_vectors, writes into the third input gls vector.
    std::array<double,3> 	crossProduct3D(std::array<double,3> u, std::array<double,3> v); ///< Helper algebraic function, calculates cross product of two arrays <double,3>
    virtual void 		setBasalNeigElementId(const std::vector<std::unique_ptr<ShapeBase>>& /*elementsList*/){ParentErrorMessage("setBasalNeigElementId");}
    bool 			isElementFlippedInPotentialNewShape(int nodeId, double newX, double newY, double newZ); ///< This function checks if the element will plip in the case that its node (nodeID) is moved to the new x,y,z coordintes specified in the input. Necessary in node collapsing.
    void 			checkForCollapsedNodes(int TissueHeightDiscretisationLayers, const std::vector<std::unique_ptr<Node>>& Nodes, const std::vector<std::unique_ptr<ShapeBase>>& Elements); ///< This function checks if any of the edges of the element is shortened to the extent that it should be collapsed.
    bool 			hasEnoughNodesOnCurve(const std::vector<std::unique_ptr<Node> > &Nodes);    ///< This function checks if the majority of the nodes of teh element reside in a curved region, to assign it to specific curvature dependent perturbations.
    void 			assignEllipseBandIdToWholeTissueColumn(size_t TissueHeightDiscretisationLayers, const std::vector<std::unique_ptr<Node>>& Nodes, const std::vector<std::unique_ptr<ShapeBase>>& Elements); ///< This function assigns the marker ID of the apical elemetn ot all its connected elements in the tissue hight (all column of the element).
    void 			assignEllipseBandId(const std::vector<std::unique_ptr<Node>>& Nodes, int selectedEllipseBandId); ///< This function assigns the marking ellipse band ID of the element depending on the definition of nodes it is consturcted of.
    void 			assignEllipseBandIdToNodes(const std::vector<std::unique_ptr<Node>>& Nodes);    ///< This function assigns the marker ID of the elemetn to all its nodes.
    void 			addToElementalElasticSystemForces(int i,int j,double value); 	/// This function is to add the input value, to the (i,j)th element of the ElementalElasticSystemForces
    void 			addToTriPointKe(int i,int j,double value); 			/// This function is to add the input value, to the (i,j)th element of the triPointKe
};

#endif
