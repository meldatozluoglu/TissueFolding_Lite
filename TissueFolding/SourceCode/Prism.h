#ifndef Prism_H
#define Prism_H


#include "ShapeBase.h"
#include <string>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>

class Prism : public ShapeBase{

protected:
    double initialApilcalEdgeLengthsSq[3];  ///< The array storing the squares of the apical edge lengths.
    double initialBasalEdgeLengthsSq[3];    ///< The array storing the squares of the apical edge lengths.

    void getCurrRelaxedShape(gsl_matrix * CurrRelaxedShape);                                        ///< This function writes the current relaxed shape of the elemetn in the gsl_matrix pointed by the provided input gsl_matrix pointer.
    void setShapeFunctionDerivatives(gsl_matrix * ShapeFuncDer,double eta, double zeta, double nu); ///< This function calculates and writes the shape function derivatives on the gsl_matrix pointed by the provided input gsl_matrix pointer, for the provided barycentric coordinates.
    void setShapeFunctionDerivativeStack(gsl_matrix* ShapeFuncDer, gsl_matrix* ShapeFuncDerStack);  ///<This function writes the input gsl_matrix ShapeFuncDer into a stack form, in the matrix ShapeFuncDerStack. Both matrices must be allocated prior to calling the function.
    void setCoeffMat();                                                                             ///< This function writes the coefficient matrix that is collating the calculated stresses in the form of the vector form of elemental stress and strain vectors.
    void calculateDVector();                                                                        ///<This function calculates the elasticity tensor for a Kirshoff material.
    void calculateD81Tensor();                                                                      ///<This function calculates the Lagrangian elasticity tensor for a neo-hookean material
    void calculateCurrk(boost::numeric::ublas::matrix<double> &currk, boost::numeric::ublas::matrix<double> &currB, boost::numeric::ublas::matrix<double>& currBE, boost::numeric::ublas::matrix<double> &currBo, double eta, double zeta, double nu); ///< This function calculates the elemental jacobian for current Gauss Point.
    void calculateCurrNodalForces(gsl_matrix *gslcurrge, gsl_matrix *gslcurrgv, gsl_matrix *gslcurrF, gsl_matrix* displacementPerDt, int pointNo); ///< This function calculates the nodal forces of teh prism for the current Gauss Point.
    void calculateCurrTriPointFForRotation(gsl_matrix *currF,int pointNo);  ///< This function calculates the deformation gradient of the prism for rigid body rotation rotation extractions, for the current Gauss Point.
    void calculateReferenceVolume();                                        ///< This function calculateds the volume of the reference element.
    std::array<std::array<double,3>,8> calculatePlaneNormals();             ///< This function calculates the plane normals for the triangular surfaces constructing the prism.
    std::array<double,3> assignNodalVector(size_t id0, size_t id1);         ///< This function calculates the 3D vector from nodes indexed at id0 to id1.
    bool checkNodePlaneConsistency(std::array<std::array<double,3>,8>& normals);    ///< This functions checks the planes of the elemetn as part of health check.
    void setInitialEdgeLenghts();                                           ///< This function sets the initial edge lengths of the element, will be essential in deciding if the elemetn should collapse its surfaces to avoid flipping.
    void checkEdgeLenghtsForBinding(std::vector<int>& masterIds, std::vector<int>& slaveIds);   ///< This function checks the edge lengths of the prism to decide if the elemetn should collapse its surfaces to avoid flipping.
    double getApicalSideLengthAverage();    ///< This function returns the average edge length of the apical surface.
    double getBasalSideLengthAverage();     ///< This function returns the average edge length of the basal surface.
    void assignExposedSurfaceAreaIndices(const std::vector<std::unique_ptr<Node>>& Nodes);  ///< This function assigns the indices of externally exposed surfaces of the prism.

public:
    Prism(const std::vector<int>& inpNodeIds, const std::vector<std::unique_ptr<Node>>& Nodes, int CurrId); ///< Constructor
    ~Prism(); ///< Destructor
    void  setElasticProperties(double EApical, double EBasal, double EMid,  double EECM, double v); ///< This function sets the elastic properties of the prism.
    std::array<double,3>  calculateBasalNormal();       ///< This function calculates the normal of the basal surface of the element.
    void  calculateApicalNormalCurrentShape();          ///< This function calculates the normal of the apical surface of the element.
    void  calculateElementShapeFunctionDerivatives();   ///< This function calculates the shape function derivatives of the prism.

    void  checkHealth(); ///< This function checks the health of the element, against fliiping.
    int getCorrecpondingApical(int currNodeId); ///< This function returns the corresponding connected apical node of a basal input node on the prism.
    bool IsThisNodeMyBasal(int currNodeId);     ///< This function checks if the input Node#Id is a basal node of the prism.
    bool IsThisNodeMyApical(int currNodeId);    ///< This function checks if the input Node#Id is an apical node of the prism.
    double getElementHeight();                  ///< This function returns the apical-basal (z) height of this prism.
    void calculateApicalArea();                 ///< This function calculates the apical area of the prism.
    void calculateBasalArea();                  ///< This function calculates the basal area of the prism.
    void updateElasticProperties();             ///< This function updates the elastic properties of the prism upon perturbation.
    void setBasalNeigElementId(const std::vector<std::unique_ptr<ShapeBase>>& elementsList); ///< This function sets the basal neighboiur of this element, which has an apical surface overlapping completely with teh basal surface of this prism.
    void constructElementStackList(const int discretisationLayers, const std::vector<std::unique_ptr<ShapeBase>>& elementsList); ///< This function constructs the apical-basal element stack list that this prosm resides in.
    void checkRotationConsistency3D();          ///< This function checks if the nodes of the prosm rotate counter-clock vise, and corrects if not.
    bool areNodesDirectlyConnected(int node0, int node1); ///< This function checks if the two nodes with input Node#Id values are directly connected to each other in topology of the prism.
};

#endif
