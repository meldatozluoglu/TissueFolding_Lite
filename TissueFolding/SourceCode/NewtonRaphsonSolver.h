#ifndef NewtonRaphsonSolver_H
#define NewtonRaphsonSolver_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <gsl/gsl_linalg.h>
#include "Node.h"
#include "ShapeBase.h"
//#include <omp.h>

/**
 *  The Newton-Raphson solver class
 *  */
class NewtonRaphsonSolver{
private:
    size_t nDim;                        ///< Dimension of the space (3D)
    size_t nNodes;                      ///< Number of nodes of the system
    double threshold;                   ///< Convergence threshold for iterations

public:
    NewtonRaphsonSolver(int nDim, int nNodes); 	///< Constructer of the N-R solver
    ~NewtonRaphsonSolver();						///< Desturctor of the N-R solver

    gsl_matrix* un;								///< The initial positions of the nodes, as calculated at the end of previous step "n"
    gsl_matrix* ge;								///< The matrix containing elastic forces on each node, size (nDim*nNodes,1). Organisation is [Node0,x ; Node0,y ; Node0,z; ... ; Noden,x ; Noden,y ; Noden,z]
    gsl_matrix* gvInternal;						///< The matrix containing internal viscous forces on each node, size (nDim*nNodes,1). Organisation is [Node0,x ; Node0,y ; Node0,z; ... ; Noden,x ; Noden,y ; Noden,z]
    gsl_matrix* gvExternal;						///< The matrix containing external viscous forces on each node, size (nDim*nNodes,1). Organisation is [Node0,x ; Node0,y ; Node0,z; ... ; Noden,x ; Noden,y ; Noden,z]
    gsl_matrix* gExt;							///< The matrix containing external forces on each node (currently includes packing forces), size (nDim*nNodes,1). Organisation is [Node0,x ; Node0,y ; Node0,z; ... ; Noden,x ; Noden,y ; Noden,z]
    gsl_vector* gSum;							///< The matrix containing sum of NewtonRaphsonSolver#ge, NewtonRaphsonSolver#gvInternal, NewtonRaphsonSolver#gvExternal, NewtonRaphsonSolver#gExt. Organisation is [Node0,x ; Node0,y ; Node0,z; ... ; Noden,x ; Noden,y ; Noden,z]
    gsl_matrix* uk;								///< The matrix storing the position of each node at iteration "k". Initiated in function NewtonRaphsonSolver#initialteUkMatrix at the beginning of each step, and updated by function NewtonRaphsonSolver#updateUkInIteration during the iteartions.
    gsl_matrix* displacementPerDt;				///< The displacement per time step of each node in current iteration "k", from its position at the end of the last time step "n"
    gsl_vector* deltaU;							///< The incremental change in positions as calculated in current iteration, resulting from the imbalance of elastic, viscous and any other external forces acting on each nodes. The solver minimises this value, convergence occurs when all incremental movements for all nodes sufficiently close to zero.
    gsl_matrix* K;								///< The Jacobian matrix, derivative of sum of forces acting on each node with respect to displacements.
    bool boundNodesWithSlaveMasterDefinition;   ///< The boolean stating if there are degrees of freedom slave to other nodes (masters).
    std::vector< std::vector<int> > slaveMasterList;    ///< The 2D integer vector storing the slave-master degrees of freedom couples, such that array [i][0] = slave to array[i][1], each i representing one dim of node position ( z of node 2 is DoF i=8);

    void setMatricesToZeroAtTheBeginningOfIteration(); 		///< The function setting the calculation matrices to zero at the beginning of each iteration.
    void setMatricesToZeroInsideIteration();												///< The function setting the relevant matrices to zero at each iteration.
    void constructUnMatrix(const std::vector<std::unique_ptr<Node> > &Nodes);               ///< This function constructs NewtonRaphsonSolver#un matrix at the beginning of the iterations.
    void initialteUkMatrix();                                                               ///< This function initiates NewtonRaphsonSolver#uk matrix at the beginning of the iterations, it is initiated to be equal to NewtonRaphsonSolver#un.
    void calculateBoundKWithSlavesMasterDoF();                                              ///< This function updates the Jacobian of the system, NewtonRaphsonSolver#K, to reflect degrees of freedom binding.
    void equateSlaveDisplacementsToMasters();                                               ///< This function moves the slaves of bound couples with a displacement equivalent to the masters'.
    void calculateDisplacementMatrix(double dt);											///< This function calculates the displacement of each node in current iteration "k", from their positions at the end of the previous step "n" (NewtonRaphsonSolver#uk - NewtonRaphsonSolver#un)
    void calcutateFixedK(const std::vector <std::unique_ptr<Node>>& Nodes);					///< This function updates the Jacobian to account for nodes  that are fixed in certain dimensions in space, as part of boundary conditions.
    void calculateForcesAndJacobianMatrixNR(const std::vector <std::unique_ptr<Node>>& Nodes, const std::vector <std::unique_ptr<ShapeBase>>&  Elements, double dt );	///< This function calculates elemental forces and Jacobians, later to be combined in NewtonRaphsonSolver#K and NewtonRaphsonSolver#gSum
    void writeForcesTogeAndgvInternal(const std::vector <std::unique_ptr<Node>>& Nodes, const std::vector <std::unique_ptr<ShapeBase>>&  Elements, std::vector<std::array<double,3>>& SystemForces);	///< This function writes the values of elemental elastic (ShapeBase#ge) and internal viscous forces (ShapeBase#gvInternal) into the system elastic and internal viscous forces, NewtonRaphsonSolver#ge, and NewtonRaphsonSolver#gvInternal, respectively.
    void writeImplicitElementalKToJacobian(const std::vector <std::unique_ptr<ShapeBase>>&  Elements);	///< This function writes the elemental values for elastic part of the Jacobian - stiffness matrix - (ShapeBase#TriPointKe) and for viscous part of Jacobian (ShapeBase#TriPointKv) into the system Jacobian NewtonRaphsonSolver#K.
    void calculateExternalViscousForcesForNR(const std::vector <std::unique_ptr<Node>>& Nodes);         ///< This function calculates the external viscous forces acting on each node, the values are sotred in NewtonRaphsonSolver#gvExternal
    void addImplicitKViscousExternalToJacobian(const std::vector <std::unique_ptr<Node>>& Nodes, double dt);	///< This function adds the external related terms of the Jacobian to the system Jacobian NewtonRaphsonSolver#K.
    void checkJacobianForAblatedNodes(std::vector <int> & AblatedNodes);                    ///< This functions checks the Jacobian to ensure the diagonal terms are non-zero for ablated nodes.
    void calculateSumOfInternalForces();                                                    ///< This function adds the ealsticity and viscosity related forces (NewtonRaphsonSolver#ge, NewtonRaphsonSolver#gvInternal, NewtonRaphsonSolver#gvExternal) to sum of forces, NewtonRaphsonSolver#gSum.
	void addExernalForces();

    void solveForDeltaU();                                                                  ///< This function solves for the displacements within the N-R step.
    //raw pointers necessary for Pardiso
	int  solveWithPardiso(double* a, double*b, int* ia, int* ja, const int n_variables);
    void constructiaForPardiso(int* ia, const int nmult, std::vector<int> &ja_vec, std::vector<double> &a_vec);
    void writeKinPardisoFormat(const int nNonzero, std::vector<int> &ja_vec, std::vector<double> &a_vec, int* ja, double* a);
	void writeginPardisoFormat(double* b, const int n);
    bool checkConvergenceViaDeltaU();                                                      ///< Check for cenvergence with the norm of displacements vector,against the threshold NewtonRaphsonSolver#threshold.
    bool checkConvergenceViaForce();                                                       ///< Check for cenvergence with the norm of forces vector,against the threshold NewtonRaphsonSolver#threshold.
    void updateUkInIteration();                                                            ///< Calulate the nodal displacemetns at the kth iteration of NR dolver.

    void displayMatrix(gsl_matrix* mat, std::string matname);
    void displayMatrix(gsl_vector* mat, std::string matname);

	bool checkIfCombinationExists(int dofSlave, int dofMaster);
    void checkMasterUpdate(int& dofMaster, int& masterId); ///< This function takes a degree of freedom number as input. This DOF is supposed to be a master. If, the dof is already a slave to another dof, then update the masetr dof. Anything that would be bound to the input dof can be bound to the already existing master of the input dof.
    bool checkIfSlaveIsAlreadyMasterOfOthers(int dofSlave, int dofMaster); ///< This function checks if the slave DOF is already master of others, if so, updates the master of said slave to the new master the current slave will be bound to.

    void updateElementPositions(const std::vector<std::unique_ptr<Node> > &Nodes, const std::vector <std::unique_ptr<ShapeBase>>&  Elements);


};

#endif
