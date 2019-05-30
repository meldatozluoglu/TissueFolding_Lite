/*
 * NewtonRaphsonSolver.cpp
 *
 *  Created on: 26 Apr 2016
 *      Author: melda
 */

#include "NewtonRaphsonSolver.h"
#include <math.h>
//#include <gsl/gsl_linalg.h>

NewtonRaphsonSolver::NewtonRaphsonSolver(int dim, int n){
    threshold = 1E-4;	/**
     *  The function initiates (memory allocates) the necessary matrices for
     *  the numerical solver.
     *
     */

	nDim = dim;
	nNodes = n;
	un = gsl_matrix_calloc(nDim*nNodes,1);
	ge = gsl_matrix_calloc(nDim*nNodes,1);
	gvInternal = gsl_matrix_calloc(nDim*nNodes,1);
	gvExternal = gsl_matrix_calloc(nDim*nNodes,1);
	gExt = gsl_matrix_calloc(nDim*nNodes,1);
	gSum = gsl_vector_calloc(nDim*nNodes);
	uk = gsl_matrix_calloc(nDim*nNodes,1);
	displacementPerDt = gsl_matrix_calloc(nDim*nNodes,1);
	deltaU = gsl_vector_calloc(nDim*nNodes);
	K = gsl_matrix_calloc(nDim*nNodes,nDim*nNodes);
	boundNodesWithSlaveMasterDefinition = false;
}

NewtonRaphsonSolver::~NewtonRaphsonSolver(){
	gsl_matrix_free(uk);
    gsl_matrix_free(un);
	gsl_matrix_free(ge);
	gsl_matrix_free(gvExternal);
	gsl_matrix_free(gvInternal);
	gsl_matrix_free(gExt);
	gsl_vector_free(gSum);
	gsl_matrix_free(displacementPerDt);
	gsl_vector_free(deltaU);
    gsl_matrix_free(K);
}

void NewtonRaphsonSolver::setMatricesToZeroAtTheBeginningOfIteration(){
	gsl_matrix_set_zero(un);
	gsl_matrix_set_zero(ge);
	gsl_matrix_set_zero(gvInternal);
	gsl_matrix_set_zero(gvExternal);
	gsl_matrix_set_zero(gExt);
	gsl_vector_set_zero(gSum);
	gsl_matrix_set_zero(uk);
	gsl_matrix_set_zero(displacementPerDt);
	gsl_vector_set_zero(deltaU);
	gsl_matrix_set_zero(K);
}

void NewtonRaphsonSolver::setMatricesToZeroInsideIteration(){
	gsl_matrix_set_zero(ge);
	gsl_matrix_set_zero(gvInternal);
	gsl_matrix_set_zero(gvExternal);
	gsl_vector_set_zero(gSum);
	gsl_matrix_set_zero(gExt);
	gsl_matrix_set_zero(K);
}

void NewtonRaphsonSolver::constructUnMatrix(const std::vector <std::unique_ptr<Node>>& Nodes){
    for (size_t i = 0; i<nNodes; ++i ){
        for (size_t j=0; j<nDim; ++j){
            gsl_matrix_set(un,3*i+j,0,Nodes[i]->Position[j]);
        }
    }
}

void NewtonRaphsonSolver::initialteUkMatrix(){
    gsl_matrix_memcpy(uk,un);
}


void NewtonRaphsonSolver::calculateDisplacementMatrix(double dt){
    /**
     * The displacement of nodes per Simulation#dt, for each iteration defined as \f$ \frac{\boldsymbol{u_k} - \boldsymbol{u_n}}{dt} \f$
     */
	gsl_matrix_memcpy(displacementPerDt,uk);
	gsl_matrix_sub(displacementPerDt,un);
	gsl_matrix_scale(displacementPerDt,1.0/dt);
}

void NewtonRaphsonSolver::calculateBoundKWithSlavesMasterDoF(){
	if (boundNodesWithSlaveMasterDefinition){
        /** The matrix calculation for node binding is in the form: \n
         * \f{eqnarray*}
          \boldsymbol{K}_{bound} & = & \boldsymbol{N^{T}} \boldsymbol{K} \boldsymbol{N} + \boldsymbol{\bar{I}}\\
          \boldsymbol{g}_{bound} & = & \boldsymbol{N^{T}}  \boldsymbol{g} \\
          \boldsymbol{N_{ij}} &= & \left\{ \begin{matrix}
                                        1 & if		& i = j \text{ and } i \text{ is not a slave}\\
                                        1 & if 	& i  \text{ is slave to } j \\
                                        0 &		& \text{elsewhere}
                                    \end{matrix} \right. \\
          \boldsymbol{\bar{I}_{ij}} & = &  \left\{ \begin{matrix}
                                             1 & if		& i = j \text{ and } i \text{ is a slave}\\
                                             0 &		& \text{elsewhere}
                                           \end{matrix} \right. \\
         * \f} \n
         * As this actual matrix calculation requires a significant memory allocation, these operations are
         *carried out on a row/column basis, rather than using full matrices. \n
         */
		int totalNumberOfDoF = K->size1; //nDim * nNodes
        for (std::vector< std::vector<int> >::iterator itSlaveMasterCouple=slaveMasterList.begin(); itSlaveMasterCouple<slaveMasterList.end(); ++itSlaveMasterCouple){
            /** First all forces on slave to master, set slave forces to zero on NewtonRaphsonSolver#gSum,
             *equivalent of \f$ \boldsymbol{g}_{bound} = \boldsymbol{N^{T}}  \boldsymbol{g} \f$.
             */
			int slaveIndex = (*itSlaveMasterCouple)[0];
			int masterIndex = (*itSlaveMasterCouple)[1];
			double slaveForce = gsl_vector_get(gSum,slaveIndex);
			double masterForce = gsl_vector_get(gSum,masterIndex);
			masterForce +=  slaveForce;
			gsl_vector_set(gSum,slaveIndex,0);
			gsl_vector_set(gSum,masterIndex,masterForce);
            /** Then start the manipulation of the Jacobian, by adding all elements of the slave degrees of freedom
             * row to master degrees of freedom row on NewtonRaphsonSolver#K, equivalent of operation \f$ \boldsymbol{N^{T}} \boldsymbol{K} \f$.
             */
			// format of view: [ matrix*, origin i, origin j, rows, columns ]
			gsl_matrix_view KSlaveRow = gsl_matrix_submatrix (K, slaveIndex, 0, 1, totalNumberOfDoF);
			gsl_matrix_view KMasterRow = gsl_matrix_submatrix (K, masterIndex, 0, 1, totalNumberOfDoF);
			gsl_matrix_add(&(KMasterRow.matrix),&(KSlaveRow.matrix));
			gsl_matrix_set_zero(&(KSlaveRow.matrix));
			//add all elements of the slave column to master column - N^T K N
            /** Followed by adding all elements of the slave degrees of freedom column to master degrees of freedom
             *  columno n NewtonRaphsonSolver#K, equivalent of operation  (cumulatively) \f$ \boldsymbol{N^{T}} \boldsymbol{K} \boldsymbol{N} \f$.
             */
			gsl_matrix_view KSlaveColumn = gsl_matrix_submatrix (K, 0, slaveIndex, totalNumberOfDoF,1);
			gsl_matrix_view KMasterColumn = gsl_matrix_submatrix (K,0, masterIndex, totalNumberOfDoF,1);
			gsl_matrix_add(&(KMasterColumn.matrix),&(KSlaveColumn.matrix));
			gsl_matrix_set_zero(&(KSlaveColumn.matrix));
            /** Finally, make the diagonal element of NewtonRaphsonSolver#K, \f$ \boldsymbol{K}_{DOFslave,DOFslave} \f$ to unity,
             * equivalent of operation  (cumulatively)
             * \f$ \boldsymbol{N^{T}} \boldsymbol{K} \boldsymbol{N} + \boldsymbol{\bar{I}}\f$.
             */
			gsl_matrix_set(K,slaveIndex,slaveIndex,1);
            //gsl_matrix views do not need to be freed, they are addresses to original matrices
		}
	}
}

void NewtonRaphsonSolver::equateSlaveDisplacementsToMasters(){
	int n = slaveMasterList.size();
	for( int slaveIterator = 0; slaveIterator <n; ++slaveIterator){
		int slaveIndex = slaveMasterList[slaveIterator][0];
		int masterIndex = slaveMasterList[slaveIterator][1];
		double masterDisplacement = gsl_vector_get(deltaU,masterIndex);
		gsl_vector_set(deltaU,slaveIndex,masterDisplacement);
    }
}

void NewtonRaphsonSolver::calcutateFixedK(const std::vector <std::unique_ptr<Node>>& Nodes){
    /** Some degrees of freedom are fixed for some nodes, as defined by the user input boundary conditions.
     * this will be recorded in the 3 dimensional boolean array of each node Node#FixedPos
     * for x,y and z coordinates. In the node has a fixed degree of freedom, then the  sum of
     * elastic and viscous forces recorded on NewtonRaphsonSolver#gSum is made zero. Then in the Jacobian, the
     * diagonal term for the degree of freedom is set to unity, all remaining terms of the column and row of
     * the degree of freedom is set to zero.
     */
    size_t dim = 3;
    size_t Ksize = K->size1;
    for(size_t i=0; i<nNodes; i++){
        for (size_t j=0; j<dim; ++j){
            if (Nodes[i]->FixedPos[j]){
                size_t index1 = i*dim+j;
                gsl_vector_set(gSum,index1,0.0); // making the forces zero
                for (size_t k =0; k<Ksize; ++k){
                    double value =0.0;
                    if (index1 == k ){value =1.0;}
                    gsl_matrix_set(K, index1, k, value);
                    gsl_matrix_set(K, k, index1, value); //K is symmetric;
                }
            }
        }
    }
}

void NewtonRaphsonSolver::calculateForcesAndJacobianMatrixNR(const std::vector <std::unique_ptr<Node>>& Nodes, const std::vector <std::unique_ptr<ShapeBase>>& Elements, double dt){
    #pragma omp parallel for
    for(std::vector<std::unique_ptr<ShapeBase>>::const_iterator itElement = Elements.begin(); itElement < Elements.end(); ++itElement){
        /** The calculation of foces and their derivatives in each element starts with
         * calculation of forces via ShapeBase#calculateForces. A series of calculations necessary for
         * Jacobian calculation are obtained at this stage. Then the Elastic and Viscous parts of the elemetnal Jacobian
         * are calculates through ShapeBase#calculateImplicitKElastic and ShapeBase#calculateImplicitKViscous ,respectively
          */
        if (!(*itElement)->IsAblated){
            (*itElement)->calculateForces(Nodes, displacementPerDt);
        }
        //std::cout<<"finished calculating forces in NR"<<std::endl;
        (*itElement)->calculateImplicitKElastic(); //This is the stiffness matrix, elastic part of the jacobian matrix
        //std::cout<<"finished calculating ImplicitK elastic in NR"<<std::endl;
        (*itElement)->calculateImplicitKViscous(displacementPerDt, dt); //This is the viscous part of jacobian matrix
        //std::cout<<"finished calculating ImplicitK viscous in NR"<<std::endl;
    }
}

void NewtonRaphsonSolver::updateElementPositions(const std::vector <std::unique_ptr<Node>>& Nodes, const std::vector <std::unique_ptr<ShapeBase>>& Elements){
    for(const auto& itElement : Elements){
        itElement->updatePositions(Nodes);
	}
}

void NewtonRaphsonSolver::writeForcesTogeAndgvInternal(const std::vector <std::unique_ptr<Node>>& Nodes, const std::vector <std::unique_ptr<ShapeBase>>&  Elements, std::vector<std::array<double,3>>& SystemForces){
    for(const auto& itElement : Elements){
        if (!itElement->IsAblated){
            itElement->writeInternalForcesTogeAndgv(ge,gvInternal,SystemForces,Nodes);
        }
    }
}

void NewtonRaphsonSolver::writeImplicitElementalKToJacobian(const std::vector <std::unique_ptr<ShapeBase>>&  Elements){
    /** All elemental elastic (ShapeBase#Ke) and internal viscous (ShapeBase#Kv) jacobians are added onto the
     * system Jacobian (NewtonRaphsonSolver#K) in this function.
     */
    for(const auto& itElement : Elements){
        itElement->writeKviscousToMainKatrix(K);
    }
    for(const auto& itElement : Elements){
        itElement->writeKelasticToMainKatrix(K);
    }
}

void NewtonRaphsonSolver::calculateExternalViscousForcesForNR(const std::vector <std::unique_ptr<Node>>& Nodes){
    for (size_t i = 0; i<nNodes; ++i ){
        for (size_t j=0; j<nDim; ++j){
            /** For all nodes of the system, the external viscous forces will depend on the
             * exposed surface associated with each node (Node#viscositySurface),
             * calculated through ShapeBase::assignViscositySurfaceAreaToNodes. This surface
             * will be muliplied by the local viscosity Node#externalViscosity and the
             * displacement of the node in current iteration k, from its position at the end of
             * previous time step n.
             */
            double surfaceAreaTimesViscosity = Nodes[i]->viscositySurface*Nodes[i]->externalViscosity[j];
            double displacementValue = gsl_matrix_get(displacementPerDt,3*i+j,0);
            gsl_matrix_set(gvExternal,3*i+j,0,surfaceAreaTimesViscosity*displacementValue);
            if (std::isnan(surfaceAreaTimesViscosity)){
                std::cerr<<" node: "<<i<<" dimention: "<<j<<" surfaceAreaTimesViscosity is nan: "<<surfaceAreaTimesViscosity<<std::endl;
            }
            if (std::isnan(displacementValue)){
                std::cerr<<" node: "<<i<<" dimention: "<<j<<" displacementValue is nan: "<<displacementValue<<std::endl;
            }
		}
	}
    /** Once all the viscous forces are calculated, the force direction will be inverted, as I am interested in the visouc
     * drag applied by the meida to the node.
     */
    gsl_matrix_scale(gvExternal,-1.0);
}

void NewtonRaphsonSolver::addImplicitKViscousExternalToJacobian(const std::vector<std::unique_ptr<Node> > &Nodes, double dt){
    /** This function will add the derivatives of external viscous drag forces with respect to
     *nodal displacement onto the system Jacobian.
     *
     */
    for (size_t i = 0; i<nNodes; ++i ){
        double surfaceAreaPerDt = Nodes[i]->viscositySurface/dt;
        for (size_t j=0; j<nDim; ++j){
            double curKValue = gsl_matrix_get(K,3*i+j,3*i+j);
            double surfaceAreaTimesViscosityPerDt = surfaceAreaPerDt*Nodes[i]->externalViscosity[j];
            gsl_matrix_set(K,3*i+j,3*i+j,surfaceAreaTimesViscosityPerDt+curKValue);
        }
	}
}

void NewtonRaphsonSolver::checkJacobianForAblatedNodes(std::vector <int> & AblatedNodes){
    /** The function clears the system Jacobian for ablated nodes, which should not have
     * any forces acting on them.
     */
	int nAblatedNode = AblatedNodes.size();
	for (int a = 0; a<nAblatedNode; ++a){
		int NodeId = AblatedNodes[a]*3;
		for (int aa= 0; aa<3; ++aa){
			double Kdiagonal = gsl_matrix_get(K,NodeId+aa,NodeId+aa);
			if (Kdiagonal == 0){
				gsl_matrix_set(K,NodeId+aa,NodeId+aa,1);
			}
		}
	}
}

void NewtonRaphsonSolver::calculateSumOfInternalForces(){
    for (size_t i=0; i<nDim*nNodes; ++i){
		gsl_vector_set(gSum,i,gsl_matrix_get(ge,i,0)+gsl_matrix_get(gvInternal,i,0)+gsl_matrix_get(gvExternal,i,0));
	}
}

void NewtonRaphsonSolver::addExernalForces(){
    /** This function adds the external forces on the system sum of forces, NewtonRaphsonSolver#gSum)
     */
    for (size_t i=0; i<nDim*nNodes; ++i){
		gsl_vector_set(gSum,i,gsl_vector_get(gSum,i)+gsl_matrix_get(gExt,i,0));
		double value = gsl_vector_get(gSum,i);
        if (std::isnan(value)){
              std::cout<<" gSUM is nan at matrix point: "<<i<<std::endl;
		}
	}
}



void NewtonRaphsonSolver::solveForDeltaU(){
    /** This funciton will solve for new displacements from the system forces and the Jacobian.
     *This requires solving a sparse system of linear equatons, and the operation is handled by Pardiso solver.
     *Please refer to the manual of Pardiso to follow the necessary steps in this function/
     */
    const int nmult  = nDim*nNodes;
    int *ia = new int[nmult+1];
    double *b = new double[nmult];
    std::vector <int> ja_vec;
    std::vector <double> a_vec;
    constructiaForPardiso(ia, nmult, ja_vec, a_vec);
    const int nNonzero = ja_vec.size();
    int* ja = new int[nNonzero];
    double* a = new double [nNonzero];
    writeKinPardisoFormat(nNonzero, ja_vec, a_vec, ja, a);
    writeginPardisoFormat(b,nmult);
    int error = solveWithPardiso(a, b, ia, ja, nmult);
    if (error != 0){std::cerr<<"Pardiso solver did not return success!!"<<std::endl;}

    if (boundNodesWithSlaveMasterDefinition){
    	equateSlaveDisplacementsToMasters();
    }
    delete[] ia;
    delete[] ja;
    delete[] a;
    delete[] b;
}





// PARDISO prototype. //
extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *, double *, int    *,    int *, int *,   int *, int *,   int *, double *, double *, int *, double *);
extern "C" void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
extern "C" void pardiso_printstats (int *, int *, double *, int *, int *, int *, double *, int *);


int NewtonRaphsonSolver::solveWithPardiso(double* a, double*b, int* ia, int* ja, const int n_variables){
    int    n = n_variables;
    int    nnz = ia[n];
    int    mtype = 11;        /* Real unsymmetric matrix */

    /* RHS and solution vectors. */
    int      nrhs = 1;          /* Number of right hand sides. */
    double   x[n_variables];//, diag[n_variables];
    /* Internal solver memory pointer pt,                  */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
    /* or void *pt[64] should be OK on both architectures  */
    void    *pt[64];

    /* Pardiso control parameters. */
    int      iparm[64];
    double   dparm[64];
    int      maxfct, mnum, phase, error, msglvl, solver;

    iparm[60] = 1; //use in-core version when there is enough memory, use out of core version when not.

    /* Number of processors. */
    int      num_procs;

    /* Auxiliary variables. */
    char    *var;
    int      i;// k;

    double   ddum;              /* Double dummy */
    int      idum;              /* Integer dummy. */


/* -------------------------------------------------------------------- */
/* ..  Setup Pardiso control parameters.                                */
/* -------------------------------------------------------------------- */

    error = 0;
    solver = 0; /* use sparse direct solver */
    pardisoinit (pt,  &mtype, &solver, iparm, dparm, &error);

    if (error != 0)
    {
        if (error == -10 )
           printf("No license file found \n");
        if (error == -11 )
           printf("License is expired \n");
        if (error == -12 )
           printf("Wrong username or hostname \n");
         return 1;
    }
    else
        //printf("[PARDISO]: License check was successful ... \n");

    /* Numbers of processors, value of OMP_NUM_THREADS */
    var = getenv("OMP_NUM_THREADS");
    if(var != NULL)
        sscanf( var, "%d", &num_procs );
    else {
        printf("Set environment OMP_NUM_THREADS to 1");
        exit(1);
    }
    iparm[2]  = num_procs;

    maxfct = 1;		    /* Maximum number of numerical factorizations.  */
    mnum   = 1;         /* Which factorization to use. */

    iparm[10] = 0; /* no scaling  */
    iparm[12] = 0; /* no matching */

    msglvl = 0;         /* Print statistical information  */
    error  = 0;         /* Initialize error flag */

/* -------------------------------------------------------------------- */
/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
/*     notation.                                                        */
/* -------------------------------------------------------------------- */
    for (i = 0; i < n+1; i++) {
        ia[i] += 1;
    }
    for (i = 0; i < nnz; i++) {
        ja[i] += 1;
    }

/* -------------------------------------------------------------------- */
/*  .. pardiso_chk_matrix(...)                                          */
/*     Checks the consistency of the given matrix.                      */
/*     Use this functionality only for debugging purposes               */
/* -------------------------------------------------------------------- */
    bool carryOutDebuggingChecks = false;
    if (carryOutDebuggingChecks){
        pardiso_chkmatrix  (&mtype, &n, a, ia, ja, &error);
        if (error != 0) {
            printf("\nERROR in consistency of matrix: %d", error);
            exit(1);
        }
    }
/* -------------------------------------------------------------------- */
/* ..  pardiso_chkvec(...)                                              */
/*     Checks the given vectors for infinite and NaN values             */
/*     Input parameters (see PARDISO user manual for a description):    */
/*     Use this functionality only for debugging purposes               */
/* -------------------------------------------------------------------- */

    if (carryOutDebuggingChecks){
        pardiso_chkvec (&n, &nrhs, b, &error);
        if (error != 0) {
            printf("\nERROR  in right hand side: %d", error);
            exit(1);
        }
    }
/* -------------------------------------------------------------------- */
/* .. pardiso_printstats(...)                                           */
/*    prints information on the matrix to STDOUT.                       */
/*    Use this functionality only for debugging purposes                */
/* -------------------------------------------------------------------- */
    if (carryOutDebuggingChecks){
        pardiso_printstats (&mtype, &n, a, ia, ja, &nrhs, b, &error);
        if (error != 0) {
            printf("\nERROR right hand side: %d", error);
            exit(1);
        }
    }
/* -------------------------------------------------------------------- */
/* ..  Reordering and Symbolic Factorization.  This step also allocates */
/*     all memory that is necessary for the factorization.              */
/* -------------------------------------------------------------------- */
    phase = 11;
    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error, dparm);
//std::cout<<"symbolic factorisation"<<std::endl;
    if (error != 0) {
        printf("\nERROR during symbolic factorization: %d", error);
        exit(1);
    }
    //printf("\nReordering completed ... ");
    //printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
    //printf("\nNumber of factorization MFLOPS = %d", iparm[18]);

/* -------------------------------------------------------------------- */
/* ..  Numerical factorization.                                         */
/* -------------------------------------------------------------------- */
    phase = 22;
//    iparm[32] = 1; /* compute determinant */

    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error,  dparm);
//std::cout<<"numerical factorisation"<<std::endl;
    if (error != 0) {
        printf("\nERROR during numerical factorization: %d", error);
        exit(2);
    }
    //printf("\nFactorization completed ...\n ");

/* -------------------------------------------------------------------- */
/* ..  Back substitution and iterative refinement.                      */
/* -------------------------------------------------------------------- */
   /* phase = 33;

    iparm[7] = 1;       // Max numbers of iterative refinement steps.

    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, b, x, &error,  dparm);

    if (error != 0) {
        printf("\nERROR during solution: %d", error);
        exit(3);
    }
    bool displayResult = false;
    if (displayResult){
        printf("\nSolve completed ... ");
        printf("\nThe solution of the system is: ");
        for (i = 0; i < n; i++) {
            printf("\n x [%d] = % f", i, x[i] );
        }
        printf ("\n");
    }
    //Write x into deltaU:
    for (int i=0; i<n_variables; ++i){
        gsl_vector_set(deltaU,i,x[i]);
    }
    */
/* -------------------------------------------------------------------- */
/* ..  Back substitution with tranposed matrix A^t x=b                  */
/* -------------------------------------------------------------------- */

	phase = 33;
	//iparm[4]  = 61;	 /*changing the precision of convergence with pre-conditioning, not sure what it does, I added as trial, but did not change anything */
	iparm[7]  = 1;       /* Max numbers of iterative refinement steps. */
	iparm[11] = 1;       /* Solving with transpose matrix. */

	pardiso (pt, &maxfct, &mnum, &mtype, &phase,
			 &n, a, ia, ja, &idum, &nrhs,
			 iparm, &msglvl, b, x, &error,  dparm);

	if (error != 0) {
		printf("\nERROR during solution: %d", error);
		exit(3);
	}

	bool displayResult = false;
	if (displayResult){
		printf("\nSolve completed ... ");
		printf("\nThe solution of the system is: ");
		for (i = 0; i < n; i++) {
			printf("\n x [%d] = % f", i, x[i] );
		}
		printf ("\n");
	}
    //Write x into deltaU:
    for (int i=0; i<n_variables; ++i){
        gsl_vector_set(deltaU,i,x[i]);
    }

/* -------------------------------------------------------------------- */
/* ..  Convert matrix back to 0-based C-notation.                       */
/* -------------------------------------------------------------------- */
    for (i = 0; i < n+1; i++) {
        ia[i] -= 1;
    }
    for (i = 0; i < nnz; i++) {
        ja[i] -= 1;
    }

/* -------------------------------------------------------------------- */
/* ..  Termination and release of memory.                               */
/* -------------------------------------------------------------------- */
    phase = -1;                 /* Release internal memory. */

    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error,  dparm);
    return 0;
}


void NewtonRaphsonSolver::constructiaForPardiso(int* ia, const int nmult, std::vector<int> &ja_vec, std::vector<double> &a_vec){
    double negThreshold = -1E-13, posThreshold = 1E-13;
    //count how many elements there are on K matrix and fill up ia:
    int counter = 0;
    for (int i =0; i<nmult; ++i){
        bool wroteiaForThisRow = false;
        for (int j=0; j<nmult; ++j){
            double Kvalue = gsl_matrix_get(K,i,j);
            if (Kvalue>posThreshold || Kvalue<negThreshold){
                ja_vec.push_back(j);
                a_vec.push_back(Kvalue);
                if (!wroteiaForThisRow){
                    //std::cout<<"writing is for row "<<i<<" column is: "<<j<<std::endl;
                    ia[i] = counter;
                    wroteiaForThisRow = true;
                }
                counter++;
            }
        }
    }
    ia[nmult] = counter;
}




void NewtonRaphsonSolver::writeKinPardisoFormat(const int nNonzero, std::vector<int> &ja_vec, std::vector<double> &a_vec, int* ja, double* a){
    //now filling up the int & double arrays for ja, a
    for (int i=0 ; i<nNonzero; ++i){
        ja[i] = ja_vec[i];
        a[i]  = a_vec [i];
    }
}

void NewtonRaphsonSolver::writeginPardisoFormat(double* b, const int n){
    for (int i=0; i<n; ++i){
        b[i] = gsl_vector_get(gSum,i);
    }
}


bool NewtonRaphsonSolver::checkConvergenceViaDeltaU(){
    /**
     * This function checks the norm of NewtonRaphsonSolver#deltaU, the change in
     *nodal positions from previous iteration to this one. If the norm is below the threshold of convergence
     *NewtonRaphsonSolver#threshold, then the solution has been achieved for the current time step.
     */
    bool converged = true;

    double d = gsl_blas_dnrm2 (deltaU);
    //displayMatrix(deltaU,"deltaUInConvergence");
    if (d>threshold){
        converged = false;
        std::cout<<" not  yet converged via du: norm "<<d<<std::endl;
    }
    else{
        std::cout<<"converged with displacement: norm"<<d<<std::endl;
    }
    return converged;
}

bool NewtonRaphsonSolver::checkConvergenceViaForce(){
    /** The system can converge with zero forces as well. This function
     *checks the norm of the sum of all nodal forces against the threshold.
     */
    bool converged = true;
    double d = gsl_blas_dnrm2 (gSum);
    if (d>threshold){
        converged = false;
        std::cout<<" not  yet converged via forces: norm "<<d<<std::endl;
    }
    else{
        std::cout<<"converged with forces: norm"<<d<<std::endl;
    }
    return converged;
}

void NewtonRaphsonSolver::updateUkInIteration(){
    /**
     * This funciton updates the positions of nodes for next iteration $ \f \boldsymbol{u}_{k+1} $ \f
     *from the positions of current iteration the $ \f \boldsymbol{u}_{k} $ \f and $ \f \boldsymbol{\delta u}_{k} $ \f.
     */
    int n = uk->size1;
    for (int i=0; i<n;++i){
    	double newValue = gsl_matrix_get(uk,i,0)+gsl_vector_get(deltaU,i);
        gsl_matrix_set(uk,i,0,newValue);
    }
}

void NewtonRaphsonSolver::displayMatrix(gsl_matrix* mat, std::string matname){
    int m = mat->size1;
    int n = mat->size2;
    std::cout<<matname<<": "<<std::endl;

    for (int i =0; i<m; i++){
        for (int j =0; j<n; j++){
            std::cout.precision(4);
            std::cout.width(6);
            std::cout<<gsl_matrix_get(mat,i,j)<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

void NewtonRaphsonSolver::displayMatrix(gsl_vector* mat, std::string matname){
    int m = mat->size;
    //int n = mat->size2;
    std::cout<<matname<<": "<<std::endl;
    for (int i =0; i<m; i++){
        std::cout.precision(4);
        std::cout.width(6);
        std::cout<<gsl_vector_get(mat,i)<<std::endl;
    }
}

bool NewtonRaphsonSolver::checkIfCombinationExists(int dofSlave, int dofMaster){
    for(auto slaveMasterCouple : slaveMasterList){
        if(slaveMasterCouple[0] == dofSlave && slaveMasterCouple[1] == dofMaster){
			return false; //continue addition? false
		}
        if(slaveMasterCouple[1] == dofSlave && slaveMasterCouple[0] == dofMaster){
			return false; //continue addition? false
		}
	}
	return true;
}

void NewtonRaphsonSolver::checkMasterUpdate(int& dofMaster, int& masterId){
    /** The function takes in the addresses for potential master and slave degrees of freedom indices.
     * If the proposed master DoF is already on the NewtonRaphsonSolver#slaveMasterList as a slave, then the
     * its master DoF should become the master  of the proposed slave DoF. The order of NewtonRaphsonSolver#slaveMasterList is
     *[slave][master].
     */
    for(auto slaveMasterCouple : slaveMasterList){
        if(slaveMasterCouple[0] == dofMaster){
            dofMaster = slaveMasterCouple[1];
			int dim = dofMaster % 3;
			masterId = (dofMaster-dim)/3;
			break;
		}
	}
}


bool NewtonRaphsonSolver::checkIfSlaveIsAlreadyMasterOfOthers(int dofSlave, int dofMaster){
    /** The function takes in the addresses for potential master and slave degrees of freedom indices.
     *If the proposed slave DoF is already on the NewtonRaphsonSolver#slaveMasterList as a master, then transfer
     *all DoF that are its slave to the new proposed master.
     */
    size_t n= slaveMasterList.size();
	bool madeChange = false;
    for(size_t i=0; i<n;++i){
		if(slaveMasterList[i][1] == dofSlave){
            std::cout<<"making change, slave was master of "<<slaveMasterList[i][0]<<std::endl;
			//proposed slave is already a master, update the master to the proposed master
			slaveMasterList[i][1] = dofMaster;
			madeChange = true;
		}
	}
	return madeChange;
}
