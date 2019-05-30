#include "Prism.h"
#include "ReferenceShapeBase.h"
#include <stdio.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

Prism::Prism(const std::vector<int>& inpNodeIds, const std::vector<std::unique_ptr<Node>>& Nodes, int CurrId){
	nNodes = 6;
	nDim = 3;
	Id = CurrId;
	ShapeDim = 3;	//3D shape
    E = 250.0;
	v = 0.3;
    internalViscosity = 0;
    lambda = E*v /(1+v)/(1-2.0*v);
    mu = E/2.0/(1+v);
    stiffnessMultiplier = 1.0;
    minimumValueOfStiffnessMultiplier = 0.00001;
    maximumValueOfStiffnessMultiplier = 100000;

    D = gsl_matrix_calloc(6,6);
    ShapeChangeRate.fill(0.0);
    GrowthRate.fill(0.0);
    relativePosInBoundingBox.fill(0.0);
    initialRelativePosInBoundingBox.fill(0.0);
	initialRelativePositionInZ = 1.0;
    apicalNormalCurrentShape.fill(0);
	columnarGrowthWeight = 1.0;
	peripodialGrowthWeight = 0.0;
	isFlipped = false;
	IsChangingShape = false;
	IsAblated = false;
	atSymetricityBoundary = false;
	capElement = false;
    rotatedGrowth = false;
	setIdentificationColour();
	setShapeType("Prism");
    ReferenceShape = std::make_unique<ReferenceShapeBase>("Prism",Id);

    readNodeIds(inpNodeIds);
	setPositionMatrix(Nodes);
	setReferencePositionMatrix();
    setInitialEdgeLenghts();
	setCoeffMat();
	calculateReferenceVolume();
	setTissuePlacement(Nodes);
	setTissueType(Nodes);

	numberOfGaussPoints =6;
    double oneOverSqrtThree = 1.0/ pow(3,0.5);
    gaussPoints[0][0]= 1.0/6.0;
    gaussPoints[0][1]= 1.0/6.0;
    gaussPoints[0][2]= oneOverSqrtThree;
    gaussPoints[1][0]= 2.0/3.0;
    gaussPoints[1][1]= 1.0/6.0;
    gaussPoints[1][2]= oneOverSqrtThree;
    gaussPoints[2][0]= 1.0/6.0;
    gaussPoints[2][1]= 2.0/3.0;
    gaussPoints[2][2]= oneOverSqrtThree;
    gaussPoints[3][0]= 1.0/6.0;
    gaussPoints[3][1]= 1.0/6.0;
    gaussPoints[3][2]= -oneOverSqrtThree;
    gaussPoints[4][0]= 2.0/3.0;
    gaussPoints[4][1]= 1.0/6.0;
    gaussPoints[4][2]= -oneOverSqrtThree;
    gaussPoints[5][0]= 1.0/6.0;
    gaussPoints[5][1]= 2.0/3.0;
    gaussPoints[5][2]= -oneOverSqrtThree;
    gaussWeights[0] = 1.0/6.0;
    gaussWeights[1] = 1.0/6.0;
    gaussWeights[2] = 1.0/6.0;
    gaussWeights[3] = 1.0/6.0;
    gaussWeights[4] = 1.0/6.0;
    gaussWeights[5] = 1.0/6.0;
    detdXdes.fill(0.0);
    detFs.fill(0.0);
    //This is a vector of 4th order tensors
    //D81 -> [nGaussPoints][3][3][3][3]
    //It is confusing enough as is, no "clever" bit semantics please.
    for (size_t i=0; i<numberOfGaussPoints; ++i){
        gsl_matrix* tmpShapeFuncDerivative = gsl_matrix_calloc(nDim, nNodes);
        ShapeFuncDerivatives.push_back(tmpShapeFuncDerivative);
        gsl_matrix* tmpShapeFuncDerStack = gsl_matrix_calloc(nDim*nDim, nDim*nNodes);
        ShapeFuncDerStacks.push_back(tmpShapeFuncDerStack);
        gsl_matrix* tmpInvdXde = gsl_matrix_calloc(nDim, nDim);
        InvdXdes.push_back(tmpInvdXde);
        gsl_matrix* tmpB = gsl_matrix_calloc(nNodes,nDim*nNodes);
        Bmatrices.push_back(tmpB);
        gsl_matrix* tmpFeMatrice = gsl_matrix_calloc(3,3);
        FeMatrices.push_back(tmpFeMatrice);
        gsl_matrix* tmpinvJShapeFuncDerStack= gsl_matrix_calloc(nDim*nDim, nDim*nNodes);
        invJShapeFuncDerStack.push_back(tmpinvJShapeFuncDerStack);
        gsl_matrix* tmpinvJShapeFuncDerStackwithFe = gsl_matrix_calloc(nDim*nDim, nDim*nNodes);
        invJShapeFuncDerStackwithFe.push_back(tmpinvJShapeFuncDerStackwithFe);
        gsl_matrix* tmpelasticStress = gsl_matrix_calloc(3,3);
        elasticStress.push_back(tmpelasticStress);
        gsl_matrix* tmpviscousStress = gsl_matrix_calloc(3,3);
        viscousStress.push_back(tmpviscousStress);
        //need to reach:  tensor[3][3][3][3], and push it back to vector D81 (81 = 3*3*3*3)
        //for each gauss point:
        std::array<std::array<std::array<std::array<double,3>,3>,3>,3> D81ForCurrentGaussPoint = {};
        for (size_t j=0; j<nDim; ++j){
             for (size_t k=0; k<nDim; ++k){
                 for (size_t l=0; l<nDim; ++l){
                     D81ForCurrentGaussPoint[j][k][l].fill(0.0);
        		 }
			}
        }
        D81.push_back(D81ForCurrentGaussPoint);
    }
    Strain = gsl_matrix_calloc(6,1);
    GrowthStrainsRotMat = gsl_matrix_alloc(3,3);
    gsl_matrix_set_identity(GrowthStrainsRotMat);
    Fg = gsl_matrix_alloc(3,3);
    gsl_matrix_set_identity(Fg);
    InvFg = gsl_matrix_calloc(3,3);
    gsl_matrix_set_identity(InvFg);
    Fsc = gsl_matrix_alloc(3,3);
    gsl_matrix_set_identity(Fsc);
    InvFsc = gsl_matrix_calloc(3,3);
    gsl_matrix_set_identity(InvFsc);

    growthIncrement = gsl_matrix_calloc(3,3);
    gsl_matrix_set_identity(growthIncrement);
    plasticDeformationIncrement= gsl_matrix_calloc(3,3);
    gsl_matrix_set_identity(plasticDeformationIncrement);
    shapeChangeIncrement = gsl_matrix_calloc(3,3);
    gsl_matrix_set_identity(shapeChangeIncrement);
    zRemodellingSoFar = 1.0;
    TriPointF = gsl_matrix_calloc(3,3);
	gsl_matrix_set_identity(TriPointF);
    TriPointKe = gsl_matrix_calloc(nDim*nNodes,nDim*nNodes);
    TriPointKv = gsl_matrix_calloc(nDim*nNodes,nDim*nNodes);
    ElementalElasticSystemForces = gsl_matrix_calloc(nNodes,nDim);
    ElementalInternalViscousSystemForces = gsl_matrix_calloc(nNodes,nDim);
	RotatedElement = false;    

	CurrShapeChangeToAdd[0] = 0;
	CurrShapeChangeToAdd[1] = 0;
	CurrShapeChangeToAdd[2] = 0;

	VolumePerNode = 0;
    ZProjectedBasalArea=0.0;
    ZProjectedApicalArea=0.0;
    BasalArea=0.0;
    ApicalArea=0.0;
    exposedLateralAreaApicalSide = 0;
    exposedLateralAreaApicalSide = 0;
    elementHasExposedApicalSurface = false;
    elementHasExposedBasalSurface = false;
    exposedApicalSurfaceNodeIds[0] = 0;
    exposedApicalSurfaceNodeIds[1] = 0;
    exposedApicalSurfaceNodeIds[2] = 0;
    exposedBasalSurfaceNodeIds[0] = 0;
    exposedBasalSurfaceNodeIds[1] = 0;
    exposedBasalSurfaceNodeIds[2] = 0;
    exposedLateralAreaApicalSideNodeIds[0] = 0;
    exposedLateralAreaApicalSideNodeIds[1] = 0;
    exposedLateralAreaApicalSideNodeIds[2] = 0;
    exposedLateralAreaApicalSideNodeIds[3] = 0;
    exposedLateralAreaBasalSideNodeIds[0] = 0;
    exposedLateralAreaBasalSideNodeIds[1] = 0;
    exposedLateralAreaBasalSideNodeIds[2] = 0;
    exposedLateralAreaBasalSideNodeIds[3] = 0;
    nLateralSurfaceAreaNodeNumber = 4;
    nSurfaceAreaNodeNumber = 3;

    remodellingPlaneRotationMatrix = gsl_matrix_calloc(3,3);
    gsl_matrix_set_identity(remodellingPlaneRotationMatrix);
    isECMMimicing = false;
    isECMMimimcingAtCircumference = false;
    atBasalBorderOfECM = false;
    isActinMimicing = false;
    atApicalBorderOfActin = false;
    insideEllipseBand = false;
    coveringEllipseBandId = -1;
    stiffnessPerturbationRateInSec = 0;
    basalNeigElementId = -1;

    isMutated = false;
    mutationGrowthRatePerSec = 0.0;

    thereIsGrowthRedistribution = false;
    growthRedistributionShrinksElement = false;
    growthRedistributionScale = 0.5;

    plasticDeformationHalfLifeMultiplier = 1.0;
}

Prism::~Prism(){
    //freeing matrices allocated
	gsl_matrix_free(growthIncrement);
	gsl_matrix_free(plasticDeformationIncrement);
    gsl_matrix_free(D);
    gsl_matrix_free(CoeffMat);
    gsl_matrix_free(Fg);
    gsl_matrix_free(InvFg);
    gsl_matrix_free(Fsc);
    gsl_matrix_free(InvFsc);
    gsl_matrix_free(TriPointF);
    gsl_matrix_free(Strain);
    gsl_matrix_free(TriPointKe);
    gsl_matrix_free(TriPointKv);
    gsl_matrix_free(GrowthStrainsRotMat);
    for (size_t i=0; i<numberOfGaussPoints; ++i){
        gsl_matrix_free (ShapeFuncDerivatives[i]);
        gsl_matrix_free (ShapeFuncDerStacks[i]);
        gsl_matrix_free (InvdXdes[i]);
        gsl_matrix_free (Bmatrices[i]);
        gsl_matrix_free (FeMatrices[i]);
        gsl_matrix_free (invJShapeFuncDerStack[i]);
        gsl_matrix_free (invJShapeFuncDerStackwithFe[i]);
        gsl_matrix_free (elasticStress[i]);
        gsl_matrix_free (viscousStress[i]);
    }
}

void Prism::setCoeffMat(){
    /** The coefficient matrix collating the calculated stresses in the form of the vector
     * form of elemental strss and strain vectors:
     *
     * \f{eqnarray*}{
        \boldsymbol{c} =
        \begin{bmatrix}
        1	& 0	& 0	& 0	& 0	& 0	& 0	& 0	& 0	\\
        0	& 0	& 0	& 0	& 1	& 0	& 0	& 0	& 0	\\
        0	& 0	& 0	& 0	& 0	& 0	& 0	& 0	& 1	\\
        0	& 1	& 0	& 1	& 0	& 0	& 0	& 0	& 0	\\
        0	& 0	& 0	& 0	& 0	& 1	& 0	& 1	& 0	\\
        0	& 0	& 1	& 0	& 0	& 0	& 1	& 0	& 0	\\
        \end{bmatrix}
        \sim
        \begin{bmatrix}
        \epsilon_{x} \\
        \epsilon_{y} \\
        \epsilon_{z}\\
        \gamma_{x,y} \\
        \gamma_{x,z} \\
        \gamma_{y,z}
        \end{bmatrix}
        \sim
        \begin{bmatrix}
        \sigma_{x} \\
        \sigma_{y} \\
        \sigma_{z}\\
        \sigma_{x,y} \\
        \sigma_{x,z} \\
        \sigma_{y,z}
        \end{bmatrix}.
        \f}
      */
    CoeffMat = gsl_matrix_calloc(6, nDim*nDim);
    gsl_matrix_set(CoeffMat,0,0,1);
    gsl_matrix_set(CoeffMat,1,4,1);
    gsl_matrix_set(CoeffMat,2,8,1);
    gsl_matrix_set(CoeffMat,3,1,1);
    gsl_matrix_set(CoeffMat,3,3,1);
    gsl_matrix_set(CoeffMat,4,5,1);
    gsl_matrix_set(CoeffMat,4,7,1);
    gsl_matrix_set(CoeffMat,5,2,1);
    gsl_matrix_set(CoeffMat,5,6,1);
}

void Prism::checkRotationConsistency3D(){
	//The nodes should be ordered in counter-clock-wise order
	//The view-vector is from apical towards basal
	//(If the system is not rotated, and is in standard coordinate system,
	//I am looking from top, view is (-)ve z;
	//If they are not, correct the order here

    std::array<double,3> vec1   = {0.0};
    std::array<double,3> vec2   = {0.0};
    std::array<double,3> view   = {0.0};
    std::array<double,3> normal = {0.0};
    for (size_t i= 0; i<nDim; ++i){
		vec1[i] = Positions[1][i] - Positions[0][i];
		vec2[i] = Positions[2][i] - Positions[0][i];
		view[i] = Positions[0][i] - Positions[3][i];
	}
    normal = crossProduct3D(vec1,vec2);
	double  dot = dotProduct3D(view,normal);
	if (dot > 0) {
        std::cerr<<"prism: "<<Id<<" nodes are ordered clockwise, correcting"<<std::endl;
        std::cout<<"Positions before swap, element: "<<Id<<std::endl;
		displayPositions();
		//swapping node ids:
		int ids[2] = { NodeIds[1], NodeIds[4]};
		NodeIds[1] = NodeIds[2];
		NodeIds[4] = NodeIds[5];
		NodeIds[2] = ids[0];
		NodeIds[5] = ids[1];
		//swapping positions:
        for (size_t i = 0; i<nDim; ++i){
			double pos[2] = {Positions[1][i],Positions[4][i]};
			double refpos[2] = {ReferenceShape->Positions[1][i],ReferenceShape->Positions[4][i]};
			Positions[1][i] = Positions[2][i];
			Positions[4][i] = Positions[5][i];
			Positions[2][i] = pos[0];
			Positions[5][i] = pos[1];
			ReferenceShape->Positions[1][i] = ReferenceShape->Positions[2][i];
			ReferenceShape->Positions[4][i] = ReferenceShape->Positions[5][i];
			ReferenceShape->Positions[2][i] = refpos[0];
			ReferenceShape->Positions[5][i] = refpos[1];
		}
        std::cout<<"Positions after swap, element: "<<Id<<std::endl;
		displayPositions();
	}
}

void  Prism::calculateApicalNormalCurrentShape(){
    apicalNormalCurrentShape.fill(0);
    std::array<double,3> u;
    std::array<double,3> v;
    for (size_t i=0; i<nDim; ++i){
		u[i] = Positions[4][i] - Positions[3][i];
		v[i] = Positions[5][i] - Positions[3][i];
	}
    apicalNormalCurrentShape = crossProduct3D(u,v);
    (void) normaliseVector3D(apicalNormalCurrentShape);
    for (size_t i=0; i<nDim; ++i){
		u[i] = Positions[0][i] - Positions[3][i];
	}
	double  dot = dotProduct3D(u,apicalNormalCurrentShape);
	if (dot<0){
        for (auto& currCoord : apicalNormalCurrentShape){
            currCoord *= -1.0;
		}
	}
}

std::array<double,3>  Prism::calculateBasalNormal(){
    std::array<double,3> u = {0.0, 0.0,0.0};
    std::array<double,3> v = {0.0, 0.0,0.0};
    for (size_t i=0; i<nDim; ++i){
		u[i] = ReferenceShape->Positions[1][i] - ReferenceShape->Positions[0][i];
		v[i] = ReferenceShape->Positions[2][i] - ReferenceShape->Positions[0][i];
	}
    std::array<double,3> normal = crossProduct3D(u,v);
    (void) normaliseVector3D(normal);
    for (size_t i=0; i<nDim; ++i){
		u[i] = ReferenceShape->Positions[3][i] - ReferenceShape->Positions[0][i];
	}
	double  dot = dotProduct3D(u,normal);
	if (dot<0){
        for (size_t i=0; i<nDim; ++i){
			normal[i] *=(-1.0);
		}
	}
    std::cerr<<"		normal after direction correction: "<<normal[0]<<" "<<normal[1]<<" "<<normal[2]<<std::endl;
    return normal;
}

void  Prism::setElasticProperties(double EApical, double EBasal, double EMid, double EECM, double v){
	this -> E = EMid;
	if (isECMMimicing){
		this -> E = EECM;
	}
	else{
		if (tissuePlacement == 0 || atBasalBorderOfECM){
			this -> E = EBasal;
		}
		else if(tissuePlacement == 1 ){
			this -> E = EApical;
		}
	}
	this -> v = v; //poisson ratio
	if (this -> v>0.5){this -> v= 0.5;}
	else if (this -> v<0.0){this -> v = 0.0;}

    lambda = E*v/(1+v)/(1-2.0*v);
    mu = E/2.0/(1+v);

    calculateDVector();
    calculateD81Tensor();
    //std::cout<<" Element: "<<Id<<" E : "<<E<<" v: "<<v<<" lambda: "<<lambda<< " mu: "<<mu<<" atBasalBorderOfECM "<<atBasalBorderOfECM<<std::endl;
}


void  Prism::updateElasticProperties(){
    lambda = stiffnessMultiplier*E*v/(1+v)/(1-2.0*v);
    mu = stiffnessMultiplier*E/2.0/(1+v);
    //These two updates are not really necessary for a neo-Hookean material, as the calculation for
    //D81 vector is carried out at each Newton-Raphson iteration. I do not know the material type now,
    //therefore I dont know if I should skip this or not. Can be eliminated for efficiency later on!
    calculateDVector();
    calculateD81Tensor();
}

void Prism::calculateDVector(){
	double multiplier = E/((1+v)*(1-2*v));
	gsl_matrix_set(D,0,0,  multiplier*(1-v));
	gsl_matrix_set(D,0,1,  multiplier*v);
	gsl_matrix_set(D,0,2,  multiplier*v);
	gsl_matrix_set(D,1,0,  multiplier*v);
	gsl_matrix_set(D,1,1,  multiplier*(1-v));
	gsl_matrix_set(D,1,2,  multiplier*v);
	gsl_matrix_set(D,2,0,  multiplier*v);
	gsl_matrix_set(D,2,1,  multiplier*v);
	gsl_matrix_set(D,2,2,  multiplier*(1-v));
	gsl_matrix_set(D,3,3,  multiplier*(1-2*v)/2);
	gsl_matrix_set(D,4,4,  multiplier*(1-2*v)/2);
	gsl_matrix_set(D,5,5,  multiplier*(1-2*v)/2);
}

void Prism::calculateD81Tensor(){
    /**
    * This function calculates the Lagrangian elasticity tensor for a neo-hookean material.
    * The Lagrangian elasticity tensor depends on the deformation gradient through the
    * Cauchy-Green Deformation tensor, and the physical material properties, namely Lame's second
    * parameter \f$ \lambda \f$ and the sheer modulus \f$ \mu \f$. \n
    * \f$\boldsymbol {\mathcal D} \f$ can be obtained in the Voigt notation:
    * \f{eqnarray*}{
            \mathcal{D}_{ijkl}  = \lambda \boldsymbol{C}^{-1}_{ij} \boldsymbol{C}^{-1}_{kl} + 2 \left( \mu - \lambda ln(J)\right) \mathcal{I}_{ijkl}
      \f}
    * where \f$ \mathcal{I}_{ijkl} \f$ is defined as:
      \f{eqnarray*}{
            \mathcal{I}_{ijkl}  = \frac{1}{2} \left[ \boldsymbol{C}^{-1}_{ik} \boldsymbol{C}^{-1}_{jl} + \boldsymbol{C}^{-1}_{il} \boldsymbol{C}^{-1}_{jk}\right]
      \f}
    * where the indices $ i,j,k,l $ go through 3 dimensions.\n
    *
    */
	// lambda is Lame s first parameter and mu is the shear modulus .
	double Idouble[3][3] = {{1.0,0.0,0.0} , {0.0,1.0,0.0}, {0.0,0.0,1.0}};
    for (size_t pointNo = 0; pointNo<numberOfGaussPoints; pointNo++){
        for (size_t I = 0; I<nDim; ++I){
            for (size_t J = 0; J<nDim; ++J){
                for (size_t K = 0; K<nDim; ++K){
                    for (size_t L = 0; L<nDim; ++L){
						D81[pointNo][I][J][K][L] = lambda*Idouble[K][L]*Idouble[I][J] + mu * ( Idouble[I][K]*Idouble[J][L] + Idouble[I][L]*Idouble[J][K] );
					}
				}
			}
		}
	}
}

void Prism::getCurrRelaxedShape(gsl_matrix* CurrRelaxedShape){
    /** This function writes the current relaxed shape of the elemetn in the gsl_matrix
     * pointed by the provided input gsl_matrix pointer. The matrix must be allocated before
     * the call to the function.
     *
     */
    for (size_t i =0; i<nNodes; ++i){
        for (size_t j=0; j<nDim; ++j){
            gsl_matrix_set(CurrRelaxedShape,i,j,ReferenceShape->Positions[i][j]);
		}
	}
}

void Prism::setShapeFunctionDerivatives(gsl_matrix* ShapeFuncDer, double eta, double zeta, double nu){
    /** This function calculates and writes the shape function derivatives on
     * the gsl_matrix pointed by the provided input gsl_matrix pointer,
     * for the provided barycentric coordinates. The matrix must be allocated prior to
     * calling the function. \n
     * For the prism, the shape functions and their derivatives are given as:
     *
     * \f{tabular}{{|c|c|c|c|c|}
            \multicolumn{5}{|c|}{ $\lambda = 1 - \xi - \eta$ }            \\
            \multicolumn{5}{|c|}{ $\alpha = \left( 1-\zeta \right) /2$  } \\
            \multicolumn{5}{|c|}{ $b = \left(1+\zeta\right)/2$ }          \\
            Node & Shape function.    & Shape function         	 & Shape function 	  	   & Shape function.         \\
                 & 			        & derivative wrt $\xi$      	 &  derivative wrt $\zeta$      & derivative wrt $\eta$ \\
                 & 		$N$ 	        & $\frac{\partial{N}}{ \partial \xi}$     & $\frac{\partial{N}}{\partial \zeta}$  &  $\frac{\partial{N}}{\partial \eta}$ \\
            1       & $\lambda \alpha$  & $- \alpha$         		 & $- \alpha$	  	   	    & $ \frac{-\lambda} {2}$        \\

            2	 & $\xi \alpha$           & $\alpha$     		 	& 0  	  	   			    & $\frac{-\xi} {2} $        \\

            3  	& $\eta \alpha$  	& 0         		 		& $\alpha$ 	  	   	    & $\frac{-\nu} {2} $        \\

            4    	& $\lambda b$  		&  $-b$      			& $-b$ 	  	   		    & $\frac{\lambda} {2} $         \\

            5	& $\xi b$  			& $b$         		 	&0 	  	   			    & $\frac{\xi} {2} $         \\

            6 	 & $\eta b$  		& 0    		 		& $b$   	  	   		    & $\frac{\eta} {2} $        \\
         \f}

     *
     * The shape function derivatives matrix has the form:
     * 	\f{eqnarray*}{
        \begin{bmatrix}
            \frac{\partial N_0}{ \partial \xi}      & ... 	& \frac{\partial N_n}{ \partial \xi}    \\
            \frac{\partial N_0}{ \partial \zeta}    & ...   & \frac{\partial N_n}{ \partial \zeta}  \\
            \frac{\partial N_0}{ \partial \eta}     & ...   & \frac{\partial N_n}{ \partial \eta}
        \end{bmatrix}
        \f}

     */

    double alpha  = (1 - zeta)/2;
	double beta = (1 + zeta)/2;
	double lambda = 1-eta-nu;

    gsl_matrix_set(ShapeFuncDer,0,0, -alpha);
    gsl_matrix_set(ShapeFuncDer,0,1,  alpha);
    gsl_matrix_set(ShapeFuncDer,0,2,  0);
    gsl_matrix_set(ShapeFuncDer,0,3, -beta);
    gsl_matrix_set(ShapeFuncDer,0,4,  beta);
    gsl_matrix_set(ShapeFuncDer,0,5, 0);

    gsl_matrix_set(ShapeFuncDer,1,0, -alpha);
    gsl_matrix_set(ShapeFuncDer,1,1,  0);
    gsl_matrix_set(ShapeFuncDer,1,2,  alpha);
    gsl_matrix_set(ShapeFuncDer,1,3, -beta);
    gsl_matrix_set(ShapeFuncDer,1,4,  0);
    gsl_matrix_set(ShapeFuncDer,1,5,  beta);

    gsl_matrix_set(ShapeFuncDer,2,0, -lambda/2);
    gsl_matrix_set(ShapeFuncDer,2,1, -eta/2);
    gsl_matrix_set(ShapeFuncDer,2,2, -nu/2);
    gsl_matrix_set(ShapeFuncDer,2,3,  lambda/2);
    gsl_matrix_set(ShapeFuncDer,2,4,  eta/2);
    gsl_matrix_set(ShapeFuncDer,2,5,  nu/2);
}

void Prism::setShapeFunctionDerivativeStack(gsl_matrix* ShapeFuncDer, gsl_matrix* ShapeFuncDerStack){
    /** This function writes the input gsl_matrix ShapeFuncDer into a stack form, in the
     * matrix ShapeFuncDerStack. Both matrices must be allocated prior to calling the function.
     * This stack is for direct multiplication with B matrix during calculation of stresses.
     * with the notation  \f$ N_{i, \xi} \f$ being the short hand for
     * \f$ \partial \boldsymbol{x_{i}} / \partial {\xi} \f$
     * for node i, the stack for of the shape function derivatives in open form is:
     *
     *
     *\f{eqnarray*}{
        \begin{bmatrix}
            N_{0,\xi} 		& 0 			& 0			& N_{1,\xi} 	& \cdots &N_{n-1,\xi} 		& 0 			& 0			\\
            N_{0,\zeta} 	& 0 			& 0 		& N_{1,\zeta}  	& \cdots &N_{n-1,\zeta}       & 0 			& 0			\\
            N_{0,\nu} 		& 0			& 0 			& N_{1,\nu}  	& \cdots &N_{n-1,\nu} 		& 0 			& 0			\\
            0 			& N_{0,\xi} 	& 0 			& 0 			& \cdots &0                 & N_{n-1,\xi} 	& 0			\\
            0 			& N_{0,\zeta} 	& 0 			& 0 	 		& \cdots &0                 & N_{n-1,\zeta} 	& 0			\\
            0 			& N_{0,\nu} 	& 0 			& 0  			& \cdots &0                 & N_{n-1,\nu} 	& 0			\\
            0 			& 0			& N_{0,\xi}         & 0  			& \cdots &0                 &0			& N_{n-1,\xi} 	\\
            0 			& 0			& N_{0,\zeta}       & 0 	 		& \cdots &0                 &0			& N_{n-1,\zeta}	\\
            0 			& 0			& N_{0,\nu}         & 0   			& \cdots &0                 &0			& N_{n-1,\nu} 	.
        \end{bmatrix}
        \f}
     *
     */
    int n = nNodes;
	int dim = nDim;
	for (int i=0; i<n;++i){
        for (int k=0; k<dim; ++k){
            gsl_matrix_set(ShapeFuncDerStack,k,i*dim, gsl_matrix_get(ShapeFuncDer,k,i));
        }
	}
    for (int k =dim; k<2*dim; ++k){
        for (int m =1; m<dim*n; ++m){
            gsl_matrix_set(ShapeFuncDerStack,k,m,gsl_matrix_get(ShapeFuncDerStack,k-dim,m-1));
        }
    }
    for (int k =2*dim; k<3*dim; ++k){
        for (int m =1; m<dim*n; ++m){
            gsl_matrix_set(ShapeFuncDerStack,k,m,gsl_matrix_get(ShapeFuncDerStack,k-dim,m-1));
        }
    }
}

void Prism::calculateElementShapeFunctionDerivatives(){
    /** This function calculates the Shape Function Derivatives for all Gauss points,
     * and then calculates the derivative of reference coordiantes with respect to barycentric coordinates
     * \f$ \frac{\partial \boldsymbol{X}}{\partial \boldsymbol{\zeta}} \f$. The determnant of htis derivative is recorded for next steps of calculation.\n
     * First the matrices for the reference shape and the \f$ \frac{\partial \boldsymbol{X}}{\partial \boldsymbol{\zeta}} \f$
     * for one gauss point are allocated.
     */
    gsl_matrix* CurrRelaxedShape = gsl_matrix_calloc(nNodes, nDim);
    gsl_matrix* dXde  = gsl_matrix_calloc(nDim, nDim);
    /** The reference shape cordinates are obtained through function ShapeBase#getCurrRelaxedShape.
     */
    getCurrRelaxedShape(CurrRelaxedShape);
    for (size_t iter =0; iter<numberOfGaussPoints;++iter){
        /**
         * Looping over all the gauss points as listed in constructor Prism#Prism,
         * the shape function derivatives are obtained in stack form via ShapeBase#setShapeFunctionDerivatives and
         * ShapeBase#setShapeFunctionDerivativeStack functions. Then \f$ \frac{\partial \boldsymbol{X}}{\partial \boldsymbol{\zeta}} \f$,
         *its determinant recorded in the array  ShapeBase#detdXdes.
         */
        double eta  = gaussPoints[iter][0];
        double nu   = gaussPoints[iter][1];
        double zeta = gaussPoints[iter][2];
        //calculating shape function derivatives:
        setShapeFunctionDerivatives(ShapeFuncDerivatives[iter],eta,zeta,nu);
        setShapeFunctionDerivativeStack(ShapeFuncDerivatives[iter],ShapeFuncDerStacks[iter]);
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, ShapeFuncDerivatives[iter], CurrRelaxedShape, 0.0, dXde);
        gsl_matrix_transpose(dXde);
        detdXdes[iter] = determinant3by3Matrix(dXde);
        bool inverted = InvertMatrix(dXde, InvdXdes[iter]);
        if (!inverted){
            std::cerr<<"dXde not inverted at point "<<iter<<"!!"<<std::endl;
        }
    }
    gsl_matrix_free(CurrRelaxedShape);
    gsl_matrix_free(dXde);
}

void Prism::calculateCurrTriPointFForRotation(gsl_matrix *currF,int pointNo){
    /**
     * This function calcultes the defrmation gradient \f$ \boldsymbol{F} \f$ for the input
     * gauss point number pointNo. The calculated deformation gradient is recorded in the in input
     * gsl_matrix pointer currF, the matrix must be allocated prior to calling the function. \n
     * The deformation gradient is the derivative of current nodal positions of an element ( \f$ \boldsymbol{x} \f$ )
     * with respect to the reference shape positions ( \f$ \boldsymbol{X} \f$ ).
     *
     * \f{eqnarray*}{
        \frac{\partial \boldsymbol{x}} {\partial \boldsymbol{X}}
        \f}
     *
     */
	const int n = nNodes;
	const int dim = nDim;
    gsl_matrix* CurrShape = gsl_matrix_alloc(n,dim);
    getPos(CurrShape);
    gsl_matrix* ShapeFuncDer = ShapeFuncDerivatives[pointNo];
	gsl_matrix* InvdXde = InvdXdes[pointNo];
	gsl_matrix* Jacobian = gsl_matrix_calloc(dim, dim);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, ShapeFuncDer, CurrShape, 0.0, Jacobian);
	gsl_matrix_transpose(Jacobian);
	//calculating F:
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Jacobian, InvdXde, 0.0, currF);
    gsl_matrix_free(CurrShape);
    gsl_matrix_free(Jacobian);
}




void Prism::calculateCurrNodalForces(gsl_matrix *currge, gsl_matrix *currgv, gsl_matrix *currF, gsl_matrix* displacementPerDt, int pointNo){
    const int n = nNodes;
    const int dim = nDim;
    gsl_matrix* currFe = gsl_matrix_alloc(dim,dim);
    gsl_matrix* CurrShape = gsl_matrix_alloc(n,dim);
    //Getting the current shape positions matrix:
    getPos(CurrShape);

    //calculating dx/de (Jacobian) and reading out dX/de, shape function derivaties:
    gsl_matrix* ShapeFuncDer = ShapeFuncDerivatives[pointNo];
    gsl_matrix* ShapeFuncDerStack = ShapeFuncDerStacks[pointNo];
    gsl_matrix* InvdXde = InvdXdes[pointNo];
    gsl_matrix* B = Bmatrices[pointNo]; //I will update and use this value, it is not a constant.
    gsl_matrix* invJShFuncDerS = invJShapeFuncDerStack[pointNo];
    gsl_matrix* invJShFuncDerSWithFe =invJShapeFuncDerStackwithFe[pointNo];
    gsl_matrix* Jacobian = gsl_matrix_calloc(dim, dim);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, ShapeFuncDer, CurrShape, 0.0, Jacobian);
    gsl_matrix_transpose(Jacobian);
    //calculating F:
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Jacobian, InvdXde, 0.0, currF);

    //calculating Fe:
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, currF, InvFg, 0.0, currFe);	///< Removing growth

    createMatrixCopy(FeMatrices[pointNo], currFe); // storing Fe for use in implicit elastic K calculation.

    //setting material type:
    bool KirshoffMaterial = false;
    bool neoHookeanMaterial = !KirshoffMaterial;

    //calculating Cauchy-Green Deformation Tensor
    gsl_matrix* E;
    gsl_matrix* S;
    gsl_matrix* C = calculateCauchyGreenDeformationTensor(currFe);

    //calculating S:
    double detFe = determinant3by3Matrix(currFe);

    double lnJ = log(detFe);
    if(std::isnan(lnJ)){
        std::cout<<"element: "<<Id<<" lnJ is nan, detFe: "<<detFe<<std::endl;
        std::cout<<" Element positions: "<<std::endl;
    	displayPositions();
    }
    if (KirshoffMaterial){
    	//calculating E (E = 1/2 *(Fe^T*Fe-I):
    	E = calculateEForNodalForcesKirshoff(C);
    	//calculating S: (S = D:E)
    	S = calculateSForNodalForcesKirshoff(E);
    }else if (neoHookeanMaterial){
    	gsl_matrix* tmpCforInversion =  gsl_matrix_calloc(nDim,nDim);
		gsl_matrix* InvC = gsl_matrix_calloc(nDim,nDim);
		createMatrixCopy(tmpCforInversion,C);
		bool inverted = InvertMatrix(tmpCforInversion, InvC);
		if (!inverted){
            std::cerr<<"C not inverted!!"<<std::endl;
		}
    	S = calculateSForNodalForcesNeoHookean(InvC,lnJ);
    	updateLagrangianElasticityTensorNeoHookean(InvC,lnJ,pointNo);

    	//I would like to keep a record of strains, therefore I am repeating this calculation here,
    	//it does not contribute to force calculation
    	E = calculateEForNodalForcesKirshoff(C);
    	gsl_matrix_set_zero(Strain);

		gsl_matrix_set(Strain,0,0, gsl_matrix_get(E,0,0));
		gsl_matrix_set(Strain,1,0, gsl_matrix_get(E,1,1));
		gsl_matrix_set(Strain,2,0, gsl_matrix_get(E,2,2));
		gsl_matrix_set(Strain,3,0, 2.0*gsl_matrix_get(E,0,1));
		gsl_matrix_set(Strain,4,0, 2.0*gsl_matrix_get(E,2,1));
		gsl_matrix_set(Strain,5,0, 2.0*gsl_matrix_get(E,0,2));
		gsl_matrix_free(tmpCforInversion);
		gsl_matrix_free(InvC);
    }
    //calculating elastic stress (elastic stress = detFe^-1 Fe S Fe^T):
    gsl_matrix_set_zero(elasticStress[pointNo]);
    gsl_matrix* compactStress  = calculateCompactStressForNodalForces(detFe, currFe,S, elasticStress[pointNo]);

    //Now from elastic stress, I will calculate nodal forces via B.
    //Calculating the inverse Jacobian stack matrix:
    gsl_matrix* InvJacobianStack = calculateInverseJacobianStackForNodalForces(Jacobian);

    //Calculating currB^T:
    detFs[pointNo] = determinant3by3Matrix(currF);
    gsl_matrix* currBT = calculateBTforNodalForces(InvJacobianStack,ShapeFuncDerStack, B, invJShFuncDerS);


    //Calculate invJShapeFuncDerStackwithFe for K calculation (using F in inverse jacobian calculation rather than Fe):
    calculateInvJShFuncDerSWithFe(currFe, InvdXde, ShapeFuncDerStack, invJShFuncDerSWithFe);

    //calculating nodal elastic forces as B^T compactStress detF
    gsl_matrix_scale(currBT,detFs[pointNo]);
    gsl_matrix_scale(currBT,detdXdes[pointNo]);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, currBT, compactStress,0.0, currge);


    //calculating viscous stress:
    gsl_matrix_set_zero(viscousStress[pointNo]);
    if (internalViscosity != 0){
    	gsl_matrix* l = calculateVelocityGradientTensor(B, displacementPerDt);
    	gsl_matrix* d = calculateRateOfDeformationTensor(l);
    	calculateViscousStress(d,viscousStress[pointNo]);
    	//The Bt I am giving into this function has already been scaled to include volume integration.
    	calculateViscousForces(currgv, currBT,viscousStress[pointNo]);
        gsl_matrix_free(l);
        gsl_matrix_free(d);
    }
    //freeing the matrices allocated in this function
    gsl_matrix_free(C);
    gsl_matrix_free(E);
    gsl_matrix_free(S);
    gsl_matrix_free(compactStress);
    gsl_matrix_free(InvJacobianStack);
    gsl_matrix_free(currBT);
    gsl_matrix_free(currFe);
    gsl_matrix_free(CurrShape);
    gsl_matrix_free(Jacobian);

    //gsl_matrix_free(viscousStress);
}



void Prism::calculateReferenceVolume(){
    size_t nTriangularFaces = 8;
    std::vector<std::array<int,3>>  triangularFaces;
    for (size_t i=0; i<nTriangularFaces; ++i){
        std::array<int,3> zeros = {0};
        triangularFaces.push_back(zeros);
    }
	//two bases:
	triangularFaces[0][0] = 0; triangularFaces[0][1] = 1; triangularFaces[0][2] = 2;
	triangularFaces[1][0] = 3; triangularFaces[1][1] = 4; triangularFaces[1][2] = 5;
	//divide the rectangular sides into trianges:
	triangularFaces[2][0] = 0; triangularFaces[2][1] = 1; triangularFaces[2][2] = 3;
	triangularFaces[3][0] = 1; triangularFaces[3][1] = 3; triangularFaces[3][2] = 4;
	triangularFaces[4][0] = 1; triangularFaces[4][1] = 2; triangularFaces[4][2] = 4;
	triangularFaces[5][0] = 2; triangularFaces[5][1] = 4; triangularFaces[5][2] = 5;
	triangularFaces[6][0] = 0; triangularFaces[6][1] = 2; triangularFaces[6][2] = 3;
	triangularFaces[7][0] = 2; triangularFaces[7][1] = 3; triangularFaces[7][2] = 5;
	//calculating the centre point of the reference shape so that I will have the
	//mid point
    std::array<double,3> centre = {0.0};
    for (size_t i = 0; i < nNodes; ++i){
        for (size_t j = 0 ;j < nDim; ++j){
			centre[j] += ReferenceShape->Positions[i][j];
		}
	}
    for (auto& centreCoord : centre){
        centreCoord /= nNodes;
	}
	ReferenceShape->Volume = calculateVolumeForInputShapeStructure(ReferenceShape->Positions, nTriangularFaces, triangularFaces,centre);
	GrownVolume = ReferenceShape->Volume;
	VolumePerNode = GrownVolume/nNodes;

	//calculating the basal area:
    std::array<double,3> vec1 = {0.0};
    std::array<double,3> vec2 = {0.0};
    for (size_t j=0 ;j < nDim; ++j){
			vec1 [j] = ReferenceShape->Positions[1][j] - ReferenceShape->Positions[0][j];
			vec2 [j] = ReferenceShape->Positions[2][j] - ReferenceShape->Positions[0][j];
	}
    std::array<double,3> baseVec = crossProduct3D(vec1,vec2);
	double normBaseVec= calculateMagnitudeVector3D (baseVec);
	double baseArea= normBaseVec/2;
	ReferenceShape->BasalArea = baseArea;
}

void Prism::checkHealth(){
    std::array<std::array<double,3>,8>  normals = calculatePlaneNormals(); //normals (8,3)
	bool elementsAreHealthy = checkNodePlaneConsistency(normals);
	if (!elementsAreHealthy){
        std::cerr<<" Element not healthy! : "<<Id<<std::endl;
	}
}

std::array<std::array<double,3>,8> Prism::calculatePlaneNormals(){
	//Calculating plane normals:
	//plane 0 -> normal 0 - > normal for nodes 0  1  2
	//plane 1 -> normal 1 - > normal for nodes 0  1  3
	//plane 2 -> normal 2 - > normal for nodes 0  2  3
	//plane 3 -> normal 3 - > normal for nodes 3  4  5
	//plane 4 -> normal 4 - > normal for nodes 3  1  4
	//plane 5 -> normal 5 - > normal for nodes 3  2  5
	//plane 6 -> normal 6 - > normal for nodes 1  2  4
	//plane 7 -> normal 7 - > normal for nodes 2  4  5
	int List[8][3]={{0,1,2},{0,3,1},{0,2,3},{3,5,4},{3,4,1},{3,2,5},{1,4,2},{2,4,5}};
    std::array<std::array<double,3>,8>  normals; //normals (8,3)
    for (auto& currentNormal : normals){
        currentNormal = {0.0, 0.0, 0.0};
    }
	for (int i=0; i<8; ++i){
        std::array<double,3> u = assignNodalVector(List[i][0],List[i][1]);
        std::array<double,3> v = assignNodalVector(List[i][0],List[i][2]);
        std::array<double,3> currNormal = crossProduct3D(u,v);
        normals[i] = {currNormal[0],currNormal[1],currNormal[2]};
    }
    return normals;
}

std::array<double,3> Prism::assignNodalVector(size_t id0, size_t id1){
	//vector from NodeId id0 to NodeId id1
    if (id0>=nNodes || id1>= nNodes){
        std::cerr<<"Error in node input in nodal vector assignment!"<<std::endl;
	}
    std::array<double,3> vec = {0.0, 0.0, 0.0};
    for (size_t i=0; i<nDim; ++i){
        vec[i] = Positions[id1][i] - Positions[id0][i];
	}
    return vec;
}

void 	Prism::setInitialEdgeLenghts(){
	int node0 = 0, node1 = 1;
	for (int i=0;i<3; ++i){
		if (i ==1 ){
			node1 = 2;
		}
		if (i ==2 ){
			node0 = 1;
		}
		double dx = Positions[node0][0]-Positions[node1][0];
		double dy = Positions[node0][1]-Positions[node1][1];
		double dz = Positions[node0][2]-Positions[node1][2];
		initialBasalEdgeLengthsSq[i] = dx*dx+dy*dy+dz*dz;
	}
	node0 = 3;
	node1 = 4;
	for (int i=0;i<3; ++i){
		if (i ==1 ){
			node1 = 5;
		}
		if (i ==2 ){
			node0 = 4;
		}
		double dx = Positions[node0][0]-Positions[node1][0];
		double dy = Positions[node0][1]-Positions[node1][1];
		double dz = Positions[node0][2]-Positions[node1][2];
		initialApilcalEdgeLengthsSq[i] = dx*dx+dy*dy+dz*dz;
	}
}

bool Prism::areNodesDirectlyConnected(int node0, int node1){
    std::vector<int>::iterator node0Placement = std::find (NodeIds.begin(), NodeIds.end(), node0);
	int node1Placement=0;
    if (node0Placement != NodeIds.end()){
		//node is on id list:
        ptrdiff_t pos = distance(NodeIds.begin(), node0Placement);
		if(pos >=3){
			node1Placement = pos-3;
		}
		else{
			node1Placement = pos+3;
		}
        //std::cout<<"node0Placement "<<pos<<" node1Placement "<<node1Placement<<std::endl;
		return NodeIds[node1Placement] == node1;
	}
	return false;
}

void Prism::checkEdgeLenghtsForBinding(std::vector <int>&masterIds, std::vector <int>&slaveIds){
    /**
     * This function returns the nodes to collapse together if any of the element surfaces
     * are shrinking below the set threshold of 10% of the original side length (set in thresholdFraction variable).\n
     */
	double thresholdFraction = 0.1; //fraction of original length the edge shrunk to
	double thresholdAngle = M_PI*8.0/9.0;// M_PI*8.0/9.0; //160 degrees, this was 100000 for most of simulations, I did not collapse based on angle

    std::array<double,3> currentBasalEdgeLengthsSq = {0.0};
    std::array<double,3> currentApicalEdgeLengthsSq = {0.0};
	int node0 = 0, node1 = 1;
    /**
     * The top side of the element will only be checked if the element top surface is at the exposed top surface of the tissue (apical
     * surface, ShapeBase#tissuePlacement = 0 or the element spans the whole tissue, ShapeBase#spansWholeTissue = true).
     */
	if (tissuePlacement == 0 || spansWholeTissue == true){
		//check angles first, and fix if necessary:
        std::array<double,3>  vec01 = {0.0};
        std::array<double,3>  vec10 = {0.0};;
        std::array<double,3>  vec02 = {0.0};
        std::array<double,3>  vec12 = {0.0};
		vec01[0] = Positions[1][0]-Positions[0][0];
		vec01[1] = Positions[1][1]-Positions[0][1];
		vec01[2] = Positions[1][2]-Positions[0][2];
		vec02[0] = Positions[2][0]-Positions[0][0];
		vec02[1] = Positions[2][1]-Positions[0][1];
		vec02[2] = Positions[2][2]-Positions[0][2];
		vec12[0] = Positions[2][0]-Positions[1][0];
		vec12[1] = Positions[2][1]-Positions[1][1];
		vec12[2] = Positions[2][2]-Positions[1][2];
		double L01 = normaliseVector3D(vec01);
		double L02 = normaliseVector3D(vec02);
		double L12 = normaliseVector3D(vec12);
		vec10[0] = vec01[0]*-1;
		vec10[1] = vec01[1]*-1;
		vec10[2] = vec01[2]*-1;
		double dot0 =dotProduct3D(vec01,vec02);
		double dot1 =dotProduct3D(vec12,vec10);
		double tet0 = acos(dot0);
		double tet1 = acos(dot1);
		if (tet0 > thresholdAngle){
			//snap to closest
			if (L01 < L02){
				masterIds.push_back(NodeIds[0]);
				slaveIds.push_back(NodeIds[1]);
			}
			else{
				masterIds.push_back(NodeIds[0]);
				slaveIds.push_back(NodeIds[2]);
			}
		}
		if (tet1 > thresholdAngle){
			if (L01 < L12){
				masterIds.push_back(NodeIds[0]);
				slaveIds.push_back(NodeIds[1]);
			}
			else{
				masterIds.push_back(NodeIds[1]);
				slaveIds.push_back(NodeIds[2]);
			}
		}
		if ( M_PI - tet0 - tet1 > thresholdAngle){
			if (L02 < L12){
				masterIds.push_back(NodeIds[0]);
				slaveIds.push_back(NodeIds[2]);
			}
			else{
				masterIds.push_back(NodeIds[1]);
				slaveIds.push_back(NodeIds[2]);
			}
		}
        //check basal side only for the most basal section, which will belong to a basal element or an element that spans the whole tissue.
		for (int i=0;i<3; ++i){
			if (i ==1 ){
				node1 = 2;
			}
			if (i ==2 ){
				node0 = 1;
			}
			double dx = Positions[node0][0]-Positions[node1][0];
			double dy = Positions[node0][1]-Positions[node1][1];
			double dz = Positions[node0][2]-Positions[node1][2];
			currentBasalEdgeLengthsSq[i] = dx*dx+dy*dy+dz*dz;
			//checking angle:
		}

	}
    /**
     * The bottom side of the element will only be checked if the element bottom surface is at the exposed bottom surface of the tissue (basal
     * surface, ShapeBase#tissuePlacement = 1 or the element spans the whole tissue, ShapeBase#spansWholeTissue = true).
     */
	if (tissuePlacement == 1 || spansWholeTissue == true){
		//check angles first, and fix if necessary:
        std::array<double,3>  vec34 = {0.0};
        std::array<double,3>  vec43 = {0.0};
        std::array<double,3>  vec35 = {0.0};
        std::array<double,3>  vec45 = {0.0};
		vec34[0] = Positions[4][0]-Positions[3][0];
		vec34[1] = Positions[4][1]-Positions[3][1];
		vec34[2] = Positions[4][2]-Positions[3][2];
		vec35[0] = Positions[5][0]-Positions[3][0];
		vec35[1] = Positions[5][1]-Positions[3][1];
		vec35[2] = Positions[5][2]-Positions[3][2];
		vec45[0] = Positions[5][0]-Positions[4][0];
		vec45[1] = Positions[5][1]-Positions[4][1];
		vec45[2] = Positions[5][2]-Positions[4][2];
		double L34 = normaliseVector3D(vec34);
		double L35 = normaliseVector3D(vec35);
		double L45 = normaliseVector3D(vec45);
		vec43[0] = vec34[0]*-1;
		vec43[1] = vec34[1]*-1;
		vec43[2] = vec34[2]*-1;
		double dot3 =dotProduct3D(vec34,vec35);
		double dot4 =dotProduct3D(vec45,vec43);
		double tet3 = acos(dot3);
		double tet4 = acos(dot4);

		if (tet3 > thresholdAngle){
			//snap to closest
			if (L34 < L35){
				masterIds.push_back(NodeIds[3]);
				slaveIds.push_back(NodeIds[4]);
			}
			else{
				masterIds.push_back(NodeIds[3]);
				slaveIds.push_back(NodeIds[5]);
			}
		}
		if (tet4 > thresholdAngle){
			if (L34 < L45){
				masterIds.push_back(NodeIds[3]);
				slaveIds.push_back(NodeIds[4]);
			}
			else{
				masterIds.push_back(NodeIds[4]);
				slaveIds.push_back(NodeIds[5]);
			}
			//checkBasalLengths = false;
		}
		if ( M_PI - tet3 - tet4 > thresholdAngle){
			if (L35 < L45){
				masterIds.push_back(NodeIds[3]);
				slaveIds.push_back(NodeIds[5]);
			}
			else{
				masterIds.push_back(NodeIds[4]);
				slaveIds.push_back(NodeIds[5]);
			}
		}
		node0 = 3;
		node1 = 4;
		for (int i=0;i<3; ++i){
			if (i ==1 ){
				node1 = 5;
			}
			if (i ==2 ){
				node0 = 4;
			}
			double dx = Positions[node0][0]-Positions[node1][0];
			double dy = Positions[node0][1]-Positions[node1][1];
			double dz = Positions[node0][2]-Positions[node1][2];
			currentApicalEdgeLengthsSq[i] = dx*dx+dy*dy+dz*dz;
		}
	}
	bool node1isSlave = false;
	bool node2isSlave = false;
	bool node4isSlave = false;
	bool node5isSlave = false;
	if (tissuePlacement == 0 || spansWholeTissue == true){
		if (currentBasalEdgeLengthsSq[0]<thresholdFraction*initialBasalEdgeLengthsSq[0]){
			//node 0 is too close to node 1, bind
			int master = NodeIds[0];
			int slave  = NodeIds[1];
			masterIds.push_back(master);
			slaveIds.push_back(slave);
			node1isSlave = true;
		}
		if (currentBasalEdgeLengthsSq[1]<thresholdFraction*initialBasalEdgeLengthsSq[1]){
			//node 0 is too close to node 2, bind
			int master = NodeIds[0];
			int slave  = NodeIds[2];
			masterIds.push_back(master);
			slaveIds.push_back(slave);
			node2isSlave = true;
		}
		if (currentBasalEdgeLengthsSq[2]<thresholdFraction*initialBasalEdgeLengthsSq[2]){
			//node 1 is too close to node 2, bind
			if (!node1isSlave || !node2isSlave){
				int master = NodeIds[1];
				if (node1isSlave){
					master = NodeIds[0];
				}
				int slave  = NodeIds[2];
				masterIds.push_back(master);
				slaveIds.push_back(slave);
			}
		}
	}
	if (tissuePlacement == 1 || spansWholeTissue == true){
		if (currentApicalEdgeLengthsSq[0]<thresholdFraction*initialApilcalEdgeLengthsSq[0]){
			//node 3 is too close to node 4, bind
			int master = NodeIds[3];
			int slave  = NodeIds[4];
			masterIds.push_back(master);
			slaveIds.push_back(slave);
			node4isSlave = true;
		}
		if (currentApicalEdgeLengthsSq[1]<thresholdFraction*initialApilcalEdgeLengthsSq[1]){
			//node 3 is too close to node 5, bind
			int master = NodeIds[3];
			int slave  = NodeIds[5];
			masterIds.push_back(master);
			slaveIds.push_back(slave);
			node5isSlave = true;
		}
		if (currentApicalEdgeLengthsSq[2]<thresholdFraction*initialApilcalEdgeLengthsSq[2]){
			if (!node4isSlave || !node5isSlave){
				int master = NodeIds[4];
				if (node1isSlave){
					master = NodeIds[3];
				}
				int slave  = NodeIds[5];
				masterIds.push_back(master);
				slaveIds.push_back(slave);
			}
		}
	}
}


bool Prism::checkNodePlaneConsistency(std::array<std::array<double,3>,8>& normals){
    //std::cout<<"inside check consistency, Id: "<<Id<<std::endl;
	//List of constricting planes for each node:
	//format is for each node (0->5): [plane1, plane2, plane3, node id for plane1, node id  for plane2&3]
    size_t List[6][5] = {{3,6,7,3,2},{3,2,5,4,2},{3,1,4,5,1},{0,6,7,0,2},{0,2,5,1,2},{0,1,4,2,1}};
	bool elementHealthy = true;
    for (size_t i =0; i<nNodes; ++i){
        std::array<double,3> u = assignNodalVector(List[i][3],i);
		double dotp[3];
        dotp[0] = dotProduct3D(u,normals[List[i][0]]);
		if (dotp[0]<0){
            std::cerr <<"The element is not consistent! - top/bottom plane, Id: "<<Id<<std::endl;
            std::cerr<<"i: "<<i<<std::endl;
            std::cerr<<"normals["<<List[i][0]<<"]: "<<normals[List[i][0]][0]<<" "<<normals[List[i][0]][1]<<" "<<normals[List[i][0]][2]<<std::endl;
            std::cerr<<"u: "<<u[0]<<" "<<u[1]<<" "<<u[2]<<std::endl;
            std::cerr<<"dotp: "<<dotp[0]<<std::endl;
			elementHealthy =  false;
		}
        u = assignNodalVector(List[i][4],i);
        dotp[1] = dotProduct3D(u,normals[List[i][1]]);
        dotp[2] = dotProduct3D(u,normals[List[i][2]]);
		if(dotp[1]<0 ||  dotp[2]<0){
            std::cerr <<"The element is not consistent! side planes, Id: "<<Id<<std::endl;
            std::cerr<<"dot 1: "<<dotp[1]<<" dot2 :"<<dotp[2]<<std::endl;
            std::cerr<<"1 : normals["<<List[i][1]<<"]: "<<normals[List[i][1]][0]<<" "<<normals[List[i][1]][1]<<" "<<normals[List[i][1]][2]<<std::endl;
            std::cerr<<"2 : normals["<<List[i][2]<<"]: "<<normals[List[i][2]][0]<<" "<<normals[List[i][2]][1]<<" "<<normals[List[i][2]][2]<<std::endl;
            for (size_t i=0; i<nNodes;++i){
                for (size_t j =0; j<nDim; ++j){
                    std::cout<<Positions[i][j]<<"  ";
				}
                std::cout<<std::endl;
			}
			elementHealthy =  false;
		}
	}
	return elementHealthy;
}

double Prism::getApicalSideLengthAverage(){
	double dx,dy,dz;
	double dsum =0.0;
	int counter = 0;
	int pairs[3][2] = {{3,4},{3,5},{4,5}};
	if (tissueType == 1){
		//peripodial
		pairs[0][0] = 0;
		pairs[0][1] = 1;
		pairs[1][0] = 0;
		pairs[1][1] = 2;
		pairs[2][0] = 1;
		pairs[2][1] = 2;
	}
	for (int i=0; i<3; ++i){
		dx = Positions[pairs[i][0]][0] - Positions[pairs[i][1]][0];
		dy = Positions[pairs[i][0]][1] - Positions[pairs[i][1]][1];
		dz = Positions[pairs[i][0]][2] - Positions[pairs[i][1]][2];
		double dsum2 = (dx*dx + dy*dy+dz*dz);
		if(dsum2 > 0){
			//the average side length can be zero for elements with collapsed nodes.
			//I do not want to add the zero length in calculation of average length.
			//I will handle the case where this function returns zero outside, in averaging all elements.
			dsum += pow((dx*dx + dy*dy+dz*dz),0.5);
			counter++;
		}
	}
	if (counter>0){
		dsum /= counter;
	}
	return dsum;
}

double Prism::getBasalSideLengthAverage(){
	double dx,dy,dz;
	double dsum =0.0;
	int counter = 0;
	int pairs[3][2] = {{0,1},{0,2},{1,2}};
	if (tissueType == 1){
		//peripodial
		pairs[0][0] = 3;
		pairs[0][1] = 4;
		pairs[1][0] = 3;
		pairs[1][1] = 5;
		pairs[2][0] = 4;
		pairs[2][1] = 5;
	}
	for (int i=0; i<3; ++i){
		dx = Positions[pairs[i][0]][0] - Positions[pairs[i][1]][0];
		dy = Positions[pairs[i][0]][1] - Positions[pairs[i][1]][1];
		dz = Positions[pairs[i][0]][2] - Positions[pairs[i][1]][2];
		double dsum2 = (dx*dx + dy*dy+dz*dz);
		if(dsum2 > 0){
			//the average side length can be zero for elements with collapsed nodes.
			//I do not want to add the zero length in calculation of average length.
			//I will handle the case where this function returns zero outside, in averaging all elements.
			dsum += pow((dx*dx + dy*dy+dz*dz),0.5);
			counter++;
		}
	}
	if (counter>0){
		dsum /= counter;
	}
	return dsum;
}

int Prism::getCorrecpondingApical(int currNodeId){
	if (NodeIds[0] == currNodeId){
		return NodeIds[3];
	}
	if (NodeIds[1] == currNodeId){
		return NodeIds[4];
	}
	if (NodeIds[2] == currNodeId){
		return NodeIds[5];
	}
	return -100;
}

bool Prism::IsThisNodeMyBasal(int currNodeId){
	if (NodeIds[0] == currNodeId || NodeIds[1] == currNodeId || NodeIds[2] == currNodeId ){
		return true;
	}
	return false;
}


bool Prism::IsThisNodeMyApical(int currNodeId){
	if (NodeIds[3] == currNodeId || NodeIds[4] == currNodeId || NodeIds[5] == currNodeId ){
		return true;
	}
	return false;
}

double Prism::getElementHeight(){
	double dx = Positions[0][0] - Positions[3][0];
	double dy = Positions[0][1] - Positions[3][1];
	double dz = Positions[0][2] - Positions[3][2];
	return pow((dx*dx + dy*dy + dz*dz),0.5);
}

void 	Prism::calculateBasalArea(){
    double Threshold = 1E-5;
	int id0 = 0, id1 = 1, id2 = 2; // this is correct for basal side, I will change it for apical calculation
	if (tissueType == 1){ //peroipodial element, the basal surface is actually on top
		id0 = 3;
		id1 = 4;
		id2 = 5;
	}
	double sideVec1[3];
	double sideVec2[3];
	double Side1 = 0.0;
	double Side2 = 0.0;
	double costet = 0.0;
	double Area = 0.0;

	for (int i = 0; i<3; ++i){
		sideVec1[i]= Positions[id1][i] - Positions[id0][i];
		sideVec2[i]= Positions[id2][i] - Positions[id0][i];
		costet += sideVec1[i] * sideVec2[i];
		Side1  += sideVec1[i] * sideVec1[i];
		Side2  += sideVec2[i] * sideVec2[i];
	}
	if (Side1 > Threshold && Side2 > Threshold){
		Side1 = pow(Side1,0.5);
		Side2 = pow(Side2,0.5);
		costet /= (Side1*Side2);
		double sinTetSq  = 1-costet*costet;
		double sintet = 0;
		if(sinTetSq>Threshold){
			sintet = pow(sinTetSq,0.5);
		}
		Area = Side1* Side2 * sintet / 2.0;
	}
	BasalArea = Area;
    //std::cout<<" Element "<<Id<<" basal area: "<<BasalArea<<std::endl;
	//ApicalArea = Area;
}

void 	Prism::calculateApicalArea(){
    double Threshold = 1E-5;
	int id0 = 3, id1 = 4, id2 = 5; // this is correct for basal side, I will change it for apical calculation
	if (tissueType == 1){ //peroipodial element, the basal surface is actually at the bottom
		id0 = 0;
		id1 = 1;
		id2 = 2;
	}
	double sideVec1[3];
	double sideVec2[3];
	double Side1 = 0.0;
	double Side2 = 0.0;
	double costet = 0.0;
	double Area = 0.0;
	for (int i = 0; i<3; ++i){
		sideVec1[i]= Positions[id1][i] - Positions[id0][i];
		sideVec2[i]= Positions[id2][i] - Positions[id0][i];
		costet += sideVec1[i] * sideVec2[i];
		Side1  += sideVec1[i] * sideVec1[i];
		Side2  += sideVec2[i] * sideVec2[i];
	}
	if (Side1 > Threshold && Side2 > Threshold){
		Side1 = pow(Side1,0.5);
		Side2 = pow(Side2,0.5);
		costet /= (Side1*Side2);
		double sinTetSq  = 1-costet*costet;
		double sintet = 0;
		if(sinTetSq>Threshold){
			sintet = pow(sinTetSq,0.5);
		}
		Area = Side1* Side2 * sintet / 2.0;
	}
	ApicalArea = Area;
}

void Prism::setBasalNeigElementId(const std::vector<std::unique_ptr<ShapeBase>>& elementsList){
    /**
     * This funciton sets the basal neighbours of all elemetns qualifying for a basal neighbour element.
     * The check is only carried out on the apical elements of the  coumnar tissue (ShapeBase#tissuePlacement = 0, ShapeBase#TissueType = 0).
     * The columnar elements bordering the basal extracellualar matrix (ShapeBase#atBasalBorderOfECM = true),
     * and the basal extracellular matrix itself (ShapeBase#isECMMimicing) are not
     * labelled, as these elements do not hold a basal neighbour that is a tissue piece.
     */
    if (tissueType== 0 //columnar layer element
			&& !tissuePlacement == 0 //not checking for basal elements, to save time
			&& !atBasalBorderOfECM  //not checking for elements bordering ECM, same as basal
			&& !isECMMimicing		 //not checking for ECM mimicking elements
			){
        /** Basal nodes of a prism are hte  nodes recorded in ShapeBase#NodeIds array indices 0-2.
         */
		int currBasalNodes[3] = {NodeIds[0], NodeIds[1],NodeIds[2]};
        for( const auto& itElement : elementsList){
            bool isApicalOwner  = itElement->IsThisNodeMyApical(currBasalNodes[0]);
			if (isApicalOwner){
				//node 0 is owned by this element, is node 1?
                isApicalOwner  = itElement->IsThisNodeMyApical(currBasalNodes[1]);
				if (isApicalOwner){
					//node 0 & 1 are owned by this element, is node 2?
                    isApicalOwner  = itElement->IsThisNodeMyApical(currBasalNodes[2]);
					if (isApicalOwner){
                        /** When all three nodes are owned as apical nodes by another element, this is the
                         * neighbour I am looking for. If the element is NOT an explicit ECM element (ShapeBase#isECMMimicing = false),
                         * then the new Id is recorded into the ShapeBase#basalNeigElementId.
                         */
                        if (!itElement-> isECMMimicing){
                            basalNeigElementId = itElement->getId();
						}
                        /** Once the element is found, break the loop regardless of the ECM status. If the loop has
                         *  reached an ECM element then this element does not have a basal bordering tissue. It should not
                         * have this case anyway, as elemetns atBasalBorderOfECM are skipped.
                         */
						break;
					}
				}
			}
		}
	}
}

void Prism::constructElementStackList(const int discretisationLayers, const std::vector<std::unique_ptr<ShapeBase>>& elementsList){
	//This is a basal element. I will take my basal nodes, and find elements that have them as apical nodes.
    //filling the array with current known apical node ids
    int currApicalNodes[3] = {NodeIds[3], NodeIds[4],NodeIds[5]};
    int currElementId = Id;
	//Starting the list with the basal element
    elementsIdsOnSameColumn.push_back(currElementId);
    for (size_t i=1; i<  size_t(discretisationLayers); ++i ){
        for( const auto& itElement : elementsList){
            int checkedElementId = itElement->getId();
			if (checkedElementId != currElementId){ //the iterator is not this element
				//I will find the elemetns that have this elements apical nodes as basal. All three ndes should match;
                bool IsBasalOwner = itElement->IsThisNodeMyBasal(currApicalNodes[0]);
				if (IsBasalOwner){
                    IsBasalOwner = itElement->IsThisNodeMyBasal(currApicalNodes[1]);
				}
				if (IsBasalOwner){
                    IsBasalOwner = itElement->IsThisNodeMyBasal(currApicalNodes[2]);
				}
				//If IsBasalOwner is true at this point, then I found the matching element.
				if (IsBasalOwner){
					//I found the element:
                    elementsIdsOnSameColumn.push_back(checkedElementId);
                    for (size_t j=0; j<3; ++j){
                        //moving apical node ids to the next elemetn, apical of this is the basal of next
                        currApicalNodes[j] = itElement->getCorrecpondingApical(currApicalNodes[j]);
					}
					currElementId = checkedElementId;
					break;
				}
			}
		}
	}

	//Now I have the list for the basal element.
	//I should construct it for at midline and apical elements:
    for (size_t i=1; i<  size_t(discretisationLayers); ++i ){
        size_t currElementOnTheSameColumnId  = elementsIdsOnSameColumn[i];
        elementsList[currElementOnTheSameColumnId]->elementsIdsOnSameColumn.insert(elementsList[currElementOnTheSameColumnId]->elementsIdsOnSameColumn.end(),elementsIdsOnSameColumn.begin(),elementsIdsOnSameColumn.end());
    }
}

void Prism::assignExposedSurfaceAreaIndices(const std::vector<std::unique_ptr<Node>>& Nodes){
	if (tissueType == 0 || tissueType == 1){//columnar or peripodial Element
		//these are written for columnar tissue, should be swapped for peorpodial tissue:
		int apicalIds[3] = {3,4,5};
		int basalIds[3] = {0,1,2};
		if (tissueType == 1){//peripodial Element
			apicalIds[0]= 0;
			apicalIds[1]= 1;
			apicalIds[2]= 2;
			basalIds[0]= 3;
			basalIds[1]= 4;
			basalIds[2]= 5;
		}
		if (tissuePlacement == 1){ //apical element
			elementHasExposedApicalSurface = true;
			//If element is apical, I will calculate the apical surface
			exposedApicalSurfaceNodeIds[0] = apicalIds[0];
			exposedApicalSurfaceNodeIds[1] = apicalIds[1];
			exposedApicalSurfaceNodeIds[2] = apicalIds[2];
		}
		else if(tissuePlacement == 0){ //basal element
			elementHasExposedBasalSurface = true;
			//If element is basal, I will calculate the basal surface
			exposedBasalSurfaceNodeIds[0] = basalIds[0];
			exposedBasalSurfaceNodeIds[1] = basalIds[1];
			exposedBasalSurfaceNodeIds[2] = basalIds[2];
		}
		else if (tissuePlacement == 2 && spansWholeTissue) {//midline element, but spans whole tissue
			elementHasExposedApicalSurface = true;
			elementHasExposedBasalSurface = true;
			exposedApicalSurfaceNodeIds[0] = apicalIds[0];
			exposedApicalSurfaceNodeIds[1] = apicalIds[1];
			exposedApicalSurfaceNodeIds[2] = apicalIds[2];
			exposedBasalSurfaceNodeIds[0] = basalIds[0];
			exposedBasalSurfaceNodeIds[1] = basalIds[1];
			exposedBasalSurfaceNodeIds[2] = basalIds[2];
		}
	}
	else if (tissueType == 2){//linker element
		int numberOfApicalNodes = 0;
		int numberOfBasalNodes = 0;
        for (size_t i=0; i<nNodes; ++i){
			if (Nodes[NodeIds[i]]->tissuePlacement == 0){ //basal Node
				exposedLateralAreaBasalSideNodeIds[numberOfBasalNodes] = i;
				numberOfBasalNodes++;
			}
			if (Nodes[NodeIds[i]]->tissuePlacement == 1){ //apical Node
				exposedLateralAreaApicalSideNodeIds[numberOfApicalNodes] = i;
				numberOfApicalNodes++;
			}
		}
	}
}
