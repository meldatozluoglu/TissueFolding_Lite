#include "ShapeBase.h"
#include "Node.h"
#include <sstream>

#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>


void 	ShapeBase::ParentErrorMessage(std::string functionName){
    std::cerr<<"You are calling the function: "<<functionName<<" from a parent here, check declaration is via pointers"<<std::endl;
}

bool 	ShapeBase::ParentErrorMessage(std::string functionName, bool returnValue){
    std::cerr<<"You are calling the function: "<<functionName<<" from a parent here, check declaration is via pointers"<<std::endl;
	return returnValue;
}

double 	ShapeBase::ParentErrorMessage(std::string functionName, double returnValue){
    std::cerr<<"You are calling the function: "<<functionName<<" from a parent here, check declaration is via pointers"<<std::endl;
	return returnValue;
}

int 	ShapeBase::ParentErrorMessage(std::string functionName, int returnValue){
    std::cerr<<"You are calling the function: "<<functionName<<" from a parent here, check declaration is via pointers"<<std::endl;
	return returnValue;
}

void	ShapeBase::setShapeType(std::string TypeName){
	/**
	 *  The function will set the shape type of the element, stored in variable ShapeBase#ShapeType.
	 *  The mapping is as follows:
	 *  	- 1 	: Prism
	 *  	- 2 	: PrismLateral
	 *  	- 3 	: Tetrahedron
	 *  	- 4 	: Triangle
	 *  	- (-100): Default number if input name is not recognised.
	 *
	 */
	if (TypeName == "Prism"){
		this->ShapeType = 1;
	}
	else if (TypeName == "PrismLateral"){
		this->ShapeType = 2;
	}
	else if (TypeName == "Tetrahedron"){
		this->ShapeType = 3;
	}
	else if (TypeName == "Triangle"){
        //std::cout<<"set shape type to triangle"<<std::endl;
		this->ShapeType = 4;
	}
	else{
		this->ShapeType= -100;
	};
    //std::cout<<"finalised set shape type"<<std::endl;
}

void	ShapeBase::setIdentificationColour(){
	/**
	 *  The function sets the unique ShapeBase#IdentifierColour, which on click, will allow the
	 *  user interface to identify this element as the selected element. ShapeBase#IdentifierColour
	 *  is an integer array of size 3, and corresponds to rgb channels, in 0-255 range. Therefore, the
	 *  selection tool works for 255^3 = 16581375 elements.
	 *  The colour storing is started from b channel moving towards r, and done with element Id
	 *
	 */
	IdentifierColour[2] = Id % 255;
	int a = (Id - IdentifierColour[2]) / 255;
	IdentifierColour[1] = ( a ) % 255;
	if (a>255){
		IdentifierColour[0] = (a - IdentifierColour[1]) / 255;
	}
	else{
		IdentifierColour[0] = 0;
	}
}

int		ShapeBase::getShapeType(){
	return ShapeType;
}

int 	ShapeBase::getId(){
	return Id;
}

size_t 	ShapeBase::getNodeNumber(){
	return nNodes;
}

const std::vector<int>&	ShapeBase::getNodeIds(){
	return NodeIds;
}

int		ShapeBase::getNodeId(int i){
	return NodeIds[i];
}

double	ShapeBase::getApicalArea(){
	return ApicalArea;
}

size_t 	ShapeBase::getDim(){
	return nDim;
}

std::string ShapeBase::getName(){
    std::string name;
	if (ShapeType == 1){
		name = "Prism";
	}
	else if (ShapeType == 2){
		name = "PrismLateral";
	}
	else if (ShapeType == 3){
		name = "Tetrahedron";
	}
	else if (ShapeType == 4){
		name = "Triangle";
	}
	else{
		name = "Unknown";
	}
    std::stringstream inter;
	inter.fill('0');
	inter.width(4);
	inter<<Id;
	name = name + inter.str();
	return name;
}

const std::vector<std::array<double,3>>& ShapeBase::getReferencePos(){
	return ReferenceShape->Positions;
}

void ShapeBase::getPos(gsl_matrix* Pos){
    for (size_t i=0; i<nNodes; ++i){
        for (size_t j =0; j<nDim; ++j){
            gsl_matrix_set (Pos, i, j, Positions[i][j]);
        }
    }
}

double 	ShapeBase::getInternalViscosity(){
	return internalViscosity;
}

double 	ShapeBase::getOriginalInternalViscosity(){
	return originalInternalViscosity;
}

double 	ShapeBase::getZRemodellingSoFar(){
	return zRemodellingSoFar;
}

double 	ShapeBase::getStiffnessMultiplier(){
	return stiffnessMultiplier;
}

void 	ShapeBase::setZRemodellingSoFar(double zRemodellingSoFar){
	this -> zRemodellingSoFar = zRemodellingSoFar;
}

double 	ShapeBase::getYoungModulus(){
	return (stiffnessMultiplier*E);
}

double 	ShapeBase::getPoissonRatio(){
	return v;
}

const std::array<double,3>& ShapeBase::getGrowthRate(){
	return GrowthRate;
}

gsl_matrix* ShapeBase::getFg(){
    gsl_matrix* tmpFg =gsl_matrix_calloc(nDim, nDim);
    createMatrixCopy(tmpFg,Fg);
    return tmpFg;
}

gsl_matrix* ShapeBase::getFe(){
    gsl_matrix* tmpFe =gsl_matrix_calloc(nDim, nDim);
    for (int iter =0; iter<3;++iter){
		gsl_matrix_add(tmpFe, FeMatrices[iter]);
	}
    gsl_matrix_scale(tmpFe,1.0/3.0);
    return tmpFe;
}

gsl_matrix* ShapeBase::getInvFg(){
    gsl_matrix* tmpInvFg =gsl_matrix_calloc(nDim, nDim);
    createMatrixCopy(tmpInvFg,InvFg);
    return tmpInvFg;
}


gsl_matrix* ShapeBase::getFsc(){
    gsl_matrix* tmpFsc =gsl_matrix_calloc(nDim, nDim);
    createMatrixCopy(tmpFsc,Fsc);
    return tmpFsc;
}

gsl_matrix* ShapeBase::getInvFsc(){
    gsl_matrix* tmpInvFsc =gsl_matrix_calloc(nDim, nDim);
    createMatrixCopy(tmpInvFsc,InvFsc);
    return tmpInvFsc;
}

void ShapeBase::createMatrixCopy(gsl_matrix* dest, gsl_matrix* src){
    int m = src->size1;
    int n = src->size2;
    gsl_matrix_set_zero(dest);
    double tmp= 0.0;
    for (int i=0; i<m; ++i){
        for (int j=0 ; j<n; ++j){
            tmp = gsl_matrix_get(src,i,j);
            gsl_matrix_set(dest,i,j,tmp);
        }
    }
}

const std::array<double,6>& ShapeBase::getShapeChangeRate(){
	return ShapeChangeRate;
}

std::array<double,3> ShapeBase::getCentre(){
    std::array<double,3> c {0.0};
    for (size_t i = 0; i<nNodes; ++i ){
        for (size_t j = 0; j< nDim; ++j){
            c[j] += Positions[i][j];
		}
	}
    c[0] /= nNodes;
    c[1] /= nNodes;
    c[2] /= nNodes;
    return c;
}

double ShapeBase::getPeripodialness(){
	return peripodialGrowthWeight;
}

double ShapeBase::getColumnarness(){
	return columnarGrowthWeight;
}


void ShapeBase::getRelativePositionInTissueInGridIndex(int nGridX, int nGridY , int& IndexX, int& IndexY, double& FracX, double& FracY){
    /**
     *  This function provides the relative position within the bounding box of the tissue, and
     * calculates which point on the growth maps should be read. The relative position is calculated through
     * the function ShapeBase#getRelativePosInBoundingBox, and converted to grid indices through
     * ShapeBase#convertRelativePosToGridIndex.
     *
     */
    std::array<double,2> reletivePos = getRelativePosInBoundingBox();
	convertRelativePosToGridIndex(reletivePos, IndexX, IndexY, FracX, FracY, nGridX, nGridY);
}

void ShapeBase::getInitialRelativePositionInTissueInGridIndex(int nGridX, int nGridY, int& IndexX, int& IndexY, double& FracX, double& FracY){
    /**
     *  This function provides the relative initial position within the bounding box of the tissue, and
     * calculates which point on the growth maps should be read. The relative position is calculated through
     * the function ShapeBase#getInitialRelativePosInBoundingBox, and converted to grid indices through
     * ShapeBase#convertRelativePosToGridIndex.
     *
     */
    std::array<double,2> reletivePos = getInitialRelativePosInBoundingBox();
	convertRelativePosToGridIndex(reletivePos, IndexX, IndexY, FracX, FracY, nGridX, nGridY);
}

bool ShapeBase::isGrowthRateApplicable( int sourceTissue, double& weight, double zmin, double zmax){
    /**
     * This function checks if the  current growth is appliable tothe element. The input sourceTissue defines which tissue
     * compartment the growth is applicable to. The weight is needed for elements that are in transition zones
     * between compartments, giving the weight of the growth that should be applide to the element, should be below
     * 1.0 by definition. \n
     *
     * The code first chacks if the relative z coordinate of teh elemetn is within the limits of growth application.
     */
    //weight is the weight of the current tissue in linker sites
	if (initialRelativePositionInZ < zmin ||  initialRelativePositionInZ > zmax){
		return false;
	}
    /**
     * Then the tissue type is checked. If compartments match, the weight is zero, if they don't it is 0 and growth is
     * not applied. If the element is in a transition compartment, the weight is assigned dependent on the
     * source tissue type.
     */
	if (sourceTissue == 0){//columnar layer growth
		if (tissueType == 0){ //columnar
			weight = 1.0;
			return true;
		}
        else if(tissueType == 2){ //linker in transition zone
			weight = columnarGrowthWeight;
			return true;
		}
	}
	else if (sourceTissue == 1){//peripodial membrane growth
		if (tissueType == 1){ //peripodial
			weight = 1.0;
			return  true;
		}
		else if ( tissueType == 2) { //linker
			weight = peripodialGrowthWeight;
			return true;
		}
	}
	return false;
}

void ShapeBase::updateGrowthWillBeScaledDueToApikobasalRedistribution(bool thisFunctionShrinksApical, std::vector<int>& ellipseBandIdsForGrowthRedistribution){
	if(tissueType == 0){ //columnar tissue
		bool insideEllipseBandWithRedistribution = false;
		if (!isECMMimicing && insideEllipseBand){
            for (size_t i=0; i<ellipseBandIdsForGrowthRedistribution.size(); ++i){
				if (coveringEllipseBandId == ellipseBandIdsForGrowthRedistribution[i]){
					insideEllipseBandWithRedistribution =  true;
					break;
				}
			}
		}
		if (insideEllipseBandWithRedistribution){
			thereIsGrowthRedistribution = true;
			if ((tissuePlacement == 1 ) || atApicalBorderOfActin) { //apical
				//this element is apically positioned, it should shrink if the function is shrinking
				//the apical layer, or expand it if not.
				growthRedistributionShrinksElement = thisFunctionShrinksApical;
			}
			else{
				//the element is basally positioned, if the function is shrinking apical, this element should
				//expand, and vice versa.
				growthRedistributionShrinksElement = !thisFunctionShrinksApical;
			}
		}
	}
}

void ShapeBase::scaleGrowthForZRedistribution( double& x, double& y, double& z){
    /** This function will modify the incremental growth deformation gradient of the element to reflect the
     * volume redistribution in the height of the tissue. First check point id for if there is such
     * redistribution in simulaiton. \n
     *
     */
    if(thereIsGrowthRedistribution){ //columnar tissue
        /**
         * IF there is distribution, hourly xy-plane growth is calculated. z growth is not affected.
         */
		double growthRedistributionTime = 24*3600; //apply in terms of the effect in 24 hours
		double shrunkRates[3] = {0,0,0}, increasedRates[3] = {0,0,0};
		double xyGrowth = exp((x+y)*growthRedistributionTime);
		if (growthRedistributionShrinksElement){
            /**
             * If the element should be loosing volume (shrunk) due to redistrtibution, the shrink rate is applied.
             */
			double scaleFactorShrinkage = log((xyGrowth-1)*growthRedistributionScale + 1)/(x+y)/growthRedistributionTime;
			shrunkRates[0] = scaleFactorShrinkage*x;
			shrunkRates[1] = scaleFactorShrinkage*y;
			shrunkRates[2] = z;
			x = shrunkRates[0];
			y = shrunkRates[1];
			z = shrunkRates[2];
		}
		else{
            /**
             * Else, the sclae is inverted for growth and applied.
             */
			double growthRedistributionScaleComplementary = 2-growthRedistributionScale;
			double scaleFactorExpansion = log((xyGrowth-1)*growthRedistributionScaleComplementary + 1)/(x+y)/growthRedistributionTime;
			increasedRates[0] = scaleFactorExpansion*x;
			increasedRates[1] = scaleFactorExpansion*y;
			increasedRates[2] = z;
			x = increasedRates[0];
			y = increasedRates[1];
			z = increasedRates[2];
		}
	}
}


double ShapeBase::getCurrentVolume(){
    /**
     * The current volume of the element is calculated from the determinant of the deformation
     * gradient and the initial volume at relaxed state. As the determinants are already recorded, this is cheaper than
     * calculation of a compex deformed shape calculation. First the average of the determinants at all Gauss points
     * (ShapeBase#numberOfGaussPoints with weights in ShapeBase#gaussWeights), as recorded in ShapeBase#detFs is calculated.
     * Then the reference shape volume is scaled.
     */
	double J = 0;
    for (size_t iter =0; iter<numberOfGaussPoints;++iter){
		J +=detFs[iter]*gaussWeights[iter];
	}
	double currentVolume =J * ReferenceShape->Volume;
	return currentVolume;
}

gsl_matrix* ShapeBase::getCurrentFe(){
    /**
     * The current elastic deformation gradient is calculated the average at all Gauss points
     * (ShapeBase#numberOfGaussPoints with weights in ShapeBase#gaussWeights), as recorded in
     * ShapeBase#FeMatrices is calculated.
     */
	gsl_matrix* FeAvr = gsl_matrix_calloc(3,3);
	gsl_matrix* Fetmp = gsl_matrix_calloc(3,3);

    for (size_t iter =0; iter<numberOfGaussPoints;++iter){
		createMatrixCopy(Fetmp, FeMatrices[iter]);
		gsl_matrix_scale(Fetmp,gaussWeights[iter]);
		gsl_matrix_add(FeAvr, Fetmp);
	}
	gsl_matrix_free(Fetmp);
	return FeAvr;
}

void ShapeBase::relaxElasticForces(){
    /**
     * The current elastic deformation gradient coped over growth deformation gradient
     * to relax all elastic forces.
     */
	gsl_matrix* tmpFgForInversion =gsl_matrix_calloc(nDim,nDim);
	createMatrixCopy(Fg,TriPointF);
	createMatrixCopy(tmpFgForInversion,Fg);
	bool inverted = InvertMatrix(tmpFgForInversion, InvFg);
	if (!inverted){
        std::cerr<<"Fg not inverted!!"<<std::endl;
	}
	double detFg = determinant3by3Matrix(Fg);
	GrownVolume = detFg*ReferenceShape->Volume;
	VolumePerNode = GrownVolume/nNodes;
	gsl_matrix_free(tmpFgForInversion);
}

void ShapeBase::calculateFgFromRates(double dt, double x, double y, double z, gsl_matrix* rotMat, gsl_matrix* increment, int sourceTissue, double zMin, double zMax){
    /**
     * The current growth deformation gradient increment is calculated from rates and rotation.\n
     * First check point is for cheking if the growth function is applicable to the element, via
     * ShapeBase#isGrowthRateApplicable. If the element is growing, the diagonal of the growth increment
     * is obtained. Then the increment is rotated with a tensor rotation (M' = R M R^T), using the input rotation matrix.
     */
    double tissueWeight;
	bool continueCalaculation = isGrowthRateApplicable(sourceTissue, tissueWeight, zMin, zMax);
	if (continueCalaculation){
        scaleGrowthForZRedistribution(x,y,z);
		double gx = exp(x*tissueWeight*dt);
		double gy = exp(y*tissueWeight*dt);
		double gz = exp(z*tissueWeight*dt);
		gsl_matrix_set(increment,0,0,gx);
		gsl_matrix_set(increment,1,1,gy);
		gsl_matrix_set(increment,2,2,gz);
		gsl_matrix* temp = gsl_matrix_calloc(3,3);
		//R * increment
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, rotMat, increment, 0.0, temp);
		//increment * R^T
		gsl_blas_dgemm (CblasNoTrans, CblasTrans,1.0, temp, rotMat, 0.0, increment);
		gsl_matrix_free(temp);
	}
	else{
		gsl_matrix_set_identity(increment);
	}
}

void ShapeBase::calculateShapeChangeIncrementFromRates(double dt, double rx, double ry, double rz, gsl_matrix* increment){
    gsl_matrix_set(increment,0,0,exp(rx*dt));
	gsl_matrix_set(increment,1,1,exp(ry*dt));
	gsl_matrix_set(increment,2,2,exp(rz*dt));
}

void ShapeBase::calculateFgFromGridCorners(int gridGrowthsInterpolationType, double dt, GrowthFunctionBase* currGF, gsl_matrix* increment, int sourceTissue,  int IndexX, int IndexY, double FracX, double FracY){
    /**
     * The current growth deformation gradient increment is calculated from the input growth magnitude and rotations map in a grid \n
     * First check point is for cheking if the growth function is applicable to the element, via
     * ShapeBase#isGrowthRateApplicable. If the element is growing, the four corners of the grid that are closest
     * to the relative position of the element are extracted from the input grid. The growth rates are obtained at four
     * the corners via the function GrowthBase#getGrowthProfileAt4Corners.
     * Depending on th einterpolation method chosen, the growth of the point can be mapped to the closest corner,
     * or can be interpolated between the points depending on the distance the relative position of the element
     * centre falls within the grid: \n
     * \verbatim
             [point 2] ------------- [point 3]
                |                        |
                |<--fracX----> (o)       |
                |               |        |
                |               |        |
                |             fracY      |
                |               |        |
             [point 0] ------------- [point 1]
      \endverbatim
     */
    double tissueWeight;
	bool continueCalaculation = isGrowthRateApplicable(sourceTissue,tissueWeight,currGF->zMin,currGF->zMax);
	if (continueCalaculation){
		//taking growth data around 4 grid points
		//
		// Grid shape:
		//    [point 2] ------------- [point 3]
		//       |                        |
		//       |<--fracX----> (o)       |
		//       |               |        |
		//       |               |        |
		//       |             fracY      |
		//       |               |        |
		//    [point 0] ------------- [point 1]
		//
        std::array<double,3> growth0{0.0};
        std::array<double,3> growth1{0.0};
        std::array<double,3> growth2{0.0};
        std::array<double,3> growth3{0.0};
        std::array<double,4> angles{0.0};
        std::array<bool,4 >  angleEliminated{false};
		currGF->getGrowthProfileAt4Corners(IndexX, IndexY, growth0, growth1, growth2, growth3, angles, angleEliminated);
        std::array<double,3> growth{0.0};
		double angle = 0;
		if (gridGrowthsInterpolationType == 0){
			//using the growth rate at the grid point:
			//if fraction is below 0,5, I will use the index availabe.
			//If it is above 0.5, it is within the range of next groid point, I will use index +1:
			if(FracX > 0.5){
				if (FracY > 0.5){
					for (int i=0; i<3; ++i){
						growth[i] = growth3[i];
					}
				}
				else{
					for (int i=0; i<3; ++i){
						growth[i] = growth1[i];
					}
				}
			}
			else{
				if (FracY > 0.5){
					for (int i=0; i<3; ++i){
						growth[i] = growth2[i];
					}
				}
				else{
					for (int i=0; i<3; ++i){
						growth[i] = growth0[i];
					}
				}
			}
		}
		else if (gridGrowthsInterpolationType == 1){
            /**  In the interpolation schenario, the angles are checked and for the corners that have an aspect ratio over
             * the threshold the angle is included in the averaging.
             */
			//calculating the angle fraction eliminated, if any:
			double FracEliminated = 0.0;
			if (angleEliminated[0]){ FracEliminated += (1.0-FracX)*(1.0-FracY);	}
			if (angleEliminated[1]){ FracEliminated += FracX*(1.0-FracY);		}
			if (angleEliminated[2]){ FracEliminated += (1.0-FracX)*FracY;		}
			if (angleEliminated[3]){ FracEliminated += FracX*FracY;				}
			//taking the linear interpolation of 4 angles at 4 grid points:
			angle = angles[0]*(1.0-FracX)*(1.0-FracY)+angles[1]*FracX*(1.0-FracY)+angles[2]*(1.0-FracX)*FracY+angles[3]*FracX*FracY;
			if (FracEliminated>0){
				if (FracEliminated >= 0.9999999){
					angle = 0.0; //if all the angles should be eliminated because all corners have low aspect ratio, then angle is arbitrary, selected as zero
				}
				else{
					angle /= (1.0-FracEliminated); //normalising the sum to the eliminated averaging
				}
			}
            /** A linera interpolation with the distances are carried out for teh growth rates and orientation angles:
             */
			//taking the linear interpolation of 4 growth rates at 4 grid points
			for (int axis =0; axis<3; axis++){
				growth[axis]  = growth0[axis]*(1.0-FracX)*(1.0-FracY)+growth1[axis]*FracX*(1.0-FracY)+growth2[axis]*(1.0-FracX)*FracY+growth3[axis]*FracX*FracY;
				growth[axis] *= tissueWeight;
			}
		}
		//write the increment from obtained growth:
		if (isActinMimicing){
            /** If the element is ShapeBase#isActinMimicing, then growth in the tissue height is ignored.
             */
			growth[2] = 0; //no z-growth in actin mimicing apical surfaces.
		}
        /** Then growth is scaled if there is any distribution of tissue volume in z axis, thourgh ShapeBase#scaleGrowthForZRedistribution.
         */
        scaleGrowthForZRedistribution(growth[0],growth[1],growth[2]);
        for (size_t axis =0; axis<3; axis++){
			double gAxis = exp(growth[axis]*dt);
			gsl_matrix_set(increment,axis,axis,gAxis);
		}
		//Rotate the growth if the angel is not zero:
		if (angle != 0.0){
			gsl_matrix* rotMat  = gsl_matrix_calloc(3,3);
			double c = cos(angle);
			double s = sin(angle);
			gsl_matrix_set(rotMat,0,0,  c );
			gsl_matrix_set(rotMat,0,1, -1.0*s);
			gsl_matrix_set(rotMat,0,2,  0.0);
			gsl_matrix_set(rotMat,1,0,  s);
			gsl_matrix_set(rotMat,1,1,  c);
			gsl_matrix_set(rotMat,1,2,  0.0);
			gsl_matrix_set(rotMat,2,0,  0.0);
			gsl_matrix_set(rotMat,2,1,  0.0);
			gsl_matrix_set(rotMat,2,2,  1.0);
			gsl_matrix* temp = gsl_matrix_calloc(nDim,nDim);
			//R * increment
			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, rotMat, increment, 0.0, temp);
			//increment * R^T
			gsl_blas_dgemm (CblasNoTrans, CblasTrans,1.0, temp, rotMat, 0.0, increment);
			gsl_matrix_free(temp);
			gsl_matrix_free(rotMat);
		}
	}
	else{
		gsl_matrix_set_identity(increment);
	}
}

gsl_matrix* ShapeBase::getGrowthIncrement(){
	return growthIncrement;
}

void ShapeBase::updateGrowthIncrement(gsl_matrix* columnar, gsl_matrix* peripodial ){
    /** This function will update the elemental growth deformation gradient from the current
     * growth deformation gradient increment. The inputs to the function provides growht for two dostonct tissue types,
     * and the growth is distributed according to the ShapeBase#tissueType.
     */
	gsl_matrix* temp = gsl_matrix_calloc(nDim,nDim);
	if (tissueType == 0){//columnar layer element, no peripodial application necessary
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, columnar, growthIncrement, 0.0, temp);
		createMatrixCopy(growthIncrement, temp);
	}
	else if (tissueType == 1){//peripodial layer element, no columnar application necessary
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, peripodial, growthIncrement, 0.0, temp);
		createMatrixCopy(growthIncrement, temp);
	}
	else if (tissueType == 2){//linker between columnar and peripodial layer element, the growths are already weighted, need to apply both
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, columnar, growthIncrement, 0.0, temp);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, peripodial, temp, 0.0, growthIncrement);
	}
	for (int i=0; i<3; ++i){
		//this is used for display purposes of the simulation. As each new value is added to growthIncrement, I an update this directly, and the result is already cumulative of multiple growth functions
		GrowthRate[i] = gsl_matrix_get(growthIncrement,i,i);
	}
	gsl_matrix_free(temp);
}

void ShapeBase::updateShapeChangeIncrement(gsl_matrix* columnarShapeChangeIncrement){
	gsl_matrix* temp = gsl_matrix_calloc(nDim,nDim);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, columnarShapeChangeIncrement, shapeChangeIncrement, 0.0, temp);
	gsl_matrix_memcpy(shapeChangeIncrement, temp);
}

void ShapeBase::setRelativePosInBoundingBox(double x, double y){
	relativePosInBoundingBox[0] = x;
	relativePosInBoundingBox[1] = y;
}

std::array<double,2> ShapeBase::getRelativePosInBoundingBox(){
     std::array<double,2>  relativePos {relativePosInBoundingBox[0],relativePosInBoundingBox[1]};
     return relativePos;
}

std::array<double,2>  ShapeBase::getInitialRelativePosInBoundingBox(){
    std::array<double,2>  relativePos {initialRelativePosInBoundingBox[0], initialRelativePosInBoundingBox[1]};
    return relativePos;
}

void ShapeBase::setInitialRelativePosInBoundingBox(){
	initialRelativePosInBoundingBox[0] = relativePosInBoundingBox[0];
	initialRelativePosInBoundingBox[1] = relativePosInBoundingBox[1];
}

void ShapeBase::setInitialZPosition(double zMin, double TissueHeight){
    std::array<double,3> c = getCentre();
	initialRelativePositionInZ = 1.0 - ( (c[2] - zMin)/TissueHeight );
	//I cannot work with maximum tissue thickness, as there may be peripodial membrane on top
	//the tissue height is applicable to columnar layer only. I need to calculate from the bottom, then convert
	//such that apical surface will have a value of 0, and basal surface will have 1.0;
	//Apical surface is 0, basal surface is 1;
	if (initialRelativePositionInZ < 0){
		initialRelativePositionInZ = 0;
	}
	if (initialRelativePositionInZ > 1.0){
		initialRelativePositionInZ = 1.0;
    }
}



void ShapeBase::convertRelativePosToGridIndex(std::array<double,2> relpos, int& indexX, int &indexY, double &fracX, double &fracY, int nGridX, int nGridY){
	relpos[0] *= (float) (nGridX-1);
	relpos[1] *= (float) (nGridY-1);
	indexX = floor(relpos[0]);
	fracX  = relpos[0] - indexX;
	indexY = floor(relpos[1]);
	fracY  = relpos[1] - indexY;
	if (indexX >= nGridX-1) { //this is for the point that is exactly the point determining the bounding box high end in X, or the side elements for columnar parameter generation (outside the bounding box by definition
		indexX = nGridX-2;
		fracX = 1.0;
	}else if (indexX<0){
		indexX = 0;
		fracX = 0.0;
	}
	if (indexY >= nGridY-1) {//this is for the point that is exactly the point determining the bounding box high end in X, or the side elements for columnar parameter generation (outside the bounding box by definition
		indexY = nGridY-2;
		fracY = 1.0;
	}else if (indexY<0){
		indexY = 0;
		fracY = 0.0;
	}
}

void 	ShapeBase::readNodeIds(const std::vector<int>& inpNodeIds){
	/**
	 *  This function will take the input of an int pointer that contains the Ids of the nodes that
	 *  construct the element, and write them into the int array ShapeBase#NodeIds of size ShapeBase#nNodes.
	 *
	 */
    for (auto idForNode : inpNodeIds){
        this->NodeIds.push_back(idForNode);
	}
}

void 	ShapeBase::displayName(){
	/**
	 *  This function will write the shape type and shape id to standard output.
	 */
    std::cout<<"Type: "<<this->ShapeType<<" Id: "<<this->Id<<std::endl;
}

void 	ShapeBase::setPositionMatrix(const std::vector<std::unique_ptr<Node>>& Nodes){
	/**
	 *  This function will take the address of a Nodes pointers vector (It should be the Simulation#Nodes vector)
	 *  It will go through the ShapeBase#NodeIds in a nested loop of ShapeBase#nNodes and ShapeBase#nDim,
	 *  read the positions of the corresponding node from the input vector storing the pointers to all the nodes,
	 *  and write the position information into the ShapeBase#Positions array. This is a double storing of the
	 *  same information, yet is practical for elasticity calculations.
	 *
	 */
    for (size_t i = 0; i<nNodes; ++i){
        std::array<double,3> currentNodePos = {0.0};
        for (size_t j = 0; j<nDim; ++j){
            currentNodePos[j] = Nodes[NodeIds[i]]->Position[j];
		}
        Positions.push_back(currentNodePos);
	}
}

void 	ShapeBase::setTissuePlacement(const std::vector<std::unique_ptr<Node>>& Nodes){
	/**
	 *  This function will take the address of a Nodes pointers vector (It should be the Simulation#Nodes vector)
	 *  Checking the nodes that construct the element, the function will decide the placement of the element
	 *  within the tissue. The placement of the element within the tissue is defined by the ShapeBase#tissuePlacement
	 *  integer such that:
	 *  	- Apical:  ShapeBase#tissuePlacement = 1
	 *  	- Basal:   ShapeBase#tissuePlacement = 0
	 *  	- Mid-line or spans the whole tissue: ShapeBase#tissuePlacement = 2
	 *  	- Lateral: ShapeBase#tissuePlacement = 3
	 *
	 *  The function will first decide which type of nodes the element is composed of.
	 *
	 */
	bool hasApicalNode = false;
	bool hasBasalNode = false;
	bool hasLateralNode = false;
	spansWholeTissue = false;
    for (auto currNodeId : NodeIds){
        if (Nodes[currNodeId]->tissuePlacement == 1){
			hasApicalNode = true;
		}
        else if (Nodes[currNodeId]->tissuePlacement == 0){
			hasBasalNode = true;
		}
        else if (Nodes[currNodeId]->tissuePlacement == 3){
			hasLateralNode = true;
		}
	}
	/**
	*  Lateral elements can be composed of any combination of lateral, mid-line, apical and basal nodes.
	*  Any element that contains at least one lateral node, must be a lateral element.
	*
	*/
	if (hasLateralNode){
		tissuePlacement = 3;
	}
	else{
		/**
		* Elements that does not contain any lateral nodes, but apical nodes must be apical.
		* The only exception to this when the whole tissue is spanned by a single layer of
		* elements. Then the element will have apical and basal nodes, and it will be defined as
		* a mid-line node, with ShapeBase#spansWholeTissue set to true.
		*
		*/
		if (hasApicalNode){
			if (hasBasalNode){
				//the element spans through the whole tissue, the mid-line value should be used
				tissuePlacement = 2;
				spansWholeTissue = true;
			}
			else{
				//the element has only apical and midline nodes, it is apical
				tissuePlacement = 1;
			}
		}
		else if (hasBasalNode){
			/**
			* Elements that does not contain any lateral or apical nodes, but
			* does contain basal nodes then must be basal.
			*
			*/
			//the element only has basal and mid-line nodes, it is basal
			tissuePlacement = 0;
		}
		else{
			/**
			* If the element does not contain any lateral, apical or basal nodes, it must be that the
			* element is composed wholly of mid-line nodes. Then, the element lies in the mid-layer of
			* tissue.
			*
			*/
			//the element has only mid-line nodes, it is mid-line
			tissuePlacement = 2;
		}
	}
}

void ShapeBase::setECMMimicing(bool IsECMMimicing){
    /** While setting the element to ECM mimicking, the Poisson's ratio is also set to zero.
     * The physical definition of the ECM material requires so.
     */
    this->isECMMimicing = IsECMMimicing;

	//setting poisson ratio to be zero, so the ECM elements will not be thinning.
	this->v = 0;
}

void ShapeBase::setActinMimicing(bool isActinMimicing){
	this->isActinMimicing = isActinMimicing;
}

void 	ShapeBase::setTissueType(const std::vector<std::unique_ptr<Node> > &Nodes){
    /**
     * The tissue can be of columnar epithelum type, peripodial type, or be connecting the two.
     * Depending on the owner nodes of the element, if the element has any linker type nodes, it will be a linker node.
     * If there are no linker elements, but there are peripodial elements, then assign peripodial, as some peripodial
     * elements can contain columnar nodes wehen the mesh does not have linkers. Only when all hte nodes are of columnar type, the element
     * is columnar.
     */
	bool hasColumnarNode = false;
	bool hasPeripodialNode = false;
	bool hasLinkerNode = false;
    for (auto currNodeId : NodeIds){
        if (Nodes[currNodeId]->tissueType == 0){
			hasColumnarNode = true;
		}
        else if (Nodes[currNodeId]->tissueType == 1){
			hasPeripodialNode = true;
		}
        else if (Nodes[currNodeId]->tissueType == 2){
			hasLinkerNode = true;
		}
	}
	if (hasLinkerNode){
		tissueType = 2;
	}
	else if (hasPeripodialNode){
		//ASK LINKER ZONE BEFORE THIS, SOME LINKER ELEMENTS CAN HAVE LINKER NODES AND OTHER TISSUE NODES, NO COLUMNAR ELEMENT OR PERIPODIAL ELEMENT SHOULD HAVE A LINKER NODE
		tissueType = 1;
		setGrowthWeightsViaTissuePlacement( 1.0);//default is set to be columnar, I will not set this for linkers, as they are set in the initiation of peripodial membrane
	}
	else if (hasColumnarNode){
		//ASK PERIPODIAL MEMBRANE BEFORE THIS, SOME PERIPODIAL ELEMENTS CAN HAVE COLUMNAR NODES, AND SOME LINKER ELEMENTS CAN HAVE COLUMNAR NODES. NO COLUMNAR ELEMENT SHOULD HAVE A PERIPODIAL NODE
		tissueType = 0;
	}
	else {
        std::cerr<<"Element is not placed into tissue correctly, Id: "<<Id<<std::endl;
	}
}

void 	ShapeBase::setGrowthWeightsViaTissuePlacement (double periWeight){
	peripodialGrowthWeight = periWeight;
	columnarGrowthWeight = 1.0 - peripodialGrowthWeight;
}

void 	ShapeBase::setReferencePositionMatrix(){
	/**
	*  This function will allocate the position array of the ReferenceShape (ReferenceShape#Positions)
	*  Then the reference will be equated to the current position of the element.
	*  It is essential this function is called at simulation initiation, and after any subsequent modifications
	*  made to the reference structure of the tissue. But it should not be called in cases where
	*  the current shape of the element has progressed to differ from its reference (such as after loading
	*  steps of a save file). In saves, the reference shapes are set at the initiation of the system,
	*  and the positions of saves steps are loaded afterwards.
	*
	*/
    for (size_t i = 0; i<nNodes; ++i){
        ReferenceShape -> Positions.push_back (Positions[i]);
	}
}

void ShapeBase::setFg(gsl_matrix* currFg){
    gsl_matrix_memcpy (Fg, currFg);
    gsl_matrix* tmpFgForInversion =gsl_matrix_calloc(nDim,nDim);
    createMatrixCopy(tmpFgForInversion, Fg);
    bool inverted = InvertMatrix(tmpFgForInversion, InvFg);
    if (!inverted){
        std::cerr<<"Fg cannot be inverted!"<<std::endl;
    }
    gsl_matrix_free(tmpFgForInversion);
}

void ShapeBase::setYoungsModulus(double E){
	this -> E = E;
}

void ShapeBase::setViscosity(double viscosity){
	this -> internalViscosity = viscosity;
}

void ShapeBase::setViscosity(double viscosityApical,double viscosityBasal, double viscosityMid){
    /** Start by setting the viscosity to the mid layer value. If the element has a specific
     * node placement (ShapeBase#tissuePlacement is apical or basal) then alter accordingly.
     */
	this -> internalViscosity = viscosityMid;
	if (tissuePlacement == 0 ){
		this -> internalViscosity = viscosityBasal;
	}
	else if(tissuePlacement == 1 ){
		this -> internalViscosity = viscosityApical;
	}
	this -> originalInternalViscosity = internalViscosity;
}

void ShapeBase::setViscosity(double viscosityApical,double viscosityBasal){
    /** Start by setting the viscosity to the mid layer value to the average of apical and basal values.
     * If the element has a specific
     * node placement (ShapeBase#tissuePlacement is apical or basal) then alter accordingly.
     */
	this -> internalViscosity = 0.5*(viscosityApical+viscosityBasal);
	if (tissuePlacement == 0 ){
		this -> internalViscosity = viscosityBasal;
	}
	else if(tissuePlacement == 1 ){
		this -> internalViscosity = viscosityApical;
	}
}

void 	ShapeBase::updateShapeFromSave(std::ifstream& file){
	file >> IsAblated;
	updateNodeIdsFromSave(file);
	updateReferencePositionMatrixFromSave(file);
}

void 	ShapeBase::updateNodeIdsFromSave(std::ifstream& file){
    for (size_t i = 0; i<nNodes; ++i){
		int savedId;
		file >> savedId;
		NodeIds[i] = savedId;
	}
}

bool 	ShapeBase::readNodeIdData(std::ifstream& file){
    for (size_t i = 0; i<nNodes; ++i){
		int savedId;
		file >> savedId;

		if (NodeIds[i] != savedId){
            std::cout<<"NodeId "<<NodeIds[i]<<" savedId "<<savedId<<"do not match"<<std::endl;
			return false;
		}
	}
	return true;
}

void 	ShapeBase::updateReferencePositionMatrixFromSave(std::ifstream& file){
    for (size_t i = 0; i<nNodes; ++i){
        for (size_t j = 0; j<nDim; ++j){
			double savedPos;
			file >> savedPos;
			ReferenceShape -> Positions[i][j] = savedPos;
		}
	}
}

void 	ShapeBase::updateReferencePositionMatrixFromInput(double** inputPos){
    for (size_t i = 0; i<nNodes; ++i){
        for (size_t j = 0; j<nDim; ++j){
			ReferenceShape -> Positions[i][j] = inputPos[i][j];
		}
	}
}

bool	ShapeBase::readReferencePositionData(std::ifstream& file){
    for (size_t i = 0; i<nNodes; ++i){
        for (size_t j = 0; j<nDim; ++j){
			double savedPos;
			file >> savedPos;
			if (ReferenceShape -> Positions[i][j] != savedPos){
				//the positions are not equal, it may be an issue of rounding, my satisfactory precision is 2%
				float percentError = (ReferenceShape -> Positions[i][j] - savedPos) / ReferenceShape -> Positions[i][j]*100.0;
				if (percentError>2.0 || percentError< -2.0){
                    std::cerr<<"ReferenceShape->Positions: "<<ReferenceShape -> Positions[i][j]<<" savedPos "<<savedPos<<" percent Error: "<<percentError<<std::endl;
					return false;
				}
			}
		}
	}
	return true;
}

void 	ShapeBase::updateReferencePositionMatrixFromMeshInput(std::ifstream& file){
	updateReferencePositionMatrixFromSave(file);
}

void ShapeBase::updateElementVolumesAndTissuePlacementsForSave(const std::vector<std::unique_ptr<Node>>& Nodes){
	calculateReferenceVolume();
	setTissuePlacement(Nodes);
	setTissueType(Nodes);
}

void 	ShapeBase::displayNodeIds(){
    for (auto currNodeId : NodeIds){
        std::cout<<currNodeId<<"  ";
        std::cout<<std::endl;
	}
}

void 	ShapeBase::displayPositions(){
    for (const auto & currNode : Positions){
        for (const auto& currDimension: currNode){
            std::cout<<currDimension<<"  ";
		}
        std::cout<<std::endl;
	}
}

void 	ShapeBase::displayReferencePositions(){
    for (const auto & currNode : ReferenceShape->Positions){
        for (const auto& currDimension: currNode){
            std::cout<<currDimension<<"  ";
        }
        std::cout<<std::endl;
    }
}

std::array<int,3>	ShapeBase::getIdentifierColour(){
	return IdentifierColour;
}

void 	ShapeBase::getStrain(int type, float &StrainMag){
	StrainMag = 0.0;
    if (type == 0){
		//DV
        StrainMag = gsl_matrix_get(Strain,0,0);
	}
    else if (type == 1){
		//AP
        StrainMag = gsl_matrix_get(Strain,1,0);
	}
    else if (type == 2){
		//AB
        StrainMag = gsl_matrix_get(Strain,2,0);
	}
    else if (type == 3){
		//xy
        StrainMag = gsl_matrix_get(Strain,3,0);
	}
    else if (type == 4){
		//yz
        StrainMag = gsl_matrix_get(Strain,4,0);
	}
    else if (type == 5){
		//xz
        StrainMag = gsl_matrix_get(Strain,5,0);
	}
	else{
		return;
	}
}

void 	ShapeBase::getNodeBasedPysProp(int type, int NodeNo, const std::vector<std::unique_ptr<Node>>& Nodes, float& PysPropMag){
	PysPropMag = 0.0;
	if (type == 0){
		PysPropMag = Nodes[NodeIds[NodeNo]] -> externalViscosity[2];
	}
}

void 	ShapeBase::getPysProp(int type, float &PysPropMag, double dt){
	if (type == 1){
		PysPropMag = getInternalViscosity();
	}
	else if (type == 2){
		PysPropMag = getYoungModulus();
	}
	else if (type == 3 ){
		PysPropMag = getPoissonRatio();
	}
	else if (type == 4){
        const std::array<double,3>& growth = getGrowthRate();
        double timescale = 24*60.0*60.0; //reporting per 24 hours
        for (size_t i =0; i< nDim ; ++i){//reporting x,y,z
			//growth is in form exp(r*dt), get r first, then adjust the time scale, and report the exponential form still:
			//And I want to rate of volume growth, that is x*y*z
			double value = exp(log(growth[i])/dt*timescale);
			PysPropMag *= value;
		}
	}
	else if (type == 5){
		PysPropMag = GrownVolume/ReferenceShape->Volume;
	}
}

double	ShapeBase::calculateEmergentShapeOrientation(){
    /** The function calculates in which direction the emergent shape is oriented.
     * We need to have the combination of growth gradient and deformation gradient.
     * It will reflect how the clones would "look":
     */
	double currEmergentVolume = calculateCurrentGrownAndEmergentVolumes();
    gsl_matrix* E;
    gsl_matrix* C = calculateCauchyGreenDeformationTensor(TriPointF);
	E = calculateEForNodalForcesKirshoff(C);
	gsl_matrix* strainBackUp = gsl_matrix_calloc(6,1);
    /** We keep a copy of the original ShapeBase#Strain to utilise the Strain matrix in eigen value
     * decomposition.
     */
	for (int i=0;i<5;++i){
		gsl_matrix_set(strainBackUp,i,0,gsl_matrix_get(Strain,i,0));
	}
	gsl_matrix_set_zero(Strain);
	gsl_matrix_set(Strain,0,0, gsl_matrix_get(E,0,0));
	gsl_matrix_set(Strain,1,0, gsl_matrix_get(E,1,1));
	gsl_matrix_set(Strain,2,0, gsl_matrix_get(E,2,2));
	gsl_matrix_set(Strain,3,0, 2.0*gsl_matrix_get(E,0,1));
	gsl_matrix_set(Strain,4,0, 2.0*gsl_matrix_get(E,2,1));
	gsl_matrix_set(Strain,5,0, 2.0*gsl_matrix_get(E,0,2));
	double e1 = 0.0, e2 = 0.0, e3 = 0.0;
	gsl_matrix* eigenVec = gsl_matrix_calloc(3,3);
	calculatePrincipalStrains2D(e1,e2,e3,eigenVec);
    /** Calculated by the funciton ShapeBase#calculatePrincipalStrains2D, the Eigen vector matrix is
     * a 3 by 3 matrix, but only stores the 2D vectors in its upper corner 2x2 terms. The vectors are
     * written in columns of the matrix. The strain I have here is Green strain, I would like to convert
     * it back to deformation gradient terms. Since  \f$ E = 1/2 *(Fe^T*Fe-I): \f$
     * \f$  F_{ii} = \sqrt{e_{i} \times 2 + 1}
     * \f$
    */
	double F11 = pow(e1*2+1,0.5);
	double F22 = pow(e2*2+1,0.5);
	double AR = F11/F22;
    /** Then with the aspect ratio, the emergent long and short axes are calculated.
     */
	if (AR < 1.0){
		AR = 1.0/AR;
		emergentShapeLongAxis[0] = AR * gsl_matrix_get(eigenVec,0,1);
		emergentShapeLongAxis[1] = AR * gsl_matrix_get(eigenVec,1,1);
		emergentShapeShortAxis[0] = gsl_matrix_get(eigenVec,0,0);
		emergentShapeShortAxis[1] = gsl_matrix_get(eigenVec,1,0);
	}
	else{
		emergentShapeLongAxis[0] = AR * gsl_matrix_get(eigenVec,0,0);
		emergentShapeLongAxis[1] = AR * gsl_matrix_get(eigenVec,1,0);
		emergentShapeShortAxis[0] = gsl_matrix_get(eigenVec,0,1);
		emergentShapeShortAxis[1] = gsl_matrix_get(eigenVec,1,1);
	}
    /** Finally, the original ShapeBase#Strain are copied over back to the Strains matrix
     */
	//copy the real strains back on the strains matrix
	for (int i=0;i<5;++i){
		gsl_matrix_set(Strain,i,0,gsl_matrix_get(strainBackUp,i,0));
	}
	gsl_matrix_free(C);
	gsl_matrix_free(E);
	gsl_matrix_free(strainBackUp);
	gsl_matrix_free(eigenVec);
    /** The function returns the volumentric strain (ratio of emergent volume to reference shape volume)
     * as autput.
     */
	return currEmergentVolume/ReferenceShape->Volume;
}

void 	ShapeBase::displayIdentifierColour(){
    std::cout <<" IdentifierColour:  "<<IdentifierColour[0]<<" "<<IdentifierColour[1]<<" "<<IdentifierColour[2]<<std::endl;
}

void 	ShapeBase::changeShapeByFsc(double dt){
    gsl_matrix* FscIncrement = gsl_matrix_calloc(nDim,nDim); ///< The increment of shape change that will be induced this step
    if (rotatedGrowth){
    	double rTemp[3] = {0.0,0.0,0.0};
		for (int i = 0; i<3; ++i){
			rTemp[i] = gsl_matrix_get(GrowthStrainsRotMat,i,0)*ShapeChangeRate[0]+gsl_matrix_get(GrowthStrainsRotMat,i,1)*ShapeChangeRate[1]+gsl_matrix_get(GrowthStrainsRotMat,i,2)*ShapeChangeRate[2];
			if ( (ShapeChangeRate[i] <0 && rTemp[i] >0.0 ) || (ShapeChangeRate[i] >0 && rTemp[i] < 0.0) ){
				rTemp[i] *= -1.0;
			}
		}
		for (int i = 0; i<3; ++i){
			ShapeChangeRate[i] = rTemp[i];
		}
	}
    for (int i=0; i<3 ;++i){
    	gsl_matrix_set(FscIncrement,i,i, exp(ShapeChangeRate[i]*dt));
    }
    gsl_matrix* temp1 = gsl_matrix_calloc(nDim,nDim);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, FscIncrement, Fsc, 0.0, temp1);
	gsl_matrix_memcpy(Fsc, temp1);
	gsl_matrix* tmpFscForInversion = gsl_matrix_calloc(nDim,nDim);
	createMatrixCopy(tmpFscForInversion,Fsc);
	bool inverted = InvertMatrix(tmpFscForInversion, InvFsc);
	if (!inverted){
        std::cerr<<"Fsc not inverted!!"<<std::endl;
	}
	gsl_matrix_free(FscIncrement);
	gsl_matrix_free(temp1);
	gsl_matrix_free(tmpFscForInversion);
}

void ShapeBase::setPlasticDeformationIncrement(double xx, double yy, double zz){
	gsl_matrix_set(plasticDeformationIncrement,0,0,xx);
	gsl_matrix_set(plasticDeformationIncrement,1,1,yy);
	gsl_matrix_set(plasticDeformationIncrement,2,2,zz);
}

void 	ShapeBase::growShapeByFg(){
     /** This function updates the current growth deformaiton gradient with the growt/shape
      * change/plastic deformation increments and their respective rotations.
      */
    if (rotatedGrowth){
        /** If the growth is rotated (ShapeBase#rotatedGrowth), the current growth increment
         *  in rotated with a tensor rotation: \f$ \mathbf{R}^{T} \mathbf{F}^{G}_{increment} \mathbf{R}  \f$.\n
         *  Where \f$ \mathbf{F}^{G}_{increment} \f$
         * is ShapeBase#growthIncrement and the rotation matrix is ShapeBase#GrowthStrainsRotMat. A similar
         * rotation is applied on shape change increment defined in ShapeBase#shapeChangeIncrement and
         * ShapeBase#GrowthStrainsRotMat.
         */
        gsl_matrix* temp = gsl_matrix_calloc(nDim,nDim);
        //R^T * growthIncrement
        gsl_blas_dgemm (CblasTrans, CblasNoTrans,1.0, GrowthStrainsRotMat, growthIncrement, 0.0, temp);
        //R^T * growthIncrement * R
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, temp, GrowthStrainsRotMat, 0.0, growthIncrement);
    	//rotate shape change increment:
    	gsl_blas_dgemm (CblasTrans, CblasNoTrans,1.0, GrowthStrainsRotMat, shapeChangeIncrement, 0.0, temp);
    	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, temp, GrowthStrainsRotMat, 0.0, shapeChangeIncrement);
    	gsl_matrix_free(temp);
    }
    //incrementing Fg with current growth rate, plastic deformation rate, and shape changes:
    gsl_matrix* temp1 = gsl_matrix_calloc(nDim,nDim);
    gsl_matrix* temp2 = gsl_matrix_calloc(nDim,nDim);
    gsl_matrix* temp3 = gsl_matrix_calloc(nDim,nDim);
    //adding plastic deformation, this increment is in already in correct orientation:
    /** The plastic deformation increment is already in the correct orientation by definiton. The increments
     * are then merged and added on the current growth deformation gradient \f$ F^{G} \f$
     */
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, plasticDeformationIncrement,growthIncrement, 0.0, temp1);
    //adding shape change, this increment is in already in correct orientation:
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, shapeChangeIncrement,temp1, 0.0, temp2);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, temp2, Fg, 0.0, temp3);
    gsl_matrix_memcpy(Fg, temp3);
    gsl_matrix* tmpFgForInversion = gsl_matrix_calloc(nDim,nDim);
    createMatrixCopy(tmpFgForInversion,Fg);
    bool inverted = InvertMatrix(tmpFgForInversion, InvFg);
    if (!inverted){
        std::cerr<<"Fg not inverted!!"<<std::endl;
    }
    /**
     * The volumentric change induced by growth is calculated via the determinant of the growth deformation
     * gradient. The current prefered volume is updated accordingly.
     */
    double detFg = determinant3by3Matrix(Fg);
    GrownVolume = detFg*ReferenceShape->Volume;
    VolumePerNode = GrownVolume/nNodes;
    //freeing matrices allocated in this function
    gsl_matrix_free(temp1);
    gsl_matrix_free(temp2);
    gsl_matrix_free(temp3);
    gsl_matrix_free(tmpFgForInversion);
}

void ShapeBase::displayDebuggingMatrices(){
    std::cout<<" Prism "<<Id<<" rotatedGrowth: "<<rotatedGrowth<<std::endl;
    displayMatrix(Fg, " Fg");
    displayMatrix(FeMatrices[0], " FeMatrices[0]");
    displayMatrix(FeMatrices[1], " FeMatrices[1]");
    displayMatrix(FeMatrices[2], " FeMatrices[2]");
    displayMatrix(GrowthStrainsRotMat, " GrowthStrainsRotMat");
}

double 	ShapeBase::calculateCurrentGrownAndEmergentVolumes(){
    /** Once the reference volume is calculated, the current prefered volume is obtained by
     * scaling the reference shape volume by the determinant of the growth deformation gradinet ShapeBase#Fg.
     * On the other hand, current emergent volume is the deformed volume and should be obtained by scaling
     * the reference volume with the determinant of the total deformation gradient, including both
     * growth and deformation.
     */
	calculateReferenceVolume();
	double detFg = determinant3by3Matrix(Fg);
    GrownVolume = detFg*ReferenceShape->Volume;
	calculateTriPointFForRatation();
	double detF =determinant3by3Matrix(TriPointF);
	double emergentVolume = detF*ReferenceShape->Volume;
    return emergentVolume;

}

bool ShapeBase::isActinStiffnessChangeAppliedToElement(bool ThereIsWholeTissueStiffnessPerturbation, bool ThereIsApicalStiffnessPerturbation, bool ThereIsBasalStiffnessPerturbation, bool ThereIsBasolateralWithApicalRelaxationStiffnessPerturbation, bool ThereIsBasolateralStiffnessPerturbation, std::vector <int> &stiffnessPerturbationEllipseBandIds, int numberOfStiffnessPerturbationAppliesEllipseBands ){
    /** The function regulating stiffness perturbations to tissue and ECM are distinct. Therefore first chack is tissue
     * specific characterisation, ShapeBase#isECMMimicing. If the the element is not, the positional coupling between
     * the input perturbation adn the element tissue properties are checked to see if the perturbation is applied to
     * this element.
      */
    if (!isECMMimicing){
		if( ThereIsWholeTissueStiffnessPerturbation  //the whole columnar tissue is perturbed
			|| ThereIsBasolateralWithApicalRelaxationStiffnessPerturbation	// there is relaxation on the apical surface and stiffenning on the rest of the tissue, further checks needed while calculating the rate
            || (ThereIsBasolateralStiffnessPerturbation && tissuePlacement != 1) //there is only basolateral stiffness perturbation, without affecting apical sides
			|| (tissuePlacement == 0 && ThereIsBasalStiffnessPerturbation    ) //the basal surface is perturbed and element is basal
			|| (tissuePlacement == 1 && ThereIsApicalStiffnessPerturbation   ) // the apical surface is perturbed and element is apical
			|| (atBasalBorderOfECM   && ThereIsBasalStiffnessPerturbation    ) //the basal surface is perturbed and element is at the layer above the "basal elements, but there is basal ECM, therefore the element is basal of the tissue (atBasalBorderOfECM)
			){
			if(insideEllipseBand){ //this covers the tissue type as well as the position, ellipses are applied to columnar and linker zones only.
				for (int stiffnessPerturbationRangeCounter =0; stiffnessPerturbationRangeCounter<numberOfStiffnessPerturbationAppliesEllipseBands; ++stiffnessPerturbationRangeCounter){
					if (coveringEllipseBandId == stiffnessPerturbationEllipseBandIds[stiffnessPerturbationRangeCounter]){
						return true;
					}
				}
			}
		}
	}
	return false;
}

bool ShapeBase::isShapeChangeAppliedToElement(std::vector<int> &ellipseBandIds, bool applyBasalECM, bool applyToLateralECM, bool applyApically, bool applyBasally, bool applyMidLayer ){
	bool checkForEllipseId = false;
	if (tissueType == 1){ //element is on the peripodial membrane, I am not applyign this to peripodial
		return false;
	}
	if (isECMMimicing){
		if  (   (applyBasalECM  && tissuePlacement == 0 )
			 || (applyToLateralECM && isECMMimimcingAtCircumference && !tissuePlacement == 0) //do not grow the basal element twice
			){
			checkForEllipseId = true;
		}
	}
	else{
		if  (   (applyBasally  && (tissuePlacement == 0 || atBasalBorderOfECM) )
				 || (applyApically && tissuePlacement == 1)
				 || (applyMidLayer && tissuePlacement == 2)
			){
				checkForEllipseId = true;
		}
	}
	int nEllipseBands = ellipseBandIds.size();
	//The element qualifies for shape change in this function, is it inside the ellipse band?
	if(checkForEllipseId && insideEllipseBand){ //this covers the tissue type as well as the position, ellipses are applied to columnar and linker zones only.
		for (int shapeChangeEllipseBandCounter =0; shapeChangeEllipseBandCounter<nEllipseBands; ++shapeChangeEllipseBandCounter){
			if (coveringEllipseBandId == ellipseBandIds[shapeChangeEllipseBandCounter]){
				return true;
			}
		}
	}
	return false;
}

bool ShapeBase::isECMChangeAppliedToElement(bool changeApicalECM, bool changeBasalECM, std::vector<int> &ECMChangeEllipseBandIds, int numberOfECMChangeEllipseBands){
	if (isECMMimicing){
		if  (    (changeApicalECM && tissuePlacement == 1 )
			  || (changeBasalECM  && tissuePlacement == 0 )
			  || (isECMMimimcingAtCircumference)
			  //|| (changeStiffnessBasalECM  && tissuePlacement == 2 )
			){
			if(insideEllipseBand){
				for (int ECMReductionRangeCounter = 0; ECMReductionRangeCounter<numberOfECMChangeEllipseBands; ++ECMReductionRangeCounter){
					if (coveringEllipseBandId == ECMChangeEllipseBandIds[ECMReductionRangeCounter]){
						return true;
					}
				}
			}

		}
	}
	return false;
}


void ShapeBase::calculateStiffnessPerturbationRate(bool ThereIsBasolateralWithApicalRelaxationStiffnessPerturbation, double stiffnessPerturbationBeginTimeInSec, double stiffnessPerturbationEndTimeInSec, double stiffnessChangedToFractionOfOriginal){
    /**
     * This function will calciulate the stiffness perturbation rate.
     */
    double totalTimePerturbationWillBeAppliedInSec = stiffnessPerturbationEndTimeInSec-stiffnessPerturbationBeginTimeInSec;
    if (totalTimePerturbationWillBeAppliedInSec <0){
        stiffnessPerturbationRateInSec = 0;
        return;
    }
    if (ThereIsBasolateralWithApicalRelaxationStiffnessPerturbation){
        /** the used rate will be different for apical elements and all the remaining elements.
         * I do not need to check for ECM, as this is called for only the elements that has
         * applied stiffness perturbations, which already excluded ECM elements.
         */
    	if(tissuePlacement == 1){ //element is apical.
            /** If the element is apical, whatever I am applying to the basal side, I will apply the
             * inverse to the apical side. If the baso-lateral side is doubling, apical surface will halve.
             * Please note this will not be feasible for elements that span the whole disc. You cannot
             * do a baso-lateral change for elements that cover the whole tissue!
             */
    		stiffnessChangedToFractionOfOriginal = 1.0/stiffnessChangedToFractionOfOriginal;
    	}
    }
    stiffnessPerturbationRateInSec =  (stiffnessChangedToFractionOfOriginal - 1.0)/totalTimePerturbationWillBeAppliedInSec;
    if (stiffnessPerturbationRateInSec<0){
    	minimumValueOfStiffnessMultiplier = stiffnessChangedToFractionOfOriginal;
    }
    else{
    	maximumValueOfStiffnessMultiplier = stiffnessChangedToFractionOfOriginal;
    }
}

void ShapeBase::updateStiffnessMultiplier(double dt){
	stiffnessMultiplier += stiffnessPerturbationRateInSec*dt;
	if (stiffnessMultiplier<minimumValueOfStiffnessMultiplier){
		stiffnessMultiplier = minimumValueOfStiffnessMultiplier;
	}
	if (stiffnessMultiplier>maximumValueOfStiffnessMultiplier){
		stiffnessMultiplier = maximumValueOfStiffnessMultiplier;
	}
}

void	ShapeBase::assignSoftHinge(double lowHingeLimit, double highHingeLimit,double softnessLevel){
	if (!isECMMimicing){
		if (relativePosInBoundingBox[0]>lowHingeLimit && relativePosInBoundingBox[0] < highHingeLimit){
			stiffnessMultiplier *= softnessLevel;
			updateElasticProperties();
		}
	}
}

bool	ShapeBase::checkZCappingInRemodelling(bool volumeConserved, double zRemodellingLowerThreshold, double zRemodellingUpperThreshold, gsl_matrix* increment, gsl_matrix* eigenVec){
    /** This function checks if the remodelling of the element in z axis have reached the specified cap.
     * First step is to figure out what the gorwth on z axis will be upon the applicatin of perturbation.
     */
	bool zCapped = false;
	gsl_matrix* tempForZCap = gsl_matrix_calloc(nDim,nDim);
	gsl_matrix* incrementToTestZCapping = gsl_matrix_calloc(nDim,nDim);
	//eigenVec * increment
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, eigenVec, increment, 0.0, tempForZCap);
	//eigenVec * increment * eigenVec^T
	gsl_blas_dgemm (CblasNoTrans, CblasTrans,1.0, tempForZCap, eigenVec, 0.0, incrementToTestZCapping);
    /** Once the gorwth increment in a non is non-volume conserved and non-scaled for z
     * approach, I will have the potential new z deformation. There are limitations to how much z axis
     * can be remodelled. I check those limits provided as inputs to the funciton and cap the z
     * deformation if necessary.
     */
	double Fzz = gsl_matrix_get(incrementToTestZCapping,2,2);
	if (volumeConserved){
		double det = determinant3by3Matrix(incrementToTestZCapping);
		double scale = 1.0/pow (det,1.0/3.0);
		double zIncrement = Fzz*scale;
		if (zRemodellingSoFar*zIncrement >= zRemodellingUpperThreshold){
			zCapped = true;
		}
		if (zRemodellingSoFar*zIncrement <= zRemodellingLowerThreshold){
			zCapped = true;
		}
	}
	else{
		if (zRemodellingSoFar*Fzz >= zRemodellingUpperThreshold){
			zCapped = true;
		}
		if (zRemodellingSoFar*Fzz <= zRemodellingLowerThreshold){
			zCapped = true;
		}
	}
	gsl_matrix_free(tempForZCap);
	gsl_matrix_free(incrementToTestZCapping);
	return zCapped;
}

void	ShapeBase::checkIfInsideEllipseBands(int nMarkerEllipseRanges, std::vector<double> markerEllipseBandXCentres, std::vector<double> markerEllipseBandR1Ranges, std::vector<double> markerEllipseBandR2Ranges, const std::vector<std::unique_ptr<Node> > &Nodes){
	for (int i=0;i<nMarkerEllipseRanges; ++i){	
		double dx  = relativePosInBoundingBox[0] - markerEllipseBandXCentres [i];
		double dy = relativePosInBoundingBox[1];
		if ( (markerEllipseBandR1Ranges[2*i]> 0 && dx <0) || (markerEllipseBandR1Ranges[2*i]< 0 && dx >0)){
			double dxOverR1 = dx/markerEllipseBandR1Ranges[2*i];
			double dyOverR2 = dy/markerEllipseBandR2Ranges[2*i];
			double d_squareLower = dxOverR1*dxOverR1 + dyOverR2*dyOverR2;
			dxOverR1 = dx/markerEllipseBandR1Ranges[2*i+1];
			dyOverR2 = dy/markerEllipseBandR2Ranges[2*i+1];
			double d_squareUpper = dxOverR1*dxOverR1 + dyOverR2*dyOverR2;
			//For element to be inside the band, calculated distance square values
			//should be larger than 1 for the smaller ellipse, and smaller than 1 for the
			//larger ellipse, hence in between two ellipses.
			if (d_squareLower> 1 && d_squareUpper<1){
				insideEllipseBand = true;
				coveringEllipseBandId = i;
                for (auto currNodeId : NodeIds){
					Nodes[currNodeId]->insideEllipseBand=true;
					Nodes[currNodeId]->coveringEllipseBandId  = i;
                    //std::cout<<"Node "<<currNodeId<<" is inside ellipse"<<i<<std::endl;
				}
			}
		}		
	}
}

void	ShapeBase::calculatePlasticDeformation3D(bool volumeConserved, double dt, double plasticDeformationHalfLife, double zRemodellingLowerThreshold, double zRemodellingUpperThreshold){
    /**
     * This function calculates the plastic deformation (remodelling) from the current elastic deformation
     * gradient. First check is against tissue specific compartments and z remodelling capping. The linker zone of the tissue and ECM mimicking
     * elements should be able to remodell in z regardless of currently active caps.
     */
    double e1 = 0.0, e2 = 0.0, e3 = 0.0;
	gsl_matrix* eigenVec = gsl_matrix_calloc(3,3);
	bool checkZCapping = true;
	if( isECMMimicing || tissueType == 2){
		checkZCapping = false;
	}
    /** The princiapl starins on the element are then calculated bie eigen vsalue decomposition
     * via ShapeBase#calculatePrincipalStrains3D.
     */
	calculatePrincipalStrains3D(e1,e2,e3,eigenVec);
    /** This gives the green strain in principal direction in the orientation of the element
     * internal coordinates. One can simply grow the element in this axis, to obtain some
     * form of plastic deformation/remodelling. Now the Green stains need to be converted
     * to deformation gradient terms. Since  \f$ E = 1/2 *(Fe^T*Fe-I): \f$
     * \f$  F_{ii} = \sqrt{e_{i} \times 2 + 1}
     * \f$
    */
	double F11 = pow(e1*2+1,0.5);
	double F22 = pow(e2*2+1,0.5);
	double F33 = pow(e3*2+1,0.5);
    /** Then remodelling with a decay half-life ShapeBase#plasticDeformationHalfLife
     * becomes: \n
     * \f$ N(t+\Delta t) = N(t) * 2 ^ (-\Delta t/\tau_{1/2})
     * \f$
     * where \f$ \tau_{1/2} = plasticDeformationHalfLifeMultiplier * plasticDeformationHalfLife/(log(2)) \f$.
     * Here the multiplier term reflects perturbations.
     */
	double tau = plasticDeformationHalfLifeMultiplier * plasticDeformationHalfLife/(log(2)); // (mean lifetime tau is half life / ln(2))
	double F11t = (F11-1)*exp(-1.0*dt/tau) + 1;
	double F22t = (F22-1)*exp(-1.0*dt/tau) + 1;
	double F33t = (F33-1)*exp(-1.0*dt/tau) + 1;
	F11 = F11/F11t;
	F22 = F22/F22t;
	F33 = F33/F33t;
    /** The obtained remodelling increment is then checked against capping in z remodelling via
     * ShapeBase#checkZCappingInRemodelling.
     */
	//writing onto an incremental matrix:
	gsl_matrix*  increment = gsl_matrix_calloc(3,3);
	gsl_matrix_set(increment,0,0,F11);
	gsl_matrix_set(increment,1,1,F22);
	gsl_matrix_set(increment,2,2,F33);
	bool zCapped = false;
	if (checkZCapping){
		zCapped = checkZCappingInRemodelling(volumeConserved, zRemodellingLowerThreshold, zRemodellingUpperThreshold, increment, eigenVec);
	}
	//If the element is ECM mimicking but not lateral, then I do not want any z remodelling:
	//I need to write the specific conditions that cover the setup where there is explicit ECM, but no
	//peripodial, therefore, the circumference elements should not be zCapped, they should be treated as lateral:
	if ((isECMMimicing) ){
		//first condition will exempt the lateral elements.
		//second condition will exempt circumferential elements that are assigned to be ECM
		if (tissueType != 2) {
			if (!isECMMimimcingAtCircumference){
				zCapped = true;
			}
		}
	}
    /** If there is z capping the procedure is repeated in 2D, with the strains calculated via
     * ShapeBase#calculatePrincipalStrains2D.
     */
	if (zCapped){
		calculatePrincipalStrains2D(e1,e2,e3,eigenVec);
		F11 = pow(e1*2+1,0.5);
		F22 = pow(e2*2+1,0.5);
		F11t = (F11-1)*exp(-1.0*dt/tau) + 1;
		F22t = (F22-1)*exp(-1.0*dt/tau) + 1;
		F33t = 1.0;
		F11 = F11/F11t;
		F22 = F22/F22t;
		F33 = 1.0;
		gsl_matrix_set(increment,0,0,F11);
		gsl_matrix_set(increment,1,1,F22);
		gsl_matrix_set(increment,2,2,F33);
	}
	//If I am conserving the volume, I need to scale:
    /** If hte volume is conserved, the incremeent is sclaed to have a determinant of unity.
     */
	if (volumeConserved){
		double det = determinant3by3Matrix(increment);
		if (zCapped){
			double scale = 1.0/pow (det,1.0/2.0); //scale the size with the square root of the determinant to keep the volume conserved. Then set z to 1 again, we do not want to affect z growth, remodelling is in x & y
			gsl_matrix_scale(increment,scale);
			//I know the eigen vector of the 2nd (in 3 dim indexed from 0) column is z axis, if z is capped
			//Then the incremental value at 2,2 position (in 3 dim indexed from 0) is z
			gsl_matrix_set(increment,2,2,1);
		}
		else{
			double scale = 1.0/pow (det,1.0/3.0); //scaling in x,y, and z
			gsl_matrix_scale(increment,scale);
		}
	}
	//The element is ECM mimicing, and it is lateral. As a temporary solution, I simply will not
	//allow any lateral elements to shrink as part of their remodelling.
	if(isECMMimicing && tissueType == 2){
		for (int i= 0; i<3; ++i){
			double value = gsl_matrix_get(increment,i,i);
			if (value < 1.0){
				gsl_matrix_set(increment,i,i,1.0);
			}
		}
	}
	//The growth I would like to apply now is written on the increment.
	//I would like to rotate the calculated incremental "growth" to be aligned with the coordinate
	//system defined by the eigen vectors.
	//In a growth setup, I calculate the rotationa matrix to be rotation by a certain angle.
	//Here, the rotation matrix is the eigen vector matrix itself. The eigen vector matrix
	//will rotate the identity matrix upon itself.
	gsl_matrix* temp = gsl_matrix_calloc(nDim,nDim);
	//eigenVec*increment
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, eigenVec, increment, 0.0, temp);
	//eigenVec*increment *eigenVec^T
	gsl_blas_dgemm (CblasNoTrans, CblasTrans,1.0, temp, eigenVec, 0.0, plasticDeformationIncrement);
	//update z remodelling so far to keep track and not update beyond limits
    /** The total z remodelling up to the current time step (ShapeBase#zRemodellingSoFar) is then updated
     * for z remodelling capping checks in following time points.
     */
	zRemodellingSoFar *= gsl_matrix_get(plasticDeformationIncrement,2,2);
	gsl_matrix_free(temp);
	gsl_matrix_free(eigenVec);
	gsl_matrix_free(increment);
}


void 	ShapeBase::CalculateGrowthRotationByF(){
    /** The rigid body rotation is extracted from the deformation gradient, and the
     * rotation around the z-axis is stored. The disassembling of the rigid body rotation is carried out via
     * ShapeBase#disassembleRotationMatrixForZ.
     */
    gsl_matrix* rotMat = gsl_matrix_alloc(3,3);
    gsl_matrix_set_identity(rotMat);
    //updating the F for the current shape positions
    //(not using leftovers from previous iteration)
    calculateTriPointFForRatation();
    rotatedGrowth = calculate3DRotMatFromF(rotMat);
    if (rotatedGrowth){
        rotatedGrowth = disassembleRotationMatrixForZ(rotMat);
        if (rotatedGrowth){
            gsl_matrix_transpose(rotMat);
            gsl_matrix_memcpy(GrowthStrainsRotMat,rotMat);
        }
    }
    gsl_matrix_free(rotMat);
}

void 	ShapeBase::calculateTriPointFForRatation(){
	gsl_matrix_set_zero(TriPointF);
    gsl_matrix* currF = gsl_matrix_calloc(nDim,nDim);
    for (size_t iter =0; iter<numberOfGaussPoints;++iter){
        calculateCurrTriPointFForRotation(currF,iter);
        gsl_matrix_scale(currF,gaussWeights[iter]);
        gsl_matrix_add(TriPointF, currF);
    }
    gsl_matrix_free(currF);
}

bool 	ShapeBase::disassembleRotationMatrixForZ(gsl_matrix* rotMat){
    /**  This function extracts the z rotation from a rotation matrix. From Extracting euler angles from
     * a rotation matrix by Mike Day of insomniac games: \n
     *
     * To extract a rotation of \f$ \mathbf{R}_{x}(\theta_{1}) \mathbf{R}_{y}(\theta_{2}) \mathbf{R}_{z}(\theta_{3})\f$ from a matrix
     * M where:
     * 	\f{eqnarray*}{
      \mathbf{M} = \begin{bmatrix}
                   m_{00} &  m_{01} &  m_{02} \\
                   m_{10} &  m_{11} &  m_{12} \\
                   m_{20} &  m_{21} &  m_{22}]
                   \end{bmatrix}
       \f}
    * and \f$ c_{1} \f$  denote \f$ cos(\theta_{1}) \f$ and \f$ s_{1} \f$ denote \f$ sin(\theta_{1}) \f$:
    * \f{eqnarray*}{
    *   tet_1 & = & atan2(m_{12},m_{22}) = atan2(s_1c_2, c_1c_2) \\
    *   c_2 & = & \sqrt{m_{00}*m_{00} + m_{01}*m_{01}}  \\
    *   tet_2 & = & atan2(-m_{02},c_2) \\
    *   tet_3 & = & atan2(m_{01},m_{00})
    * \f}
    */
    double tethaZ = atan2(gsl_matrix_get(rotMat,0,1),gsl_matrix_get(rotMat,0,0));
    if (tethaZ > 0.017 || tethaZ < -0.017){ //0.017 rad is ~1 degrees
        //rotation is more than 1 degrees, element incremental growth should be rotated
        double c = cos(tethaZ);
        double s = sin(tethaZ);
        gsl_matrix_set(rotMat,0,0,  c );
        gsl_matrix_set(rotMat,0,1, -1.0*s);
        gsl_matrix_set(rotMat,0,2,  0.0);
        gsl_matrix_set(rotMat,1,0,  s);
        gsl_matrix_set(rotMat,1,1,  c);
        gsl_matrix_set(rotMat,1,2,  0.0);
        gsl_matrix_set(rotMat,2,0,  0.0);
        gsl_matrix_set(rotMat,2,1,  0.0);
        gsl_matrix_set(rotMat,2,2,  1.0);
        return true;
    }
    else{
        return false;   //rotation is less than 1 degrees;
    }
}

void 	ShapeBase::calculateEmergentRotationAngles(){
	calculateTriPointFForRatation();
	gsl_matrix* rotMat = gsl_matrix_calloc(3,3);
	calculate3DRotMatFromF(rotMat);
	gsl_matrix_free (rotMat);
}

bool 	ShapeBase::calculate3DRotMatFromF(gsl_matrix* rotMat){
    /** The rigid body rotation is extracted via single value decomposition.
     */
    gsl_matrix* Sgsl = gsl_matrix_alloc (3, 3);
    gsl_matrix* V = gsl_matrix_alloc (3, 3);
    gsl_matrix* R = gsl_matrix_alloc (3, 3);
    gsl_vector* Sig = gsl_vector_alloc (3);
    gsl_vector* workspace = gsl_vector_alloc (3);

    //Singular Value Decomposition
    createMatrixCopy (Sgsl,TriPointF);

    (void) gsl_linalg_SV_decomp (Sgsl, V, Sig, workspace); //This function does have a return value, but I do not need to use it, hence cast it to void!
    /** The decomposition isdone by the gsl routine gsl_linalg_SV_decomp, and here, the output gives
     * Sgsl as \f$ \mathbf{U} \f$ , the rotation matrix \f$ \mathbf{U} \mathbf{V}^{T}  \f$ in the
     * decomposition \f$ \mathbf{A} = \mathbf{U} \mathbf{S} \mathbf{V}^{T} \f$
    */
    //Sgsl ended up as U, I need the rotation matrix U V^T in the decomposition A = U S V^T (jose writes as F = V S U^T in emails, so be careful)
    gsl_blas_dgemm (CblasNoTrans, CblasTrans,1.0, Sgsl, V,0.0, rotMat);

    double det = determinant3by3Matrix(rotMat);
    if (det<0){
        std::cout<<"Error! Flipped element, Id: "<<Id<<std::endl;
        isFlipped = true;
    }
    gsl_matrix_free (Sgsl);
    gsl_matrix_free (V);
    gsl_matrix_free (R);
    gsl_vector_free (Sig);
    gsl_vector_free (workspace);
    //Now I need to check if there is only numerical error accumulationg on rotMat, or there is an actual rotation (above 1 degrees):
    double threshold = 0.017; //this is sine 1 degrees
    for (int i=0;i<3;++i){
        for (int j=0;j<3;++j){
            if(i != j){
                if (gsl_matrix_get(rotMat,i,j)>threshold || gsl_matrix_get(rotMat,i,j)< (-1.0*threshold)) {
                    return true;
                }
            }
        }
    }
    return false; //none of the off - diagonal terms of the matrix are above the threshold, the current rotation is only numerical error.
}

void ShapeBase::mutateElement(double growthFold, double growthRatePerHour){
	isMutated = true;
	//growth will be uniform in x and y
	mutationGrowthRatePerSec= growthRatePerHour/2.0/3600.0;
	mutationGrowthFold= growthFold;
}

void ShapeBase::updateGrowthByMutation(double dt){
    /** Mutation can be defined in a scaled growth increase, stored in ShapeBase#mutationGrowthFold, or a
     * direct overwriting of hte growth, with ShapeBase#mutationGrowthRatePerSec. If the growth fold in non-zero,
     * growth in xy plane is scaled accordingly, and set in the rates ShapeBase#setGrowthRate.
     * If the growth is set directly, then the increment is updated directly to the given value.
     * Then the incremnt is updated via ShapeBase#updateGrowthIncrementFromRate.
     */
	if (mutationGrowthFold>0){
		//the fold increase is not zero, which means I should be inducing fold increase in growth:
		//I will take the curretn increment, calculate the determinant (absolute growth).
		//Then I will redistribute this growth to x & y:
		double growthPerDt = determinant3by3Matrix(growthIncrement);
		double ratePerSec = log(growthPerDt)/dt;
		double growthPer24hrs = exp(ratePerSec*3600.0*24.0);
		double newGrowthPer24hrs = mutationGrowthFold*growthPer24hrs;
		double newGrowthRatePerSec = log(newGrowthPer24hrs)/24.0/3600.0/2.0;
		setGrowthRate(dt,newGrowthRatePerSec,newGrowthRatePerSec,0.0);
		updateGrowthIncrementFromRate();
	}
	else{
		//overwriting up any growth that might be there, with uniform growth in x & y:
        setGrowthRate(dt,mutationGrowthRatePerSec,mutationGrowthRatePerSec,0.0);
		updateGrowthIncrementFromRate();
	}
}

void 	ShapeBase::calculateRelativePosInBoundingBox(double boundingBoxXMin, double boundingBoxYMin, double boundingBoxLength, double boundingBoxWidth){
	relativePosInBoundingBox = getCentre();
	relativePosInBoundingBox[0] = (relativePosInBoundingBox[0] -boundingBoxXMin) / boundingBoxLength;
	relativePosInBoundingBox[1] = (relativePosInBoundingBox[1] - boundingBoxYMin) / boundingBoxWidth;

}

bool ShapeBase::isElementFlippedInPotentialNewShape(int nodeId, double newX, double newY, double newZ){
    /**
     * The current positions are taken on a temporary matrix, and the position update is carried out to the new
     * position specified in the funciton input. Then the new elastic deformation gradient is calculated
     * through the standard procedure. If the determinant of the deformation gradient  is negative, then the
     * element is flipped and the node position update the the new position should be avoided.
     */
    const int n = nNodes;
    const int dim = nDim;
    gsl_matrix* currF = gsl_matrix_alloc(dim,dim);
    gsl_matrix* currFe = gsl_matrix_alloc(dim,dim);
    gsl_matrix* CurrShape = gsl_matrix_alloc(n,dim);
    getPos(CurrShape);
    for (int i=0;i<n;++i){
    	if (NodeIds[i] == nodeId){
    		gsl_matrix_set(CurrShape,i,0,newX);
    		gsl_matrix_set(CurrShape,i,1,newY);
    		gsl_matrix_set(CurrShape,i,2,newZ);
    	}
    }
    bool elementWillFlip = false;
    for (size_t pointNo =0; pointNo<numberOfGaussPoints;++pointNo){
		//calculating dx/de (Jacobian) and reading out dX/de, shape function derivaties:
		gsl_matrix* ShapeFuncDer = ShapeFuncDerivatives[pointNo];
		gsl_matrix* InvdXde = InvdXdes[pointNo];
		gsl_matrix* Jacobian = gsl_matrix_calloc(dim, dim);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, ShapeFuncDer, CurrShape, 0.0, Jacobian);
		gsl_matrix_transpose(Jacobian);
		//calculating F:
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Jacobian, InvdXde, 0.0, currF);
		//calculating Fe:
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, currF, InvFg, 0.0, currFe);	///< Removing growth
		double detFe = determinant3by3Matrix(currFe);
		if (detFe<0){
		  //element will be flipped if I collapse its node
			elementWillFlip = true;
		}
		gsl_matrix_free(Jacobian);
		if (elementWillFlip){
			break;
		}
    }
	gsl_matrix_free(CurrShape);
	gsl_matrix_free(currFe);
	gsl_matrix_free(currF);
	return elementWillFlip;
}


void ShapeBase::checkForCollapsedNodes(int TissueHeightDiscretisationLayers, const std::vector<std::unique_ptr<Node>>& Nodes, const std::vector<std::unique_ptr<ShapeBase>>& Elements){
    /** If the element has any nodes that are collapsed with any other node, then the element is
     * considered collapsed. This definition does not make any distinction betwenn elements collapsing on
     * themselves to avoid flips, or the nodal collapse due to adhesion of two elements.
     */
    bool elementCollapsed = false;
    for (size_t j =0 ; j<nNodes; ++j){
		int nodeId = NodeIds[j];
		int collapsedNodeNumber = Nodes[nodeId]->collapsedWith.size();
		if (collapsedNodeNumber>0){elementCollapsed = true;break;}
	}
	if (elementCollapsed){
		insideEllipseBand = true;
		int selectedEllipseBandId = 100;
		if (tissuePlacement == 1){ //apical collapse: ECM relaxation, cell shortening, volume redistribution to shrink top
			selectedEllipseBandId = 100;
		}
		else{ //basal collapse, volume redistribution to shrink bottom
			selectedEllipseBandId = 101;
		}
		coveringEllipseBandId = selectedEllipseBandId;
		assignEllipseBandIdToWholeTissueColumn(TissueHeightDiscretisationLayers,Nodes,Elements);
	}
}

bool ShapeBase::hasEnoughNodesOnCurve(const std::vector<std::unique_ptr<Node>>& Nodes){
	int threshold = 3;
	int nodesOnCurveCounter = 0;
    for (size_t i=0;i<nNodes; i++){
		if (Nodes[NodeIds[i]]->onFoldInitiation){
			nodesOnCurveCounter++;
		}
		if (nodesOnCurveCounter>=threshold){
			return true;
		}
	}
	return false;
}

void ShapeBase::assignEllipseBandIdToWholeTissueColumn(size_t TissueHeightDiscretisationLayers, const std::vector<std::unique_ptr<Node>>& Nodes, const std::vector<std::unique_ptr<ShapeBase>>& Elements){
    for (size_t i=0; i < TissueHeightDiscretisationLayers;++i){
		int idOfElementOnSameColumn = elementsIdsOnSameColumn[i];
		Elements[idOfElementOnSameColumn]->assignEllipseBandId(Nodes,coveringEllipseBandId);
	}
}

void ShapeBase::assignEllipseBandId(const std::vector<std::unique_ptr<Node>>& Nodes, int selectedEllipseBandId){
	insideEllipseBand = true;
	coveringEllipseBandId = selectedEllipseBandId;
	assignEllipseBandIdToNodes(Nodes);
}

void ShapeBase::assignEllipseBandIdToNodes(const std::vector<std::unique_ptr<Node>>& Nodes){
    for (size_t j =0 ; j<nNodes; ++j){
		int nodeId = NodeIds[j];
		Nodes[nodeId]->insideEllipseBand = true;
		Nodes[nodeId]->coveringEllipseBandId = coveringEllipseBandId;
		Nodes[nodeId]->checkOwnersforEllipseAsignment = true;
	}
}

double ShapeBase::getElementalElasticForce(int nodeIndex, int dimIndex){
	return gsl_matrix_get(ElementalElasticSystemForces,nodeIndex,dimIndex);
}

void ShapeBase::setElementalElasticForce(int nodeIndex, int dimIndex, double value){
	 gsl_matrix_set(ElementalElasticSystemForces,nodeIndex,dimIndex,value);
}

void ShapeBase::addToElementalElasticSystemForces(int i,int j,double value){
	double baseValue = gsl_matrix_get(ElementalElasticSystemForces,i,j);
	baseValue +=value;
	gsl_matrix_set(ElementalElasticSystemForces,i,j,baseValue);
}

void ShapeBase::addToTriPointKe(int i,int j,double value){
    /** This funciton adds the input value to the elemental stiffness matrix, elastic part of the elemental
     * Jacobian, ShapeBase#TriPointKe, at the input indices (i,j).
     */
	double baseValue = gsl_matrix_get(TriPointKe,i,j);
	baseValue +=value;
	gsl_matrix_set(TriPointKe,i,j,baseValue);
}


void 	ShapeBase::displayRelativePosInBoundingBox(){
        std::cout<<"Element: "<<Id<<"  relative position in the tissue bounding box: "<<relativePosInBoundingBox[0]<<" "<<relativePosInBoundingBox[1]<<std::endl;
}

bool 	ShapeBase::DoesPointBelogToMe(int IdNode){
    for (size_t i = 0; i<nNodes; ++i){
		if (NodeIds[i] == IdNode){
			return true;
		}
	}
	return false;
}
/*
double 	ShapeBase::determinant3by3Matrix(double* rotMat){
	double det =0.0;
	det  =  rotMat[0]*(rotMat[4]*rotMat[8]-rotMat[5]*rotMat[7]);
	det -= rotMat[1]*(rotMat[3]*rotMat[8]-rotMat[5]*rotMat[6]);
	det += rotMat[2]*(rotMat[3]*rotMat[7]-rotMat[4]*rotMat[6]);
	return det;
}*/

double 	ShapeBase::determinant3by3Matrix(boost::numeric::ublas::matrix<double>& Mat){
	double det =0.0;
	det  =  Mat(0,0)*(Mat(1,1)*Mat(2,2)-Mat(1,2)*Mat(2,1));
	det -= Mat(0,1)*(Mat(1,0)*Mat(2,2)-Mat(1,2)*Mat(2,0));
	det += Mat(0,2)*(Mat(1,0)*Mat(2,1)-Mat(1,1)*Mat(2,0));
	return det;
}

double 	ShapeBase::determinant3by3Matrix(gsl_matrix* Mat){
    double det =0.0;
    det  = gsl_matrix_get(Mat,0,0)* ( gsl_matrix_get(Mat,1,1)*gsl_matrix_get(Mat,2,2)-gsl_matrix_get(Mat,1,2)*gsl_matrix_get(Mat,2,1) );
    det -= gsl_matrix_get(Mat,0,1)* ( gsl_matrix_get(Mat,1,0)*gsl_matrix_get(Mat,2,2)-gsl_matrix_get(Mat,1,2)*gsl_matrix_get(Mat,2,0) );
    det += gsl_matrix_get(Mat,0,2)* ( gsl_matrix_get(Mat,1,0)*gsl_matrix_get(Mat,2,1)-gsl_matrix_get(Mat,1,1)*gsl_matrix_get(Mat,2,0) );
    return det;
}
double 	ShapeBase::determinant2by2Matrix(boost::numeric::ublas::matrix<double>& Mat){
	double det = Mat(0,0) * Mat(1,1) - Mat(0,1) * Mat(1,0);
	return det;
}

void	ShapeBase::calculateRotationAngleSinCos(std::array<double,3>& u, std::array<double,3>& v, double& c, double& s){
	//aligning u onto v:
	c = dotProduct3D(u,v);
	if (c > 1.0){
		c = 1.0;
		s = 0.0;
	}
	else if( c<-1.0){
		c = -1.0;
		s = 0.0;
	}
	else{
		double tet = acos(c);
		s = sin(tet);
	}
}

void	ShapeBase::calculateRotationAxis(const std::array<double,3>& u, const std::array<double,3> &v, std::array<double,3>& rotAx, double c){
	//aligning u onto v:
	if (c>-0.99998){
        rotAx = crossProduct3D(u,v);
        (void) normaliseVector3D(rotAx);
	}
	else{
		//the angle is 180 degree, the standard rotation axis calculation will be wrong, I am rotating over x axis at all times;
		rotAx[0]= 1;rotAx[1]= 0;rotAx[2]= 0;
	}
}

void	ShapeBase::constructRotationMatrix(double c, double s, double* rotAx, double* rotMat){
	rotMat[0] = c + rotAx[0]*rotAx[0]*(1 - c);
	rotMat[1] = rotAx[0]*rotAx[1]*(1 - c) - rotAx[2]*s;
	rotMat[2] = rotAx[0]*rotAx[2]*(1 - c) + rotAx[1]*s;

	rotMat[3] = rotAx[1]*rotAx[0]*(1 - c) + rotAx[2]*s;
	rotMat[4] = c + rotAx[1]*rotAx[1]*(1 - c);
	rotMat[5] = rotAx[1]*rotAx[2]*(1 - c) - rotAx[0]*s;

	rotMat[6] = rotAx[2]*rotAx[0]*(1 - c) - rotAx[1]*s;
	rotMat[7] = rotAx[2]*rotAx[1]*(1 - c) + rotAx[0]*s;
	rotMat[8] = c + rotAx[2]*rotAx[2]*(1 - c);
}

void	ShapeBase::constructRotationMatrix(double c, double s, std::array<double,3>& rotAx, std::array<double,9>& rotMat){
    rotMat[0] = c + rotAx[0]*rotAx[0]*(1 - c);
    rotMat[1] = rotAx[0]*rotAx[1]*(1 - c) - rotAx[2]*s;
    rotMat[2] = rotAx[0]*rotAx[2]*(1 - c) + rotAx[1]*s;

    rotMat[3] = rotAx[1]*rotAx[0]*(1 - c) + rotAx[2]*s;
    rotMat[4] = c + rotAx[1]*rotAx[1]*(1 - c);
    rotMat[5] = rotAx[1]*rotAx[2]*(1 - c) - rotAx[0]*s;

    rotMat[6] = rotAx[2]*rotAx[0]*(1 - c) - rotAx[1]*s;
    rotMat[7] = rotAx[2]*rotAx[1]*(1 - c) + rotAx[0]*s;
    rotMat[8] = c + rotAx[2]*rotAx[2]*(1 - c);
}

void	ShapeBase::rotateVectorByRotationMatrix(double* u,double* rotMat){
	double x = rotMat[0]*u[0]+rotMat[1]*u[1]+rotMat[2]*u[2];
	double y = rotMat[3]*u[0]+rotMat[4]*u[1]+rotMat[5]*u[2];
	double z = rotMat[6]*u[0]+rotMat[7]*u[1]+rotMat[8]*u[2];
	u[0] = x;
	u[1] = y;
	u[2] = z;
}

void	ShapeBase::rotateVectorByRotationMatrix(std::array<double,3>& u, std::array<double,9> rotMat){
    double x = rotMat[0]*u[0]+rotMat[1]*u[1]+rotMat[2]*u[2];
    double y = rotMat[3]*u[0]+rotMat[4]*u[1]+rotMat[5]*u[2];
    double z = rotMat[6]*u[0]+rotMat[7]*u[1]+rotMat[8]*u[2];
    u[0] = x;
    u[1] = y;
    u[2] = z;
}

void	ShapeBase::rotateVectorByRotationMatrix(double* u,gsl_matrix* rotMat){
	double x = gsl_matrix_get(rotMat,0,0)*u[0]+gsl_matrix_get(rotMat,0,1)*u[1]+gsl_matrix_get(rotMat,0,2)*u[2];
	double y = gsl_matrix_get(rotMat,1,0)*u[0]+gsl_matrix_get(rotMat,1,1)*u[1]+gsl_matrix_get(rotMat,1,2)*u[2];
	double z = gsl_matrix_get(rotMat,2,0)*u[0]+gsl_matrix_get(rotMat,2,1)*u[1]+gsl_matrix_get(rotMat,2,2)*u[2];
	u[0] = x;
	u[1] = y;
	u[2] = z;
}

void  ShapeBase::rotateReferenceElementByRotationMatrix(std::array<double,9> rotMat){
    for (size_t i=0; i<nNodes; ++i){
        std::array<double,3> u = {0.0, 0.0, 0.0};
        for (size_t j=0; j<nDim; ++j){
            u[j] = ReferenceShape->Positions[i][j];
        }
        rotateVectorByRotationMatrix(u,rotMat);
        for (size_t j=0; j<nDim; ++j){
            ReferenceShape->Positions[i][j] = u[j];
        }
    }
}
void	ShapeBase::calculateForces(const std::vector<std::unique_ptr<Node>>& Nodes, gsl_matrix* displacementPerDt){
	if (ShapeDim == 3){		//3D element
        calculateForces3D(Nodes, displacementPerDt);
    }
}

void ShapeBase::writeInternalForcesTogeAndgv(gsl_matrix* ge, gsl_matrix* gvInternal, std::vector<std::array<double,3>>& SystemForces, const std::vector<std::unique_ptr<Node>>& Nodes){
	//now all the forces are written on SysyemForces
    //Now I will add the forces into ge, this step can be made faster by separating calculate forces function into two,
    //and filling up either ge or System forces depending on the solution method:
    for (size_t i = 0; i< nNodes; ++i){
        for ( size_t j=0; j<nDim; j++){
            size_t indexI = nDim*NodeIds[i]+j;
        	double elementalvalue = gsl_matrix_get(ElementalElasticSystemForces,i,j);
        	double matrixValue = gsl_matrix_get(ge,indexI,0);
            gsl_matrix_set(ge, indexI,0,matrixValue + elementalvalue);
        	elementalvalue = gsl_matrix_get(ElementalInternalViscousSystemForces,i,j);
        	matrixValue = gsl_matrix_get(gvInternal,indexI,0);
            gsl_matrix_set(gvInternal, indexI,0,matrixValue + elementalvalue);
        }
    }
    int counter = 0;
    for (size_t i = 0; i<nNodes; ++i){
        for (size_t j = 0; j<nDim; ++j){
            if (!Nodes[NodeIds[i]]->FixedPos[j]){
               // std::cout<<"element: "<<Id<<" writing to system forces for NodeID: "<<NodeIds[i]<<" dim: "<<j<<" initial value "<<SystemForces[NodeIds[i]][j]<<" Felastic: "<<gsl_matrix_get(ElementalElasticSystemForces,i,j)<<" Fvisc: "<<gsl_matrix_get(ElementalInternalViscousSystemForces,i,j)<<std::endl;
                SystemForces[NodeIds[i]][j] = SystemForces[NodeIds[i]][j] + gsl_matrix_get(ElementalElasticSystemForces,i,j) + gsl_matrix_get(ElementalInternalViscousSystemForces,i,j);

            }
            counter++;
        }
    }
}

void	ShapeBase::calculateForces3D(const std::vector<std::unique_ptr<Node>>& Nodes,  gsl_matrix* displacementPerDt){
    size_t dim = nDim;
    size_t n = nNodes;
    //calculating F and B in a 3 point gaussian:
    gsl_matrix* TriPointge  = gsl_matrix_calloc(dim*n,1);
    gsl_matrix* TriPointgv  = gsl_matrix_calloc(dim*n,1);
    gsl_matrix_set_zero(TriPointF);
    gsl_matrix_set_zero(ElementalElasticSystemForces);
    gsl_matrix_set_zero(ElementalInternalViscousSystemForces);
    gsl_matrix* currge = gsl_matrix_calloc(dim*n,1);
    gsl_matrix* currgv = gsl_matrix_calloc(dim*n,1);
    gsl_matrix* currF = gsl_matrix_calloc(dim,dim);
    //The point order is established in shape function derivative calculation!
    /** The nodal forces are calculated as an average of all Gauss Points. The number of Gauss Points
     * and their weights are stored in ShapeBase#numberOfGaussPoints and ShapeBase#gaussWeights, respectively.
     * The order of points weights should be consistent with point definition order in shape function derivative
     * calculation. The nodal forces for the selected Gauss point are calculated via ShapeBase#calculateCurrNodalForces.
     * Then the output nodal elastic forces, viscous forces, and the deformation gradient are scaled with weights.
     */
    for (size_t iter =0; iter<numberOfGaussPoints;++iter){
    	calculateCurrNodalForces(currge, currgv, currF, displacementPerDt, iter);
        gsl_matrix_scale(currge,gaussWeights[iter]);
        gsl_matrix_add(TriPointge, currge);
        gsl_matrix_scale(currgv,gaussWeights[iter]);
        gsl_matrix_add(TriPointgv, currgv);
        gsl_matrix_scale(currF,gaussWeights[iter]);
        gsl_matrix_add(TriPointF, currF);
    }
    /**
     * Then the elemental forces are written on the system forces.
     */
    int counter = 0;
    for (size_t i = 0; i<nNodes; ++i){
            for (size_t j = 0; j<nDim; ++j){
            	if (!Nodes[NodeIds[i]]->FixedPos[j]){
            		double value = gsl_matrix_get(ElementalElasticSystemForces,i,j);
            		value -= gsl_matrix_get(TriPointge,counter,0);
            		gsl_matrix_set(ElementalElasticSystemForces,i,j,value);
                    if(std::isnan(value)){
                        std::cout<<"force from ElementalElasticSystemForces for element "<<Id<<" dim: "<<j<<" is nan"<<std::endl;
					}
            		value = gsl_matrix_get(ElementalInternalViscousSystemForces,i,j);
            		value -= gsl_matrix_get(TriPointgv,counter,0);
            		gsl_matrix_set(ElementalInternalViscousSystemForces,i,j,value);
                    if(std::isnan(value)){
                        std::cout<<"force from ElementalInternalViscousSystemForces for element "<<Id<<" dim: "<<j<<" is nan"<<std::endl;
            		}
				}
				counter++;
            }
    }
    //freeing matrices allocated in this function
    gsl_matrix_free(TriPointge);
    gsl_matrix_free(TriPointgv);
    gsl_matrix_free(currge);
    gsl_matrix_free(currgv);
    gsl_matrix_free(currF);
}


gsl_matrix* ShapeBase::calculateCauchyGreenDeformationTensor(gsl_matrix* Fe){
    /** The Cauchy Green deformation Tensor \f$ \mathbf{C} \f$ is defned as: \n
     * \f$ \mathbf{C}  = \left( \mathbf{F}^{eT} \mathbf{F}^{e}  \right) \f$
     * where \f$ \mathbf{F}^{e} \f$  is the elastic part of the deformation gradient.
     */
	//calculating C (C = (Fe^T*Fe):
	gsl_matrix* C =  gsl_matrix_alloc(nDim, nDim);
	gsl_blas_dgemm (CblasTrans, CblasNoTrans,1.0, Fe, Fe,0.0, C);
	return C;
}

gsl_matrix* ShapeBase::calculateEForNodalForcesKirshoff(gsl_matrix* C){
    /** \f$ \mathbf{E} \f$ for Kirshoff material is: \n
     * \f$ \mathbf{E}  = \frac{1}{2} \times \left( \mathbf{C}- \mathbf{I}  \right) \f$
     * where \f$ \mathbf{C} \f$  is the Cauchy Green deformation Tensor calculated in
     * ShapeBase#calculateCauchyGreenDeformationTensor and \f$ \mathbf{I} \f$ is identity.
     */
    //calculating E ( E = 1/2 *(Fe^T*Fe-I) ; E = 1/2 *(C-I):):
    gsl_matrix* E =  gsl_matrix_alloc(nDim, nDim);
	createMatrixCopy(E,C);
    gsl_matrix* I = gsl_matrix_alloc(nDim, nDim);
    gsl_matrix_set_identity(I);
    gsl_matrix_sub(E,I);
    gsl_matrix_scale(E, 0.5);
    gsl_matrix_free(I);
    return E;
}

gsl_matrix* ShapeBase::calculateSForNodalForcesKirshoff(gsl_matrix* E){
    //calculating S: (S = D:E)
    /** The Second order Piola- Kirshoff Stress Tensor \f$ \mathbf{S} \f$ for Kirshoff material is: \n
     * \f$ \mathbf{S}  =  \mathbf{D}: \mathbf{E} \f$
     * where \f$ \mathbf{E} \f$  is calculated in
     * ShapeBase#calculateEForNodalForcesKirshoff and \f$ \mathbf{D} \f$ is ShapeBase#D.
     */
    gsl_matrix_set_zero(Strain);
    gsl_matrix* compactS = gsl_matrix_calloc(6,1);
    gsl_matrix_set(Strain,0,0, gsl_matrix_get(E,0,0));
    gsl_matrix_set(Strain,1,0, gsl_matrix_get(E,1,1));
    gsl_matrix_set(Strain,2,0, gsl_matrix_get(E,2,2));
    gsl_matrix_set(Strain,3,0, 2.0*gsl_matrix_get(E,0,1));
    gsl_matrix_set(Strain,4,0, 2.0*gsl_matrix_get(E,2,1));
    gsl_matrix_set(Strain,5,0, 2.0*gsl_matrix_get(E,0,2));
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, D, Strain,0.0, compactS);

    gsl_matrix* S =  gsl_matrix_alloc(nDim, nDim);
    gsl_matrix_set(S,0,0,gsl_matrix_get(compactS,0,0));
    gsl_matrix_set(S,1,1,gsl_matrix_get(compactS,1,0));
    gsl_matrix_set(S,2,2,gsl_matrix_get(compactS,2,0));
    gsl_matrix_set(S,1,0,gsl_matrix_get(compactS,3,0));
    gsl_matrix_set(S,1,2,gsl_matrix_get(compactS,4,0));
    gsl_matrix_set(S,2,0,gsl_matrix_get(compactS,5,0));
    gsl_matrix_set(S,0,1,gsl_matrix_get(compactS,3,0));
    gsl_matrix_set(S,2,1,gsl_matrix_get(compactS,4,0));
    gsl_matrix_set(S,0,2,gsl_matrix_get(compactS,5,0));

    gsl_matrix_free(compactS);
    return S;
}

gsl_matrix* ShapeBase::calculateSForNodalForcesNeoHookean(gsl_matrix* invC, double lnJ){
    /** The Second order Piola- Kirshoff Stress Tensor for a Neo-Hookean material is:
     * \f$
       \mathbf{S}^e  = \mu (\mathbf {I} - \mathbf {C^{-1}}) + \lambda (ln J^e) \mathbf {C^{-1}}
     * \f$
     */
	//S = mu (I - C^-1) + lambda (lnJ) C^-1
	gsl_matrix* S =  gsl_matrix_alloc(nDim, nDim);
	createMatrixCopy(S,invC);
	gsl_matrix* I = gsl_matrix_alloc(nDim, nDim);
	gsl_matrix_set_identity(I);
    gsl_matrix_sub(I,invC);  //(I - C^-1)
    gsl_matrix_scale(I, mu); // mu (I - C^-1)
    gsl_matrix_scale(S, lambda*lnJ); //lambda (lnJ) C^-1
    gsl_matrix_add(S,I); // mu (I - C^-1) + lambda (lnJ) C^-1
    gsl_matrix_free(I);
	return S;
}

void ShapeBase::updateLagrangianElasticityTensorNeoHookean(gsl_matrix* invC, double lnJ, int pointNo){
    /** The voigt notation of the transformed elasticity tensor \f$ \mathcal D \f$ is calculated  via
     * \f$ \mathcal{D}_{ijkl}  = \lambda \boldsymbol{C}^{-1}_{ij} \boldsymbol{C}^{-1}_{kl} + 2 \left( \mu - \lambda ln(J)\right) \mathcal{I}_{ijkl} \f$ \n
     * where \f$ \lambda \f$ Lame's first parameter and \f$ \mu \f$ is the shear modulus. Each value at index ijkl is
     * calculated for each Gauss Point and the information is stored in an array of
     * D[gauss points] [system dimension (3D)] [system dimension (3D)] [system dimension (3D)] [system dimension (3D)].
     */
	double multiplier = 2*(mu - lambda*lnJ);
    for (size_t I = 0; I<nDim; ++I){
        for (size_t J = 0; J<nDim; ++J){
            for (size_t K = 0; K<nDim; ++K){
                for (size_t L = 0; L<nDim; ++L){
                    double Iijkl = 0.5* (gsl_matrix_get(invC,I,K)*gsl_matrix_get(invC,J,L) + gsl_matrix_get(invC,I,L)*gsl_matrix_get(invC,J,K));
                    D81[pointNo][I][J][K][L] = lambda*gsl_matrix_get(invC,I,J)*gsl_matrix_get(invC,K,L) + multiplier * Iijkl;
                }
            }
        }
    }
}

gsl_matrix* ShapeBase::calculateCompactStressForNodalForces(double detFe, gsl_matrix* Fe, gsl_matrix* S, gsl_matrix* Stress){
    /** Stress is calculated via: \n
     * \f$ \mathbf{\sigma}^e  = J^{e-1} \mathbf{F}^{e} \mathbf{S}^{e} \mathbf{F}^{eT}  \f$
     */
    gsl_matrix* tmpMat1 =  gsl_matrix_calloc(nDim, nDim);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Fe, S,0.0, tmpMat1);
    gsl_blas_dgemm (CblasNoTrans, CblasTrans,1.0, tmpMat1, Fe,0.0, Stress);
    gsl_matrix_scale(Stress, 1.0/detFe);
    gsl_matrix* compactStress =  gsl_matrix_calloc(6,1);
    gsl_matrix_set(compactStress,0,0,gsl_matrix_get(Stress,0,0));
    gsl_matrix_set(compactStress,1,0,gsl_matrix_get(Stress,1,1));
    gsl_matrix_set(compactStress,2,0,gsl_matrix_get(Stress,2,2));
    gsl_matrix_set(compactStress,3,0,gsl_matrix_get(Stress,0,1));
    gsl_matrix_set(compactStress,4,0,gsl_matrix_get(Stress,2,1));
    gsl_matrix_set(compactStress,5,0,gsl_matrix_get(Stress,0,2));

    gsl_matrix_free(tmpMat1);
    return compactStress;
}

void ShapeBase::calculateInvJShFuncDerSWithFe(gsl_matrix* currFe, gsl_matrix* InvDXde, gsl_matrix* ShapeFuncDerStack, gsl_matrix *invJShFuncDerSWithFe){
	//I want InvJe, normally J InvDXde = F, I can get Je from
	// Je InvDXde = Fe
	// but I can also get InvJe directly from:
	// InvJe Je InvdXde = InvJe Fe => I InvdXde = InvJe Fe => InvdXde InvFe = InvJe I => InvJe = InvdXde InvFe
	gsl_matrix* tmpFeforInversion =  gsl_matrix_calloc(nDim,nDim);
	gsl_matrix* InvFe = gsl_matrix_calloc(nDim,nDim);
	gsl_matrix* InvJe = gsl_matrix_calloc(nDim,nDim);
	createMatrixCopy(tmpFeforInversion,currFe);
	bool inverted = InvertMatrix(tmpFeforInversion, InvFe);
	if (!inverted){
        std::cerr<<"Fe not inverted!!"<<std::endl;
	}
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, InvDXde, InvFe,0.0, InvJe);

	int dim2 = nDim*nDim;
	//Generating the inverse Jacobian(elastic) stack:
	gsl_matrix* InvJacobianElasticStack =  gsl_matrix_calloc(dim2,dim2);
    for (size_t i =0; i<nDim; i++){
        for (size_t m=0; m<nDim; ++m){
            for (size_t n=0; n<3; ++n){
				gsl_matrix_set(InvJacobianElasticStack,i*nDim+m,i*nDim+n,gsl_matrix_get(InvJe,n,m));
			}
		}
	}

	//I am calculating this for k calculation, in case there is growth. Under conditions that there is no growth, this function is not necessary,
	//the values of invJShFuncDerSWithF and  invJShFuncDerS will be equal
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, InvJacobianElasticStack, ShapeFuncDerStack,0.0, invJShFuncDerSWithFe);
	gsl_matrix_free(tmpFeforInversion);
	gsl_matrix_free(InvFe);
	gsl_matrix_free(InvJe);
	gsl_matrix_free(InvJacobianElasticStack);
}

gsl_matrix* ShapeBase::calculateBTforNodalForces(gsl_matrix* InvJacobianStack, gsl_matrix* ShapeFuncDerStack, gsl_matrix* B, gsl_matrix *invJShFuncDerS){
     //calculating B:
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, InvJacobianStack, ShapeFuncDerStack,0.0, invJShFuncDerS);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, CoeffMat, invJShFuncDerS,0.0, B);
    //generating B^T:
    gsl_matrix* BT = gsl_matrix_alloc(nNodes*nDim,6);
    gsl_matrix_transpose_memcpy(BT,B);
    return BT;
}

gsl_matrix* ShapeBase::calculateInverseJacobianStackForNodalForces(gsl_matrix* Jacobian){
    int dim2 = nDim*nDim;
    //invrting the Jacobian:
    gsl_matrix* tmpJacobianForInversion =  gsl_matrix_calloc(nDim,nDim);
    gsl_matrix* InvJacobian = gsl_matrix_calloc(nDim,nDim);
    createMatrixCopy(tmpJacobianForInversion,Jacobian);
    bool inverted = InvertMatrix(tmpJacobianForInversion, InvJacobian);
    if (!inverted){
        std::cerr<<"Jacobian not inverted!!"<<std::endl;
    }
    //displayMatrix(Jacobian,"Jacobian");
    //Generating the inverse Jacobian stack:
    gsl_matrix* InvJacobianStack =  gsl_matrix_calloc(dim2,dim2);
    for (size_t i =0; i<nDim; i++){
        for (size_t m=0; m<nDim; ++m){
            for (size_t n=0; n<3; ++n){
                gsl_matrix_set(InvJacobianStack,i*nDim+m,i*nDim+n,gsl_matrix_get(InvJacobian,n,m));
            }
        }
    }
    gsl_matrix_free(tmpJacobianForInversion);
    gsl_matrix_free(InvJacobian);
    return InvJacobianStack;
}

gsl_matrix* ShapeBase::calculateVelocityGradientTensor(gsl_matrix* B, gsl_matrix* displacementPerDt){
	/**
	 * Inputs:
	 * -# The elemental B matrix (6 , ShapeBase#nDim x ShapeBase#nNodes).
	 * -# The displacement of all nodes of the system, divided by the time
	 * step (ShapeBase#nDim x Simulation#nNodes, 1).
	 *
	 * Output:
	 * -# Velocity gradient tensor in Voigt notation (6, 1).
	 *
	 * This function calculates the velocity gradient tensor from elemental B matrix and elemental displacement.
	 * The elemental B matrix is composed of a stack of B matrices for each node of the element:
	 * 	\f{eqnarray*}{
        	\textbf{B}  &=& \left[ \left[ \textbf{B}_{0} \right] \left[ \textbf{B}_{1} \right] ... \left[ \textbf{B}_{nNode}\right] \right]\\
						&=& \left[
		\begin{bmatrix}
			\partial_x N_0 	& 0 				& 0	\\
			0 				& \partial_y N_0  	& 0	\\
			0 				& 0 				& \partial_z N_0 \\
			\partial_y N_0 	& \partial_x N_0 	& 0\\
			\partial_z N_0 	& 0 				& \partial_x N_0 \\
			0 				& \partial_z N_0 	& \partial_y N_0
		\end{bmatrix}

		\begin{bmatrix}
			\partial_x N_1 	& 0 				& 0	\\
			0 				& \partial_y N_1  	& 0	\\
			0 				& 0 				& \partial_z N_1 \\
			\partial_y N_1 	& \partial_x N_1 	& 0\\
			\partial_z N_1 	& 0 				& \partial_x N_1 \\
			0 				& \partial_z N_1 	& \partial_y N_1
		\end{bmatrix}

		...

		\begin{bmatrix}
			\partial_x N_{nNode} 	& 0 					& 0	\\
			0 						& \partial_y N_{nNode}  & 0	\\
			0 						& 0 					& \partial_z N_{nNode} \\
			\partial_y N_{nNode} 	& \partial_x N_{nNode} 	& 0\\
			\partial_z N_{nNode} 	& 0 					& \partial_x N_{nNode} \\
			0 						& \partial_z N_{nNode} 	& \partial_y N_{nNode}
		\end{bmatrix}
		\right]
		\f}
	 * The elemental displacement matrix is extracted from the system displacement matrix via
	 * the function ShapeBase#constructElementalDisplacementMatrix. The displacement is calculated
	 * as the displacement of a node from its position at the end of last time step, \f$ u_{n}\f$ to the position
	 * at the current Newton-Raphson iteration \f$ u_{k}\f$. With the velocities (displacement per time step),
	 * and the \f$\textbf{B}\f$ matrix, velocity gradient tensor can be calculated through:
	 * \f{eqnarray*}{
	 	 	 	 \boldsymbol{l} & = \boldsymbol{B}  \boldsymbol{v_{n+1}}\nonumber \\
								& = \boldsymbol{B} \frac{{u_{n+1}^{k} - u_{n}}} {\delta t}.
		\f}
	 *
	 * Procedure:
	 * - construct the ElementalDisplacementMatrix.
	 */
	gsl_matrix* elementalDisplacementMatrix = constructElementalDisplacementMatrix(displacementPerDt);
	/**
	 * - Allocate the velocity gradient tensor in Voigt notation.
	 */
	gsl_matrix* l =  gsl_matrix_calloc(6,1);
	/**
	 * - calculate velocity gradient tensor.
	 */
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, B, elementalDisplacementMatrix,0.0, l);
    /**
	 * - free allocated memory.
	 */
	gsl_matrix_free(elementalDisplacementMatrix);
	//displayMatrix(l,"l_forForceCalc");
    /**
	 * - return velocity gradient tensor.
	 */
	return l;
}

gsl_matrix* ShapeBase::calculateRateOfDeformationTensor(gsl_matrix* l){
	/**
	 * Inputs:
	 * -# The velocity gradient tensor given in Voigt notation (6 x 1).
	 *
	 * Output:
	 * -# Rate of deformation tensor (3 x 3).
	 *
	 * This function primarily rearranges the elements of the
	 * velocity gradient tensor given in Voigt notation, to from the rate of deformation tensor
	 *
	 * Procedure:
	 * - Allocate the memory for rate of deformation tensor
	 */
	gsl_matrix* d =  gsl_matrix_calloc(3,3);
	/**
	 * - Write the terms of velocity gradient tensor into rate of deformation tensor
	 */
	gsl_matrix_set(d,0,0,gsl_matrix_get(l,0,0) );
	gsl_matrix_set(d,1,1,gsl_matrix_get(l,1,0) );
	gsl_matrix_set(d,2,2,gsl_matrix_get(l,2,0) );
	gsl_matrix_set(d,0,1,0.5*gsl_matrix_get(l,3,0));
	gsl_matrix_set(d,2,1,0.5*gsl_matrix_get(l,4,0));
	gsl_matrix_set(d,0,2,0.5*gsl_matrix_get(l,5,0));
	gsl_matrix_set(d,1,0,0.5*gsl_matrix_get(l,3,0));
	gsl_matrix_set(d,1,2,0.5*gsl_matrix_get(l,4,0));
	gsl_matrix_set(d,2,0,0.5*gsl_matrix_get(l,5,0));
	/**
	 * - Return rate of deformation tensor
	 */
	return d;
}

void ShapeBase::calculateViscousStress(gsl_matrix* d, gsl_matrix* viscousStress){
	/**
	 * Inputs:
	 * -# The rate of deformation matrix (ShapeBase#nDim x ShapeBase#nDim).
	 * -# the viscous stress of current Gauss point, the result will be written on this matrix
	 *
	 * This function will calculate the internal viscous stress of the element using rate of deformation matrix
	 * and ShapeBase#internalViscosity, \f$\eta\f$ , via:
	 * \f[\sigma^{v} = \eta  \textbf{d} \f]
	 *
	 * Procedure:
	 * - Copy rate of deformation tensor over to viscous stress tensor.
	 *
	 */
	createMatrixCopy(viscousStress, d);
	/**
	 *  - Scale with the internal viscosity to obtain viscous stress tensor
	 *
	 */
	gsl_matrix_scale(viscousStress,internalViscosity);
}

gsl_matrix* ShapeBase::constructElementalDisplacementMatrix(gsl_matrix* displacement){
	/**
	 * Inputs:
	 * -# The displacement matrix (ShapeBase#nDim x Simulation#nNodes, 1).
	 *
	 * This function calculates the elemental displacement matrix from the
	 * displacement matrix of the whole system, given as input. In
	 * current usage, under normal circumstances, the input matrix is displacement
	 * divided by time step. The displacement is calculated by Simulation#calculateDisplacementMatrix
	 * Both matrices, the displacement matrix of the whole system and the
	 * elemental displacement matrix are in vector form:
	 *
   	    \f$ displacement =
   	    \begin{bmatrix}
			\Delta x_{0}\\
			\Delta y_{0}\\
			\Delta z_{0}\\
			... ,\\
			\Delta x_{N}\\
			\Delta y_{N}\\
			\Delta z_{N}
			\end{bmatrix}
   	    \f$

	 */
	gsl_matrix* elementalDisplacementMatrix = gsl_matrix_calloc(nDim*nNodes,1);
    for (size_t i=0; i<nNodes; ++i){
		int index = NodeIds[i];
        for (size_t j=0; j<nDim; ++j){
			double value = gsl_matrix_get(displacement,index*nDim+j,0);
			gsl_matrix_set(elementalDisplacementMatrix,i*nDim+j,0,value);
		}
	}
	return elementalDisplacementMatrix;
}

void ShapeBase::calculateViscousForces(gsl_matrix*  gv, gsl_matrix*  BTdetFdetdXde, gsl_matrix* viscousStress){
	/**
	 * Inputs:
	 * -# Elemental matrix for internal viscous forces (ShapeBase#nDim x ShapeBase#nNodes, 1)
	 * the resulting forces will be written on this matrix
	 * -# Transpose of elemental B matrix, multiplied by the determinant
	 * of the deformation gradient, \f$\textbf{F}\f$, and the determinant of \f$ \delta \textbf{X}/\delta \boldsymbol{\xi}\f$
	 * -#  The viscous stresses calculated in ShapeBase#calculateViscousStress
	 *
	 * This function will calculate the elemental viscous forces from viscous stress, via:
	 * 	\f{eqnarray*}{
        	\textbf{g}^v &=& \int_{V} \textbf{B}^{T} \sigma^{v} dV \\
          	  	  	  	 &=& det(\textbf{F})det\left( \frac{\delta \textbf{X} }{\delta \boldsymbol{\xi}} \right) \textbf{B}^{T} \sigma^{v}
		\f}
	 * Procedure:
	 * - Allocate the memory for stress in Voigt notation
	 * */
	gsl_matrix* compactStress =  gsl_matrix_calloc(6,1);
	/**
	 * - Write stress in Voigt notation
	 */
    gsl_matrix_set(compactStress,0,0,gsl_matrix_get(viscousStress,0,0));
    gsl_matrix_set(compactStress,1,0,gsl_matrix_get(viscousStress,1,1));
    gsl_matrix_set(compactStress,2,0,gsl_matrix_get(viscousStress,2,2));
    gsl_matrix_set(compactStress,3,0,gsl_matrix_get(viscousStress,0,1));
    gsl_matrix_set(compactStress,4,0,gsl_matrix_get(viscousStress,2,1));
    gsl_matrix_set(compactStress,5,0,gsl_matrix_get(viscousStress,0,2));
	/**
	 * - Calculate nodal forces
	 */
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, BTdetFdetdXde, compactStress,0.0, gv);
	/**
	 * - Free memory
	 */
    gsl_matrix_free(compactStress);
}

void ShapeBase::calculateImplicitKElastic(){
    //std::cout<<"calculating implicit K elastic for element: "<<Id<<std::endl;
    int dim = nDim;
    int n = nNodes;
    if (IsAblated){
    	gsl_matrix_set_zero(TriPointKe);
    }
    else{
    	//calculating Kelastic in a 6 point gaussian:
		gsl_matrix* currK = gsl_matrix_calloc(dim*n,dim*n);
		gsl_matrix_set_zero(TriPointKe);
        for (size_t iter =0; iter<numberOfGaussPoints;++iter){
			gsl_matrix_set_zero(currK);
			calculateElasticKIntegral1(currK,iter);
			calculateElasticKIntegral2(currK,iter);
			gsl_matrix_scale(currK,gaussWeights[iter]);
			gsl_matrix_add(TriPointKe, currK);
		}
	    gsl_matrix_free(currK);
    }
}

void ShapeBase::calculateImplicitKViscous(gsl_matrix* displacementPerDt, double dt){
    /**	The viscous part of the elemental Jacobian will be calculated over the sum of two integrals, via
     * ShapeBase#calculateViscousKIntegral1 and ShapeBase#calculateViscousKIntegral2.
     * The inputs are: \n
     * 1) the displacement between the current N-R iteration positions \f$ \mathbf{u}_k^{t+\Delta t} \f$ and
     * the nodal positions at teh end of previous step \f$ \mathbf{u}_n^{t} \f$. \n
     * 2) the time step, \f$ \Delta t \f$. \n
     * The calculation is carried out over all Gauss Points.
     */
	//the object has access to all the necessary matrices to carry out the calculation.
	int dim = nDim;
	int n = nNodes;
	if (IsAblated || internalViscosity == 0){
	    //This is for efficiency, I do not calculate if the
	    //current element is laser ablated
		gsl_matrix_set_zero(TriPointKv);
	}
	else{
        //calculating Kviscous in an average of all Gauss Points.
		//assign a temporary matix, 18 x 18 for a prism (6 nodes, 3D).
		gsl_matrix* currK = gsl_matrix_calloc(dim*n,dim*n);
		//set the temporary matrix to zero.
		gsl_matrix_set_zero(TriPointKv);
		//the weights of the gauss points

		//define the matrices for velocity gradient, and the
		//term in parentheses in calculation of the first integral.
		//These will be calculated once per Gauss point for the element.
		gsl_matrix* velocityGradient = gsl_matrix_calloc(dim,dim);
		gsl_matrix* paranthesisTermForKv1 =  gsl_matrix_calloc(dim,dim);
	    //loop over Gauss points
        for (size_t iter =0; iter<numberOfGaussPoints;++iter){
			//set the temporary matrix to zero.
			gsl_matrix_set_zero(currK);
			//set the velocity gradient to zero.
			gsl_matrix_set_zero(velocityGradient);
			//set the parentheses term for first integral to zero.
			gsl_matrix_set_zero(paranthesisTermForKv1);
			//calculate the velocity gradient:
			calculateVelocityGradient(velocityGradient, displacementPerDt, iter);
			//calculate ( I / dt - velocityGradient) for first term:
			//set the term to identity:
			gsl_matrix_set_identity(paranthesisTermForKv1);
			//divide the term by dt to obtain I/dt:
			gsl_matrix_scale(paranthesisTermForKv1,1.0/dt);
			//substract velocity gradient to obtain ( I / dt - velocityGradient):
			gsl_matrix_sub(paranthesisTermForKv1,velocityGradient);
			//calculate the first integral:
			calculateViscousKIntegral1(currK, paranthesisTermForKv1, iter);
			//calculate the second integral:
			calculateViscousKIntegral2(currK, iter);
		    //scaling the resulting temporary matrix with Gauss point weight.
			gsl_matrix_scale(currK,gaussWeights[iter]);
		    //Adding the temporary matrix to the elemental Kviscous.
			gsl_matrix_add(TriPointKv, currK);
		}
		//free the memory allocated in this function
		gsl_matrix_free(currK);
		gsl_matrix_free(velocityGradient);
		gsl_matrix_free(paranthesisTermForKv1);
	}
}

void ShapeBase::writeKelasticToMainKatrix(gsl_matrix* K){
    for (size_t a=0; a<nNodes; ++a){
        for (size_t b=0; b<nNodes; ++b){
            int NodeId1 = NodeIds[a];
            int NodeId2 = NodeIds[b];
            NodeId1 *= nDim;
            NodeId2 *= nDim;
            for (size_t i=0; i<nDim; ++i){
                for (size_t j=0; j<nDim; ++j){
                    double valueij = gsl_matrix_get(K,NodeId1+i,NodeId2+j);
					valueij	+= gsl_matrix_get(TriPointKe,a*nDim+i,b*nDim+j);
                    gsl_matrix_set(K,NodeId1+i,NodeId2+j,valueij);
                    if (std::isnan(valueij)){
                		double tripointvalue = gsl_matrix_get(TriPointKe,a*nDim+i,b*nDim+j);
                        std::cout<<" element: "<<Id<<" K elastic dimention: "<<i<<" "<<j<<" is NaN after addition: "<<valueij<<" tri point value: "<<tripointvalue<<std::endl;
                	}
                }
            }
        }
    }
}

void ShapeBase::writeKviscousToMainKatrix(gsl_matrix* K){
	if (internalViscosity != 0){
        for (size_t a=0; a<nNodes; ++a){
            for (size_t b=0; b<nNodes; ++b){
				int NodeId1 = NodeIds[a];
				int NodeId2 = NodeIds[b];
				NodeId1 *= nDim;
				NodeId2 *= nDim;
                for (size_t i=0; i<nDim; ++i){
                    for (size_t j=0; j<nDim; ++j){
						double valueij = gsl_matrix_get(K,NodeId1+i,NodeId2+j);
						valueij	+= gsl_matrix_get(TriPointKv,a*nDim+i,b*nDim+j);
						gsl_matrix_set(K,NodeId1+i,NodeId2+j,valueij);
                        if (std::isnan(valueij)){
	                		double tripointvalue = gsl_matrix_get(TriPointKv,a*nDim+i,b*nDim+j);
                            std::cout<<" element: "<<Id<<" K viscous dimention: "<<i<<" "<<j<<" is NaN after addition: "<<valueij<<" tri point value: "<<tripointvalue<<std::endl;
	                	}
					}
				}
			}
		}
    }
}

void	ShapeBase::calculateForceFromStress(int nodeId, gsl_matrix* Externalstress, gsl_matrix *ExternalNodalForces){
    gsl_matrix_set_zero(ExternalNodalForces);
    int nodeNo = 0;
    for (size_t i=0; i<nNodes; i++){
        if (NodeIds[i] == nodeId){
            nodeNo = i;
            break;
        }
    }
    for (size_t pointNo = 0; pointNo<3; pointNo++){
        gsl_matrix* BaT = gsl_matrix_calloc(nDim,6);
        gsl_matrix* Bb = gsl_matrix_calloc(6,nDim);
        gsl_matrix* B = Bmatrices[pointNo];
        consturctBaTBb(B, BaT,Bb,nodeNo,0);
        gsl_matrix* NodeForces = gsl_matrix_calloc(3,1);
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, BaT, Externalstress,0.0, NodeForces);
        gsl_matrix_scale(NodeForces,1.0/3.0);
        gsl_matrix_scale(NodeForces,detFs[pointNo]);
        gsl_matrix_scale(NodeForces,detdXdes[pointNo]);
        gsl_matrix_add(ExternalNodalForces,NodeForces);
        gsl_matrix_free(BaT);
        gsl_matrix_free(Bb);
        gsl_matrix_free(NodeForces);
    }
}
void	ShapeBase::calculateElasticKIntegral1(gsl_matrix* currElementalK,int pointNo){
    gsl_matrix* invJShFuncDerS = invJShapeFuncDerStack[pointNo];
    double detF = detFs[pointNo];
    double detdXde = detdXdes[pointNo];
    gsl_matrix* Fe = FeMatrices[pointNo];
	double detFe = determinant3by3Matrix(Fe);
    for (size_t a =0; a<nNodes; ++a){
        for (size_t b=0; b<nNodes; ++b){
            gsl_matrix* Keab = gsl_matrix_calloc(3,3);
            double DNa[3] = {0.0,0.0,0.0};
            double DNb[3] = {0.0,0.0,0.0};

            for (size_t i=0;i<nDim;++i){
                // original version: DNa[i] = gsl_matrix_get(invJShFuncDerS,i,nDim*a);
                // original version: DNb[i] = gsl_matrix_get(invJShFuncDerS,i,nDim*b);
            	DNa[i] = gsl_matrix_get(invJShFuncDerS,i,nDim*a);
                DNb[i] = gsl_matrix_get(invJShFuncDerS,i,nDim*b);
            }
            //writing Kab:
            for (size_t i = 0 ; i<nDim; ++i){
                for (size_t k=0; k<nDim; ++k){
                    double value = 0;
                    //the sum over j,l,I,J,K,L, to get Kab(i,k):
                    for (size_t j = 0; j<nDim; ++j){
                        for (size_t l=0; l<nDim; ++l){
                            for (size_t I=0; I<nDim; ++I){
                                for (size_t J=0; J<nDim; ++J){
                                    for (size_t K=0; K<nDim; ++K){
                                        for (size_t L=0; L<nDim; ++L){
                                            value += (gsl_matrix_get(Fe,i,I)*gsl_matrix_get(Fe,j,J)*gsl_matrix_get(Fe,k,K)*gsl_matrix_get(Fe,l,L)*D81[pointNo][I][J][K][L]*DNb[l]*DNa[j]);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    value *= detF*detdXde;
                    value /= detFe;
                    value += gsl_matrix_get(Keab,i,k);
                    gsl_matrix_set(Keab,i,k,value);
                }
            }
            //now I have Kab for current gauss point, I need to write in into currK:
            for (size_t i=0; i<nDim; ++i){
                for (size_t j=0; j<nDim; ++j){
                    double value = gsl_matrix_get(currElementalK,a*nDim+i, b*nDim+j);
                    value += gsl_matrix_get(Keab,i, j);
                    gsl_matrix_set(currElementalK,a*nDim+i, b*nDim+j,value);
                }
            }
            gsl_matrix_free(Keab);
        }
    }
}


void	ShapeBase::calculateElasticKIntegral2(gsl_matrix* currElementalK,int pointNo){
    gsl_matrix* invJShFuncDerS = invJShapeFuncDerStack[pointNo];
    gsl_matrix* Stress = elasticStress[pointNo];
    double detF = detFs[pointNo];
    double detdXde = detdXdes[pointNo];

    gsl_matrix* DNaT = gsl_matrix_calloc(1,nDim);
    gsl_matrix* DNb = gsl_matrix_calloc(nDim,1);
    gsl_matrix* Keab2 = gsl_matrix_calloc(1,1);
    for (size_t a =0; a<nNodes; ++a){
        for (size_t b=0; b<nNodes; ++b){
            for (size_t i=0;i<nDim;++i){
                gsl_matrix_set(DNaT,0,i,gsl_matrix_get(invJShFuncDerS,i,nDim*a));
                gsl_matrix_set(DNb,i,0,gsl_matrix_get(invJShFuncDerS,i,nDim*b));
            }
            gsl_matrix* tmp1 = gsl_matrix_calloc(1,nDim);
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, DNaT, Stress,0.0, tmp1);
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, tmp1, DNb,0.0, Keab2);
            double value = gsl_matrix_get(Keab2,0,0)*detF*detdXde;
            for (size_t i=0; i<nDim; ++i){
                int index1 = a*nDim+i;
                int index2 = b*nDim+i;
                double addedValue = gsl_matrix_get(currElementalK,index1,index2) + value;//adding the calculated value to current K matirx
                gsl_matrix_set(currElementalK,index1,index2,addedValue);
            }
            gsl_matrix_free(tmp1);
        }
    }
    gsl_matrix_free(DNaT);
    gsl_matrix_free(DNb);
    gsl_matrix_free(Keab2);
}

void ShapeBase::calculateVelocityGradient( gsl_matrix* velocityGradient, gsl_matrix* displacementPerDt, int pointNo){
	gsl_matrix* invJShFuncDerS = invJShapeFuncDerStack[pointNo];
    for (size_t c=0; c<nNodes; ++c){
		gsl_matrix* delVc = gsl_matrix_calloc(nDim,nDim);
		//get \DelNc^T
		gsl_matrix* DNc = gsl_matrix_calloc(nDim,1);
        for (size_t i=0;i<nDim;++i){
			gsl_matrix_set(DNc,i,0,gsl_matrix_get(invJShFuncDerS,i,nDim*c));
		}
		//calculate velocity of node c:
		gsl_matrix* vc = gsl_matrix_calloc(nDim,1);
		int id = NodeIds[c];
		gsl_matrix_set(vc,0,0,gsl_matrix_get(displacementPerDt,id*nDim,0));
		gsl_matrix_set(vc,1,0,gsl_matrix_get(displacementPerDt,id*nDim+1,0));
		gsl_matrix_set(vc,2,0,gsl_matrix_get(displacementPerDt,id*nDim+2,0));
		//calculate vc *DNc^T:
		calculateOuterProduct(vc,DNc,delVc);
		//add the nodal calculation to velocity gradient:
		gsl_matrix_add(velocityGradient, delVc);
		//free memory allocated in this loop:
		gsl_matrix_free(delVc);
		gsl_matrix_free(DNc);
		gsl_matrix_free(vc);
	}
}

void ShapeBase::calculateViscousKIntegral1(gsl_matrix* currElementalK, gsl_matrix* paranthesisTermForKv1, int pointNo){
	gsl_matrix* BaT = gsl_matrix_calloc(nDim,6);
    gsl_matrix* Bb = gsl_matrix_calloc(6,nDim);
    gsl_matrix* BaTBb = gsl_matrix_calloc(nDim,nDim);
    gsl_matrix* KvabInt1 = gsl_matrix_calloc(nDim,nDim); //First integral in calculation of Kv for nodes a and b
    gsl_matrix* B = Bmatrices[pointNo];
    for (size_t a=0;a<nNodes;++a){
        for (size_t b=0; b<nNodes; ++b){
    		consturctBaTBb(B, BaT,Bb,a,b);
    		//Bb matrix should have the last three rows as 0.5, this matrix stems from the
    		//"d" matrix, as opposed to viscous stresses, therefore the definition
    		//should have 0.5 on off diagonal terms, therefore the last three rows of Bb.
            for (size_t i=3; i<6; ++i){
                for (size_t j=0; j<nDim; ++j){
    				double value = gsl_matrix_get(Bb, i,j);
    				gsl_matrix_set(Bb,i,j,0.5*value);
    			}
    		}
    		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, BaT, Bb,0.0, BaTBb);
    		//the paranthesis term is: ( I / dt - velocityGradient)
    		//calculate all multiplication: BaT*Bb* (I/dt - \Del v):
			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, BaTBb, paranthesisTermForKv1,0.0, KvabInt1);
			//volume integration:
    	    gsl_matrix_scale(KvabInt1,detFs[pointNo]);
    	    gsl_matrix_scale(KvabInt1,detdXdes[pointNo]);
    	    //scaling by viscosity:
    	    gsl_matrix_scale(KvabInt1,internalViscosity);

            for (size_t i=0; i<nDim; ++i){
                for (size_t j=0; j<nDim; ++j){
        	    	int index2 = a*nDim+i;
        	    	int index1 = b*nDim+j;
        	    	double value = gsl_matrix_get(KvabInt1,i,j);
        	    	double addedValue = gsl_matrix_get(currElementalK,index1,index2) + value; //adding the calculated value to current K matrix
        	    	gsl_matrix_set(currElementalK,index1,index2,addedValue);
        	    }
			}
    	}
    }
    gsl_matrix_free(BaT);
    gsl_matrix_free(Bb);
    gsl_matrix_free(BaTBb);
    gsl_matrix_free(KvabInt1);
}

void ShapeBase::calculateViscousKIntegral2(gsl_matrix* currElementalK,int pointNo){
    gsl_matrix* invJShFuncDerS = invJShapeFuncDerStack[pointNo];
    gsl_matrix* Stress = viscousStress[pointNo];
    gsl_matrix* DNa = gsl_matrix_calloc(nDim,1);
    gsl_matrix* DNb = gsl_matrix_calloc(nDim,1);
    gsl_matrix* KvabInt2 = gsl_matrix_calloc(nDim,nDim);
    for (size_t a =0; a<nNodes; ++a){
        for (size_t b=0; b<nNodes; ++b){
            for (size_t i=0;i<nDim;++i){
                gsl_matrix_set(DNa,i,0,gsl_matrix_get(invJShFuncDerS,i,nDim*a));
                gsl_matrix_set(DNb,i,0,gsl_matrix_get(invJShFuncDerS,i,nDim*b));
            }
            gsl_matrix* DNaDNbOuterProduct = gsl_matrix_calloc(nDim, nDim);
            gsl_matrix* DNbDNaOuterProduct = gsl_matrix_calloc(nDim, nDim);
            calculateOuterProduct(DNa, DNb, DNaDNbOuterProduct);
            calculateOuterProduct(DNb, DNa, DNbDNaOuterProduct);
            gsl_matrix* paranthesisTerm = gsl_matrix_calloc(nDim, nDim);
            gsl_matrix_memcpy(paranthesisTerm,DNaDNbOuterProduct);
            gsl_matrix_sub(paranthesisTerm,DNbDNaOuterProduct);
            //gsl_matrix_add(paranthesisTerm,DNbDNaOuterProduct);
            //gsl_matrix_scale(paranthesisTerm, 0.5);
			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Stress, paranthesisTerm,0.0, KvabInt2);
			/*if (Id == 0){
				//displayMatrix(paranthesisTerm,"paranthesisTerm");
				displayMatrix(Stress,"Stress");
				//displayMatrix(KvabInt2,"KvabInt2-beforeVolumeIntegration");
			}*/
			gsl_matrix_scale(KvabInt2, detFs[pointNo]);
			gsl_matrix_scale(KvabInt2, detdXdes[pointNo]);
            for (size_t i=0; i<nDim; ++i){
                for (size_t j=0; j<nDim; ++j){
					int index2 = a*nDim+i;
					int index1 = b*nDim+j;
					double value = gsl_matrix_get(KvabInt2,i,j);
					double addedValue = gsl_matrix_get(currElementalK,index1,index2) + value; //adding the calculated value to current K matrix
					gsl_matrix_set(currElementalK,index1,index2,addedValue);
				}
			}
            gsl_matrix_free(DNaDNbOuterProduct);
            gsl_matrix_free(DNbDNaOuterProduct);
            gsl_matrix_free(paranthesisTerm);
        }
    }
    gsl_matrix_free(DNa);
    gsl_matrix_free(DNb);
    gsl_matrix_free(KvabInt2);
}

void	ShapeBase::calculateOuterProduct(gsl_matrix* a, gsl_matrix* b, gsl_matrix* outerProduct){
    /** The outer product is defined as: \n
     * \f$ \mathbf{a} outer \mathbf{b} = \mathbf{a} \mathbf{b}^{T}
     */
    size_t size2 = a->size2;
    if (b->size2 != size2){
        std::cerr<<"matrix dimension mismatch in outer product calculation"<<std::endl;
	}
	gsl_blas_dgemm (CblasNoTrans, CblasTrans,1.0, a, b,0.0, outerProduct);
}

gsl_matrix*	ShapeBase::calculateSymmetricisedTensorProduct(gsl_matrix* a, gsl_matrix* b){
	int size1 = a->size1;
	int size2 = a->size2;
	if ((int) b->size1 != size1){
        std::cerr<<"matrix dimension mismatch in symmetricised outer product calculation"<<std::endl;
	}
	if ((int) b->size2 != size2){
        std::cerr<<"matrix dimension mismatch in symmetricised outer product calculation"<<std::endl;
	}
	//calculating individual outer products a x b = a bT
	gsl_matrix* abOuterProduct = gsl_matrix_calloc(size1,size1);
	calculateOuterProduct(a,b,abOuterProduct);
	gsl_matrix* baOuterProduct = gsl_matrix_calloc(size1,size1);
	calculateOuterProduct(b,a,abOuterProduct);
	//calculating the averaged contraction term:
	gsl_matrix* averagedContraction = gsl_matrix_calloc(size1,size1);
	gsl_matrix_add(averagedContraction,abOuterProduct);
	gsl_matrix_add(averagedContraction,baOuterProduct);
	gsl_matrix_scale(averagedContraction,0.5);
	gsl_matrix_free(abOuterProduct);
	gsl_matrix_free(baOuterProduct);
	return averagedContraction;
}

void ShapeBase::consturctBaTBb(gsl_matrix* B, gsl_matrix* BaT, gsl_matrix* Bb, int a, int b){
    for (size_t i=0; i<6; ++i){
        for (size_t j=0; j<nDim; ++j){
            gsl_matrix_set(BaT,j,i,gsl_matrix_get(B,i,a*nDim+j)); //transpose of Ba
            gsl_matrix_set(Bb,i,j,gsl_matrix_get(B,i,b*nDim+j)); //Bb
        }
    }
}

void	ShapeBase::fillNodeNeighbourhood(const std::vector<std::unique_ptr<Node>>& Nodes){
    for (size_t i = 0; i<nNodes; ++i){
        for (size_t j = 0; j<nNodes; ++j){
			if ( i !=j ){
                size_t n = Nodes[NodeIds[i]]->immediateNeigs.size();
				bool alreadyOnList = false;
                for (size_t k=0; k<n; ++k){
					if (NodeIds[j] == Nodes[NodeIds[i]]->immediateNeigs[k]){
						alreadyOnList = true;
						break;
					}
				}
				if (!alreadyOnList){
					Nodes[NodeIds[i]]->immediateNeigs.push_back(NodeIds[j]);
				}
			}
		}
	}
}

void	ShapeBase::updatePositions(const std::vector<std::unique_ptr<Node> > &Nodes){
    for (size_t i = 0; i<nNodes; ++i){
        for (size_t j = 0; j<nDim; ++j){
			Positions[i][j] = Nodes[NodeIds[i]]->Position[j];
		}
	}
}

void 	ShapeBase::setGrowthRate(double dt, double rx, double ry, double rz){
	GrowthRate[0] = exp(rx*dt);
	GrowthRate[1] = exp(ry*dt);
	GrowthRate[2] = exp(rz*dt);
}

void 	ShapeBase::setGrowthRateViaInputTimeMultipliedMagnitude(double x, double y, double z){
    GrowthRate[0] = x;
    GrowthRate[1] = y;
    GrowthRate[2] = z;
}

void 	ShapeBase::updateGrowthIncrementFromRate(){
	gsl_matrix_set_identity(growthIncrement);
	gsl_matrix_set(growthIncrement,0,0,GrowthRate[0]);
	gsl_matrix_set(growthIncrement,1,1,GrowthRate[1]);
	gsl_matrix_set(growthIncrement,2,2,GrowthRate[2]);
}

void 	ShapeBase::calculatePrincipalStrains2D(double& e1, double &e2,  double &e3, gsl_matrix* eigenVec){
	gsl_matrix* Strain2D = gsl_matrix_calloc(2,2);
	gsl_matrix_set(Strain2D,0,0, gsl_matrix_get(Strain,0,0));
	gsl_matrix_set(Strain2D,1,1, gsl_matrix_get(Strain,1,0));
	gsl_matrix_set(Strain2D,0,1, 0.5 * gsl_matrix_get(Strain,3,0));
	gsl_matrix_set(Strain2D,1,0, 0.5 * gsl_matrix_get(Strain,3,0));
	gsl_vector* eigenValues = gsl_vector_calloc(2);
	gsl_matrix* eigenVec2D = gsl_matrix_calloc(2,2);
	gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(2);
	gsl_eigen_symmv(Strain2D, eigenValues, eigenVec2D, w);
	gsl_eigen_symmv_free(w);
	gsl_eigen_symmv_sort(eigenValues, eigenVec2D, GSL_EIGEN_SORT_ABS_ASC);
	e1 = gsl_vector_get(eigenValues,0);
	e2 = gsl_vector_get(eigenValues,1);
	e3 = 0;
	gsl_matrix_set_identity(eigenVec);
    for (size_t i=0; i<2; ++i){
        for (size_t j=0; j<2; ++j){
			gsl_matrix_set(eigenVec,i,j,gsl_matrix_get(eigenVec2D,i,j));
		}
	}
	gsl_vector_free(eigenValues);
	gsl_matrix_free(eigenVec2D);
	gsl_matrix_free(Strain2D);
}

void 	ShapeBase::calculatePrincipalStrains3D(double& e1, double &e2,  double &e3, gsl_matrix* eigenVec){
	gsl_matrix* Strain3D = gsl_matrix_calloc(3,3);
	gsl_matrix_set(Strain3D,0,0, gsl_matrix_get(Strain,0,0));
	gsl_matrix_set(Strain3D,1,1, gsl_matrix_get(Strain,1,0));
	gsl_matrix_set(Strain3D,0,1, 0.5 * gsl_matrix_get(Strain,3,0));
	gsl_matrix_set(Strain3D,1,0, 0.5 * gsl_matrix_get(Strain,3,0));
	gsl_matrix_set(Strain3D,2,2, gsl_matrix_get(Strain,2,0));
	gsl_matrix_set(Strain3D,2,1, 0.5 * gsl_matrix_get(Strain,4,0));
	gsl_matrix_set(Strain3D,1,2, 0.5 * gsl_matrix_get(Strain,4,0));
	gsl_matrix_set(Strain3D,0,2, 0.5 * gsl_matrix_get(Strain,5,0));
	gsl_matrix_set(Strain3D,2,0, 0.5 * gsl_matrix_get(Strain,5,0));
	gsl_vector* eigenValues = gsl_vector_calloc(3);
	gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(3);
	gsl_eigen_symmv(Strain3D, eigenValues, eigenVec, w);
	gsl_eigen_symmv_free(w);
	gsl_eigen_symmv_sort(eigenValues, eigenVec, GSL_EIGEN_SORT_ABS_ASC);
	e1 = gsl_vector_get(eigenValues,0);
	e2 = gsl_vector_get(eigenValues,1);
	e3 = gsl_vector_get(eigenValues,2);
	gsl_vector_free(eigenValues);
	gsl_matrix_free(Strain3D);
}

double  ShapeBase::calculateVolumeForInputShapeStructure(std::vector<std::array<double,3>> shapePositions, size_t nTriangularFaces, std::vector<std::array<int,3>>  triangularFaces, std::array<double,3> midPoint ){
	double totalVolume = 0;
    for (size_t i=0; i< nTriangularFaces; ++i){
		//calculateTetrahedraVolume:
        int node0 = triangularFaces[i][0];
        int node1 = triangularFaces[i][1];
        int node2 = triangularFaces[i][2];
        std::array<double,3> vec1 = {0.0};
        std::array<double,3> vec2 = {0.0};
        std::array<double,3> vecMid = {0.0};
        for (size_t j=0 ;j < nDim; ++j){
				vec1   [j] = shapePositions[node1][j] - shapePositions[node0][j];
				vec2   [j] = shapePositions[node2][j] - shapePositions[node0][j];
				vecMid [j] = midPoint[j] - shapePositions[node0][j];
		}
        std::array<double,3> baseVec = crossProduct3D(vec1,vec2);
        double normBaseVec = calculateMagnitudeVector3D (baseVec);
		double baseArea= normBaseVec/2;
		if (i == 0){
			ReferenceShape->BasalArea = baseArea;
		}
		double height = dotProduct3D(vecMid,baseVec) / normBaseVec;
		if (height <0){
			height *= (-1.0);
		}
		ReferenceShape->height = height;
		double currVolume = 1.0/3.0 *(height *baseArea);
		totalVolume += currVolume;
	}
	return totalVolume;
}

void 	ShapeBase::setShapeChangeInrementToIdentity(){
    gsl_matrix_set_identity(shapeChangeIncrement);
}

void 	ShapeBase::setShapeChangeRate(double x, double y, double z, double xy, double yz, double xz){
    ShapeChangeRate[0] = x;
    ShapeChangeRate[1] = y;
    ShapeChangeRate[2] = z;
    ShapeChangeRate[3] = xy;
    ShapeChangeRate[4] = yz;
    ShapeChangeRate[5] = xz;
}

//R Fginc R^T version:
void ShapeBase::calculateCurrentGrowthIncrement(gsl_matrix* resultingGrowthIncrement, double dt, double growthx, double growthy, double growthz, gsl_matrix* ShearAngleRotationMatrix){
    gsl_matrix_set_identity(resultingGrowthIncrement);
    gsl_matrix_set(resultingGrowthIncrement,0,0,exp(growthx*dt));
    gsl_matrix_set(resultingGrowthIncrement,1,1,exp(growthy*dt));
    gsl_matrix_set(resultingGrowthIncrement,2,2,exp(growthz*dt));
    gsl_matrix* temp = gsl_matrix_calloc(nDim,nDim);;
    //R * resultingGrowthIncrement
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, ShearAngleRotationMatrix, resultingGrowthIncrement, 0.0, temp);
    //R * resultingGrowthIncrement * R^T
    gsl_blas_dgemm (CblasNoTrans, CblasTrans,1.0, temp, ShearAngleRotationMatrix, 0.0, resultingGrowthIncrement);
    gsl_matrix_free(temp);

}

void 	ShapeBase::updateShapeChangeRate(double x, double y, double z, double xy, double yz, double xz){
    ShapeChangeRate[0] += x;
    ShapeChangeRate[1] += y;
    ShapeChangeRate[2] += z;
    ShapeChangeRate[3] += xy;
    ShapeChangeRate[4] += yz;
    ShapeChangeRate[5] += xz;
}

bool 	ShapeBase::InvertMatrix(gsl_matrix* input, gsl_matrix* inverse){
    // Define the dimension n of the matrix
    // and the signum s (for LU decomposition)
    int s;

    // Define all the used matrices
    gsl_permutation * perm = gsl_permutation_alloc (input->size1);

    // Make LU decomposition of matrix m
    gsl_linalg_LU_decomp (input, perm, &s);

    // Invert the matrix m
    gsl_linalg_LU_invert (input, perm, inverse);

    return true;
}

bool 	ShapeBase::InvertMatrix(boost::numeric::ublas::matrix<double>& input, boost::numeric::ublas::matrix<double>& inverse){
    //Matrix inversion routine.
    //Uses lu_factorize and lu_substitute in uBLAS to invert a matrix
    using namespace boost::numeric::ublas;
    typedef permutation_matrix<std::size_t> pmatrix;

    // create a working copy of the input
    matrix<double> A(input);

    // create a permutation matrix for the LU-factorization
    pmatrix pm(A.size1());

    // perform LU-factorization
    int res = lu_factorize(A, pm);
    if (res != 0)
        return false;

    // create identity matrix of "inverse"
    inverse.assign(identity_matrix<double> (A.size1()));

    // backsubstitute to get the inverse
    lu_substitute(A, pm, inverse);

    return true;
}

void	ShapeBase::crossProduct3D(gsl_vector* u, gsl_vector* v, gsl_vector* cross){
    gsl_vector_set(cross,0, ( gsl_vector_get(u,1)*gsl_vector_get(v,2) - gsl_vector_get(u,2)*gsl_vector_get(v,1) ) );
    gsl_vector_set(cross,1, ( gsl_vector_get(u,2)*gsl_vector_get(v,0) - gsl_vector_get(u,0)*gsl_vector_get(v,2) ) );
    gsl_vector_set(cross,2, ( gsl_vector_get(u,0)*gsl_vector_get(v,1) - gsl_vector_get(u,1)*gsl_vector_get(v,0) ) );
}

std::array<double,3> 	ShapeBase::crossProduct3D(std::array<double,3> u, std::array<double,3> v){
    std::array<double,3> cross = {0.0, 0.0, 0.0};
    cross[0] = u[1]*v[2] - u[2]*v[1];
    cross[1] = u[2]*v[0] - u[0]*v[2];
    cross[2] = u[0]*v[1] - u[1]*v[0];
    return cross;
}

double	ShapeBase::calculateMagnitudeVector3D(std::array<double,3> v){
    double mag = v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
    mag = pow(mag,0.5);
    return mag;
}

void	ShapeBase::normaliseVector3D(gsl_vector* v){
    double x = gsl_vector_get(v,0);
    double y = gsl_vector_get(v,1);
    double z = gsl_vector_get(v,2);
    double mag2 = x*x + y*y + z*z;
    if (fabs(mag2) > 1E-14 && fabs(mag2 - 1.0f) > 1E-14) {
        double mag = pow(mag2,0.5);
        gsl_vector_scale(v,1.0/mag);
    }
}

double ShapeBase::normaliseVector3D(std::array<double,3>& v){
    double x = v[0];
    double y = v[1];
    double z = v[2];
    double mag2 = x*x + y*y + z*z;
    if (fabs(mag2) > 1E-14 && fabs(mag2 - 1.0f) > 1E-14) {
        double mag = pow(mag2,0.5);
        v[0] /= mag;
        v[1] /= mag;
        v[2] /= mag;
        return mag;
    }
    else{
        return 0.0;
    }
}

double	ShapeBase::getNormVector3D(gsl_vector* v){
    double x = gsl_vector_get(v,0);
    double y = gsl_vector_get(v,1);
    double z = gsl_vector_get(v,2);
    double mag2 = x*x + y*y + z*z;
    return pow(mag2,0.5);
}

double 	ShapeBase::dotProduct3D(std::array<double,3>& u, std::array<double,3>& v){
    double dot = 0;
    dot = u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
    return dot;
}

void 	ShapeBase::displayMatrix(boost::numeric::ublas::matrix<double>& mat, std::string matname){
    size_t m = mat.size1();
    size_t n = mat.size2();
    std::cout<<matname<<": "<<std::endl;

    for (size_t i =0; i<m; i++){
        for (size_t j =0; j<n; j++){
            std::cout.precision(4);
            std::cout.width(6);
            std::cout<<mat(i,j)<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

void 	ShapeBase::displayMatrix(gsl_matrix* mat, std::string matname){
    size_t m = mat->size1;
    size_t n = mat->size2;
    std::cout<<matname<<": "<<std::endl;

    for (size_t i =0; i<m; i++){
        for (size_t j =0; j<n; j++){
            std::cout.precision(4);
            std::cout.width(6);
            std::cout<<gsl_matrix_get(mat,i,j)<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

void 	ShapeBase::displayMatrix(gsl_vector* mat, std::string matname){
    size_t m = mat->size;
    std::cout<<matname<<": "<<std::endl;

    for (size_t i =0; i<m; i++){
        std::cout.precision(4);
        std::cout.width(6);
        std::cout<<gsl_vector_get(mat,i)<<std::endl;
    }
    std::cout<<std::endl;
}

void 	ShapeBase::displayMatrix(boost::numeric::ublas::matrix<int>& mat, std::string matname){
    size_t m = mat.size1();
    size_t n = mat.size2();
    std::cout<<matname<<": "<<std::endl;

    for (size_t i =0; i<m; i++){
        for (size_t j =0; j<n; j++){
            std::cout.precision(4);
            std::cout.width(6);
            std::cout<<mat(i,j)<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

void	ShapeBase::displayMatrix(boost::numeric::ublas::vector<double>& vec, std::string matname){
    size_t m = vec.size();
    std::cout<<matname<<": "<<std::endl;
    for (size_t i =0; i<m; i++){
        std::cout.precision(4);
        std::cout.width(6);
        std::cout<<vec(i)<<" ";
    }
    std::cout<<std::endl;
}

void 	ShapeBase::assignVolumesToNodes(const std::vector<std::unique_ptr<Node>>& Nodes){
    for (size_t i=0; i<nNodes; i++){
        Nodes[NodeIds[i]]->mass += VolumePerNode;
    }
}
void 	ShapeBase::calculateViscositySurfaces(){
    if (elementHasExposedApicalSurface){
        calculateApicalArea();
    }
    if (elementHasExposedBasalSurface){
        calculateBasalArea();
    }
}

void ShapeBase::assignViscositySurfaceAreaToNodes(const std::vector<std::unique_ptr<Node> > &Nodes){
    if (elementHasExposedApicalSurface){
        for(size_t i=0;i<nSurfaceAreaNodeNumber; ++i){
            Nodes[NodeIds[exposedApicalSurfaceNodeIds[i]]]->viscositySurface+=ApicalArea/nSurfaceAreaNodeNumber;
        }
    }
    if (elementHasExposedBasalSurface){
        for(size_t i=0;i<nSurfaceAreaNodeNumber; ++i){
            Nodes[NodeIds[exposedBasalSurfaceNodeIds[i]]]->viscositySurface+=BasalArea/nSurfaceAreaNodeNumber;
        }
    }
}

void 	ShapeBase::calculateZProjectedAreas(){
    double Threshold = 1E-5;
    int id0 = 0, id1 = 1, id2 = 2; // this is correct for basal side, I will change it for apical calculation
    for (size_t tissueSide = 0; tissueSide<2; tissueSide++){
        if ( tissueSide == 1){
            //I am calculating basal area,
            id0 = 3;
            id1 = 4;
            id2 = 5;
        }
        double sideVec1[2];
        double sideVec2[2];
        double Side1 = 0.0;
        double Side2 = 0.0;
        double costet = 0.0;
        double Area = 0.0;
        for (size_t i = 0; i<2; ++i){
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
            double sintet = pow((1-costet*costet),0.5);
            Area = Side1* Side2 * sintet / 2.0;
        }
        if(tissueSide == 0){
            ZProjectedBasalArea = Area;
        }
        else{
            ZProjectedApicalArea = Area;
        }
    }
}


void 	ShapeBase::assignZProjectedAreas(const std::vector<std::unique_ptr<Node>>& Nodes){
    if (ShapeType == 1 ){ //only written for prisms
        for (size_t i=0; i<3; i++){
            Nodes[NodeIds[i]]->zProjectedArea +=ZProjectedBasalArea/3.0;
        }
        for (size_t i=3; i<6; i++){
            Nodes[NodeIds[i]]->zProjectedArea +=ZProjectedApicalArea/3.0;
        }
    }
}

void 	ShapeBase:: assignElementToConnectedNodes(const std::vector<std::unique_ptr<Node>>& Nodes){
    for (size_t i=0; i<nNodes; i++){
        Nodes[NodeIds[i]]->connectedElementIds.push_back(Id);
        double weightfFraction = (ReferenceShape->Volume/nNodes)/Nodes[NodeIds[i]]->mass;
        Nodes[NodeIds[i]]->connectedElementWeights.push_back(weightfFraction);
        if(tissueType == 2){
            Nodes[NodeIds[i]]->hasLateralElementOwner = true;
        }
    }
}

void ShapeBase::scaleGrowthIncrement(double multiplier){
    gsl_matrix_scale(growthIncrement,multiplier);
}

bool ShapeBase::areanyOfMyNodesAtCircumference(const std::vector<std::unique_ptr<Node>>& Nodes){
    bool thereIsNodeAtCircumference = false;
    for (size_t i=0; i< nNodes; ++i){
        if (Nodes[NodeIds[i]]->atCircumference){
            thereIsNodeAtCircumference = true;
            break;
        }
    }
    return thereIsNodeAtCircumference;
}
