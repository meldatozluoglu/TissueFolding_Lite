/*
 * GrowthFunctionTypes.h
 *
 *  Created on: 24 Apr 2015
 *      Author: melda
 */

#ifndef GROWTHFUNCTIONTYPES_H
#define GROWTHFUNCTIONTYPES_H

#include "GrowthFunctionBase.h"

class UniformGrowthFunction : public GrowthFunctionBase{
private:

public:
    std::array<double,3> GrowthRate; ///< The double array stating the uniform growth rate throughout the tissue. Growth rate is in 1/sec, format: [ DV axis (x), AP axis (y), and AB axis (z)]
	gsl_matrix* ShearAngleRotationMatrix; ///< The rotation  matrix for the orientation of the growth on x-y plane. This matrix is constructed through  UniformGrowthFunction#angle
	double angle;	///< The rotation angle for the orientation of the growth on x-y plane.

	UniformGrowthFunction(int id, int type, float initTime, float endTime, bool applyToColumnarLayer, bool applyToPeripodialMembrane, bool applyToBasalECM, bool applyToLateralECM, double DVGrowth, double APGrowth, double ABGrowth, double angle) : GrowthFunctionBase(id, type, initTime, endTime, applyToColumnarLayer, applyToPeripodialMembrane, applyToBasalECM, applyToLateralECM){
		/**
		 *  The first six parameters will be directed to the parent constructor, GrowthFunctionBase#GrowthFunctionBase. \n
		 *  doubles DVGrowth, APGrowth and ABGrowth will set the  UniformGrowthFunction#GrowthRate, in the given order.
		 */
		this -> angle = angle;
		ShearAngleRotationMatrix = gsl_matrix_calloc(3,3);
		gsl_matrix_set_identity(ShearAngleRotationMatrix);
		double c = cos(angle);
		double s = sin(angle);
		gsl_matrix_set(ShearAngleRotationMatrix,0,0,  c );
		gsl_matrix_set(ShearAngleRotationMatrix,0,1, -1.0*s);
		gsl_matrix_set(ShearAngleRotationMatrix,0,2,  0.0);
		gsl_matrix_set(ShearAngleRotationMatrix,1,0,  s);
		gsl_matrix_set(ShearAngleRotationMatrix,1,1,  c);
		gsl_matrix_set(ShearAngleRotationMatrix,1,2,  0.0);
		gsl_matrix_set(ShearAngleRotationMatrix,2,0,  0.0);
		gsl_matrix_set(ShearAngleRotationMatrix,2,1,  0.0);
		gsl_matrix_set(ShearAngleRotationMatrix,2,2,  1.0);
		GrowthRate[0] = DVGrowth;	//x
		GrowthRate[1] = APGrowth;	//y
		GrowthRate[2] = ABGrowth;	//z
	}///< The constructor of UniformGrowthFunction

    ~UniformGrowthFunction(){}

    std::array<double,3> getGrowthRate(){
		/**
         *  This function will return the UniformGrowthFunction#GrowthRate of the current growth function.
		 */
        std::array<double,3> maxValues {0};
        maxValues[0] = GrowthRate[0];
		maxValues[1] = GrowthRate[1];
		maxValues[2] = GrowthRate[2];
        return maxValues;
	} ///< The function is to get the 3D growth rate of the current growth function.

    void setGrowtRate(double ex, double ey, double ez){
		/**
		 *  This function will set the UniformGrowthFunction#GrowthRate of the current growth function to the input values
		 *  The parameters are in the order [ DV axis (x), AP axis (y), and AB axis (z)].
		 */
		GrowthRate[0] = ex;
		GrowthRate[1] = ey;
		GrowthRate[2] = ez;
	}///< The function is to set the 3D growth rate of the current growth function.

	gsl_matrix* getShearAngleRotationMatrix(){
		return ShearAngleRotationMatrix;
	}///< The function return the matrix pointer to the oriented growth rotation matrix.

	double getShearAngle(){
		return angle;
	}///< The function returns the oriented growth angle in radians (double).

    void writeSummary(std::ofstream &saveFileSimulationSummary,double dt){
		/**
		 *  This function will write the UniformGrowthFunction details into the simulation summary file, provided as the first input.
		 *  Time step (dt) of the simulation is provided as second input, to report the growth rates per hour.
		 *  The output should look like: \n
		 *			Growth Type:  Uniform (1)
		 *			Initial time(sec): UniformGrowthFunction#initTime	FinalTime time(sec): UniformGrowthFunction#endTime	GrowthRate(fraction/hr): DVGrowth(in 1/hr)  APGrowth(in 1/hr) ABGrowth(in 1/hr)
		 */
        saveFileSimulationSummary<<"Growth Type:  Uniform (1)"<<std::endl;
		saveFileSimulationSummary<<"	Initial time(sec): ";
		saveFileSimulationSummary<<initTime;
		saveFileSimulationSummary<<"	FinalTime time(sec): ";
		saveFileSimulationSummary<<endTime;
		saveFileSimulationSummary<<"	GrowthRate(fraction/hr): ";
		saveFileSimulationSummary<<GrowthRate[0]/dt*3600.0;
		saveFileSimulationSummary<<"  ";
		saveFileSimulationSummary<<GrowthRate[1]/dt*3600.0;
		saveFileSimulationSummary<<"  ";
		saveFileSimulationSummary<<GrowthRate[2]/dt*3600.0;
		saveFileSimulationSummary<<"	Growth Angle(degrees): ";
		saveFileSimulationSummary<<angle/M_PI*180.0;
        saveFileSimulationSummary<<std::endl;
	}///< The function is to write the growth function summary to simulation summary file
};

class UniformShapeChangeFunction : public UniformGrowthFunction{
private:

public:
	int ShapeChangeType;
    UniformShapeChangeFunction(int id, int type, float initTime, float endTime, bool applyToColumnarLayer, bool applyToPeripodialMembrane, bool applyToBasalECM, bool applyToLateralECM, int ShapeChangeType, double ShapeChangeRate) : UniformGrowthFunction( id,  type,  initTime,  endTime,  applyToColumnarLayer,  applyToPeripodialMembrane, applyToBasalECM,  applyToLateralECM, 0.0,  0.0,  0.0, 0.0){
		/**
		 *  Forst six parameters will be directed to the parent constructor, UniformGrowthFunction#UniformGrowthFunction.
		 *  The growth rates in Dv, AB and AP will be fed as 0 to the parent constructor. \n
		 *  The parameter ShapeChangeType will define the type of the shape change as 1: ColumnarToSquamous change, 2: Apical shrinkage, and 2: basal shrinkage.
		 *  The paremeter ShapeChangeRate will define the rate of defined shape change. The constructor will calculate the corrected rate, as to conserve the volume
		 *  while shape change is occuring.
		 */
		this->ShapeChangeType = ShapeChangeType;
		if (ShapeChangeType == 1){
			//Thisis change form columnat to cuboidal
			GrowthRate[0] =   0.5 * ShapeChangeRate;  //xx
			GrowthRate[1] =   0.5 * ShapeChangeRate;  //yy
			GrowthRate[2] =  -1.0 *ShapeChangeRate;		  //zz
		}
	}///< The constructor of UniformShapeChangeFunction

    ~UniformShapeChangeFunction(){}

    void writeSummary(std::ofstream &saveFileSimulationSummary,double dt){
		/**
		 *  This function will write the UniformGrowthFunction details into the simulation summary file, provided as the first input.
		 *  Time step (dt) of the simulation is provided as second input, to report the growth rates per hour.
		 *  The output should look like: \n
		 *			Growth Type:  Uniform (1)
		 *			Initial time(sec): UniformGrowthFunction#initTime	FinalTime time(sec): UniformGrowthFunction#endTime	GrowthRate(fraction/hr): DVGrowth(in 1/hr)  APGrowth(in 1/hr) ABGrowth(in 1/hr)
		 */
        saveFileSimulationSummary<<"Shape Change Type:  Uniform (1)"<<std::endl;
		saveFileSimulationSummary<<"	Shape change type: ";
		saveFileSimulationSummary<<ShapeChangeType;
		saveFileSimulationSummary<<"	Initial time(sec): ";
		saveFileSimulationSummary<<initTime;
		saveFileSimulationSummary<<"	FinalTime time(sec): ";
		saveFileSimulationSummary<<endTime;
		saveFileSimulationSummary<<"	ShapeChangeRate(fraction/hr): ";
		saveFileSimulationSummary<<GrowthRate[0]/dt*3600.0;
		saveFileSimulationSummary<<"  ";
		saveFileSimulationSummary<<GrowthRate[1]/dt*3600.0;
		saveFileSimulationSummary<<"  ";
		saveFileSimulationSummary<<GrowthRate[2]/dt*3600.0;
        saveFileSimulationSummary<<std::endl;
	}///< The function is to write the growth function summary to simulation summary file
};


class RingGrowthFunction : public GrowthFunctionBase{
private:

public:
	double centre[2];		///< The double array of 2, giving the centre of the ring in micro-meters, format [x, y].
	double innerRadius;		///< The inner radius of the ring, inner boundary of the growth region in micro-meters. The growth will be zero at and inside the inner radius. This value can be set to zero to have circular growth.
	double outerRadius; 	///< The outer radius of the ring, outer boundary of the growth region in micro-meters. The growth will be at the maximum value set by RingGrowthFunction#GrowthRate at the outer radius.
	double GrowthRate[3]; 	///< The maximum growth rate at the RingGrowthFunction#outerRadius, in (1/sec), format: [ DV axis (x), AP axis (y), and AB axis (z)]
	gsl_matrix* ShearAngleRotationMatrix; ///< The rotation  matrix for the orientation of the growth on x-y plane. This matrix is constructed through  UniformGrowthFunction#angle
	double angle;	///< The rotation angle for the orientation of the growth on x-y plane.
	RingGrowthFunction(int id, int type, float initTime, float endTime, bool applyToColumnarLayer, bool applyToPeripodialMembrane, bool applyToBasalECM, bool applyToLateralECM, double Cx, double Cy, double innerR, double outerR, double DVGrowth, double APGrowth, double ABGrowth, double angle) : GrowthFunctionBase(id, type, initTime, endTime, applyToColumnarLayer, applyToPeripodialMembrane, applyToBasalECM, applyToLateralECM){
		/**
		 *  The first six parameters will be directed to the parent constructor, GrowthFunctionBase#GrowthFunctionBase. \n
		 *  doubles Cx and Cy will set RingGrowthFunction#centre[0] and RingGrowthFunction#centre[1], respectively. \n
		 *  doubles innerR and outerR will set RingGrowthFunction#innerRadius and RingGrowthFunction#outerRadius, respectively. \n
		 *  doubles DVGrowth, APGrowth and ABGrowth will set the  RingGrowthFunction#GrowthRate, in the given order.
		 */
		this->innerRadius = innerR;
		this->outerRadius = outerR;
		this->centre[0] = Cx;
		this->centre[1] = Cy;
		this -> angle = angle;
		ShearAngleRotationMatrix = gsl_matrix_calloc(3,3);
		gsl_matrix_set_identity(ShearAngleRotationMatrix);
		double c = cos(angle);
		double s = sin(angle);
		gsl_matrix_set(ShearAngleRotationMatrix,0,0,  c );
		gsl_matrix_set(ShearAngleRotationMatrix,0,1, -1.0*s);
		gsl_matrix_set(ShearAngleRotationMatrix,0,2,  0.0);
		gsl_matrix_set(ShearAngleRotationMatrix,1,0,  s);
		gsl_matrix_set(ShearAngleRotationMatrix,1,1,  c);
		gsl_matrix_set(ShearAngleRotationMatrix,1,2,  0.0);
		gsl_matrix_set(ShearAngleRotationMatrix,2,0,  0.0);
		gsl_matrix_set(ShearAngleRotationMatrix,2,1,  0.0);
		gsl_matrix_set(ShearAngleRotationMatrix,2,2,  1.0);

		GrowthRate[0] = DVGrowth;	//x
		GrowthRate[1] = APGrowth;	//y
		GrowthRate[2] = ABGrowth;	//z
	}///< The constructor of RingGrowthFunction
    ~RingGrowthFunction(){}

	void 	getCentre(float &centreX, float &centreY){
		centreX = centre[0];
		centreY = centre[1];
	} ///< This function writes RingGrowthFunction#centre into the input float addresses centreX and centreY, in this order.

	float 	getInnerRadius(){
		return innerRadius;
	}///< This function returns RingGrowthFunction#innerRadius.

	float 	getOuterRadius(){
		return outerRadius;
	}///< This function returns RingGrowthFunction#outerRadius.

    std::array<double,3> getGrowthRate(){
		/**
         *  This function will return the RingGrowthFunction#GrowthRate of the current growth function.
        */
        std::array<double,3> maxValues {0};
		maxValues[0] = GrowthRate[0];
		maxValues[1] = GrowthRate[1];
		maxValues[2] = GrowthRate[2];
        return maxValues;
	}///< The function is to set the 3D maximum growth rate of the current ring growth function.

	void setGrowtRate(double ex, double ey, double ez){
		/**
		 *  This function will set the RingGrowthFunction#GrowthRate of the current growth function to the input values
		 *  The parameters are in the order [ DV axis (x), AP axis (y), and AB axis (z)].
		 */
		GrowthRate[0] = ex;
		GrowthRate[1] = ey;
		GrowthRate[2] = ez;
	}///< The function is to set the 3D maximum growth rate of the current ring growth function.

	gsl_matrix* getShearAngleRotationMatrix(){
		return ShearAngleRotationMatrix;
	}///< The function return the matrix pointer to the oriented growth rotation matrix.

	double getShearAngle(){
		return angle;
	}///< The function returns the oriented growth angle in radians (double).

    void writeSummary(std::ofstream &saveFileSimulationSummary,double dt){
		/**
		 *  This function will write the RingGrowthFunction details into the simulation summary file, provided as the first input.
		 *  Time step (dt) of the simulation is provided as second input, to report the growth rates per hour.
		 *  The output should look like: \n
		 *			Growth Type:  Ring (2)
		 *			Initial time(sec): RingGrowthFunction#initTime	FinalTime time(sec): RingGrowthFunction#endTime	Centre(micron): RingGrowthFunction#centre[0]  RingGrowthFunction#centre[1]	Inner radius(micron): RingGrowthFunction#innerRadius	Outer radius(micron): RingGrowthFunction#outerRadius	GrowthRate(fraction/hr):DVGrowth(in 1/hr)  APGrowth(in 1/hr) ABGrowth(in 1/hr)
		 */
        saveFileSimulationSummary<<"Growth Type:  Ring (2)"<<std::endl;
		saveFileSimulationSummary<<"	Initial time(sec): ";
		saveFileSimulationSummary<<initTime;
		saveFileSimulationSummary<<"	FinalTime time(sec): ";
		saveFileSimulationSummary<<endTime;
		saveFileSimulationSummary<<"	Centre(micron): ";
		saveFileSimulationSummary<<centre[0];
		saveFileSimulationSummary<<"  ";
		saveFileSimulationSummary<<centre[1];
		saveFileSimulationSummary<<"	Inner radius(micron): ";
		saveFileSimulationSummary<<innerRadius;
		saveFileSimulationSummary<<"	Outer radius(micron): ";
		saveFileSimulationSummary<<outerRadius;
		saveFileSimulationSummary<<"	GrowthRate(fraction/hr): ";
		saveFileSimulationSummary<<GrowthRate[0]/dt*3600.0;
		saveFileSimulationSummary<<"  ";
		saveFileSimulationSummary<<GrowthRate[1]/dt*3600.0;
		saveFileSimulationSummary<<"  ";
		saveFileSimulationSummary<<GrowthRate[2]/dt*3600.0;
		saveFileSimulationSummary<<"	Growth Angle(degrees): ";
		saveFileSimulationSummary<<angle/M_PI*180.0;
        saveFileSimulationSummary<<std::endl;
	}///< The function is to write the growth function summary to simulation summary file
};


class GridBasedGrowthFunction : public GrowthFunctionBase{
private:

public:
    size_t nGridX;	///< The number of grid points that discretise the tissue in x
    size_t nGridY;	///< The number of grid points that discretise the tissue in y
    //GrowthMatrix : [grid_i] [grid_j][x,y,z]
    std::vector<std::vector<std::array<double,3>>> GrowthMatrix;	///<The matrix of growth rates in (1/sec). It is a matrix of double triplets for growth rate at each grid point. The dimensions of the matrix are equal to (GridBasedGrowthFunction::nGridX, GridBasedGrowthFunction::nGridY), and set in constructor of the GridBasedGrowthFunction. The triplets store the growth rate in [ DV axis (x), AP axis (y), and AB axis (z)].

    //xyShearAngleMatrix : [grid_i] [grid_j][orientation angle]
    std::vector<std::vector<double>> xyShearAngleMatrix; ///<The matrix of xy shear rate (rad/sec). It is a matrix of doubles at each grid point. he dimensions of the matrix are equal to (GridBasedGrowthFunction::nGridX, GridBasedGrowthFunction::nGridY), and set in constructor of the GridBasedGrowthFunction.
    //xyShearRotationsMatrix : [grid_i] [grid_j][orientation rotation mat ]
    gsl_matrix*** xyShearRotationsMatrix;
    //aspectRatioOverThresoldMatrix : [grid_i] [grid_j][AR over threshold bool]
    std::vector<std::vector<bool>> aspectRatioOverThresoldMatrix;

    //[grid_i] [grid_j][4-corners][x,y,z]
    std::vector<std::vector<std::array<std::array<double,3>,4>>> compatibleGrowths;
    //[grid_i] [grid_j][4-corners][(bool)]
    std::vector<std::vector<std::array<bool,4>>> compatibleAngleEliminated;
    //[grid_i] [grid_j][4-corners][(double) angles ]
    std::vector<std::vector<std::array<double,4>>> compatibleAngles;

    GridBasedGrowthFunction(int id, int type, float initTime, float endTime, bool applyToColumnarLayer, bool applyToPeripodialMembrane, bool applyToBasalECM, bool applyToLateralECM, int nX, int nY, std::vector<std::vector<std::array<double,3>>>& GrowthMat, std::vector<std::vector<double>>& AngleMat) : GrowthFunctionBase(id, type, initTime, endTime, applyToColumnarLayer, applyToPeripodialMembrane, applyToBasalECM,  applyToLateralECM){
		/**
		 *  The first six parameters will be directed to the parent constructor, GrowthFunctionBase#GrowthFunctionBase. \n
		 *  integers nX and nY will set GridBasedGrowthFunction#nGridX and GridBasedGrowthFunction#nGridY, respectively.  GridBasedGrowthFunction::GrowthMatrix will be initiated to point at a 2 dimensional matrix of double triplets the size(nX, nY). \n
		 *  double*** GrowthMat is the pointer to the 2-dimensional matrix of double triplets, holding the 3D growth rates at each grid point. Values stored in GrowthMat will set the values in GridBasedGrowthFunction::GrowthMatrix.
		 *  The matrix storing the growth rates have been read from an input file through ModelInputObject#readGrowthOptions and related functions therein \n
		 */
		double aspectRatioThreshold = 1.1;
		this ->nGridX = nX;
		this ->nGridY = nY;
        //GrowthMatrix : [grid_i] [grid_j][x,y,z]
        for (size_t i= 0; i<nGridX; ++i){
            std::vector<std::array<double,3>> tmpGridMatrixY(nGridY,std::array<double,3>{0.0});
            GrowthMatrix.push_back(tmpGridMatrixY);
            std::vector<double> tmpShearY(nGridY,0.0);
            xyShearAngleMatrix.push_back(tmpShearY);
            std::vector<bool> tmpAROverThresY (nGridY,true); //set AR to be significant, and assign it to be too small if necesarry
            aspectRatioOverThresoldMatrix.push_back(tmpAROverThresY);

            //matrices with 4 corners,for compatable averaging:
            std::vector<std::array<std::array<double,3>,4>> tmpGrid4CornerMatrixY(nGridY,std::array<std::array<double,3>,4>{std::array<double,3>{0.0}});
            compatibleGrowths.push_back(tmpGrid4CornerMatrixY);
            std::vector<std::array<bool,4>> tmpGrid4CornerAngleEliminatedY(nGridY,std::array<bool,4>{false});
            compatibleAngleEliminated.push_back(tmpGrid4CornerAngleEliminatedY);
            std::vector<std::array<double,4>> tmpGrid4CornerAngleY(nGridY,std::array<double,4>{0.0});
            compatibleAngles.push_back(tmpGrid4CornerAngleY);
        }
        for (size_t i=0; i<nGridX; ++i){
            for (size_t j=0; j<nGridY; ++j){
                xyShearAngleMatrix[i][j] = AngleMat[i][j]*M_PI/180.0; //converting to radians
                for (size_t k=0; k<3; ++k){
                    GrowthMatrix[i][j][k] = GrowthMat[i][j][k];
                }
                double AR = GrowthMatrix[i][j][0] / GrowthMatrix[i][j][1];
				if (AR < 0 ){
					//some of the growth rates have higher aspect ratio than their nuclei count, for instance, an aspect ratio of 5 with nuclei count 4
					//then the tissue actually needs to actively shring in one direction while growing in the other.
					//this will lead to a negative aspect ratio, and this extremely high aspect ratio would be eliminated by the threshold, as it would be less than the positive threshold.
					AR *=-1.0;
				}
				if (AR < 1.0){
					AR = 1/AR;
				}
				if(AR < aspectRatioThreshold){
					aspectRatioOverThresoldMatrix[i][j] = false;
                }
			}
		}

		xyShearRotationsMatrix = new gsl_matrix**[(const int) nGridX];
        for (size_t i=0; i<nGridX; ++i){
			xyShearRotationsMatrix[i] = new gsl_matrix*[(const int) nGridY];
            for (size_t j=0; j<nGridY; ++j){
				xyShearRotationsMatrix[i][j] = new gsl_matrix;
				xyShearRotationsMatrix[i][j] = gsl_matrix_calloc(3,3);
				double c = cos(xyShearAngleMatrix[i][j]);
				double s = sin(xyShearAngleMatrix[i][j]);
                //cout<<"xyShearRotationsMatrix ["<<i<<"]["<<j<<"], angle: "<<xyShearAngleMatrix[i][j]<<" cos: "<<c<<" sin "<<s<<std::endl;
				gsl_matrix_set(xyShearRotationsMatrix[i][j],0,0,  c );
				gsl_matrix_set(xyShearRotationsMatrix[i][j],0,1, -1.0*s);
				gsl_matrix_set(xyShearRotationsMatrix[i][j],0,2,  0.0);
				gsl_matrix_set(xyShearRotationsMatrix[i][j],1,0,  s);
				gsl_matrix_set(xyShearRotationsMatrix[i][j],1,1,  c);
				gsl_matrix_set(xyShearRotationsMatrix[i][j],1,2,  0.0);
				gsl_matrix_set(xyShearRotationsMatrix[i][j],2,0,  0.0);
				gsl_matrix_set(xyShearRotationsMatrix[i][j],2,1,  0.0);
				gsl_matrix_set(xyShearRotationsMatrix[i][j],2,2,  1.0);
			}
		}
		//now I want to fill up growth matrix for corners:
		preCalculateAnglesForCompatibleAveraging();
	} ///< The constructor of GridBasedGrowthFunction

	~GridBasedGrowthFunction(){
        for (size_t i=0; i<nGridX; ++i){
            for (size_t j=0; j<nGridY; ++j){
                gsl_matrix_free(xyShearRotationsMatrix[i][j]);
			}
		}
	}

    size_t getGridX(){
		return nGridX;
	}///< This function returns GridBasedGrowthFunction#nGridX.

    size_t getGridY(){
		return nGridY;
	}///< This function returns GridBasedGrowthFunction#nGridY.

	double getGrowthMatrixElement(int i, int j, int k){
		return GrowthMatrix[i][j][k];
	}///< This function returns the growth rate at grid point [i]\[j\] (in dimensions GridBasedGrowthFunction#nGridX, GridBasedGrowthFunction#nGridY), for the growth dimension [k] (as in  [ DV axis (x), AP axis (y), and AB axis (z)] ).

	double getXyShearAngleMatrixElement(int i, int j){
		return xyShearAngleMatrix[i][j];
	}///< This function returns the xyShear angle at grid point [i]\[j\] (in dimensions GridBasedGrowthFunction#nGridX, GridBasedGrowthFunction#nGridY).

	gsl_matrix* getXyShearRotationsMatrixElement(int i, int j){
		return xyShearRotationsMatrix[i][j];
	}///< This function returns the xyShear rotation matrix at grid point [i]\[j\] (in dimensions GridBasedGrowthFunction#nGridX, GridBasedGrowthFunction#nGridY).

	bool isAspectRatioOverOne(int i, int j){
		return aspectRatioOverThresoldMatrix[i][j];
	}

	void setGrowthMatrixElement(double ex, double ey, double ez, int i, int j){
		GrowthMatrix[i][j][0] = ex;
		GrowthMatrix[i][j][1] = ey;
		GrowthMatrix[i][j][2] = ez;
	}///< This function sets the growth rate at grid point [i]\[j\] (in dimensions GridBasedGrowthFunction#nGridX, GridBasedGrowthFunction#nGridY), to the growth rate [ex, ey, ez] in the format [ DV axis (x), AP axis (y), and AB axis (z)].

	void preCalculateAnglesForCompatibleAveraging(){
        for (size_t IndexX = 0; IndexX<nGridX-1; ++IndexX){
            for (size_t IndexY = 0; IndexY<nGridY-1; ++IndexY){
				//read the growth rates at four corners for this grid point:
				double growth0[3] = {0.0,0.0,0.0};
				double growth1[3] = {0.0,0.0,0.0};
				double growth2[3] = {0.0,0.0,0.0};
				double growth3[3] = {0.0,0.0,0.0};
				for (int axis =0; axis<3; axis++){
					growth0[axis] = getGrowthMatrixElement(IndexX,IndexY,axis);
					growth1[axis] = getGrowthMatrixElement(IndexX+1,IndexY,axis);
					growth2[axis] = getGrowthMatrixElement(IndexX,IndexY+1,axis);
					growth3[axis] = getGrowthMatrixElement(IndexX+1,IndexY+1,axis);
				}
				//now reading the angles:
				double angles[4] = {0.0,0.0,0.0,0.0};
				//Not summing aspect ratios close to one:
				//and recording the eliminated angles for ease of use later on:
				bool AngleEliminated[4] = {false, false, false, false};
				if (isAspectRatioOverOne( IndexX,   IndexY))   {angles[0] = getXyShearAngleMatrixElement(IndexX  , IndexY);  }
				else {AngleEliminated[0] = true;}
				if (isAspectRatioOverOne( IndexX+1, IndexY))   {angles[1] = getXyShearAngleMatrixElement(IndexX+1, IndexY);  }
				else {AngleEliminated[1] = true;}
				if (isAspectRatioOverOne( IndexX,   IndexY+1)) {angles[2] = getXyShearAngleMatrixElement(IndexX,   IndexY+1);}
				else {AngleEliminated[2] = true;}
				if (isAspectRatioOverOne( IndexX+1, IndexY+1)) {angles[3] = getXyShearAngleMatrixElement(IndexX+1, IndexY+1);}
				else {AngleEliminated[3] = true;}
				//bring all angles in the interval 0 - 180:
				for (int i=0; i<4; ++i){
					while (angles[i]<0){angles[i] += M_PI;}
					while (angles[i]>M_PI){angles[i] -= M_PI;}
				}
				//now I am bringing the angles to the vicinity of each other:
				//selecting the base angle, the first one that was not eliminated will do:
				bool angleSelected = false;
				int selectedIndice = 3;
				double tet = 0;
				for (int i=0; i<4; ++i){
					if(!AngleEliminated[i] && !angleSelected){
						tet = angles[i];
						angleSelected = true;
					}
					if (angleSelected){
						selectedIndice = i;
						break;
					}
				}
				//selected the first angle, now bring all the rest near this one:
				//I keep a record if I need to flip the growths:
				bool flipGrowth[4] = {false, false, false, false};
				for (int i=selectedIndice; i<4; ++i){
					double dtet = angles[i] - tet; // angle from selected angle (tet) to current angle (angles[i]).
					if (dtet>0){
						//selected tetha is larger than this angle:
						if (dtet >= M_PI/4.0 && dtet < 3.0*M_PI/4.0 ){
							//the difference is between 45-135 degrees, I flip the growth, and rotate angle -90 degrees.
							//The difference comes to the interval (-45) - 45
							angles[i] = angles[i] - M_PI/2.0;
							flipGrowth[i] = true;
						}
						else if (dtet >= 3.0*M_PI/4.0 && dtet< 5.0*M_PI/4.0 ){
							//the difference is between 135-225 degrees, I only need to rotate -180 degrees,
							//no flip necessary, the angle will fall into the interval (-45) - 45
							angles[i] = angles[i] - M_PI;
						}
						else if (dtet >= 5.0*M_PI/4.0 && dtet< 7.0*M_PI/4.0){
							//the difference is between 225 and 315 degrees, I need to rotate and flip as in the first case.
							//I need to rotate by -270 degrees to bring to the interval (-45) - 45
							angles[i] = angles[i] - 6.0*M_PI/4.0;
							flipGrowth[i] = true;
						}
						else if (dtet >= 7.0*M_PI/4.0){
							//the difference is higher than 315 degrees. I need to rotate the angle by -360 degrees, and
							//the difference will come into (-45) - 0 (as it cannot be higher than 360. No need to flip
							angles[i] = angles[i] - 2.0*M_PI;
						}
					}
					else if (dtet<0){
						//selected tetha is smaller than this angle:
						if (dtet <= -M_PI/4.0 && dtet > -3.0*M_PI/4.0 ){
							//the difference is between (-45)- (-135) degrees, I flip the growth, and rotate angle +90 degrees.
							//The difference comes to the interval (-45) - 45
							angles[i] = angles[i] + M_PI/2.0;
							flipGrowth[i] = true;
						}
						else if (dtet <= -3.0*M_PI/4.0 && dtet > -5.0*M_PI/4.0 ){
							//the difference is between (-135) - (-225) degrees, I only need to rotate +180 degrees,
							//no flip necessary, the angle will fall into the interval (-45) - 45
							angles[i] = angles[i] + M_PI;
						}
						else if (dtet <= -5.0*M_PI/4.0 && dtet > -7.0*M_PI/4.0){
							//the difference is between (-225) and (-315) degrees, I need to rotate and flip as in the first case.
							//I need to rotate by +270 degrees to bring to the interval (-45) - 45
							angles[i] = angles[i] + 6.0*M_PI/4.0;
							flipGrowth[i] = true;
						}
						else if (dtet <= -7.0*M_PI/4.0){
							//the difference is less than (-315) degrees. I need to rotate the angle by +360 degrees, and
							//the difference will come into 0 - 45 (as it cannot be lower than (-360). No need to flip
							angles[i] = angles[i] + 2.0*M_PI;
						}
					}
				}
				//all angles are brought to the vicinity of the selected base angle;
				//flipping the growth where necessary:
				if (flipGrowth[0]){
					double x = growth0[0];
					growth0[0] = growth0[1];
					growth0[1] = x;
				}
				if (flipGrowth[1]){
					double x = growth1[0];
					growth1[0] = growth1[1];
					growth1[1] = x;
				}
				if (flipGrowth[2]){
					double x = growth2[0];
					growth2[0] = growth2[1];
					growth2[1] = x;
				}
				if (flipGrowth[3]){
					double x = growth3[0];
					growth3[0] = growth3[1];
					growth3[1] = x;
				}
                for (size_t axis=0; axis<3; ++axis){
					compatibleGrowths[IndexX][IndexY][0][axis]= growth0[axis];
					compatibleGrowths[IndexX][IndexY][1][axis]= growth1[axis];
					compatibleGrowths[IndexX][IndexY][2][axis]= growth2[axis];
					compatibleGrowths[IndexX][IndexY][3][axis]= growth3[axis];
				}
                for (size_t i=0; i<4; ++i){
					compatibleAngles[IndexX][IndexY][i] =  angles[i];
					compatibleAngleEliminated[IndexX][IndexY][i] = AngleEliminated[i];
				}
			} // end of loop over IndexY
		} //end of loop over IndexX
	}

    void getGrowthProfileAt4Corners(int IndexX, int IndexY, std::array<double,3>& growth0, std::array<double,3>& growth1, std::array<double,3>& growth2, std::array<double,3>& growth3, std::array<double,4>& angles, std::array<bool,4>& anglesEliminated){
        for (size_t i=0; i<3; ++i){
			growth0[i] = compatibleGrowths[IndexX][IndexY][0][i];
			growth1[i] = compatibleGrowths[IndexX][IndexY][1][i];
			growth2[i] = compatibleGrowths[IndexX][IndexY][2][i];
			growth3[i] = compatibleGrowths[IndexX][IndexY][3][i];
		}
        for (size_t i=0; i<4; ++i){
			angles[i]  = compatibleAngles[IndexX][IndexY][i];
			anglesEliminated[i] = compatibleAngleEliminated[IndexX][IndexY][i];
		}
	}

    void writeSummary(std::ofstream &saveFileSimulationSummary){
		/**
		 *  This function will write the GridBasedGrowthFunction details into the simulation summary file, provided as the first input.
		 *  The output should look like: \n
		 *			Growth Type:  growth From File (3)
		 *			Initial time(sec): GridBasedGrowthFunction#initTime	FinalTime time(sec): GridBasedGrowthFunction#endTime	Growth matrix mesh size: GridBasedGrowthFunction#nGridX GridBasedGrowthFunction#nGridY
		 */
        saveFileSimulationSummary<<"Growth Type:  growth From File (3)"<<std::endl;
		saveFileSimulationSummary<<"	Initial time(sec): ";
		saveFileSimulationSummary<<initTime;
		saveFileSimulationSummary<<"	FinalTime time(sec): ";
		saveFileSimulationSummary<<endTime;
		saveFileSimulationSummary<<"	Growth matrix mesh size: ";
        saveFileSimulationSummary<<nGridX <<" "<<nGridY<<std::endl;
	}///< The function is to write the growth function summary to simulation summary file
};

class markerEllipseBasedShapeChangeFunction : public GrowthFunctionBase{
private:

public:
	int nEllipseBands;
	double GrowthRate[3];
	int ShapeChangeType;
	bool conserveVolume;
	double ShapeChangeECMLimit;
    markerEllipseBasedShapeChangeFunction(int id, int type, float initTime, float endTime, bool applyTissueApical, bool applyTissueBasal, bool applyTissueMidline, bool applyToBasalECM, bool applyToLateralECM, int ShapeChangeType, double ShapeChangeFractionPerHr, std::vector <int> markerEllipseIds, bool conserveVolume) : GrowthFunctionBase(id, type, initTime, endTime, true, false, applyToBasalECM,  applyToLateralECM){
		/**
		 *  First six parameters will be directed to the parent constructor, UniformGrowthFunction#UniformGrowthFunction.
		 *  The growth rates in Dv, AB and AP will be fed as 0 to the parent constructor. \n
		 *  The parameter ShapeChangeType will define the type of the shape change as 1: ColumnarToSquamous change, 2: Apical shrinkage, and 2: basal shrinkage.
		 *  The parameter ShapeChangeRate will define the rate of defined shape change. The constructor will calculate the corrected rate, as to conserve the volume
		 *  while shape change is occurring.
		 */
		this->ShapeChangeType = ShapeChangeType;
		this->applyTissueBasal = applyTissueBasal;
		this->applyTissueApical = applyTissueApical;
		this->applyTissueMidLine = applyTissueMidline;
		this->conserveVolume = conserveVolume;
        //this->ShapeChangeECMLimit = ShapeChangeECMLimit;
		//This is change form columnar to cuboidal
		//I have how much I want the shape to change x+y in a second.
		//Then exp(rx*1)*exp(ry*1) = ShapeChangeFractionPerSec
		//Then rx+ry = ln(ShapeChangeFractionPerSec)
		//At the same time, I would like to compansate the change in x+y with
		//the change in z, therefore:
		//exp(rz*1) = 1/ShapeChangeFractionPerSec
		//rz = ln(1/ShapeChangeFractionPerSec);
		GrowthRate[0] =  0.5 * log(1.0+ShapeChangeFractionPerHr)/3600;  //xx
		GrowthRate[1] =  0.5 * log(1.0+ShapeChangeFractionPerHr)/3600;  //yy
		if (conserveVolume){
			GrowthRate[2] =  log(1.0/(1.0+ShapeChangeFractionPerHr))/3600;		 //zz
		}
		else{
			GrowthRate[2] = 0;
		}
		nEllipseBands = markerEllipseIds.size();
		for (int i=0;i<nEllipseBands; ++i){
			appliedEllipseBandIds.push_back(markerEllipseIds[i]);
		}
        std::cout<<" Shape Fraction change per h: "<<ShapeChangeFractionPerHr<<" rates of shape change: "<<GrowthRate[0]<<" "<<GrowthRate[1]<<" "<<GrowthRate[2]<<std::endl;
	}///< The constructor of UniformShapeChangeFunction

	~markerEllipseBasedShapeChangeFunction(){};

    std::array<double,3> getShapeChangeRateRate(){
			/**
			 *  This function will write the markerEllipseBasedShapeChangeFunction#GrowthRate of the current growth function to the input double array
			 *  pointer. The double array pointer should be set to point at a double array of size 3 (or higher) before calling the function.
			 */
            std::array<double,3> rates {0.0};
			rates[0] = GrowthRate[0];
			rates[1] = GrowthRate[1];
			rates[2] = GrowthRate[2];
            return rates;
		} ///< The function is to get the 3D growth rate of the current shape change function.


    void writeSummary(std::ofstream &saveFileSimulationSummary,double dt){
		/**
		 *  This function will write the UniformGrowthFunction details into the simulation summary file, provided as the first input.
		 *  Time step (dt) of the simulation is provided as second input, to report the growth rates per hour.
		 *  The output should look like: \n
		 *			Growth Type:  Uniform (1)
		 *			Initial time(sec): UniformGrowthFunction#initTime	FinalTime time(sec): UniformGrowthFunction#endTime	GrowthRate(fraction/hr): DVGrowth(in 1/hr)  APGrowth(in 1/hr) ABGrowth(in 1/hr)
		 */
        saveFileSimulationSummary<<"Shape Change Type:  MArkerEllipseBased (2)"<<std::endl;
		saveFileSimulationSummary<<"	Shape change type: ";
		saveFileSimulationSummary<<ShapeChangeType;
		saveFileSimulationSummary<<"	Initial time(sec): ";
		saveFileSimulationSummary<<initTime;
		saveFileSimulationSummary<<"	FinalTime time(sec): ";
		saveFileSimulationSummary<<endTime;
		saveFileSimulationSummary<<"	ShapeChangeRate(fraction/hr): ";
		saveFileSimulationSummary<<GrowthRate[0]/dt*3600.0;
		saveFileSimulationSummary<<"  ";
		saveFileSimulationSummary<<GrowthRate[1]/dt*3600.0;
		saveFileSimulationSummary<<"  ";
		saveFileSimulationSummary<<GrowthRate[2]/dt*3600.0;
		saveFileSimulationSummary<<"	Applied ellipses: ";
		for (int i=0;i<nEllipseBands; ++i){
			saveFileSimulationSummary<<appliedEllipseBandIds[i]<<" ";
		}
        saveFileSimulationSummary<<std::endl;
	}///< The function is to write the growth function summary to simulation summary file
};


#endif /* GROWTHFUNCTIONTYPES_H */
