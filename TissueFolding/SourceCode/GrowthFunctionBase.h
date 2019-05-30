#ifndef GrowthFunctionBase_H
#define GrowthFunctionBase_H

#include <stdio.h>
#include <iostream>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <array>
#include <vector>

class GrowthFunctionBase{
private:

public:
	GrowthFunctionBase(int id, int type, float initTime, float endTime, bool applyToColumnarLayer, bool applyToPeripodialMembrane, bool applyToBasalECM, bool applyToLateralECM){
		/**
		 *  integer id will set GrowthFunctionBase#Id. \n
		 *  integer type will set GrowthFunctionBase#Type. \n
		 *  floats initTime and endTime will set GrowthFunctionBase#initTime and GrowthFunctionBase#endTime respectively. \n
		 *	booleans applyToColumnarLayer and applyToPeripodialMembrane will set GrowthFunctionBase#applyToColumnarLayer and GrowthFunctionBase#applyToPeripodialMembrane, respectively. \n
		 *   */
		this->Id = id;
		this->Type = type;
		this->initTime = initTime;
		this->endTime = endTime;
		this->applyToColumnarLayer = applyToColumnarLayer;
		this->applyToPeripodialMembrane = applyToPeripodialMembrane;
		this->applyToBasalECM = applyToBasalECM;
		this->applyToLateralECM = applyToLateralECM;
		zMin = 0;
		zMax = 1.0;
	} ///< The constructor of GrowthFunctionBase. Different growth functions will be derived from this class
    ~GrowthFunctionBase(){}

	int Type; 						///< The type of the growth function, 1: uniform growth, 2: Ring shaped growth, 3: Grid based growth, where growth rates read from a separate input file
	int Id;							///< The unique identification number of the growth function
	float initTime;					///< The initiation time of the growth, in seconds
	float endTime;					///< The end time of the growth, in seconds.
	bool applyToColumnarLayer;		///< Boolean stating if the growth should be applied to columnar layer
	bool applyToPeripodialMembrane; ///< Boolean stating if the growth should be applied to peripodial membrane
	bool applyToBasalECM;			///< Boolean stating if the growth should be applied to basal explicit ECM layer
	bool applyToLateralECM ;		///< Boolean stating if the growth should be applied to lateral explicit ECM layer
    bool applyTissueApical;         ///< Boolean stating if the growth should be applied to apical side of the tissue
    bool applyTissueBasal;          ///< Boolean stating if the growth should be applied to basal side of the tissue
    bool applyTissueMidLine;        ///< Boolean stating if the growth should be applied to mid-line side of the tissue
    double zMin;                    ///< Double giving the minimum z fraction of the tissue height that the growth is applied to
    double zMax;                    ///< Double giving the maximum z fraction of the tissue height that the growth is applied to
    std::vector <int> appliedEllipseBandIds;    ///< The vector of integers holding the iIDs of the marker ellipse bands that the growth isapplied to

    void 	ParentErrorMessage(std::string functionName){
        std::cerr<<"You are calling the function: "<<functionName<<" from a parent here, check declaration is via pointers"<<std::endl;
	}
    double 	ParentErrorMessage(std::string functionName, double returnValue){
        std::cerr<<"You are calling the function: "<<functionName<<" from a parent here, check declaration is via pointers"<<std::endl;
		return returnValue;
	}
    int 	ParentErrorMessage(std::string functionName, int returnValue){
        std::cerr<<"You are calling the function: "<<functionName<<" from a parent here, check declaration is via pointers"<<std::endl;
		return returnValue;
	}
    virtual void		writeSummary(std::ofstream&/*saveFileSimulationSummary*/){ParentErrorMessage("writeSummary_no_dt");} /// Virtual function of the parent to write the summary of growthFunction.
    virtual void		writeSummary(std::ofstream&/*saveFileSimulationSummary*/, double /*dt*/){ParentErrorMessage("writeSummary_with_dt");} /// Virtual function of the parent to write the summary of growthFunction.
    virtual void 		getCentre(float &/*centreX*/, float &/*centreY*/){ParentErrorMessage("getCentre");} /// Virtual function of the parent to get centre of the growth function definition.
    virtual float 		getInnerRadius(){return ParentErrorMessage("getInnerRadius",0.0);} /// Virtual function of the parent to get inner radius of the growth function definition.
    virtual float 		getOuterRadius(){return ParentErrorMessage("getOuterRadius",0.0);} /// Virtual function of the parent to get outer radius of the growth function definition.
    virtual std::array<double,3> 		getGrowthRate(){ParentErrorMessage("getGrowthRate");return std::array<double,3>{0.0};} /// Virtual function of the parent to get growth rate of the growth function definition.
    virtual std::array<double,3> 		getShapeChangeRateRate(){ParentErrorMessage("getShapeChangeRateRate");return std::array<double,3>{0.0};} /// Virtual function of the parent to get shape change rate of the growth function definition.
    virtual gsl_matrix* getShearAngleRotationMatrix(){ParentErrorMessage("getShearAngleRotationMatrix"); gsl_matrix* dummy = gsl_matrix_calloc(0,0); return dummy;}  /// Virtual function of the parent to get rotation angle matrix of the growth function definition.
    virtual double 		getShearAngle(){ParentErrorMessage("getShearAngle");return 0.0;} /// Virtual function of the parent to get rotation angle of the growth function definition.
    virtual size_t      getGridX(){return ParentErrorMessage("getGridX",0);} /// Virtual function of the parent to get gridX of the growth function definition.
    virtual size_t      getGridY(){return ParentErrorMessage("getGridY",0);} /// Virtual function of the parent to get gridY of the growth function definition.
    virtual	double 		getGrowthMatrixElement(int /*i*/, int /*j*/, int /*k*/){return ParentErrorMessage("getGrowthMatrixElement",0.0);} /// Virtual function of the parent to get grid matrix element of the growth function definition.
    virtual	double 		getXyShearAngleMatrixElement(int /*i*/, int /*j*/){return ParentErrorMessage("getXyShearhMatrixElement",0.0);} /// Virtual function of the parent to get grid matrix element of the growth function definition.
    virtual bool 		isAspectRatioOverOne(int /*i*/, int /*j*/){return ParentErrorMessage("isAspectRatioOverOne",0);}/// Virtual function of the parent to check aspect ratioof the growth function definition.
    virtual gsl_matrix* getXyShearRotationsMatrixElement(int /*i*/, int /*j*/){ParentErrorMessage("getShearAngleRotationMatrixElement");gsl_matrix* dummy = gsl_matrix_calloc(0,0); return dummy;} //this is used by grid based growth
    virtual void 		getGrowthProfileAt4Corners(int /*IndexX*/, int /*IndexY*/, std::array<double,3>& /*growth0*/, std::array<double,3>&/*growth1*/, std::array<double,3>&/*growth2*/, std::array<double,3>&/*growth3*/, std::array<double,4>&/*angles*/, std::array<bool,4>&/*anglesEliminated*/){ParentErrorMessage("getGrowthProfileAt4Corners");}

    virtual void		setGrowtRate(double /*ex*/, double /*ey*/, double /*ez*/){ParentErrorMessage("setGrowtRate");}/// Virtual function of the parent to set growth rateof the growth function definition.
    virtual void		setGrowthMatrixElement(double /*ex*/, double /*ey*/, double /*ez*/, int /*i*/, int /*j*/){ParentErrorMessage("setGrowtRate");} /// Virtual function of the parent to set growth rateof the growth function definition.
};
#endif

