/*
 * CellMigration.cpp
 *
 *  Created on: 14 Sep 2016
 *      Author: melda
 */


#include "CellMigration.h"
//#include "Node.h"
//#include <gsl/gsl_linalg.h>

using namespace std;

CellMigration::CellMigration(int nElements, double rate){
	this->r1 = 70;
	this->r2 = 45;
	this->numberOfElements = nElements;
	this->migrationRate = rate;
	//this->migrationOrigin[0] = centerX;
	//this->migrationOrigin[1] = centerY;
	elementAnglesList = new double[(const int) numberOfElements];
	listOfLeavingVolumePerStep = new double[(const int) numberOfElements];
	listOfGainedVolumePerStep = new double[(const int) numberOfElements];
	listOfNumberOfNeigsEachElementSendsMaterialTo = new int[(const int) numberOfElements];
	listOfNumberOfNeigsEachElementGetsMaterialFrom = new int[(const int) numberOfElements];
	listOfRateFractions = new double[(const int) numberOfElements];
	for (int i=0; i<nElements; ++i){
		elementConnectivityMap.push_back(vector<int> ());
	}
	centreInZAtOriginOfMigration = 0;
	//this->migrationRadius = r;
}


CellMigration::~CellMigration(){
	delete[] elementAnglesList;
	delete[] listOfLeavingVolumePerStep;
	delete[] listOfGainedVolumePerStep;
	delete[] listOfNumberOfNeigsEachElementSendsMaterialTo;
	delete[] listOfNumberOfNeigsEachElementGetsMaterialFrom;
	delete[] listOfRateFractions;
}

void CellMigration::assignElementRadialVectors(vector <ShapeBase*>& Elements){
	cout<<"inside element angle  asignment"<<endl;
	int counter = 0;
	for( vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		elementAnglesList[counter] = 5;
		if ((*itElement)->tissueType == 0 || (*itElement)->tissueType == 1){
			// For columnar and peripodial elements, the angle will scale with distance from centre;
			double *c = (*itElement)->getCentre();
			double f = calculateXYPositionFractionScaledToEllipse(c[0],c[1]);
			//If f is 1, then the element is at the border of linkers, if 0, it is at tissue centre.
			//For peripodial elements, I want the angle to be 0 at centre, 1 at border, use f as is.
			if ((*itElement)->tissueType == 1){
				elementAnglesList[counter]  = f;
			}
			//For columnar elements, I want the angle to be 2 at the border, and 3 at centre, modify f accordingly.
			if ((*itElement)->tissueType == 0){
				elementAnglesList[counter]  = 3-f;
			}
		}
		if ((*itElement)->tissueType == 2){
			//For Linker elements, the value will go from 1 to 2. It will be scaled with peripodialness.
			//for peripodialness = 1, angle = 1, for peripodialness = 0. angle = 2
			elementAnglesList[counter] = 2 - (*itElement)->getPeripodialness();
		}
		counter++;
	}
}

void CellMigration::assignOriginOfMigration(vector <ShapeBase*>& Elements, double originAngle){
	for( vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if ((*itElement)->tissueType == 2){
			(*itElement)->setCellMigration(true);
		}
	}
}

double CellMigration::calculateAngleOfVector(double dx, double dy, double dz){
	double r = pow((dx*dx + dy*dy + dz*dz),0.5);
	if (r > 1E-6){
		dz = dz/r;
	}
	else {
		dz = 1.0;
	}
	double tet = acos(dz);
	return tet;
}

double CellMigration::getRateFractionForElement(int i){
	if (i < numberOfElements){
		return listOfRateFractions[i];
	}
	else{
		return 0;
	}
}

double CellMigration::getMigrationAngleForElement(int i){
	if (i < numberOfElements){
		return elementAnglesList[i];
	}
	else{
		return 10;
	}
}


double CellMigration::calculateXYPositionFractionScaledToEllipse(double posX, double posY){
	//For columnar elements, I calculate the distance from the edge. First 20% gets to migrate with dampening.
	//xy-planar radius from centre to point:
	double dx = posX;
	double dy = posY;
	double r = pow(dx*dx + dy*dy,0.5);
	//xy-planar angle of the point:
	double t = acos(dx/r);
	//the actual radius of the tissue at this angle, using tissue dimensions as r1 = 70, r2= 45;
	double x = r1*cos(t);
	double y = r2*sin(t);
	double r0 = pow(x*x + y*y,0.5);
	//the relative position of point, 1 will be at the border of linker zone, 0 will be tissue centre
	return r/r0;

}
void CellMigration::generateListOfRateFractions(vector <ShapeBase*>& Elements){
	double threshold = 0.7; //the threshold at which migration will stop, 0.8 means 20% of the border zone will migrate
	for (int i=0; i<numberOfElements; ++i){
		listOfRateFractions[i] = 0.0;
		if (Elements[i]->tissueType == 1 || Elements[i]->tissueType == 2 ){
			//All linker and peripodial elements migrate at max rate;
			listOfRateFractions[i] = 1.0;
		}
		if (Elements[i]->tissueType == 0 ){
			//For columnar elements, I calculate the distance from the edge. First 20% gets to migrate with dampening.
			//xy-planar radius from centre to point:
			double *c = Elements[i]->getCentre();
			double f = calculateXYPositionFractionScaledToEllipse(c[0],c[1]);
			//scale the fraction to be 1 at f= 1, and 0 at f = 0.8;
			// if (f < 0.8) then listOfRateFractions[i] = 0.0, as done above
			if (f > 1){
				//this can only be a numerical error:
				listOfRateFractions[i] = 1;
			}
			else if (f > threshold){
				listOfRateFractions[i] = (f - threshold ) / ( 1.0 - threshold);
			}
		}
	}
}

void CellMigration::assignElementConnectivity(vector <Node*>& Nodes, vector <ShapeBase*>& Elements){
	for (int i=0; i< numberOfElements ;++i){
		Elements[i]->fillLateralNeighbours(Nodes,elementConnectivityMap[i]);
		listOfNumberOfNeigsEachElementSendsMaterialTo[i] = 0;
		listOfNumberOfNeigsEachElementGetsMaterialFrom[i] = 0;
		int nNeigs = elementConnectivityMap[i].size();
		for (int j=0; j< nNeigs; ++j){
			int currNeig = elementConnectivityMap[i][j];
			if (elementAnglesList[i] < elementAnglesList[currNeig]){
				listOfNumberOfNeigsEachElementSendsMaterialTo[i]++;
			}
			else{
				listOfNumberOfNeigsEachElementGetsMaterialFrom[i]++;
			}
		}
	}
	//This element have been assigned to a migrating element. I will generate its connectivity map now:
	//For non lateral elements:
	//bibap
	//currElement->getLateralNeighbours();
}

void CellMigration::updateMigrationLists(vector <ShapeBase*>& Elements, double dt){
	int counter = 0;
	for( vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if (listOfRateFractions[counter]>0 && (*itElement)->getCellMigration()){
			double totalmassToSendOut = listOfRateFractions[counter] * migrationRate/3600*dt * (*itElement)->GrownVolume;
			double massToSendOutPerElement = totalmassToSendOut /listOfNumberOfNeigsEachElementSendsMaterialTo[counter];
			listOfLeavingVolumePerStep[counter] += totalmassToSendOut;
			//the element is migrating, I should check the which of neigs it will deliver material to
			int nNeigs = elementConnectivityMap[counter].size();
			if ((*itElement)->Id == 1076 || (*itElement)->Id == 1929 ){
				cout<<"Element: "<<(*itElement)->Id<<"angle: "<<elementAnglesList[counter]<<" nNeigs: "<<nNeigs<<" grown volume: "<<(*itElement)->GrownVolume<<endl;
			}
			//Element will only migrate out material to angles higher than itself:
			for (int i=0; i< nNeigs; ++i){
				int currNeig = elementConnectivityMap[counter][i];
				if ((*itElement)->Id == 1076 || (*itElement)->Id == 1929 ){
					cout<<" Element: "<<(*itElement)->Id<<" currNeig: "<<currNeig<<" neigAngle: "<<elementAnglesList[currNeig]<<" neig grown volume : "<<Elements[currNeig]->GrownVolume<<endl;
				}
				if (elementAnglesList[counter] < elementAnglesList[currNeig]){
					//Current element is sending material to the element on the neighbour list:
					listOfGainedVolumePerStep[currNeig]+=massToSendOutPerElement;
				}
			}
		}
		counter++;
	}

	int idList[6] = {1076, 1141, 1073, 1929, 1933, 2009};
	for(int j=0;j<6;++j){
		int i = idList[j];
		cout<<"Element "<<Elements[i]->Id<<" gainedVolume: "<<listOfGainedVolumePerStep[i]<<" lostVolume: "<<listOfLeavingVolumePerStep[i]<<endl;
	}
}

void CellMigration::updateVolumesWithMigration(vector <ShapeBase*>& Elements){
	int counter = 0;
	for( vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		double netVolumeChange = listOfGainedVolumePerStep[counter] - listOfLeavingVolumePerStep[counter];
		//I should apply a multiplicative growth rate such that the net volume change of this element will be equal to "netVolumeChange"
		double determinantToBe = ( (*itElement)->GrownVolume + netVolumeChange)/(*itElement)->GrownVolume;


		//Now increase the growth increment for this step:
		gsl_matrix* migrationIncrement = gsl_matrix_calloc(3,3);
		if((*itElement)->tissueType == 2){
			//linker elements should grow in 3D:
			double incrementX = pow(determinantToBe,1.0/3.0);
			gsl_matrix_set(migrationIncrement,0,0,incrementX);
			gsl_matrix_set(migrationIncrement,1,1,incrementX);
			gsl_matrix_set(migrationIncrement,2,2,incrementX);
		}
		else{
			//The element is columnar of peripodial,
			//Now I distribute this among x & y only:
			double incrementX = pow(determinantToBe,0.5);
			gsl_matrix_set(migrationIncrement,0,0,incrementX);
			gsl_matrix_set(migrationIncrement,1,1,incrementX);
			gsl_matrix_set(migrationIncrement,2,2,1);
		}
		if ((*itElement)->Id == 1076 || (*itElement)->Id == 1929 ){
			cout<<" Element: "<<(*itElement)->Id<<" netVolumeChange: "<<netVolumeChange<<" det to be: "<<determinantToBe<<" grown volume before update: "<<(*itElement)->GrownVolume<<endl;
			(*itElement)->displayMatrix(migrationIncrement,"migrationIncrement");
		}
		(*itElement)->addMigrationIncrementToGrowthIncrement(migrationIncrement);
		//reset volume change due to migration
		listOfGainedVolumePerStep[counter] = 0.0;
		listOfLeavingVolumePerStep[counter] = 0.0;
		//free the matrix:
		gsl_matrix_free(migrationIncrement);
		counter++;
	}
}

void CellMigration::updateMigratingElements(vector <ShapeBase*>& Elements){
	int counter = 0;
	vector<int> elementsToUpdateThisStep;
	for( vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if (listOfRateFractions[counter]>0 && (*itElement)->getCellMigration()){
			int nNeigs = elementConnectivityMap[counter].size();
			for (int i=0; i< nNeigs; ++i){
				int currNeig = elementConnectivityMap[counter][i];
				if (!Elements[currNeig]->getCellMigration()){
					elementsToUpdateThisStep.push_back(currNeig);
				}
			}
		}
		counter++;
	}
	int n = elementsToUpdateThisStep.size();
	for (int i=0; i<n; ++i){
		Elements[elementsToUpdateThisStep[i]]->setCellMigration(true);
	}
}

