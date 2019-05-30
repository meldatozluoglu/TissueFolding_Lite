using namespace std;

#include "Simulation.h"
#include <vector>

int main(int argc, char * argv[])
{
	std::shared_ptr<Simulation> Sim01 = std::make_shared<Simulation>();
	Sim01->displayIsOn = false;
	if (argc<2){
		Sim01->DisplaySave = false;
		cerr<<"Using default settings"<<endl;
	}
	else{
		bool Success = Sim01->readExecutableInputs(argc, argv);
		if (!Success){
			cerr<<"Error in input to executable"<<endl;
			return 1;
		}
	}
	if (Sim01->DisplaySave){
		cerr<<"This is the executable for running the model without display"<<endl;
		return true;
	}
	else{
		Sim01->initiateSystem();
		for (const auto& itElement : Sim01->Elements){
           		itElement->updatePositions(Sim01->Nodes);
		}
		cout<<"Initiating simulation in the background"<<endl;
		while (Sim01->currSimTimeSec <= Sim01->SimLength){
			//cout<<"running step: "<<Sim01.timestep<<", this is time: "<<Sim01.timestep*Sim01.dt<<" sec"<<endl;
			bool Success = Sim01->runOneStep();
			if (Success == false ){
				break;
			}
		}
		Sim01->wrapUpAtTheEndOfSimulation();
	}
	cout<<"Finished Simulation"<<endl;
}

