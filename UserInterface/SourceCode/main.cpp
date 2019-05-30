#include <QGraphicsWidget>
#include <QtGui>
#include <QtWidgets>
#include "MainWindow.h"
#include "GLWidget.h"

#include "../TissueFolding/SourceCode/Simulation.h"
#include <vector>

class MainWindow;
class SurfaceBase;
class GLWidget;

//Simulation* Sim01;
int main(int argc, char **argv)
{
	bool Success = false;
    std::shared_ptr<Simulation> Sim01 = std::make_shared<Simulation>(); //  Sim01 is a unique_ptr that owns a new simulation

    //Sim01 = new Simulation();
	Sim01->displayIsOn = true;
	if (argc<2){
		Sim01->DisplaySave = false;
        cerr<<"Using default settings"<<std::endl;
		Success = true;
	}
	else{
		Success = Sim01->readExecutableInputs(argc, argv);
	}
	if (Success == 0 ){
        cout<<"Error in input to executable"<<std::endl;
		return true;
	}

    QGLFormat glFormat;
    glFormat.setVersion( 3, 2 );
    glFormat.setProfile( QGLFormat::CoreProfile );
    QGLFormat::setDefaultFormat(glFormat);

    QApplication app(argc, argv);

	if (Sim01->DisplaySave){
        cout<<"Initiating simulation display"<<std::endl;
		Success = Sim01->initiateSavedSystem();
	}
	else{
		Success = Sim01->initiateSystem();
		if (Success == 0 ){
            cout<<"System is not initiated successfully, terminating"<<std::endl;
			return true;
		}
		else{
            cout<<"system initiated"<<std::endl;
		}
        for (const auto& itElement : Sim01->Elements){
            itElement->updatePositions(Sim01->Nodes);
		}
	}

	if (Success == 0 ){
        cout<<"System is not initiated successfully, terminating"<<std::endl;
		return true;
	}

    MainWindow mw(Sim01);
	mw.show();
	mw.MainGLWidget->show();
	mw.raise();
    mw.setGeometry(50, 50, 900, 550); //mac
    mw.MainScene->update();

	return app.exec();
}




