#include <iostream>
#include <QtGui>
#include <QtWidgets>
#include <QGraphicsWidget>
#include <QPushButton>
#include <QLabel>
#include <QLineEdit>
#include "MainWindow.h"
#include "GLWidget.h"
#include <sstream>

class MainWindow;

MainWindow::MainWindow(std::shared_ptr<Simulation> Sim01)
 {

	interatorForPressure = 0;

	MainScene = new QGraphicsScene;
    MainScene->setSceneRect(0, 0, 1000, 480);
    MainGrid = new QGridLayout;
    CentralWidget = new QWidget;
    this->Sim01 = Sim01;
    Sim01->calculateBoundingBox();
    nCoordBox = 6;

    setWindowTitle(tr("Tissue Origami"));
    generateControlPanel();
    setUpGLWidget();
    setUpCentralWidget();

    simulationStartClock = std::clock();
    simulationStartTime = time(0);
    displayedSimulationLength = false;

    MainScene->update();

    timer = new QTimer(this);
    connect(timer, SIGNAL(timeout()), this, SLOT(timerSimulationStep()));
	timer->start(0);
 }

MainWindow::~MainWindow(){
	delete MainGLWidget;
	delete MainGrid;
	delete MainScene;
}

void MainWindow::generateControlPanel(){
    /**
     * This function will generate the layout of the whole button control panel.
     * All layouts will be contained in a vertical box layout, MainWindow#ControlPanelMainHBox.
     */

	ControlPanelMainHBox = new QVBoxLayout();
	ControlPanelMainHBox->setSpacing(2);

	//Generating Selection Display Panel:
    /** The grid of displays to contain all the selected item information is contained in
     * MainWindow#SelectionDisplayGrid, and is organised in function MainWindow#setUpSelectionDisplayGrid.
     * This grid is at the top of the MainWindow#ControlPanelMainHBox.
     */
	QGridLayout *SelectionDisplayGrid = new QGridLayout;
	setUpSelectionDisplayGrid(SelectionDisplayGrid);
	ControlPanelMainHBox->addLayout(SelectionDisplayGrid,Qt::AlignTop);

	//Generating project display options panel:
    /** The grid of control tools to contain all the display options for the specifications
     * of current simulaiton are contained in MainWindow#ProjectDisplayOptionsGrid, and are organised
     * in funciton MainWindow#setUpProjectDisplayOptionGrid. This grid is the secod grid from the top of the
     * MainWindow#ControlPanelMainHBox.
     */
	QGridLayout *ProjectDisplayOptionsGrid = new QGridLayout;
	setUpProjectDisplayOptionGrid(ProjectDisplayOptionsGrid);
    ControlPanelMainHBox->addLayout(ProjectDisplayOptionsGrid,Qt::AlignTop);

	//Generating view options Panel:
    /** The grid of display output options are contained in MainWindow#ViewOptionsGrid and are set up in
     * function MainWindow#VsetUpViewOptionsGrid.
     */
	QGridLayout *ViewOptionsGrid = new QGridLayout;
	setUpViewOptionsGrid(ViewOptionsGrid);
	ControlPanelMainHBox->addLayout(ViewOptionsGrid,Qt::AlignBottom);

	//Generating the quit button:
    /** The very bottom of the MainWindow#ControlPanelMainHBox contains the Quit button of the user interaface,
     * in its own horizontal layout box, BottomLineBox. The button is connected to close the main window, which will call
     * the destructor and tidy up the memory followed by the terminaiton of the user interface.
     */
	QHBoxLayout *BottomLineBox = new QHBoxLayout; // the bottom line will include quit button only for now
	QPushButton	*QuitButton = new QPushButton("Quit",this);
	QuitButton->setFixedWidth(100);
    connect(QuitButton, SIGNAL(clicked()), this,SLOT(close()));
	BottomLineBox->addWidget(QuitButton,Qt::AlignBottom| Qt::AlignRight);
	ControlPanelMainHBox->addLayout(BottomLineBox,Qt::AlignBottom);

    /** The MainWindow#ControlPanelMainHBox is itself added onto the MainWindow#MainGrid, which will contain the
     * OpenGL window (MainWindow#MainGLWidget, set up in MainWindow#setUpGLWidget) and the control panel.
     */
	MainGrid->addLayout(ControlPanelMainHBox,0,1,Qt::AlignLeft);
	MainGrid->setColumnStretch(1,-2);
}

void MainWindow::setUpGLWidget(){
    /**
     * This funciton will generate the OpenGL window. The initail values from the user controlled display variables
     * such as the physical propoerty and strain display ranges, are initited. Their updates are linked to the parametrs of the
     * MainWindow#MainGLWidget through signals.
     */

	MainGLWidget = new GLWidget();
	MainGLWidget->Sim01 = Sim01;
	MainGLWidget->DisplayStrainRange[0] = StrainSpinBoxes[0]->value();
	MainGLWidget->DisplayStrainRange[1] = StrainSpinBoxes[1]->value();
    for (int i =0; i<4; ++i){
    	MainGLWidget->DisplayPysPropRange[i][0] = PysPropSpinBoxes[0]->value();
    	MainGLWidget->DisplayPysPropRange[i][0] = PysPropSpinBoxes[1]->value();
    }
    /** Once generated, the MainWindow#MainGLWidget is contained in the MainWindow#MainGrid. and the stretch of the first column
     * containing the OpenGL window is set high, to ensure the display window is as large as possible next to the control panel.
     */
    MainGrid->addWidget(MainGLWidget,0,0,Qt::AlignLeft);
	MainGrid->setColumnStretch(0,10);
}

void MainWindow::setUpCentralWidget(){
    /**
     * All the contents of the MainWindow#MainGrid are containded in the MainWindow#CentralWidget. This widget controls the signal
     * coupling to the parameter updated of the MainWindow#MainGLWidget.
     */
    setCentralWidget(CentralWidget);
    CentralWidget->setParent(this);
    CentralWidget->setLayout(MainGrid);
    connect(MainGLWidget, SIGNAL(SelectedItemChanged()), this, SLOT(SelectedItemChange()));
    connect(MainGLWidget, SIGNAL(NeedToClearManualElementSelection()), this, SLOT(ManualElementSelectionReset()));
}

void MainWindow::setUpSelectionDisplayGrid(QGridLayout *SelectionDisplayGrid){
    /**
     * The item seleciton display grid is set up in this function, first the text titles of the
     * display boxes are assigned in MainWinfow#setItemSelectionTitles. Then the gird of boxes
     * to display the coordinantes of the selected element are sey up in MainWindow#setCoordBoxes. Finally,
     * the box to idsplay the selected item Id, which is used as both an input source for selection and display
     * box for selection of elements with mouse clicks, is set up in MainWindow#setSelectionByIdSelection.
     */
	QFont boldFont("SansSerif", 10, QFont::Bold,true);
	QFont font("SansSerif", 10);
	setItemSelectionTitles(font, boldFont, SelectionDisplayGrid);
	setCoordBoxes(font, boldFont, SelectionDisplayGrid);
    setSelectionByIdSelection(font, SelectionDisplayGrid);
}

void MainWindow::setSelectionByIdSelection(QFont font1, QGridLayout *SelectionDisplayGrid){
    /**
     * The box will diplay the selected element id, hence its placeholder text will contain the range of values it can have, from
     * the number of elements of the simulation to be displayed. The manual input value is also validated with same criteria. Any changes
     * of the id text is connected to the MainWindow#manualElementSelection function to update selection information.
     */
    QLabel *ElementSelectTitle = new QLabel("Select <br> Element:");
    ElementSelectTitle->setFont(font1);
    ElementSelectBox = new QLineEdit();
    ElementSelectBox->setPlaceholderText ( QString("# 0-%1").arg(Sim01->Elements.size()-1) );
    ElementSelectBox->setFont(font1);
    ElementSelectBox->setStyleSheet("background-color: white");
    ElementSelectBox->setFixedWidth(50);
    ElementSelectBox->setValidator( new QIntValidator(0, Sim01->Elements.size()-1, this) );
    SelectionDisplayGrid->addWidget(ElementSelectTitle,0,4,1,1,Qt::AlignLeft);
    SelectionDisplayGrid->addWidget(ElementSelectBox,1,4,1,1,Qt::AlignLeft);
    connect(ElementSelectBox, SIGNAL(textChanged(const QString &)), this, SLOT(manualElementSelection(const QString &)));
}
void MainWindow::setUpProjectDisplayOptionGrid(QGridLayout *ProjectDisplayOptionsGrid){
    /**
     * The checkboxes, combo boxes and spinners of system property display selections are set up via
     * MainWindow#setStrainDisplayMenu and MainWindow#setPysPropDisplayMenu. The size preferences of the grid is set up such that
     * to push the upper grid of selection display to the top of the window:
     */
    QFont boldFont("SansSerif", 10, QFont::Bold,true);
	QFont font("SansSerif", 10);
	setStrainDisplayMenu(ProjectDisplayOptionsGrid);
	setPysPropDisplayMenu(ProjectDisplayOptionsGrid);
	setDisplayPreferences(ProjectDisplayOptionsGrid);
	int widthPhysProp = PysPropComboBox->minimumSizeHint().width();
	int widthStrain = StrainComboBox->minimumSizeHint().width();
	if (widthPhysProp>widthStrain){
        StrainComboBox->setMinimumWidth(widthPhysProp);
	}
	else{
        PysPropComboBox->setMinimumWidth(widthStrain);
	}
	ProjectDisplayOptionsGrid->setRowStretch(10,10);
}

void MainWindow::setUpViewOptionsGrid(QGridLayout *ViewOptionsGrid){
    /**
     * This grid currently contains the symmetricity display button, whihc controls the symettricity display booleans of the
     * MainWindow#MainGLWidget through the slot GLWidget#pdateDrawSymmetricityViewToggle.
     */
	ViewOptionsGrid->setRowStretch(0,10);
	ViewOptionsGrid->setColumnStretch(10,10);
	SymmetricityDisplayButton = new QPushButton("Hide \n Symmetric",this);
	connect(SymmetricityDisplayButton, SIGNAL(clicked()), this,SLOT(updateDrawSymmetricityViewToggle()));
    SymmetricityDisplayButton->setFixedWidth(100);
    ViewOptionsGrid->addWidget(SymmetricityDisplayButton,4,4,1,1);
}

void MainWindow::setStrainDisplayMenu(QGridLayout *ProjectDisplayOptionsGrid){
    /** The strain display options will contain first the checkbox to decide if the starins are
     * displayed. The checkbox state change is connected to MainWindow#updateStrainCheckBox, which will
     * enable the control panel for the strain option inputs once checked, as well as unchecking the pysical property discplay menu.
     * These tow display options are mutually exclusive as they both colourcode the tissue with the parameter
     * of interest.
     */
    DisplayCheckBoxes[0] = new QCheckBox("Strain");
	DisplayCheckBoxes[0]->setChecked(false);
	connect(DisplayCheckBoxes[0] , SIGNAL(stateChanged(int)),this,SLOT(updateStrainCheckBox(int)));

	StrainComboBox = new QComboBox();
    /** The combo box lists the available strain displays, and the updates of the combobox change which strain is read from the
     *  simulation data through the funciton MainWindow#updateStrain.
     */
	StrainComboBox->addItem("Strain in DV");
	StrainComboBox->addItem("Strain in AP");
	StrainComboBox->addItem("Strain in AB");
	StrainComboBox->addItem("Shear in xy");
	StrainComboBox->addItem("Shear in yz");
	StrainComboBox->addItem("Shear in xz");
	StrainComboBox->setEnabled(false);
	connect(StrainComboBox , SIGNAL(currentIndexChanged(int)),this,SLOT(updateStrain(int)));

    StrainSpinBoxes[0] = new  QDoubleSpinBox();
    StrainSpinBoxes[1] = new  QDoubleSpinBox();
    StrainSpinBoxes[0]->setRange( -1.0, -0.00 );
    StrainSpinBoxes[0]->setSingleStep( 0.005 );
    StrainSpinBoxes[0]->setValue ( -0.35 );
    StrainSpinBoxes[0]->setDecimals(5);
    StrainSpinBoxes[0]->setEnabled(false);
    StrainSpinBoxes[1]->setRange( 0.0, 1.0 );
    StrainSpinBoxes[1]->setSingleStep( 0.005 );
    StrainSpinBoxes[1]->setValue( 0.35 );
    StrainSpinBoxes[1]->setDecimals(3);
    StrainSpinBoxes[1]->setEnabled(false);
	connect(StrainSpinBoxes[0] , SIGNAL(valueChanged(double)),this,SLOT(updateStrainSpinBoxes()));
    connect(StrainSpinBoxes[1], SIGNAL(valueChanged(double)), this, SLOT(updateStrainSpinBoxes()));

    ProjectDisplayOptionsGrid->addWidget(DisplayCheckBoxes[0],0,0,1,2,Qt::AlignLeft);
    ProjectDisplayOptionsGrid->addWidget(StrainComboBox,1,0,1,2,Qt::AlignLeft);
    ProjectDisplayOptionsGrid->addWidget(StrainSpinBoxes[0],2,0,1,1,Qt::AlignLeft);
    ProjectDisplayOptionsGrid->addWidget(StrainSpinBoxes[1],2,1,1,1,Qt::AlignLeft);
}

void MainWindow::setPysPropDisplayMenu(QGridLayout *ProjectDisplayOptionsGrid){
    /** The physical property display options will contain first the checkbox to decide if the starins are
     * displayed. The checkbox state change is connected to MainWindow#updatePysCheckBox, which will
     * enable the control panel for the physical property option inputs once checked, as well as unchecking the strain discplay menu.
     * These tow display options are mutually exclusive as they both colourcode the tissue with the parameter
     * of interest.
     */
	DisplayCheckBoxes[1] = new QCheckBox("Physical Properties");
	DisplayCheckBoxes[1]->setChecked(false);
	connect(DisplayCheckBoxes[1] , SIGNAL(stateChanged(int)),this,SLOT(updatePysCheckBox(int)));
    /** The combo box lists the available physical property displays, and the updates of the combobox change which physical property
     *  is read from the simulation data through the funciton MainWindow#updatePysProp.
     */
	PysPropComboBox = new QComboBox();
	PysPropComboBox->addItem("External Viscosity");
	PysPropComboBox->addItem("Internal Viscosity");
	PysPropComboBox->addItem("Young Modulus");
	PysPropComboBox->addItem("Poisson Ratio");
    PysPropComboBox->addItem("Growth Rate");
	PysPropComboBox->setEnabled(false);
	connect(PysPropComboBox , SIGNAL(currentIndexChanged(int)),this,SLOT(updatePysProp(int)));

	PysPropSpinBoxes[0] = new  QDoubleSpinBox();
	PysPropSpinBoxes[1] = new  QDoubleSpinBox();
	PysPropSpinBoxes[0]->setRange ( 0, 10 );
	PysPropSpinBoxes[0]->setSingleStep( 0.1 );
	PysPropSpinBoxes[0]->setValue ( 3 );
	PysPropSpinBoxes[0]->setEnabled(false);
	PysPropSpinBoxes[1]->setRange( 0, 10.0 );
	PysPropSpinBoxes[1]->setSingleStep( 0.1 );
	PysPropSpinBoxes[1]->setValue( 5.0 );
	PysPropSpinBoxes[1]->setEnabled(false);
    connect(PysPropSpinBoxes[0], SIGNAL(valueChanged (double)), this, SLOT(updatePysPropSpinBoxes()));
    connect(PysPropSpinBoxes[1], SIGNAL(valueChanged (double)), this, SLOT(updatePysPropSpinBoxes()));

    ProjectDisplayOptionsGrid->addWidget(DisplayCheckBoxes[1],0,2,1,2,Qt::AlignLeft);
    ProjectDisplayOptionsGrid->addWidget(PysPropComboBox,1,2,1,2,Qt::AlignLeft);
    ProjectDisplayOptionsGrid->addWidget(PysPropSpinBoxes[0],2,2,1,1,Qt::AlignLeft);
    ProjectDisplayOptionsGrid->addWidget(PysPropSpinBoxes[1],2,3,1,1,Qt::AlignLeft);
}

void MainWindow::setItemSelectionTitles(QFont font, QFont boldFont, QGridLayout *SelectionDisplayGrid){
	QLabel *PanelTitle = new QLabel("Selected Item Properties");
	PanelTitle->setFont(boldFont);
	QLabel *NameTitle = new QLabel("Name:");
	NameTitle->setFont(boldFont);
	NameBox = new QLineEdit();
	NameBox->setPlaceholderText ( "No input" );
	NameBox->setReadOnly(true);
	NameBox->setFont(font);

	QLabel *NodeIdTitle = new QLabel("id");
	NodeIdTitle ->setFont(boldFont);
	QLabel *CoordTitlex = new QLabel("x");
	CoordTitlex ->setFont(boldFont);
	QLabel *CoordTitley = new QLabel("y");
	CoordTitley ->setFont(boldFont);
	QLabel *CoordTitlez = new QLabel("z");
	CoordTitlez ->setFont(boldFont);

	SelectionDisplayGrid->addWidget(PanelTitle,0,0,1,2,Qt::AlignLeft);
	SelectionDisplayGrid->addWidget(NameTitle,1,0,1,1,Qt::AlignLeft);
	SelectionDisplayGrid->addWidget(NameBox,1,1,1,2,Qt::AlignLeft);

	SelectionDisplayGrid->addWidget(NodeIdTitle,2,1,1,1,Qt::AlignHCenter);
	SelectionDisplayGrid->addWidget(CoordTitlex,2,2,1,1,Qt::AlignHCenter);
	SelectionDisplayGrid->addWidget(CoordTitley,2,3,1,1,Qt::AlignHCenter);
	SelectionDisplayGrid->addWidget(CoordTitlez,2,4,1,1,Qt::AlignHCenter);
}

void MainWindow::setCoordBoxes(QFont font, QFont boldFont, QGridLayout *SelectionDisplayGrid){
	CoordLabel_n[0] = new QLabel("Node 1");
	CoordLabel_n[1] = new QLabel("Node 2");
	CoordLabel_n[2] = new QLabel("Node 3");
	CoordLabel_n[3] = new QLabel("Node 4");
	CoordLabel_n[4] = new QLabel("Node 5");
	CoordLabel_n[5] = new QLabel("Node 6");
	for (int i = 0 ;i<nCoordBox; ++i){
		CoordLabel_n[i] ->setFont(boldFont);
		CoordBox_id[i] = new QLineEdit();
		CoordBox_x[i] = new QLineEdit();
		CoordBox_y[i] = new QLineEdit();
		CoordBox_z[i] = new QLineEdit();
		CoordBox_id[i]->setPlaceholderText( "No input" );
		CoordBox_x[i]->setPlaceholderText( "No input" );
		CoordBox_y[i]->setPlaceholderText( "No input" );
		CoordBox_z[i]->setPlaceholderText( "No input" );
		CoordBox_id[i]->setReadOnly(true);
		CoordBox_x[i]->setReadOnly(true);
		CoordBox_y[i]->setReadOnly(true);
		CoordBox_z[i]->setReadOnly(true);
		CoordBox_id[i]->setFont(font);
		CoordBox_x[i]->setFont(font);
		CoordBox_y[i]->setFont(font);
		CoordBox_z[i]->setFont(font);
        CoordBox_id[i]->setFixedWidth(50);
        CoordBox_x[i]->setFixedWidth(50);
        CoordBox_y[i]->setFixedWidth(50);
        CoordBox_z[i]->setFixedWidth(50);
		SelectionDisplayGrid->addWidget(CoordLabel_n[i],i+3,0,1,1,Qt::AlignLeft);
		SelectionDisplayGrid->addWidget(CoordBox_id[i],i+3,1,1,1,Qt::AlignLeft);
		SelectionDisplayGrid->addWidget(CoordBox_x[i],i+3,2,1,1,Qt::AlignLeft);
		SelectionDisplayGrid->addWidget(CoordBox_y[i],i+3,3,1,1,Qt::AlignLeft);
		SelectionDisplayGrid->addWidget(CoordBox_z[i],i+3,4,1,1,Qt::AlignLeft);
	}
}

void MainWindow::setDisplayPreferences(QGridLayout *ProjectDisplayOptionsGrid){
    //draw pipette aspiration CheckBox
    DisplayPreferencesCheckBoxes[0] = new QCheckBox("Display Pipette");
	DisplayPreferencesCheckBoxes[0]->setChecked(false);
    ProjectDisplayOptionsGrid->addWidget(DisplayPreferencesCheckBoxes[0],3,0,1,2,Qt::AlignLeft);  // display pipette
}

void  MainWindow::updateDrawSymmetricityViewToggle(){
	if (MainGLWidget->drawSymmetricity){
		MainGLWidget->drawSymmetricity = false;
		SymmetricityDisplayButton->setText("Show \n Symmetric");
	}
	else{
		MainGLWidget->drawSymmetricity = true;
		SymmetricityDisplayButton->setText("Hide \n Symmetric");
	}
}

void  MainWindow::updateDisplayPipette(int s){
    if ( s == 2 )
        MainGLWidget->displayPipette = true;
    else
        MainGLWidget->displayPipette = false;
}

void MainWindow::updateStrain(int s){
    MainGLWidget->updateVBO = true;
	MainGLWidget->StrainToDisplay = s;
	MainGLWidget->update();
}

void MainWindow::updatePysProp(int s){
    MainGLWidget->updateVBO = true;
	MainGLWidget->PysPropToDisplay = s;
	float low = MainGLWidget->DisplayPysPropRange[MainGLWidget->PysPropToDisplay][0];
	float high = MainGLWidget->DisplayPysPropRange[MainGLWidget->PysPropToDisplay][1];
	float min[2] = {MainGLWidget->DisplayPysPropBounds[MainGLWidget->PysPropToDisplay][0],MainGLWidget->DisplayPysPropBounds[MainGLWidget->PysPropToDisplay][2]};
	float max[2] = {MainGLWidget->DisplayPysPropBounds[MainGLWidget->PysPropToDisplay][1],MainGLWidget->DisplayPysPropBounds[MainGLWidget->PysPropToDisplay][3]};
	int   decimals = MainGLWidget->DisplayPysPropDecimals[MainGLWidget->PysPropToDisplay];
	float step = MainGLWidget->DisplayPysPropSteps[MainGLWidget->PysPropToDisplay];
	PysPropSpinBoxes[0]->setRange( min[0], max[0] );
	PysPropSpinBoxes[0]->setValue (low);
	PysPropSpinBoxes[0]->setDecimals(decimals);
	PysPropSpinBoxes[0]->setSingleStep( step );
	PysPropSpinBoxes[1]->setRange( min[1], max[1] );
	PysPropSpinBoxes[1]->setValue (high);
	PysPropSpinBoxes[1]->setDecimals(decimals);
	PysPropSpinBoxes[1]->setSingleStep( step );
}
void MainWindow::updateStrainCheckBox(int s){

    MainGLWidget->updateVBO = true;
	if (s == 0){
        /** If the check lead to disabling of the strain display, the control panel items
         * and disabled.
         */
		MainGLWidget->StrainToDisplay = -1;
        StrainComboBox->setEnabled(false);
		StrainSpinBoxes[0]->setEnabled(false);
		StrainSpinBoxes[1]->setEnabled(false);
		MainGLWidget->DisplayStrains = false;
	}
	else{
        /** If the check lead to enabling of the strain display, the control panel items
         * and enabled, and the physical property display checkbox is unchecked. This will
         * emit its own change signal and make necessary changes on physical property control panel.
         */
		DisplayCheckBoxes[1]->setChecked(false);
		MainGLWidget->StrainToDisplay = StrainComboBox->currentIndex ();
		StrainComboBox->setEnabled(true);
		StrainSpinBoxes[0]->setEnabled(true);
		StrainSpinBoxes[1]->setEnabled(true);
		MainGLWidget->DisplayStrains = true;
	}
}

void MainWindow::updatePysCheckBox(int s){
    MainGLWidget->updateVBO = true;
	if (s == 0){
        /** If the check lead to disabling of the physical property display, the control panel items
         * and disabled.
         */
        MainGLWidget->PysPropToDisplay = -1;
		PysPropComboBox->setEnabled(false);
		PysPropSpinBoxes[0]->setEnabled(false);
		PysPropSpinBoxes[1]->setEnabled(false);
		MainGLWidget->DisplayPysProp = false;
	}
	else{
        /** If the check lead to enabling of the physical property display, the control panel items
         * and enabled, and the strain display checkbox is unchecked. This will
         * emit its own change signal and make necessary changes on strain control panel.
         */
        DisplayCheckBoxes[0]->setChecked(false);
		MainGLWidget->PysPropToDisplay = PysPropComboBox->currentIndex();
		float low = MainGLWidget->DisplayPysPropRange[MainGLWidget->PysPropToDisplay][0];
		float high = MainGLWidget->DisplayPysPropRange[MainGLWidget->PysPropToDisplay][1];
		float min[2] = {MainGLWidget->DisplayPysPropBounds[MainGLWidget->PysPropToDisplay][0],MainGLWidget->DisplayPysPropBounds[MainGLWidget->PysPropToDisplay][2]};
		float max[2] = {MainGLWidget->DisplayPysPropBounds[MainGLWidget->PysPropToDisplay][1],MainGLWidget->DisplayPysPropBounds[MainGLWidget->PysPropToDisplay][3]};
		int decimals = MainGLWidget->DisplayPysPropDecimals[MainGLWidget->PysPropToDisplay];
		float step =  MainGLWidget->DisplayPysPropSteps[MainGLWidget->PysPropToDisplay];

		PysPropSpinBoxes[0]->setRange( min[0], max[0] );
		PysPropSpinBoxes[0]->setValue (low);
		PysPropSpinBoxes[0]->setDecimals(decimals);
		PysPropSpinBoxes[0]->setSingleStep( step );
		PysPropSpinBoxes[1]->setRange( min[1], max[1] );
		PysPropSpinBoxes[1]->setValue (high);
		PysPropSpinBoxes[1]->setDecimals(decimals);
		PysPropSpinBoxes[1]->setSingleStep( step );
		PysPropComboBox->setEnabled(true);
		PysPropSpinBoxes[0]->setEnabled(true);
		PysPropSpinBoxes[1]->setEnabled(true);
		MainGLWidget->DisplayPysProp = true;
	}
}

void MainWindow::updateStrainSpinBoxes(){
    MainGLWidget->updateVBO = true;
	MainGLWidget->DisplayStrainRange[0] = StrainSpinBoxes[0]->value();
	MainGLWidget->DisplayStrainRange[1] = StrainSpinBoxes[1]->value();
}

void MainWindow::updatePysPropSpinBoxes(){
    MainGLWidget->updateVBO = true;
	MainGLWidget->DisplayPysPropRange[MainGLWidget->PysPropToDisplay][0] = PysPropSpinBoxes[0]->value();
	MainGLWidget->DisplayPysPropRange[MainGLWidget->PysPropToDisplay][1] = PysPropSpinBoxes[1]->value();
}

void MainWindow::SelectedItemChange(){
    QString tmpstring = QString::fromStdString(MainGLWidget->SelectedItemName);
    NameBox->setText(tmpstring);
    for (int i = 0 ;i<nCoordBox; ++i){
    	 CoordBox_id[i]->setText ( "" );
    	 CoordBox_x[i]->setText ( "" );
    	 CoordBox_y[i]->setText ( "" );
    	 CoordBox_z[i]->setText ( "" );
    	 CoordBox_id[i]->setEnabled(false);
    	 CoordBox_x[i]->setEnabled(false);
    	 CoordBox_y[i]->setEnabled(false);
    	 CoordBox_z[i]->setEnabled(false);
    	 if ((signed int)MainGLWidget->SelectedPos.size()>i*3){
    		 CoordBox_id[i]->setText ( MainGLWidget->SelectedId[i] );
			 CoordBox_x[i]->setText ( MainGLWidget->SelectedPos[i*3] );
			 CoordBox_y[i]->setText ( MainGLWidget->SelectedPos[i*3+1] );
			 CoordBox_z[i]->setText ( MainGLWidget->SelectedPos[i*3+2] );
			 CoordBox_id[i]->setEnabled(true);
			 CoordBox_x[i]->setEnabled(true);
			 CoordBox_y[i]->setEnabled(true);
			 CoordBox_z[i]->setEnabled(true);
		}
    }
 }

void MainWindow::manualElementSelection(const QString &newValue){
	MainGLWidget->manualElementSelection(newValue.toInt());
}

void MainWindow::ManualElementSelectionReset(){
	ElementSelectBox->blockSignals(true);
	ElementSelectBox->setValidator( new QIntValidator(0, Sim01->Elements.size()-1, this) );
	ElementSelectBox->setPlaceholderText ( QString("# 0 to %1").arg(Sim01->Elements.size()-1) );
	ElementSelectBox->setText("");
	ElementSelectBox->blockSignals(false);
}

void MainWindow::timerSimulationStep(){
	bool 	automatedSave =false ;
    bool 	slowstepsOnDisplay = false;
	bool 	slowstepsOnRun = false;
	int 	slowWaitTime = 10;


	if (Sim01->DisplaySave){
        /** If the user interface is called with the option to display a previously saved simulation, the data will be read
         * from the save files until the flag Simulation#reachedEndOfSaveFile is set to true.
         */
		if (Sim01->timestep == 0){
			if( automatedSave ){
				Sim01->assignTips();
            }
		}
		if(!Sim01->reachedEndOfSaveFile){
            if (takeScreenShotAfterUpdate){
                bool filesaved = takeScreenshot();
                if (filesaved == false){
                    cerr<< "Cannot create image file, will not save images! "<<std::endl;
                    Sim01->saveImages = false;
                }
                takeScreenShotAfterUpdate=false;
            }
            /** The data will be read via the function Simulation#updateOneStepFromSave. Upon each update of the simulation data,
             * vertex buffer will be set to update.
             */
			Sim01->updateOneStepFromSave();
            MainGLWidget->updateVBO = true;
            MainGLWidget->update();
            QTime dieTime= QTime::currentTime().addSecs(1);
            while( QTime::currentTime() < dieTime ){
			    QCoreApplication::processEvents(QEventLoop::AllEvents, 1);
            }
            /** If the display is saving images, then the toggle to inform taking a screenshot of the GLWidget after this step is displayed
             * will be set to on after each read step.
             */
			if (Sim01->saveImages){
                takeScreenShotAfterUpdate = true;
			}
            /** If the display is set to process slowly, then the process is waited accordingly.
             */
			if (slowstepsOnDisplay){
				QTime dieTime= QTime::currentTime().addSecs(slowWaitTime);
				while( QTime::currentTime() < dieTime ){
					QCoreApplication::processEvents(QEventLoop::AllEvents, slowWaitTime);
				}
			}
		}
		else{
            /** Once the end of the save file is reached, if the toggle for automatedSave is set true, the window will close itself.
             * This toggle is exists to enable users to process multiple simulations sequantially and save thier screen putputs within a script.\n
             */
			if (automatedSave){
				MainGLWidget->close();
				close();
			}
		}
	}
	else{
        /** If the user interface is called with the option to run a new simulation, either from scratch or continuing form an existing save file,
         * then simulation will be run for each time step via Simulation#runOneStep, until the desired simulation length is reached, such that
         * Simulation#currSimTimeSec is larger or equal to Simulation#SimLength.
         */

		if (Sim01->currSimTimeSec <= Sim01->SimLength){
            if (takeScreenShotAfterUpdate){
                bool filesaved = takeScreenshot();
                if (filesaved == false){
                    cerr<< "Cannot create image file, will not save images! "<<std::endl;
                    Sim01->saveImages = false;
                }
                takeScreenShotAfterUpdate=false;
            }

            cout<<"started step: "<<Sim01->currSimTimeSec<<" length: "<<Sim01->SimLength<<std::endl;
			bool Success = Sim01->runOneStep();
            /** The simulation data will be updated at the end of each ran step, and
             * vertex buffer will be set to update accordingly.
             */
            MainGLWidget->updateVBO = true;
            MainGLWidget->update();
			if (slowstepsOnRun){
				QTime dieTime= QTime::currentTime().addSecs(slowWaitTime);
				while( QTime::currentTime() < dieTime ){
					QCoreApplication::processEvents(QEventLoop::AllEvents, slowWaitTime);
				}
			}
			if (Sim01->saveImages && Sim01->timestep%Sim01->imageSaveInterval == 0){
                /** If the simulation options require saving screenshots of the GLWidget, and the current simulation time is
                 * coincidin gwith the desired image save interval (Simulation#imageSaveInterval), then the toggle to save
                 * the screenshot once the display is updated is toggled on.
                 */
                takeScreenShotAfterUpdate = true;
			}
			if (Success == false ){
                /** If Simulation#runOneStep did not retun success, then ther must be a flipped element in the system, and the user is notified.
                 * The simulation is forced to finish by setting thecurrent time higher than desired dimulation length,
                 * to avoid unidentifiable crashes in matrix opertion failure in the next step.
                 */
                std::cerr<<"there is a flipped element, I am not continuing simulation"<<std::endl;
				Sim01->timestep = 2*Sim01->SimLength;
			}
		}
        else if(!displayedSimulationLength){
            /** If the simulation has reached its final length, the first screen update after the finalisation of the simulation
             * tidies up the final simulation parameters, takes one last screenshot, and reports simulation process time.
             * The toggle MainWindow#displayedSimulationLength is toggled to ensure this simulation end reporting is carried out once.
             */

        	displayedSimulationLength = true;
            Sim01->wrapUpAtTheEndOfSimulation();
            if (Sim01->saveImages){
                bool filesaved = takeScreenshot();
                if (filesaved == false){
                    cerr<< "Cannot create image file, will not save images! "<<std::endl;
                    Sim01->saveImages = false;
                }
            }
            double durationClock = ( std::clock() - simulationStartClock ) / (double) CLOCKS_PER_SEC;
            double durationTime = std::difftime(std::time(0), simulationStartTime);
            std::cout<<"Simulation time: "<<durationTime<<" sec, Simulation clock: "<<durationClock<<" sec"<<std::endl;
        }
    }
    MainGLWidget->update();
}

bool MainWindow::takeScreenshot(){
    /** The screenshots directory should be created before the executable is started. Ideally, the simulations are
     * managed through scripts therefore the file management is left to the user at that stage.
     */
    std::cout<<"taking screenshot"<<std::endl;
	QPixmap originalPixmap;
	QScreen *screen = qApp->primaryScreen();
	int x = MainGLWidget->geometry().x();
	int y = MainGLWidget->geometry().y();
	int w = MainGLWidget->width();
	int h = MainGLWidget->height();
	originalPixmap = screen->grabWindow(this->winId(),x,y,w,h);
	QString timepoint = QString("%1").arg(Sim01->timestep,6,10,QChar('0'));
	QString timeinsec = QString("%1").arg(Sim01->timestep*Sim01->dt);
    std::string saveScreenshotsDirectory = Sim01->saveDirectory + Sim01->saveScreenshotsDirectory;
    QString fileName = QString(saveScreenshotsDirectory.c_str())+"frame"+ timepoint +"-"+timeinsec+"sec.png";
    if (originalPixmap.save(fileName, "png") != true ){
        std::cerr<<"Cannot save to "<<fileName.toStdString()<<std::endl;
        std::cerr<<"Have you created the directory:  "<<saveScreenshotsDirectory<<std::endl;
        std::cerr<<" Directory creation is not automated to avoid accidental overwrites"<<std::endl;
        return false;
    }
    return true;
}
