/*
 * MainWindow.h
 *
 *  Created on: 18 Mar 2014
 *      Author: melda
 */

#ifndef MAINWINDOW_H_
#define MAINWINDOW_H_

#include <QtGui>
#include <QtWidgets>
#include <QGraphicsWidget>
#include <QMainWindow>
#include <QLineEdit>
#include <QTimer>
#include <QSlider>
#include <ctime>
class GLWidget;

#include "../TissueFolding/SourceCode/Simulation.h"

class MainWindow : public QMainWindow
 {
     Q_OBJECT

 public:
    MainWindow(std::shared_ptr<Simulation> Sim01);
    ~MainWindow();
    QGraphicsScene	*MainScene;
    QVBoxLayout		*ControlPanelMainHBox;
    QGridLayout		*MainGrid;
    GLWidget 		*MainGLWidget;
    std::shared_ptr<Simulation> Sim01;
	int interatorForPressure;


public slots:
    void 	SelectedItemChange();                                                               ///< When the selected item is changed, this function updates on screen displayed information.
    void 	manualElementSelection(const QString &);                                            ///< When the selected element is changed, this function updates on screen displayed information.
    void 	ManualElementSelectionReset();                                                      ///< This function resets the element selection informaton when no element is selected upon selection update.
    void 	timerSimulationStep();                                                              ///< This function iterates the simulation by one time step.
    void 	updateStrain(int);                                                                  ///< This function will update the strain to display, between DV, AP and AB axes (corresponding to x,y,z in the current coordinate system)
    void 	updateStrainCheckBox(int);                                                          ///< This function will update the strain checkbox, only one of Strain or Pysical property can be displayed at a time.
    void 	updateStrainSpinBoxes();                                                            ///< This function will update the strain spin box to select which strain is displayed.
    void 	updatePysProp(int s);                                                               ///< This function will update the physical property to display.
    void 	updatePysCheckBox(int);                                                             ///< This function will update the physical property checkbox, only one of Strain or Pysical property can be displayed at a time.
    void 	updatePysPropSpinBoxes();                                                           ///< This function will update the physical property spin box to select which strain is displayed.
    void  	updateDrawSymmetricityViewToggle();                                                 ///< This function will toggle the boolean to draw the symmetricity or not.
    void    updateDisplayPipette(int);


private:
    void generateControlPanel();                                                                ///< Generate the control panel with view options.
    void setUpView();                                                                           ///< Set up the view.
    void setUpGLWidget();                                                                       ///< Set up the openGL window.
    void setUpCentralWidget();                                                                  ///< Set up the central widget defining the whole window structure.
    void setUpSelectionDisplayGrid(QGridLayout *SelectionDisplayGrid);                          ///< Set up the control panel grid for selected item options.
    void setUpProjectDisplayOptionGrid(QGridLayout *ProjectDisplayOptionsGrid);                 ///< Set up the control panel grid for display options.
    void setUpViewOptionsGrid(QGridLayout *ViewOptionsGrid);                                    ///< The panel for control buttons (Quit)
    void setSelectionByIdSelection(QFont font, QGridLayout *SelectionDisplayGrid);              ///< The function sets the element selection from manual input.
    void setCoordBoxes(QFont font, QFont boldFont, QGridLayout *SelectionDisplayGrid);          ///< The function sets up the coordinate display grid.
    void setItemSelectionTitles(QFont font, QFont boldFont, QGridLayout *SelectionDisplayGrid); ///< The function sets up the selected item name & Id displays.
    void setStrainDisplayMenu(QGridLayout *DisplayOptionsGrid);                                 ///< The function sets up strain displays.
    void setPysPropDisplayMenu(QGridLayout *DisplayOptionsGrid);                                ///< The function sets up physical property displays.
    void setDisplayPreferences(QGridLayout *DisplayOptionsGrid);
    bool takeScreenshot();                                                                      ///< This function saves the screenshots.
    QWidget         *CentralWidget;                                                             ///< The main widget controlling the window setup.
    QLineEdit       *NameBox;                                                                   ///< The text containing the name of the current selected element.
    QLineEdit       *ElementSelectBox;                                                          ///< The text box containing the manually selected ShapeBase#Id
    QTimer          *timer;                                                                     ///< The timer controlles drawing frequency, in cases where we would like to have slow display.
    int             nCoordBox;                                                                  ///< The number of coordinate boxes on display
    QLineEdit 		*CoordBox_id[6];                                                            ///< The text array containing the Node#Id values of the selected element.
    QLineEdit 		*CoordBox_x[6];                                                             ///< The text array containing the x-coordinates of the nodes for the selected element.
    QLineEdit 		*CoordBox_y[6];                                                             ///< The text array containing the y-coordinates of the nodes for the selected element.
    QLineEdit 		*CoordBox_z[6];                                                             ///< The text array containing the z-coordinates of the nodes for the selected element.
    QLabel  		*CoordLabel_n[6];                                                           ///< The labels for the nodes
    QCheckBox		*DisplayCheckBoxes[2];                                                      ///< The array containing the check boxes of display options, [Strain][Physical Properties]
    QComboBox   	*StrainComboBox;                                                            ///< The combo box storing strain display options.
    QDoubleSpinBox	*StrainSpinBoxes[2];                                                        ///< The spin boxes listing the minimum and maximum of displayed strain.
    QComboBox   	*PysPropComboBox;                                                           ///< The combo box storing physical property display options.
    QDoubleSpinBox 	*PysPropSpinBoxes[2];                                                       ///< The spin boxes listing the minimum and maximum of displayed physical property.
    QGroupBox   	*ColourCodingBox;
    QCheckBox		*DisplayPreferencesCheckBoxes[1];
    QLabel          *SimTime;                                                                   ///< Simulation time in seconds.
    QPushButton		*SymmetricityDisplayButton;                                                 ///< The button to toggle symmetric display
    std::clock_t    simulationStartClock;	                                                    ///< simulation clock as in cpu clock time (will be the same regardless of using single or multi-processors. sleep during process etc.
    std::time_t 	simulationStartTime;	                                                    ///< simulation time in real time, given in seconds
    bool            displayedSimulationLength = false;                                          ///< The boolean stating if the display tool showed all the length of the simulation.
    bool takeScreenShotAfterUpdate=false;
};


#endif /* MAINWINDOW_H_ */
