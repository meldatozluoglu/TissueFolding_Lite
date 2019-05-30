/*
 * GLWidget.h
 *
 *  Created on: 20 Mar 2014
 *      Author: melda
 */

#ifndef GLWIDGET_H_
#define GLWIDGET_H_

#include <QGLWidget>
#include <QtOpenGL>
#include <vector>
#include <array>
#include <string.h>
#include <math.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtx/string_cast.hpp>

#include "../TissueFolding/SourceCode/Simulation.h"

/*! GLWidget class */
 class GLWidget : public QGLWidget
 {
    Q_OBJECT
    GLuint vao;                         ///< Vertex array object
    GLuint vbo;                         ///< Vertex buffer
    GLuint ebo;                         ///< GLElements buffer
    GLuint vertexShader;                ///< Vertex Shader for tissue
    GLuint fragmentShader;              ///< The fragment Shader for tissue
    GLuint shaderProgram;               ///< The shader program to manage all the shaders.
    GLint posAttrib;                    ///< Position attribute id for shaders.
    GLint colAttrib;                    ///< Colour attribute id for shaders.
    GLint normalAttrib;                 ///< Normal attribute id for shaders.
    GLint pickingColourAttrib;          ///< Element picking colour attribute id for shaders.
    vector<float> GLVertices;           ///< Vector of vertices
    vector<int> GLElements;             ///< Vector of elements
    GLint uniProj;                      ///< Projection uniform
    GLint uniTrans;                     ///< Projection uniform
    GLint uniTransToOrigin;             ///< Translation to coordinate system origin uniform
    GLint uniTransBackToPosition;       ///< Translation back to original position uniform
    GLint uniRot;                       ///< Rotation uniform
    GLint uniPicking;                   ///< element picking uniform
    GLint uniXSymmetric;                ///< symmetricity display - X uniform
    GLint uniYSymmetric;                ///< symmetricity display - Y uniform
    GLint uniXYSymmetric;               ///< symmetricity display - XY uniform
    glm::mat4 proj;                     ///< Projection matrix
    glm::mat4 modelRotate;              ///< Model rotation matrix
    glm::mat4 modelTranslate;           ///< Model translation matrix
    glm::mat4 modelTranslateToOrigin;   ///< Model translate to cooridinate system origin before rotation matrix
    glm::mat4 modelTranslateBackToPosition;    ///< Model translate back to original position after rotation matrix
    GLint isDrawingForPicking;          ///< toggle for drawing picking
    GLint isDrawingXSymmetric;          ///< toggle for symmetricity display - X
    GLint isDrawingYSymmetric;          ///< toggle for symmetricity display - Y
    GLint isDrawingXYSymmetric;         ///< toggle for symmetricity display - XY
    float orthoViewLimits[6];           ///< The ortagonal view limits

 public:
    GLWidget(QWidget *parent = 0); ///< Constructor
    ~GLWidget();                   ///< Destructor

     std::shared_ptr<Simulation> Sim01; ///< Shared pointer to the simulation object.
     bool       updateVBO;              ///< The boolean flagging the GLWidget#vbo should be updated.
     QSize		minimumSizeHint() const;
     QSize 		sizeHint() const;
     void 		manualElementSelection(int i);  ///< This funciton will make screen and vbo updates for selecting an element from id input menu bar manually.
     bool 		ItemSelected;                   ///< Boolean stating if an element is selected.
     string 	SelectedItemName;               ///< The name of the selected element.
     int 		SelectedItemIndex;              ///< The index of the selected element.
     bool 		DisplayStrains;                 ///< Boolean stating if the strains are visualised.
     float 		DisplayStrainRange[2];          ///< The double array giving the scale of minimum and maximum strain for colour mapping.
     int 		StrainToDisplay;                ///< The option for the selected strain.
     bool 		DisplayPysProp;                 ///< Boolean stating if the physical propoerties are visualised.
     int 		PysPropToDisplay;               ///< The option for the selected physical property.
     float 		DisplayPysPropRange[5][2];      ///< The current ranges for physical property, for each possible physical property: [min][max].
     float 		DisplayPysPropBounds[5][4];     ///< The minimum and maximum of the range for selected physical property, for each possible physical property: [min-lower bound][min-upper bound][max-lower bound][max-upper bound].
     int		DisplayPysPropDecimals[5];      ///< The decimal points for spin boxes for physical property display.
     float 		DisplayPysPropSteps[5];         ///< current step of the spinbox for physical property display.
     vector <QString> SelectedPos;              ///< The vector of selected positions in string form for display
     vector <QString> SelectedId;               ///< The vector of selected ids in string form for display
     bool       displayPipette;
     bool 		drawSymmetricity;
     float 		obj_pos[3];                     ///< The array in 3D, stating the traslation of the current view.

 signals:
     void SelectedItemChanged();                ///< Signal emmitted when the selected item has changed, either with manu input, or as mouse click.
     void NeedToClearManualElementSelection();  ///< Signal emmitted to clear existing element selection

 protected:
     void   initializeGL();                     ///< This function sets up the GL object
     void 	paintGL();                          ///< This function paints the screen at each update
     void 	resizeGL(int width, int height);    ///< This function ensures aspect ratio preserving update of the view.
     void 	mousePressEvent (QMouseEvent *event);   ///< Mouse press event for view manipulation.
     void 	mouseReleaseEvent (QMouseEvent *event); ///< Mouse release event for view manipulation.
     void 	mouseMoveEvent (QMouseEvent *event);    ///< Mouse move event for view manipulation.
     void 	wheelEvent(QWheelEvent *event);         ///< Mouse wheel event for view manipulation.
     void 	ObjectSelection(QPoint LastPos);        ///< This function checks for selection of an element under the last mouse click, updates model as needed.
     void   updateVBOData();                        ///< This funciton updates all the GLWidget#vbo information from the simulation.
     void   updateVBOColors();                      ///< This funciton updates the GLWidget#vbo node colour information from the simulation.
     void 	resetItemSelectionInfo(int source);     ///< This function reseats the element selection information, as a new selection identified no elements.
     void 	findElement();                          ///< This function finds the element under mouse click from picked unique identifier colour.
     bool 	findElement(int i);                     ///< This function finds the element with input index.
     void 	getColourOfPoint(QPoint LastPos);       ///< This funciton returns the screen colour of the input point.
     void 	drawForPicking ();                      ///< This function draws the elements on screen buffer with their unique picking identifier colours.
     void 	initialiseNodeColourList();             ///< This funciton initialises the nodal colours from simulation information and display preferences.
     void 	fillItemSelectionInfo(int i);           ///< This function fills in the item selection information for the element with input indice i.

 private:
     QPoint lastPos;                                ///< Last clicked position for short clinks, relaese position for dragged clicks.
     QPoint InitialClickPos;                        ///< Initial clicked position for dragged clicks.
     int 	MouseButton;                            ///< The mouse button identifier.
     
     float 	aspectratio;                            ///< The current aspect ratio of the openGL window.
     int 	PickedColour[4];                        ///< The array giving the picked colour from screen buffer [r][g][b][alpha]
     vector<std::array<float,3>> NodeColourList;    ///< The vector of 3D colour descriptions for all nodes.
     double 	Qcurr[4];                           ///< The current quaternian defining the rotation of the world coordinate system.
     double     Qlast[4];                           ///< The last quaternian defining the rotation of the world coordinate system.
     float 		MatRot[16];                         ///< 4x4 rotation curernt matrix of the world written in vector form.
     bool 		checkPickedColour(std::array<int,3> ElementColour);                 ///< This functon checks if the unique identifier colour of the element (as returned by ShapeBase#getIdentifierColour) is the same as the current picked colour (GLWidget#PickedColour).
     void 		rotateByQuaternians(std::array<double,4>& Qrot);                    ///< This function rotates the system rotation definition quaternians by the input quaternian, and updates the vbo accordingly.
     void 		normaliseCurrentRotationAngle (std::array<double,4>& Qrot);         ///< This function normalises the input quaternian.
     void 		rotateCurrentRotationQuaternian(std::array<double,4>& Qrot);        ///< This function rotates the quaternian defining the current system rotation by the input quaternian.
     void 		rotateMatrix();                                                     ///< This function rotates the world rotation matrix by the input quaternian.
     std::array<float,3> getDisplayColour(float Data);                              ///< This function calculated the display colour of a strain or physical property, from simulation information and the display range specified at used interface.
     void 		constructNodeColourList();                                          ///< This function constructs the node colour list of all nodes in the simulation, from user interface preferences and simulation data.
     std::array<std::array<float,3>,6> 	getElementColourList(int i);                ///< This function constructs the element colour list of all nodes in the simulation, from user interface preferences and simulation data.
 };


#endif /* GLWIDGET_H_ */
