#include <iostream>
#include <math.h>
#include <vector>
#include <string.h>

#include <QtGui>
#include <QtOpenGL>
#include <QWheelEvent>

#include "GLWidget.h"
#include "../TissueFolding/SourceCode/ShapeBase.h"


const char* vertexSource = R"glsl(
    #version 150 core
    in vec3 position;
    in vec3 color;
    in vec3 pickingColour;
    in vec3 normal;
    out vec3 Color;

    uniform mat4 modelRotate;
    uniform mat4 modelTranslate;
    uniform mat4 proj;
    uniform mat4 modelTranslateToOrigin;
    uniform mat4 modelTranslateBackToPosition;
    uniform int isDrawingForPicking;
    uniform int isDrawingXSymmetric;
    uniform int isDrawingYSymmetric;
    uniform int isDrawingXYSymmetric;

    void main()
    {
        if (isDrawingForPicking == 1){ //if I am drawing the borders, the colour is black, else, use vertex colour
            Color = pickingColour;
        }
        else{
            float ambientFrac = 0.6;
            vec4 rotatedNormal =  modelRotate* vec4(normal, 1.0);
            float cosTheta = clamp( rotatedNormal[2], 0,1 );
            Color = color*vec3(ambientFrac)+vec3(1.0f-ambientFrac)*color*cosTheta;
        }
        if (isDrawingXSymmetric == 1){
            gl_Position =  proj * modelTranslate * modelTranslateBackToPosition* modelRotate * modelTranslateToOrigin* vec4(-1.0*position[0],position[1],position[2], 1.0);
        }
        else if (isDrawingYSymmetric == 1){
            gl_Position =  proj * modelTranslate * modelTranslateBackToPosition* modelRotate * modelTranslateToOrigin* vec4(position[0], -1.0f*position[1], position[2], 1.0);
        }
        else if(isDrawingXYSymmetric == 1){
            gl_Position =  proj * modelTranslate * modelTranslateBackToPosition* modelRotate * modelTranslateToOrigin* vec4(-1.0f*position[0],-1.0f*position[1],position[2], 1.0);
        }
        else{
            gl_Position =  proj * modelTranslate * modelTranslateBackToPosition* modelRotate * modelTranslateToOrigin* vec4(position, 1.0);
        }
    }
)glsl";

const char* fragmentSource = R"glsl(
    #version 150 core

    in vec3 Color;
    out vec4 outColor;

    uniform int isDrawingForPicking;
    void main()
    {
        outColor = vec4(Color, 1.0);
    }
)glsl";


using namespace std;

 GLWidget::GLWidget(QWidget *parent) : QGLWidget(parent)
 {
     updateVBO = true;
     obj_pos[0] =  0.0f;
     obj_pos[1] =  0.0f;
     obj_pos[2] =  200.f;
	 MatRot[0]  = 1.0; MatRot[1]  = 0.0; MatRot[2]  = 0.0; MatRot[3]  = 0.0;
	 MatRot[4]  = 0.0; MatRot[5]  = 1.0; MatRot[6]  = 0.0; MatRot[7]  = 0.0;
	 MatRot[8]  = 0.0; MatRot[9]  = 0.0; MatRot[10] = 1.0; MatRot[11] = 0.0;
	 MatRot[12] = 0.0; MatRot[13] = 0.0; MatRot[14] = 0.0; MatRot[15] = 1.0;

     Qcurr[0] =  0.82262; Qcurr[1] =  -0.529113; Qcurr[2] =  0.00377568; Qcurr[3] =   -0.208138;
     Qlast[0] =  0.82262; Qlast[1] =  -0.529113; Qlast[2] =  0.00377568; Qlast[3] =   -0.208138;
	 PickedColour[0] = 255;
	 PickedColour[1] = 255;
	 PickedColour[2] = 255;
	 PickedColour[3] = 1;
	 ItemSelected = false;
	 SelectedItemIndex = -1;
	 SelectedItemName = "N/A";
     aspectratio =1.0;
     DisplayStrains = false;
     DisplayPysProp = false;
     drawSymmetricity = true;
     //current ranges:
     DisplayPysPropRange[0][0] = 1000.0; DisplayPysPropRange[0][1] = 20000.0; 	//External Viscosity
     DisplayPysPropRange[1][0] = 1000.0; DisplayPysPropRange[1][1] = 50000.0; 	//Internal Viscosity
     DisplayPysPropRange[2][0] = 10.0; DisplayPysPropRange[2][1] = 1600.0;     //Young's modulus
     DisplayPysPropRange[3][0] = 0.05; DisplayPysPropRange[3][1] = 0.49;        //Poisson's ratio
     DisplayPysPropRange[4][0] = 0.0; DisplayPysPropRange[4][1] = 12;           //volumetric (xyz) growth rate

     //the minimum and maximum they can get:
     DisplayPysPropBounds[0][0] = 0.0;   DisplayPysPropBounds[0][1] = 10000.0;      //External Viscosity - lower limit, min max range
     DisplayPysPropBounds[0][2] = 100.0; DisplayPysPropBounds[0][3] = 100000.0;     //External Viscosity - upper limit, min max range
     DisplayPysPropBounds[1][0] = 0.0;   DisplayPysPropBounds[1][1] = 5000.0;       //Internal Viscosity - lower limit, min max range
     DisplayPysPropBounds[1][2] = 100.0; DisplayPysPropBounds[1][3] = 100000.0;     //Internal Viscosity - upper limit, min max range
     DisplayPysPropBounds[2][0] = 1.0;   DisplayPysPropBounds[2][1] = 1000.0;       //Young's modulus - lower limit, min max range
     DisplayPysPropBounds[2][2] = 51.0;  DisplayPysPropBounds[2][3] = 100000.0;  	//Young's modulus - upper limit, min max range
     DisplayPysPropBounds[3][0] = 0.0;   DisplayPysPropBounds[3][1] = 0.1;          //Poisson's ratio - lower limit, min max range
     DisplayPysPropBounds[3][2] = 0.11;  DisplayPysPropBounds[3][3] = 0.5;          //Poisson's ratio - upper limit, min max range
     DisplayPysPropBounds[4][0] = 0.0;   DisplayPysPropBounds[4][1] = 10;           //xy-planar growth rate - lower limit, min max range
     DisplayPysPropBounds[4][2] = 1;     DisplayPysPropBounds[4][3] = 150; 			//xy-planar growth rate - upper limit, min max range
     //the decimals to display:
     DisplayPysPropDecimals[0] = 0;
     DisplayPysPropDecimals[1] = 0;
     DisplayPysPropDecimals[2] = 0;
	 DisplayPysPropDecimals[3] = 2;
	 DisplayPysPropDecimals[4] = 0;
     DisplayPysPropSteps[0] = 1;
	 DisplayPysPropSteps[1] = 1;
	 DisplayPysPropSteps[2] = 10;
	 DisplayPysPropSteps[3] = 0.05;
	 DisplayPysPropSteps[4] = 1;
     setSizePolicy(QSizePolicy ::Expanding , QSizePolicy ::Expanding );
 }

 GLWidget::~GLWidget()
 {   
     glDeleteProgram(shaderProgram);
     glDeleteShader(fragmentShader);
     glDeleteShader(vertexShader);

     glDeleteBuffers(1, &ebo);
     glDeleteBuffers(1, &vbo);

     glDeleteVertexArrays(1, &vao);
 }

 QSize GLWidget::minimumSizeHint() const
 {
     return QSize(50, 50);
 }

 QSize GLWidget::sizeHint() const
 {
     return QSize(6000, 6000);
 }

 void GLWidget::initializeGL()
 {
     qglClearColor(QColor::fromRgbF(1, 1, 1, 1));
     glEnable(GL_DEPTH_TEST);
     initialiseNodeColourList();
     DisplayStrains = false;

     //objects and shaders:
     for (const auto& simElement :Sim01->Elements){
         /** In initialisation the vertex buffers will be constructed.
          *Each prism of the system is actually made up of 6 nodes, the ids
          *rangeing from 0-5. For correct lighting, same nodes should hold different normals at the edges.
          *Therefore, for each prism of 6 nodes, a total of 18 vertices are stored: 3 upper surface, 3 lower surface,
          *4 for each of 3 side rectangles. The indice order is stored in the array vertexOrderToAddForPrism.
          */
         int addedVertexNumberPerPrism = 18;
         size_t vertexOrderToAddForPrism[18] ={0, 1, 2,
                                               3, 4, 5,
                                               0, 1, 3, 4,
                                               1, 2, 4, 5,
                                               0, 2, 3, 5
                                              };

         for(size_t i : vertexOrderToAddForPrism){
             /** For each of the stored vertices, first the position array is inserted into the GLWidget#GLVertices.
              */
             float pos[3] = {(float) simElement->Positions[i][0],
                            (float) simElement->Positions[i][1],
                            (float) simElement->Positions[i][2]};
             GLVertices.insert(GLVertices.end(), std::begin(pos), std::end(pos));
             /** then the current colour of the vertex,
              */
             //r,g,b:
             GLVertices.push_back(1.0);
             GLVertices.push_back(0.0);
             GLVertices.push_back(0.0);
             /** then the current normal pointing out of the vertex,
              */
             //normal:
             GLVertices.push_back(0.0);
             GLVertices.push_back(0.0);
             GLVertices.push_back(1.0);
             /** and lastly, the unique identifier colour of ElementColour are appended to the vertex buffer.
              */
             //picking colour:
             std::array<int,3> ElementColour = simElement->getIdentifierColour();
             GLVertices.push_back(((float) ElementColour[0])/255.0);
             GLVertices.push_back(((float) ElementColour[1])/255.0);
             GLVertices.push_back(((float) ElementColour[2])/255.0);
         }
         /** Each prism will be defined by size trienagles (upper/lower surfaces and two triangles each for each of the three side rectangles).
          * These triangles will be defined by the 18 vertices added to the vectex buffer for this triangle. Their
          * order in terms of the elemental vertex index is stored in trianglesForNodeOrderAddedByElement.
          */

         int trianglesForNodeOrderAddedByElement[24] = {
                      0, 1, 2,
                      3, 4, 5,
                      7, 6, 8,
                      7, 9, 8,
                      11, 10, 12,
                      11, 12, 13,
                      14, 15, 17,
                      14, 16, 17
                     };
         /** As these indices are in the elemental definition, they are be positioned in the vertex
          *array with the offset: (element index from the simulation) *18. With the offset correction,
          *the tirangle information is appended to the GLWidget#GLElements
          */
         int offsetForNodeIndex = simElement->Id*addedVertexNumberPerPrism;
         for (auto &index : trianglesForNodeOrderAddedByElement){
             index += offsetForNodeIndex;
         }
         GLElements.insert(GLElements.end(), std::begin(trianglesForNodeOrderAddedByElement), std::end(trianglesForNodeOrderAddedByElement));
     }
     /** Then the vertex (GLWidget#vbo) and element (GLWidget#veo) buffers are generated with the size of
      * the vectors GLWidget#GLVertices and GLWidget#GLElements.
      */
     glGenVertexArrays(1, &vao);
     glBindVertexArray(vao);

     glGenBuffers(1, &vbo);
     glBindBuffer(GL_ARRAY_BUFFER, vbo);

     glGenBuffers(1, &ebo);
     glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);

     glBufferData(GL_ELEMENT_ARRAY_BUFFER,GLElements.size() * sizeof(int), &GLElements.front(), GL_DYNAMIC_DRAW);
     glBufferData(GL_ARRAY_BUFFER,GLVertices.size() * sizeof(float), &GLVertices.front(), GL_DYNAMIC_DRAW);
     /** The shaders are set up.
      */
     vertexShader = glCreateShader(GL_VERTEX_SHADER);
     glShaderSource(vertexShader, 1, &vertexSource, NULL);
     glCompileShader(vertexShader);

     fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
     glShaderSource(fragmentShader, 1, &fragmentSource, NULL);
     glCompileShader(fragmentShader);

     /** In development, the shader status can be checked in a hard coded loop within this function.
      *The loop is not active in release, but can be activated with setting checkShaderStatus = true insisde the code.
      * This loop will also print the log buffer.
      */
     bool checkShaderStatus = false; //hard coded debugging option to check if the shaders compiled
     if (checkShaderStatus){
         //check if shaders compiled:
         GLint vertexStatus;
         glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &vertexStatus);
         cerr<<"vertex shader status: "<<vertexStatus<<std::endl;

         char vertesBuffer[512];
         glGetShaderInfoLog(vertexShader, 512, NULL, vertesBuffer);
         cerr<<"vertex shader buffer: "<<vertesBuffer<<std::endl;

         GLint fragmentStatus;
         glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &fragmentStatus);
         cerr<<"fragment shader status: "<<fragmentStatus<<std::endl;

         char fragmentBuffer[512];
         glGetShaderInfoLog(fragmentShader, 512, NULL, fragmentBuffer);
         cerr<<"vertex shader buffer: "<<fragmentBuffer<<std::endl;
     }
     shaderProgram = glCreateProgram();
     glAttachShader(shaderProgram, vertexShader);
     glAttachShader(shaderProgram, fragmentShader);

     glBindFragDataLocation(shaderProgram, 0, "outColor");
     glLinkProgram(shaderProgram);
     glUseProgram(shaderProgram);
     /** The position attribute is assiged for the first 3 items on the stride of 12 floats in the vertex buffer.
      */
     posAttrib = glGetAttribLocation(shaderProgram, "position");
     glVertexAttribPointer(posAttrib, 3, GL_FLOAT, GL_FALSE, 12*sizeof(float), 0);
     glEnableVertexAttribArray(posAttrib);
     /** The colour attribute is assiged after an offset of 3 floats (for position) to the next
      *  3 items on the stride of 12 floats in the vertex buffer.
      */
     colAttrib = glGetAttribLocation(shaderProgram, "color");
     glEnableVertexAttribArray(colAttrib);
     glVertexAttribPointer(colAttrib, 3, GL_FLOAT, GL_FALSE, 12*sizeof(float), (void*)(3*sizeof(float)));

     /** The normal attribute is assiged after an offset of 6 floats (for position + colour) to the next
      *  3 items on the stride of 12 floats in the vertex buffer.
      */
     normalAttrib = glGetAttribLocation(shaderProgram, "normal");
     glEnableVertexAttribArray(normalAttrib);
     glVertexAttribPointer(normalAttrib, 3, GL_FLOAT, GL_FALSE, 12*sizeof(float), (void*)(6*sizeof(float)));
     /** Fianllt, the unique colour identifier of the elementm picking colour attribute is assiged after
      * an offset of 9 floats (for position + colour + normal) to the next
      *  3 items on the stride of 12 floats in the vertex buffer.
      */
     pickingColourAttrib = glGetAttribLocation(shaderProgram, "pickingColour");
     glEnableVertexAttribArray(pickingColourAttrib);
     glVertexAttribPointer(pickingColourAttrib, 3, GL_FLOAT, GL_FALSE, 12*sizeof(float), (void*)(9*sizeof(float)));

     /** Then the orthagonal view and rotation matrices are set followed by the drawing options.
      */

     orthoViewLimits[2] = -53.5f;
     orthoViewLimits[3] =  53.5f;
     orthoViewLimits[0] = orthoViewLimits[2]*aspectratio;
     orthoViewLimits[1] = orthoViewLimits[3]*aspectratio;
     orthoViewLimits[4] = -1000.f;
     orthoViewLimits[5] = 1000.f;

     proj = glm::ortho(orthoViewLimits[0],orthoViewLimits[1],orthoViewLimits[2],orthoViewLimits[3],  orthoViewLimits[4],  orthoViewLimits[5]);
     uniProj = glGetUniformLocation(shaderProgram, "proj");
     glUniformMatrix4fv(uniProj, 1, GL_FALSE, glm::value_ptr(proj));

     glm::quat newQuaternion = glm::quat(Qcurr[0],Qcurr[1],Qcurr[2],Qcurr[3]);
     modelRotate = glm::mat4_cast(newQuaternion);
     uniRot = glGetUniformLocation(shaderProgram, "modelRotate");
     glUniformMatrix4fv(uniRot, 1, GL_FALSE, glm::value_ptr(modelRotate));

     modelTranslate = glm::mat4(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1);
     uniTrans = glGetUniformLocation(shaderProgram, "modelTranslate");
     glUniformMatrix4fv(uniTrans, 1, GL_FALSE, glm::value_ptr(modelTranslate));

     modelTranslateToOrigin = glm::mat4(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1);
     uniTransToOrigin = glGetUniformLocation(shaderProgram, "modelTranslateToOrigin");
     glUniformMatrix4fv(uniTransToOrigin, 1, GL_FALSE, glm::value_ptr(modelTranslateToOrigin));

     modelTranslateBackToPosition = glm::mat4(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1);
     uniTransBackToPosition = glGetUniformLocation(shaderProgram, "modelTranslateBackToPosition");
     glUniformMatrix4fv(uniTransToOrigin, 1, GL_FALSE, glm::value_ptr(modelTranslateBackToPosition));

     isDrawingForPicking = 0; //start with faces, then raw borders
     uniPicking = glGetUniformLocation(shaderProgram, "isDrawingForPicking");
     glUniform1i(uniPicking, isDrawingForPicking);

     isDrawingXSymmetric = 0; //start without symmetricity
     uniXSymmetric = glGetUniformLocation(shaderProgram, "isDrawingXSymmetric");
     glUniform1i(uniXSymmetric, isDrawingXSymmetric);

     isDrawingYSymmetric = 0; //start without symmetricity
     uniYSymmetric = glGetUniformLocation(shaderProgram, "isDrawingYSymmetric");
     glUniform1i(uniYSymmetric, isDrawingYSymmetric);

     isDrawingXYSymmetric = 0; //start without symmetricity
     uniXYSymmetric = glGetUniformLocation(shaderProgram, "isDrawingXYSymmetric");
     glUniform1i(uniXYSymmetric, isDrawingXYSymmetric);
 }

 void GLWidget::initialiseNodeColourList(){
     size_t n = Sim01->Nodes.size();
     std::array<float,3> dummyColourBlack = {0.0,0.0,0.0};
     NodeColourList.insert(NodeColourList.end(), n, dummyColourBlack);
 }

 void GLWidget::updateVBOColors(){
     constructNodeColourList();
     std::array<std::array<float,3>,6> NodeColours; //r, g, b colours for 6 nodes of element
     size_t addedVertexNumberPerPrism = 18;
     int numberOfAttributesOnVertex = 12; //x,y,z, r,g,b, normal_x,normal_y,normal_z, pickColour_r, pickColour_g, pickColour_b
     for (const auto& simElement : Sim01->Elements){
         NodeColours = getElementColourList(simElement->Id);
         size_t vertexOrderToAddForPrism[18] ={0, 1, 2,
                                               3, 4, 5,
                                               0, 1, 3, 4,
                                               1, 2, 4, 5,
                                               0, 2, 3, 5
                                              };
         int offsetForNodeIndex = simElement->Id*addedVertexNumberPerPrism*numberOfAttributesOnVertex;
         for(size_t i = 0; i<addedVertexNumberPerPrism; ++i){
            size_t curentNodeIndexOnPrism = vertexOrderToAddForPrism[i];
            int indexOnGLVertices = offsetForNodeIndex +i*numberOfAttributesOnVertex;
            if (ItemSelected && SelectedItemIndex == simElement->Id){
                GLVertices[indexOnGLVertices + 3] = 0.1; //r //picked alamet coulur is dark grey
                GLVertices[indexOnGLVertices + 4] = 0.1; //g
                GLVertices[indexOnGLVertices + 5] = 0.1; //b
            }
            else{
                GLVertices[indexOnGLVertices + 3] = NodeColours[curentNodeIndexOnPrism][0]; //r
                GLVertices[indexOnGLVertices + 4] = NodeColours[curentNodeIndexOnPrism][1]; //g
                GLVertices[indexOnGLVertices + 5] = NodeColours[curentNodeIndexOnPrism][2]; //b
            }
         }
    }
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER,GLVertices.size() * sizeof(float), &GLVertices.front(), GL_DYNAMIC_DRAW);
 }

 void GLWidget::updateVBOData(){
     updateVBO = false;
     constructNodeColourList();
     //update colours & positions:
     int numberOfAttributesOnVertex = 12; //x,y,z, r,g,b, normal_x,normal_y,normal_z, pickColour_r, pickColour_g, pickColour_b
     for (const auto& simElement : Sim01->Elements){
         //calculate nodal normals:
         glm::vec3 vecnodes[6][6];
         glm::vec3 triNormals[8];
         glm::vec3 nodeNormals[18];
         size_t nodeToTriMap[18] = {0, 0, 0,
                                    1, 1, 1,
                                    2, 2, 2, 3,
                                    4, 4, 4, 5,
                                    6, 6, 6, 7};
         size_t trianglesForNodesOfElement[24] = {0, 2, 1,
                       3, 4, 5,
                       1, 3, 0,
                       1, 4, 3,
                       2, 4, 1,
                       2, 5, 4,
                       0, 5, 2,
                       0, 3, 5
                      };
         size_t vectorIndices[11][2]={ {0,1},
                                       {0,2},
                                       {0,3},
                                       {0,5},
                                       {1,3},
                                       {1,4},
                                       {2,1},
                                       {2,4},
                                       {2,5},
                                       {3,4},
                                       {3,5}};
         //each of the 6 nodes have 4 tirangles attributed to it:
         for (const auto& indexPair : vectorIndices){
             size_t index0 =indexPair[0];
             size_t index1 =indexPair[1];
             float dx = (float) (simElement->Positions[index1][0] - simElement->Positions[index0][0]);
             float dy = (float) (simElement->Positions[index1][1] - simElement->Positions[index0][1]);
             float dz = (float) (simElement->Positions[index1][2] - simElement->Positions[index0][2]);
             vecnodes[index0][index1] = glm::vec3(dx, dy, dz);
             vecnodes[index1][index0] = -1.0f*vecnodes[index0][index1];
         }

         for (size_t i=0; i<8; ++i){
             size_t index0 = trianglesForNodesOfElement[3*i];
             size_t index1 = trianglesForNodesOfElement[3*i+1];
             size_t index2 = trianglesForNodesOfElement[3*i+2];
             triNormals[i] = glm::normalize(glm::cross(vecnodes[index0][index1], vecnodes[index0][index2]));
         }
         for (size_t i=0; i<18; ++i){
             nodeNormals[i] = triNormals[nodeToTriMap[i]];
         }
         //continue;
         std::array<std::array<float,3>,6> NodeColours;
         NodeColours = getElementColourList(simElement->Id);
         size_t addedVertexNumberPerPrism = 18;
         size_t vertexOrderToAddForPrism[18] ={0, 1, 2,
                                               3, 4, 5,
                                               0, 1, 3, 4,
                                               1, 2, 4, 5,
                                               0, 2, 3, 5
                                              };
         size_t offsetForNodeIndex = simElement->Id*addedVertexNumberPerPrism*numberOfAttributesOnVertex;
         for(size_t i = 0; i<addedVertexNumberPerPrism; ++i){
            size_t curentNodeIndexOnPrism = vertexOrderToAddForPrism[i];
            size_t indexOnGLVertices = offsetForNodeIndex +i*numberOfAttributesOnVertex;
            GLVertices[indexOnGLVertices + 0] = (float) simElement->Positions[curentNodeIndexOnPrism][0];
            GLVertices[indexOnGLVertices + 1] = (float) simElement->Positions[curentNodeIndexOnPrism][1];
            GLVertices[indexOnGLVertices + 2] = (float) simElement->Positions[curentNodeIndexOnPrism][2];
            if (ItemSelected && SelectedItemIndex == simElement->Id){
                GLVertices[indexOnGLVertices + 3] = 0.1; //r //picked alamet coulur is dark grey
                GLVertices[indexOnGLVertices + 4] = 0.1; //g
                GLVertices[indexOnGLVertices + 5] = 0.1; //b
            }
            else{
                GLVertices[indexOnGLVertices + 3] = NodeColours[curentNodeIndexOnPrism][0]; //r
                GLVertices[indexOnGLVertices + 4] = NodeColours[curentNodeIndexOnPrism][1]; //g
                GLVertices[indexOnGLVertices + 5] = NodeColours[curentNodeIndexOnPrism][2]; //b
            }
            GLVertices[indexOnGLVertices + 6] = nodeNormals[i][0]; //normal_x
            GLVertices[indexOnGLVertices + 7] = nodeNormals[i][1]; //normal-y
            GLVertices[indexOnGLVertices + 8] = nodeNormals[i][2]; //normal_z

         }
     }
     glBindBuffer(GL_ARRAY_BUFFER, vbo);
     glBufferData(GL_ARRAY_BUFFER,GLVertices.size() * sizeof(float), &GLVertices.front(), GL_DYNAMIC_DRAW);
 }

 void GLWidget::paintGL()
 {
     QSize viewport_size = size();
     glViewport(0, 0, viewport_size.width()*this->devicePixelRatio(), viewport_size.height()*this->devicePixelRatio());
     glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
     if (updateVBO){
         //make changes on the vertex buffer and bind it again here:
        updateVBOData();
     }
     modelTranslate =  glm::translate(glm::mat4(1.0f), glm::vec3(obj_pos[0], obj_pos[1], -obj_pos[2]));
     glUniformMatrix4fv(uniTrans, 1, GL_FALSE, glm::value_ptr(modelTranslate));
     modelTranslateToOrigin = glm::translate(glm::mat4(1.0f), glm::vec3(Sim01->SystemCentre[0], Sim01->SystemCentre[1], -Sim01->SystemCentre[2]));
     glUniformMatrix4fv(uniTransToOrigin, 1, GL_FALSE, glm::value_ptr(modelTranslateToOrigin));
     modelTranslateBackToPosition = glm::translate(glm::mat4(1.0f), glm::vec3(-Sim01->SystemCentre[0], -Sim01->SystemCentre[1], Sim01->SystemCentre[2]));
     glUniformMatrix4fv(uniTransBackToPosition, 1, GL_FALSE, glm::value_ptr(modelTranslateBackToPosition));

     //always strat drawing from non-symmetric version:
     isDrawingXSymmetric = 0; //start without symmetricity
     glUniform1i(uniXSymmetric, isDrawingXSymmetric);
     isDrawingYSymmetric = 0; //start without symmetricity
     glUniform1i(uniYSymmetric, isDrawingYSymmetric);
     isDrawingXYSymmetric = 0; //start without symmetricity
     glUniform1i(uniXYSymmetric, isDrawingXYSymmetric);

     glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
     glBufferData(GL_ELEMENT_ARRAY_BUFFER,GLElements.size() * sizeof(int), &GLElements.front(), GL_DYNAMIC_DRAW);
     glDrawElements(GL_TRIANGLES, GLElements.size(), GL_UNSIGNED_INT, 0);

     if(drawSymmetricity){
         if (Sim01->symmetricX && Sim01->symmetricY){
             isDrawingXYSymmetric = 1; //draw x &y symmetric, reflected on both axes
             glUniform1i(uniXYSymmetric, isDrawingXYSymmetric);
             glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
             glBufferData(GL_ELEMENT_ARRAY_BUFFER,GLElements.size() * sizeof(int), &GLElements.front(), GL_DYNAMIC_DRAW);
             glDrawElements(GL_TRIANGLES, GLElements.size(), GL_UNSIGNED_INT, 0);
         }
         if(Sim01->symmetricX){
             isDrawingXSymmetric = 1; //draw x symmetric, reflected on x axis
             glUniform1i(uniXSymmetric, isDrawingXSymmetric);
             glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
             glBufferData(GL_ELEMENT_ARRAY_BUFFER,GLElements.size() * sizeof(int), &GLElements.front(), GL_DYNAMIC_DRAW);
             glDrawElements(GL_TRIANGLES, GLElements.size(), GL_UNSIGNED_INT, 0);
         }
         if(Sim01->symmetricY){
             isDrawingYSymmetric = 1; //draw y symmetric, reflected on y axis
             glUniform1i(uniYSymmetric, isDrawingYSymmetric);
             glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
             glBufferData(GL_ELEMENT_ARRAY_BUFFER,GLElements.size() * sizeof(int), &GLElements.front(), GL_DYNAMIC_DRAW);
             glDrawElements(GL_TRIANGLES, GLElements.size(), GL_UNSIGNED_INT, 0);
         }
    }
 }

 void GLWidget::constructNodeColourList(){
    for (const auto& node : Sim01->Nodes){
        std::array<float,3> currColour;
        if(!DisplayStrains && !DisplayPysProp){
			//I am not displaying ant data on the colours, therefore I do not need any calculaitons on element basis, the colour is constant
            if (node->tissueType == 0){ // columnar layer
                NodeColourList[node->Id][0]=0.75;
                NodeColourList[node->Id][1]=1.0;
                NodeColourList[node->Id][2]=1.0;
			}
            else if (node->tissueType == 1){ // Peripodial membrane
                NodeColourList[node->Id][0]=1.0;
                NodeColourList[node->Id][1]=1.0;
                NodeColourList[node->Id][2]=0.75;
			}
            else if (node->tissueType == 2){ // Linker Zone
                NodeColourList[node->Id][0]=0.87;
                NodeColourList[node->Id][1]=1.0;
                NodeColourList[node->Id][2]=0.87;
			}
		}
		else{
			if(DisplayStrains){
				float StrainMag = 0.0;
                int nConnectedElements = node->connectedElementIds.size();
                for (int i=0;i<nConnectedElements; ++i){
                    float TmpStrainMag =0.0;
                    Sim01->Elements[node->connectedElementIds[i]]->getStrain(StrainToDisplay, TmpStrainMag);
                    StrainMag += (TmpStrainMag)*node->connectedElementWeights[i];
                }
                currColour = getDisplayColour(StrainMag);
			}
			else if (DisplayPysProp){
				float PysPropMag = 0.0;
                //If the physical property is external viscosity, then get the colour directly
				if (PysPropToDisplay == 0){
                    PysPropMag = (float) node->externalViscosity[0];
				}
				else{
                    int nConnectedElements = node->connectedElementIds.size();
					for (int i=0;i<nConnectedElements; ++i){
						float TmpPysPropMag = 0.0;
						if ( PysPropToDisplay == 4 || PysPropToDisplay == 5 || PysPropToDisplay == 6){
							//Growth is multiplicative, base should be 1.0:
							TmpPysPropMag = 1.0;
						}
                        Sim01->Elements[node->connectedElementIds[i]]->getPysProp(PysPropToDisplay, TmpPysPropMag, Sim01->dt);
                        PysPropMag += TmpPysPropMag*node->connectedElementWeights[i];
					}
				}
                currColour = getDisplayColour(PysPropMag);
			}
            NodeColourList[node->Id][0]=currColour[0];
            NodeColourList[node->Id][1]=currColour[1];
            NodeColourList[node->Id][2]=currColour[2];
		}
	}
}

std::array<std::array<float,3>,6> GLWidget::getElementColourList(int i){
    size_t n = Sim01->Elements[i]->getNodeNumber();
    const std::vector<int>& NodeIds = Sim01->Elements[i]->getNodeIds();
    std::array<std::array<float,3>,6> NodeColours;
    for (size_t j = 0; j<n; j++){
		if (DisplayPysProp && PysPropToDisplay == 4){
			//Displaying growth rate, I want this on an elemental basis
            std::array<float,3> currColour;
			//Growth is multiplicative, base should be 1.0:
			float PysPropMag = 1.0;
			Sim01->Elements[i]->getPysProp(PysPropToDisplay, PysPropMag, Sim01->dt);
            currColour = getDisplayColour(PysPropMag);
			NodeColours[j][0]=currColour[0];
			NodeColours[j][1]=currColour[1];
			NodeColours[j][2]=currColour[2];
		}
		else if (Sim01->thereIsExplicitECM && Sim01->Elements[i]->isECMMimicing){
            if(!DisplayStrains && !DisplayPysProp){
                NodeColours[j][0]=0.6;
                NodeColours[j][1]=0.6;
                NodeColours[j][2]=0.0;
			}
			else{
				NodeColours[j][0]=NodeColourList[NodeIds[j]][0];
				NodeColours[j][1]=NodeColourList[NodeIds[j]][1];
				NodeColours[j][2]=NodeColourList[NodeIds[j]][2];
			}
		}
		else if (Sim01->thereIsExplicitActin && Sim01->Elements[i]->isActinMimicing){
            if(!DisplayStrains && !DisplayPysProp){
				NodeColours[j][0]=0.0;
				NodeColours[j][1]=0.6;
				NodeColours[j][2]=0.0;
			}
			else{
				NodeColours[j][0]=NodeColourList[NodeIds[j]][0];
				NodeColours[j][1]=NodeColourList[NodeIds[j]][1];
				NodeColours[j][2]=NodeColourList[NodeIds[j]][2];
			}
		}
		else{
			NodeColours[j][0]=NodeColourList[NodeIds[j]][0];
			NodeColours[j][1]=NodeColourList[NodeIds[j]][1];
			NodeColours[j][2]=NodeColourList[NodeIds[j]][2];
		}
	}
	return NodeColours;
 }

 std::array<float,3> GLWidget::getDisplayColour(float Data){
     std::array<float,3> OutputColour;
     float DataMin=0, DataMax=0;
	 if(DisplayStrains){
		 DataMin = DisplayStrainRange[0];
		 DataMax = DisplayStrainRange[1];
	 }
	 else if (DisplayPysProp){
		 DataMin = DisplayPysPropRange[PysPropToDisplay][0];
		 DataMax = DisplayPysPropRange[PysPropToDisplay][1];
	 }
	 float segment = (DataMax - DataMin)/5.0;
	 float r,g,b;
	 float minR = 0.4;
	 float minB = 0.4;
	 OutputColour[0] = 0.0;
	 OutputColour[1] = 0.0;
	 OutputColour[2] = 0.0;
	 if ((Data - DataMin) < segment/2){
		 float d = (Data - DataMin);
		 r = (1.0 - minR)/ (segment/2.0) * d + minR;
		 g = 0;
		 b = 0;
	 }
	 else if ((Data - DataMin) <segment*1.5){
		 float d = (Data - DataMin) - segment/2.0;
		 r = 1.0;
		 g = d/segment;
		 b = 0;
	 }
	 else if ((Data - DataMin)< segment*2.5){
		 float d = (Data - DataMin) - 1.5*segment;
		 r = 1.0 - d/segment;
		 g = 1.0;
		 b = 0;
	 }
	 else if ((Data - DataMin)< segment*3.5){
		 float d = (Data - DataMin) - 2.5* segment;
		 r = 0;
		 g = 1.0;
		 b = d/segment;
	 }
	 else if ((Data - DataMin)< segment*4.5){
		 float d = (Data - DataMin) - 3.5* segment;
		 r = 0;
		 g = 1.0 - d/segment;
		 b = 1;
	 }
	 else{
		 float d = (Data - DataMin) - 4.5* segment;
		 r = 0;
		 g = 0;
		 b = 1.0 - d*(1-minB)/(segment/2.0);
	 }
	 //invert the display from blue to red for growth (rate and total growth)!:
	 if (DisplayPysProp && (PysPropToDisplay==4 ||  PysPropToDisplay==5 ||  PysPropToDisplay==6) ){
		 double tmpred = r;
		 r=b;
		 b=tmpred;
	 }
	 OutputColour[0] = r;
	 OutputColour[1] = g;
	 OutputColour[2] = b;
     return OutputColour;
 }

 void GLWidget::resizeGL(int width, int height)
 {
     /** Upon resizing the window, the width and height are scaled with the aspect ratio to keep the
      * aspect ratio of the drawn simulation results constant.
      */
     width *= this->devicePixelRatio();
     height *= this->devicePixelRatio();
     int  side = qMin(width, height);
     glViewport((width - side) / 2, (height - side) / 2, width, height);
     aspectratio = float(width) /float (height);
     proj = glm::ortho(orthoViewLimits[0],orthoViewLimits[1],orthoViewLimits[2],orthoViewLimits[3],orthoViewLimits[4],orthoViewLimits[5]);
     glUniformMatrix4fv(uniProj, 1, GL_FALSE, glm::value_ptr(proj));
 }

 void GLWidget::mousePressEvent(QMouseEvent *event)
 {
     /** The mouse press event will update the clicked button identifier and record the initial position of click.
      */
     lastPos = event->pos();
     if(event->button()==Qt::LeftButton){
         MouseButton = 0;
         InitialClickPos = event->pos();
     }
     else if(event->button()==Qt::MidButton){
    	 MouseButton = 1;
     }
     else if(event->button()==Qt::RightButton){
         MouseButton = 2;
     }
 }

 void GLWidget::mouseReleaseEvent (QMouseEvent *event)
  {
      /** If a left click is released in the immediate neighbourhood of the initial click (allow for 5 pixel
       *error), then this is identified as an attempt at object selection, and the relevant function GLWidget#ObjectSelection is called.
        */
      lastPos = event->pos();
      if(event->button()==Qt::LeftButton){
          MouseButton = 0;
          int dx = InitialClickPos.x() - lastPos.x();
          int dy = InitialClickPos.x() - lastPos.x();
          dx *=this->devicePixelRatio();
          dy *=this->devicePixelRatio();
          if (dx<0) {
              dx *= -1;
          }
          if (dy<0) {
              dy *= -1;
          }
          if (dx<5 && dy<5){
        	  ObjectSelection(lastPos);
          }
          else{
              /** For hardware setups that do not have an inherent middle click,
               * there is ian alternative view movement option added. If the left button is dragged, this will also let
               * traslation of the coordinates.
               */
              //replica of MouseButton == 1 event:
              int dx = lastPos.x() - InitialClickPos.x();
              int dy = lastPos.y() - InitialClickPos.y();
              dx *=this->devicePixelRatio();
              dy *=this->devicePixelRatio();
              //cout<<"dx: "<<dx<<" dy: "<<dy<<std::endl;
              float speed = 0.1;
              obj_pos[0] +=  dx*speed;
              obj_pos[1] += -dy*speed;
              lastPos = event->pos();
              updateGL();
          }
      }
  }

 void GLWidget::mouseMoveEvent(QMouseEvent *event)
 {
     /** If a left click is being moved, then the last position is updated until relase.
      */
	 if (MouseButton==0){
		 lastPos = event->pos();
	 }
	 else if(MouseButton==1){
         /** If a middle click is being moved, then the last position is updated and the
          * view is transposed accordingly.
          */
         int dx = event->x() - lastPos.x();
    	 int dy = event->y() - lastPos.y();
         dx *=this->devicePixelRatio();
         dy *=this->devicePixelRatio();
    	 float speed = 0.1;
    	 obj_pos[0] +=  dx*speed;
    	 obj_pos[1] += -dy*speed;
    	 lastPos = event->pos();
         updateGL();
     }
	 else if(MouseButton==2){
         /** If a right click is being moved, then the system is rotated with the projected movement of the mouse.
          *Rotation is done with quaternians.
          */
         float speed = 50;
		 QSize viewport_size = size();
         float width  = viewport_size.width()/this->devicePixelRatio();
         float height = viewport_size.height()/this->devicePixelRatio();
		 double initialPos[3] = {lastPos.x()/(width*2.0) - 1.0, lastPos.y()/(height*2.0) - 1.0,  0.0};
		 double finalPos[3]   = {event->x()/(width*2.0) - 1.0, event->y()/(height*2.0)  - 1.0,  0.0};
         initialPos[1] = (-1.0) * initialPos[1];
         finalPos[1]   = (-1.0) * finalPos[1];
         float r = 20;
		 float r2 = r*r;
		 //Projecting event location:
		 double lengthSq = finalPos[0]*finalPos[0]+finalPos[1]*finalPos[1];
		 if (lengthSq  <= r2){
			 finalPos[2] = pow((r2 - lengthSq),0.5);
		 }
		 else
		 {
			 double length = pow(lengthSq,0.5);
			 finalPos[0] = finalPos[0]/length*r;
			 finalPos[1] = finalPos[1]/length*r;
			 finalPos[2] = finalPos[2]/length*r;
		 }
		 double mag = pow((finalPos[0]*finalPos[0] + finalPos[1]*finalPos[1] + finalPos[2]*finalPos[2]),0.5);
		 if (fabs(mag) > 0.001 && fabs(mag - 1.0f) > 0.001) {
			 finalPos[0] = finalPos[0]/mag;
			 finalPos[1] = finalPos[1]/mag;
			 finalPos[2] = finalPos[2]/mag;
		 }
		 //projecting last pos:
		 lengthSq = initialPos[0]*initialPos[0]+initialPos[1]*initialPos[1];
		 if (lengthSq  <= r2)
			 initialPos[2] = pow((r2 - lengthSq),0.5);
		 else
		 {
			 double  length = pow(lengthSq ,0.5);
			 initialPos[0] = initialPos[0]/length*r;
			 initialPos[1] = initialPos[1]/length*r;
			 initialPos[2] = initialPos[2]/length*r;
		 }
		 mag = pow((initialPos[0]*initialPos[0] + initialPos[1]*initialPos[1] + initialPos[2]*initialPos[2]),0.5);
		 if (fabs(mag) > 0.001 && fabs(mag - 1.0f) > 0.001) {
			 initialPos[0] = initialPos[0]/mag;
			 initialPos[1] = initialPos[1]/mag;
			 initialPos[2] = initialPos[2]/mag;
		 }

		 double cross[3];
		 cross[0] = initialPos[1]*finalPos[2] - initialPos[2]*finalPos[1];
		 cross[1] = initialPos[2]*finalPos[0] - initialPos[0]*finalPos[2];
		 cross[2] = initialPos[0]*finalPos[1] - initialPos[1]*finalPos[0];
		 mag = pow((cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]),0.5);
		 cross[0] /= mag;
		 cross[1] /= mag;
		 cross[2] /= mag;

		 double dot;
		 dot = initialPos[0]*finalPos[0] + initialPos[1]*finalPos[1] + initialPos[2]*finalPos[2];
         float angle = acosf(min(1.0,dot));
		 angle *= speed;

         if (angle>1E-4 || angle < -1E-4){
             std::array<double,4> Qrot = {cosf(angle), sinf(angle)* cross[0], sinf(angle)*cross[1], sinf(angle)*cross[2]};
             rotateByQuaternians(Qrot);
         }
		 lastPos = event->pos();
		 updateGL();
     }
 }

 void GLWidget::rotateByQuaternians (std::array<double,4>& Qrot){
	//ArcBall rotation
     /** Once the qaternian to roate the system is defined in GLWidget#mouseMoveEvent,
      * this function carries out the actual rotation. The current rotation increment quaternian is first normalised
      * via GLWidget#normaliseCurrentRotationAngle, then the existing rotation quaterneian is rotated by the
      * increment in GLWidget#rotateCurrentRotationQuaternian. Once the model quaternian is updated, the
      * GLWidget#modelRotate model rotation matrix is recast from the quaternian, then the model rotation
      * uniform is updated.
      */
	normaliseCurrentRotationAngle(Qrot);
    rotateCurrentRotationQuaternian(Qrot);

    glm::quat newQuaternion = glm::quat(Qcurr[0],Qcurr[1],Qcurr[2],Qcurr[3]);
    modelRotate = glm::mat4_cast(newQuaternion);
    glUniformMatrix4fv(uniRot, 1, GL_FALSE, glm::value_ptr(modelRotate));
 }

 void GLWidget::normaliseCurrentRotationAngle (std::array<double,4>& Qrot){
     /** Normalising the input quaternian, the input is modified as it is provided as an address.
      */
	double mag2 = Qrot[0] * Qrot[0] + Qrot[1] * Qrot[1] + Qrot[2] * Qrot[2] + Qrot[3] * Qrot[3];
	if (fabs(mag2) > 0.00000001 && fabs(mag2 - 1.0f) > 0.00000001) {
		double mag = pow(mag2,0.5);
		Qrot[0] /= mag;
		Qrot[1] /= mag;
		Qrot[2] /= mag;
		Qrot[3] /= mag;
	}
 }

 void GLWidget::rotateCurrentRotationQuaternian(std::array<double,4>& Qrot){
	//Multiplying current rotation quaternian with the rotation I want:
	Qcurr[0] = 	Qlast[0] * Qrot[0] - Qlast[1] * Qrot[1] - Qlast[2] * Qrot[2] - Qlast[3] * Qrot[3];
	Qcurr[1] =  Qlast[0] * Qrot[1] + Qlast[1] * Qrot[0] - Qlast[2] * Qrot[3] + Qlast[3] * Qrot[2];
	Qcurr[2] =  Qlast[0] * Qrot[2] + Qlast[1] * Qrot[3] + Qlast[2] * Qrot[0] - Qlast[3] * Qrot[1];
	Qcurr[3] =  Qlast[0] * Qrot[3] - Qlast[1] * Qrot[2] + Qlast[2] * Qrot[1] + Qlast[3] * Qrot[0];

	Qlast[0] = Qcurr[0];
	Qlast[1] = Qcurr[1];
	Qlast[2] = Qcurr[2];
	Qlast[3] = Qcurr[3];
 }

 void GLWidget::rotateMatrix(){
	double x2 = Qcurr[1] * Qcurr[1];  double y2 = Qcurr[2] * Qcurr[2];  double z2 = Qcurr[3] * Qcurr[3];
	double xy = Qcurr[1] * Qcurr[2];  double xz = Qcurr[1] * Qcurr[3];  double yz = Qcurr[2] * Qcurr[3];
	double wx = Qcurr[0] * Qcurr[1];  double wy = Qcurr[0] * Qcurr[2];  double wz = Qcurr[0] * Qcurr[3];

	MatRot[0]  = 1.0f - 2.0f * (y2 + z2);
	MatRot[1]  = 2.0f * (xy + wz);
	MatRot[2]  = 2.0f * (xz - wy);
	MatRot[3]  = 0.0f;

	MatRot[4]  = 2.0f * (xy - wz);
	MatRot[5]  = 1.0f - 2.0f * (x2 + z2);
	MatRot[6]  = 2.0f * (yz + wx);
	MatRot[7]  = 0.0f;

	MatRot[8]  = 2.0f * (xz + wy);
	MatRot[9]  = 2.0f * (yz - wx);
	MatRot[10] = 1.0f - 2.0f * (x2 + y2);
	MatRot[11] = 0.0f;

	MatRot[12] = 0.0f;
	MatRot[13] = 0.0f;
	MatRot[14] = 0.0f;
	MatRot[15] = 1.0f;
 }

 void GLWidget::wheelEvent(QWheelEvent *event)
  {
     /** If the mouse wheel is activated, the system is zooming in or out.
      */
	 float numDegrees = event->delta() / 8;
	 float numSteps = numDegrees / 30;
     float speed = 10.0;
     obj_pos[2] += -speed*numSteps;
     orthoViewLimits[2] += speed*numSteps;
     orthoViewLimits[3] -= speed*numSteps;
     orthoViewLimits[0] = orthoViewLimits[2]*aspectratio;
     orthoViewLimits[1] = orthoViewLimits[3]*aspectratio;
     proj = glm::ortho(orthoViewLimits[0],orthoViewLimits[1],orthoViewLimits[2],orthoViewLimits[3], orthoViewLimits[4], orthoViewLimits[5]);
     glUniformMatrix4fv(uniProj, 1, GL_FALSE, glm::value_ptr(proj));
	 updateGL();
  }

 void GLWidget::ObjectSelection(QPoint LastPos){
     /** This function starts by drawing the schene into the buffer with the GLWidget#drawForPicking function. Here,
      * all visible elements of the schene are drawn as usual, but instead of colourd defiend by teh user inputs,
      * they are drawn with their unique drawing colours. Please bear in mind this is a transition phase for screen picking,
      * and thse colours are not displatyed at any point.
      */
	 drawForPicking();
     /**
      * The picking attempt is registered bu the function GLWidget#mouseReleaseEvent, and the last pointwhere the mouse
      * was release is recorded in LastPos, an input of this function. Picked colour is read into GLWidget#PickedColour via function
      * GLWidget#getColourOfPoint. Then the owner of the picked unique color is identified in GLWidget#findElement.
      *
      * On the display, there may be item labels left, and those are cleaned by  GLWidget#resetItemSelectionInfo. Then the element is
      * found from the picked colour via GLWidget#findElement.
      */
	 getColourOfPoint(LastPos);
	 resetItemSelectionInfo(1);
	 findElement();
     isDrawingForPicking=0;
     glUniform1i(uniPicking, isDrawingForPicking);
     updateVBOColors();
	 emit SelectedItemChanged();
}


 void GLWidget::manualElementSelection(int i){
     /** This function is selected to update elemetn selection when the id of the element is typed in.
      */
     resetItemSelectionInfo(2);
	 bool validElementId = findElement(i);
	 if (validElementId){
        updateVBOColors();
        emit SelectedItemChanged();
	 }
 }

 void GLWidget::resetItemSelectionInfo(int source){
	 if (source == 1){
		 //The source is item selection via screen click, I will clear up the changes of both manual inputs
		 emit NeedToClearManualElementSelection();
	 }
     //if source == 2, source is a manual element selection
	 ItemSelected = false;
	 SelectedItemName = "";
	 SelectedItemIndex = -1;
	 while ( SelectedPos.size()>0){
		 SelectedPos.pop_back();
 	 }
	 while ( SelectedId.size()>0){
		 SelectedId.pop_back();
 	 }
     updateVBOColors();
 }

 void GLWidget::getColourOfPoint(QPoint LastPos){
	 QSize viewport_size = size();
	 unsigned char pixels[4];
     glReadPixels(LastPos.x()*this->devicePixelRatio(), viewport_size.height()*this->devicePixelRatio() - LastPos.y()*this->devicePixelRatio(), 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, &pixels);
	 PickedColour[0] = (int)pixels[0];
	 PickedColour[1] = (int)pixels[1];
	 PickedColour[2] = (int)pixels[2];
	 PickedColour[3] = (int)pixels[3];
 }

 void GLWidget::findElement(){
	 int n = Sim01->Elements.size();
	 for (int i =0; i<n;i++){
         std::array<int,3> ElementColour = Sim01->Elements[i]->getIdentifierColour();
		 ItemSelected = checkPickedColour(ElementColour);
		 if (ItemSelected){
			fillItemSelectionInfo(i);
			SelectedItemIndex = i;
		    update();
			break;
		}
	}
 }

bool GLWidget::findElement(int i){
	int n = Sim01->Elements.size();
	if (i<n){
        ItemSelected = true;
        fillItemSelectionInfo(i);
        SelectedItemIndex = i;
        return true;
	}
	return false;
 }

 void GLWidget::fillItemSelectionInfo(int i){
	SelectedItemName = Sim01->Elements[i]->getName();
    size_t nNodes = Sim01->Elements[i]->getNodeNumber();
    size_t nDim = Sim01->Elements[i]->getDim();
    for (size_t j=0;j<nNodes;j++){
        for (size_t k =0 ;k<nDim; k++){
            QString tmpstring = QString::number(Sim01->Elements[i]->Positions[j][k], 'f', 2);
            SelectedPos.push_back(tmpstring);
		}
		QString tmpstring = QString::number(Sim01->Elements[i]->NodeIds[j], 'f', 0);
		SelectedId.push_back(tmpstring);
	}
 }

 bool GLWidget::checkPickedColour(std::array<int,3> ElementColour){
	 for (int i=0; i<3; ++i){
		 if (ElementColour[i] != PickedColour[i]){
			  return false;
		 }
	 }
	 return true;

 }

 void GLWidget::drawForPicking(){
     isDrawingForPicking=1;
     glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
     glUniform1i(uniPicking, isDrawingForPicking);
     glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
     glBufferData(GL_ELEMENT_ARRAY_BUFFER,GLElements.size() * sizeof(int), &GLElements.front(), GL_DYNAMIC_DRAW);
     glDrawElements(GL_TRIANGLES, GLElements.size(), GL_UNSIGNED_INT, 0);
 }
