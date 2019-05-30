CurrPath = {MainClonePath}/UserInterface/

TARGET = $$CurrPath/Release/TissueFoldingUI

QMAKE_CFLAGS_RELEASE += -fopenmp
QMAKE_CFLAGS_DEBUG += -fopenmp
QMAKE_LFLAGS +=  -fopenmp
QMAKE_CXXFLAGS += -std=c++17 -fopenmp -openmp

#Most straightfirward -but ugly- way to make openGL VBOs work on generic linux. Modify for your curren OS
DEFINES += "GL_GLEXT_PROTOTYPES"

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
QT += opengl

# Input
OBJECTS_DIR += $$CurrPath/Obj/

HEADERS       +=  	$$CurrPath/SourceCode/MainWindow.h \
                        $$CurrPath/SourceCode/GLWidget.h \
                        $$CurrPath/../TissueFolding/SourceCode/*.h

SOURCES += $$CurrPath/SourceCode/main.cpp \
        $$CurrPath/SourceCode/MainWindow.cpp \
        $$CurrPath/SourceCode/GLWidget.cpp \
        $$CurrPath/../TissueFolding/SourceCode/Prism.cpp \
        $$CurrPath/../TissueFolding/SourceCode/ReferenceShapeBase.cpp \
        $$CurrPath/../TissueFolding/SourceCode/ShapeBase.cpp \
        $$CurrPath/../TissueFolding/SourceCode/Simulation.cpp \
        $$CurrPath/../TissueFolding/SourceCode/Node.cpp \
        $$CurrPath/../TissueFolding/SourceCode/ModelInputObject.cpp \
        $$CurrPath/../TissueFolding/SourceCode/RandomGenerator.cpp \
        $$CurrPath/../TissueFolding/SourceCode/NewtonRaphsonSolver.cpp

#libs and includes for linux for independent license pardiso:
LIBS += -L/usr/include -lgsl -lgslcblas -lpardiso600-GNU720-X86-64  -fopenmp -llapack -lgomp -lpthread -lgfortran -lm 

# install
target.path   =  $$CurrPath/
sources.files =  $$SOURCES $$HEADERS $$RESOURCES TissueFoldingUI.pro
sources.path  =  $$CurrPath
INSTALLS += target sources

