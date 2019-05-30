CurrPath = {MainClonePath}/TissueFolding/

TARGET = $$CurrPath/Release/TissueFolding

QMAKE_CFLAGS_RELEASE += -fopenmp
QMAKE_CFLAGS_DEBUG += -fopenmp
QMAKE_CXXFLAGS += -fopenmp -openmp -std=c++17
QMAKE_LFLAGS +=  -fopenmp

# Input
OBJECTS_DIR += $$CurrPath/Obj/

HEADERS       +=  	$$CurrPath/SourceCode/*.h

SOURCES += $$CurrPath/SourceCode/main.cpp \
	$$CurrPath/SourceCode/Prism.cpp \
	$$CurrPath/SourceCode/ReferenceShapeBase.cpp \
 	$$CurrPath/SourceCode/ShapeBase.cpp \
        $$CurrPath/SourceCode/Simulation.cpp \
        $$CurrPath/SourceCode/Node.cpp \
        $$CurrPath/SourceCode/ModelInputObject.cpp \
	$$CurrPath/SourceCode/RandomGenerator.cpp \
	$$CurrPath/SourceCode/NewtonRaphsonSolver.cpp 

LIBS += -L/usr/include -lgsl -lgslcblas -lgomp -fopenmp -lpardiso600-GNU720-X86-64 -llapack 

# install
target.path   =  $$CurrPath/
sources.files =  $$SOURCES $$HEADERS $$RESOURCES TissueFolding.pro
sources.path  =  $$CurrPath
INSTALLS += target sources
