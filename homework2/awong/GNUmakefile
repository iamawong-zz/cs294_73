#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.

RESOURCE_DIR=../../resources
VW=$(RESOURCE_DIR)/VisitWriter
CFLAGS = -g -Wall
#CFLAGS = -O3
CXX = g++

all: ShowMesh TestMatrix a.out
 
a.out: FEMain.o FEGrid.o FEPoissonOperator.o Element.o Node.o VisitWriter.o ReInsert.o JacobiSolver.o SparseMatrix.o
	  $(CXX) FEMain.o FEGrid.o Element.o Node.o VisitWriter.o ReInsert.o FEPoissonOperator.o JacobiSolver.o SparseMatrix.o

FEMain.o: FEMain.cpp FEGrid.H SparseMatrix.H FEPoissonOperator.H ReInsert.H
	$(CXX) $(CFLAGS) -c -o FEMain.o FEMain.cpp

TestMatrix: TestMatrix.cpp SparseMatrix.o
	$(CXX) $(CFLAGS) -o TestMatrix TestMatrix.cpp SparseMatrix.o

vars:
	echo $(VW)

ShowMesh: VisitWriter.o ShowMesh.cpp
	$(CXX) -I$(VW) $(CFLAGS) -o ShowMesh ShowMesh.cpp VisitWriter.o

VisitWriter.o: $(VW)/VisitWriter.cpp $(VW)/VisitWriter.H
	echo $(VW)
	$(CXX) -I$(VW) -c -O2 -o VisitWriter.o $(VW)/VisitWriter.cpp

SparseMatrix.o: SparseMatrix.H SparseMatrix.cpp
	$(CXX)  $(CFLAGS) -c -o SparseMatrix.o SparseMatrix.cpp 

ReInsert.o: ReInsert.H ReInsert.cpp
	$(CXX)  $(CFLAGS) -c -o ReInsert.o ReInsert.cpp 

JacobiSolver.o: JacobiSolver.H JacobiSolver.cpp
	$(CXX)  $(CFLAGS) -c -o JacobiSolver.o JacobiSolver.cpp 


Node.o: Node.cpp Node.H
	$(CXX) $(CFLAGS) -c -o Node.o Node.cpp

Element.o: Element.H Element.cpp 
	$(CXX) $(CFLAGS) -c -o Element.o Element.cpp

FEGrid.o: FEGrid.H FEGrid.cpp Element.H Node.H $(VW)/VisitWriter.H
	$(CXX) -I$(VW) $(CFLAGS) -c -o FEGrid.o FEGrid.cpp

FEPoissonOperator.o: FEPoissonOperator.H FEPoissonOperator.cpp FEGrid.H Node.H Element.H SparseMatrix.H
	$(CXX) $(CFLAGS) -c -o FEPoissonOperator.o FEPoissonOperator.cpp

testFEGrid: FEGrid.o Element.o Node.o testFEGrid.cpp VisitWriter.o
	$(CXX) $(CFLAGS) -o testFEGrid testFEGrid.cpp FEGrid.o Element.o Node.o VisitWriter.o
