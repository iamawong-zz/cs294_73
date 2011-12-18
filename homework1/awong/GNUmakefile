

a.out: MDArray.H Box.H MDArray.o Box.o mdarrayMain.o WriteMDArray.o VisitWriter.o
	g++ MDArray.o Box.o mdarrayMain.o  WriteMDArray.o VisitWriter.o

MDArray.o: MDArray.H Box.H MDArray.cpp
	g++ -c -g -o MDArray.o MDArray.cpp

Box.o: Box.H Box.cpp
	g++ -c -g -o Box.o Box.cpp

mdarrayMain.o: mdarrayMain.cpp Box.H MDArray.H WriteMDArray.H ../../resources/VisitWriter/VisitWriter.H
	g++ -c -g -o mdarrayMain.o mdarrayMain.cpp

VisitWriter.o: ../../resources/VisitWriter/VisitWriter.cpp ../../resources/VisitWriter/VisitWriter.H
	g++ -I ../../resources/VisitWriter -c -g  -o VisitWriter.o ../../resources/VisitWriter/VisitWriter.cpp

WriteMDArray.o: WriteMDArray.cpp WriteMDArray.H ../../resources/VisitWriter/VisitWriter.H MDArray.H Box.H
	g++ -I ../../resources/VisitWriter -c -g -o WriteMDArray.o WriteMDArray.cpp
