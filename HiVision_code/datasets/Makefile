HEADERS =-I./ -I/usr/local/include/
LIBS =-L./ -L/usr/local/lib

OBO:  pointindex.cpp linestringindex.cpp polygonindex.cpp
	g++ $(HEADERS) $(LIBS) -std=c++11 -fopenmp -Wall -O3 -o Pointindex pointindex.cpp -lgdal -lhiredis 
	g++ $(HEADERS) $(LIBS) -std=c++11 -fopenmp -Wall -O3 -o Linestringindex linestringindex.cpp -lgdal -lhiredis 
	g++ $(HEADERS) $(LIBS) -std=c++11 -fopenmp -Wall -O3 -o Polygonindex polygonindex.cpp -lgdal -lhiredis 
clean:
	rm  Pointindex Linestringindex Polygonindex
	
