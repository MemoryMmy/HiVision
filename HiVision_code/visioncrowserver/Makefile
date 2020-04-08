HEADERS =-I./ -I/usr/local/include/ -I/home/dbrg/APPS/boost_1_64_0/include 
LIBS =-L./ -L/usr/local/lib -L/home/dbrg/APPS/boost_1_64_0/lib

OBO:  server.cpp
	g++ $(HEADERS) $(LIBS) -std=c++11 -fopenmp -Wall -O3 -o Crowserver server.cpp -lboost_system -lhiredis -lpng -lgdal
clean:
	rm  Server
	
