all:
	g++ -fopenmp -g sg_metropolis.cpp -o sg_metropolis -lboost_system -lboost_filesystem
clean:
	rm sg_metropolis
