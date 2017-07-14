CXX := g++
CXXFLAGS := -O3 -march=native
all: mcsp mcsp_tighter_bounding

mcsp: mcsp.c graph.c graph.h
	$(CXX) $(CXXFLAGS) -Wall -std=c++11 -o mcsp graph.c mcsp.c -pthread

mcsp_tighter_bounding: mcsp.c graph.c graph.h
	$(CXX) $(CXXFLAGS) -Wall -std=c++11 -o mcsp_tighter_bounding graph.c mcsp.c -pthread -DTIGHTER_BOUNDING

clean:
	rm mcsp mcsp_tighter_bounding
