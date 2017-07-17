CXX := g++
CXXFLAGS := -O3 -march=native
all: mcsp mcsp_tighter_bounding mcsp_path_len mcsp_path_len_tighter_bounding mcsp_2paths mcsp_2paths_tighter_bounding

mcsp: mcsp.c graph.c graph.h
	$(CXX) $(CXXFLAGS) -Wall -std=c++11 -o mcsp graph.c mcsp.c -pthread

mcsp_tighter_bounding: mcsp.c graph.c graph.h
	$(CXX) $(CXXFLAGS) -Wall -std=c++11 -o mcsp_tighter_bounding graph.c mcsp.c -pthread -DTIGHTER_BOUNDING

mcsp_path_len: mcsp_path_len.c graph.c graph.h
	$(CXX) $(CXXFLAGS) -Wall -std=c++11 -o mcsp_path_len graph.c mcsp_path_len.c -pthread

mcsp_path_len_tighter_bounding: mcsp_path_len.c graph.c graph.h
	$(CXX) $(CXXFLAGS) -Wall -std=c++11 -o mcsp_path_len_tighter_bounding graph.c mcsp_path_len.c -pthread -DTIGHTER_BOUNDING

mcsp_2paths: mcsp_2paths.c graph.c graph.h
	$(CXX) $(CXXFLAGS) -Wall -std=c++11 -o mcsp_2paths graph.c mcsp_2paths.c -pthread

mcsp_2paths_tighter_bounding: mcsp_2paths.c graph.c graph.h
	$(CXX) $(CXXFLAGS) -Wall -std=c++11 -o mcsp_2paths_tighter_bounding graph.c mcsp_2paths.c -pthread -DTIGHTER_BOUNDING

clean:
	rm mcsp mcsp_tighter_bounding
