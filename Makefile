CXX := g++
#CXXFLAGS := -O0 -g -ggdb -fsanitize=address
CXXFLAGS := -O3 -march=native -g
PROGRAMS := mcsp mcsp_path_len mcsp_2paths mcsp_lazy non_ind_cp_2paths non_ind_mcsp_2paths non_ind_mcsp_2paths_with_bells_and_whistles non_ind_mcsp_2paths_restarts non_ind_mcsp_2paths_restarts_simple_nogoods mcsp_optimised mcsp_tmp
TIGHTER_BOUNDING_PROGRAMS := mcsp_tighter_bounding mcsp_path_len_tighter_bounding mcsp_2paths_tighter_bounding

all: $(PROGRAMS) $(TIGHTER_BOUNDING_PROGRAMS) mcsp_sparse

mcsp_sparse: mcsp_sparse.c sparse_graph.c sparse_graph.h
	$(CXX) $(CXXFLAGS) -Wall -std=c++11 -o mcsp_sparse sparse_graph.c mcsp_sparse.c -pthread

mcsp_bitsets: mcsp_bitsets.c sparse_graph.c sparse_graph.h
	$(CXX) $(CXXFLAGS) -Wall -std=c++11 -o mcsp_bitsets sparse_graph.c mcsp_bitsets.c -pthread

mcsp_non_ind_bitsets: mcsp_non_ind_bitsets.c sparse_graph.c sparse_graph.h
	$(CXX) $(CXXFLAGS) -Wall -std=c++11 -o mcsp_non_ind_bitsets sparse_graph.c mcsp_non_ind_bitsets.c -pthread

mcsp_non_ind_bitsets_supps: mcsp_non_ind_bitsets_supps.c sparse_graph.c sparse_graph.h
	$(CXX) $(CXXFLAGS) -Wall -std=c++11 -o mcsp_non_ind_bitsets_supps sparse_graph.c mcsp_non_ind_bitsets_supps.c -pthread

define PROGRAM_template =
$(1): $(1).c graph.c graph.h
	$(CXX) $(CXXFLAGS) -Wall -std=c++11 -o $(1) graph.c $(1).c -pthread
endef

$(foreach program,$(PROGRAMS),$(eval $(call PROGRAM_template,$(program))))

define TIGHTER_BOUNDING_PROGRAM_template =
$(1): $(subst _tighter_bounding,,$(1)).c graph.c graph.h
	$(CXX) $(CXXFLAGS) -Wall -std=c++11 -o $(1) graph.c $(subst _tighter_bounding,,$(1)).c -pthread -DTIGHTER_BOUNDING
endef

$(foreach program,$(TIGHTER_BOUNDING_PROGRAMS),$(eval $(call TIGHTER_BOUNDING_PROGRAM_template,$(program))))

clean:
	rm $(PROGRAMS)
