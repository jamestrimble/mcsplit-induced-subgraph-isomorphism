CXX := g++
CXXFLAGS := -O0 -march=native -g -ggdb -fsanitize=address
#CXXFLAGS := -O3 -march=native -g
PROGRAMS := mcsp_sparse mcsp mcsp_path_len mcsp_2paths mcsp_lazy non_ind_cp_2paths non_ind_mcsp_2paths non_ind_mcsp_2paths_with_bells_and_whistles non_ind_mcsp_2paths_restarts non_ind_mcsp_2paths_restarts_simple_nogoods
TIGHTER_BOUNDING_PROGRAMS := mcsp_tighter_bounding mcsp_path_len_tighter_bounding mcsp_2paths_tighter_bounding

all: $(PROGRAMS) $(TIGHTER_BOUNDING_PROGRAMS)

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
