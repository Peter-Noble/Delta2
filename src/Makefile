CC = icx

# Excluded: $(wildcard ./timestepping/*.cpp)
strategy_src = $(wildcard ./strategies/broad_phase/*.cpp) $(wildcard ./strategies/contact_detection/*.cpp) $(wildcard ./strategies/contact_detection_continuous/*.cpp) $(wildcard ./strategies/contact_force/*.cpp) $(wildcard ./strategies/friction/*.cpp) $(wildcard ./strategies/PDE/*.cpp) $(wildcard ./strategies/time_step_size/*.cpp)
all_src = $(wildcard *.cpp) $(wildcard ./model/*.cpp) $(wildcard ./common/*.cpp) $(wildcard ./collision_detection/*.cpp) $(shell ls ./embree_common/**/*.cpp | grep -v "cmake" | grep -v "taskschedulerppl" | grep -v "taskschedulerinternal") $(wildcard ./quickhull/*.cpp) $(wildcard ./quickhull/Structs/*.cpp) $(wildcard ./scenarios/*.cpp) $(strategy_src)
test_src = $(wildcard test_*.cpp) $(wildcard ./model/test_*.cpp) $(wildcard ./common/test_*.cpp) $(wildcard ./collision_detection/test_*.cpp)
src = $(filter-out $(test_src), $(all_src))
clean_src = $(all_src) runner/main.cpp runner/unit_test_runner.cpp
obj = $(src:.cpp=.o)
release_obj = $(addsuffix .release, $(obj))
test_obj = $(test_src:.cpp=.o)
all_obj = $(obj) $(test_obj)
clean_obj = $(clean_src:.cpp=.o)
dep = $(obj:.o=.d)  # one dependency file for each source
test_dep = $(test_obj:.o=.d)
all_dep = $(dep) $(test_dep)
clean_dep = $(clean_obj:.o=.d)

is_hamilton = 0

ifeq ($(is_hamilton), 1)
	ifeq ($(CC), $(filter $(CC),icx dpcpp))
		VINCLUDE = -I/home/wxmv95/tetgen
		VLIBNAMES = /home/wxmv95/tetgen/libtet.a
		VCFLAGS = -qopenmp
		VFLAGS = -qopenmp
	else
		VINCLUDE = -I/home/wxmv95/tetgen
		VLIBNAMES = /home/wxmv95/tetgen/libtet.a
		VCFLAGS = -fopenmp
		VFLAGS = -fopenmp
	endif
	DEFINES = -DTASKING_TBB -DNOGL

	INCLUDE = -I/home/wxmv95/eigen -I/home/wxmv95/simde -I/home/wxmv95/embree-3.13.2.x86_64.linux/include/embree3 -I/home/wxmv95/embree-3.13.2.x86_64.linux/include -I/home/wxmv95/Catch2/single_include -I/home/wxmv95/libigl/include -I/usr/include -I/opt/intel/oneapi/tbb/latest/include -I/home/wxmv95/CGAL-5.5.1/include $(VINCLUDE)
	LIBS = -L/home/wxmv95/embree-3.13.2.x86_64.linux/lib -lgmp -lstdc++ -lm -ldl -ltbb
	LIBNAMES = /home/wxmv95/embree-3.13.2.x86_64.linux/lib/libembree3.so.3 $(VLIBNAMES)
else
	ifeq ($(CC), $(filter $(CC),icx dpcpp))
		VINCLUDE = -I/home/peter/Documents/tetgen1.6.0 -I/home/peter/Documents/glad/include
		VLIBNAMES = /home/peter/Documents/tetgen1.6.0/libtet.a /home/peter/Documents/glad/src/libglad.a
		VCFLAGS = -qopenmp
		VFLAGS = -qopenmp
	else
		VINCLUDE = -I/home/peter/Documents/tetgen1.6.0-clang -I/home/peter/Documents/glad-clang/include
		VLIBNAMES = /home/peter/Documents/tetgen1.6.0-clang/libtet.a /home/peter/Documents/glad-clang/src/libglad.a
		VCFLAGS = -fopenmp
		VFLAGS = -fopenmp
	endif
	DEFINES = -DTASKING_TBB

	INCLUDE = -I/home/peter/Documents/eigen -I/opt/intel/oneapi/vtune/latest/include -I/home/peter/Documents/simde -I/home/peter/Documents/embree-3.13.2.x86_64.linux/include/embree3 -I/home/peter/Documents/embree-3.13.2.x86_64.linux/include -I/home/peter/Documents/Catch2/single_include -I/home/peter/Documents/libigl/include -I/usr/include -I/opt/intel/oneapi/tbb/latest/include $(VINCLUDE)
	LIBS = -L/home/peter/Documents/embree-3.13.2.x86_64.linux/lib -lgmp -lstdc++ -lm -ldl -ltbb
	LIBNAMES = /opt/intel/oneapi/vtune/latest/sdk/lib64/libittnotify.a /home/peter/Documents/embree-3.13.2.x86_64.linux/lib/libembree3.so.3 /usr/lib/x86_64-linux-gnu/libglfw.so.3.3 $(VLIBNAMES)
endif

CFLAGS = -std=c++17 -DCGAL_DISABLE_ROUNDING_MATH_CHECK $(VCFLAGS) # -qopt-report
FLAGS = -march=native $(VFLAGS)

#  -fsanitize=memory -fno-optimize-sibling-calls

# DEFINES = -DEXCLUDE_GOOGLE_BENCH -DLIKWID_PERFMON

#-I/ddn/home/wxmv95/catch2 -I/ddn/apps/Cluster-Apps/likwid/5.0.1/include/ 
# -I/home/wxmv95/likwid-5.2.0/src/includes
# -L/home/wxmv95/likwid-5.2.0

release: mainrelease
debug: maindebug
alldebug: maindebug test
debugnolink: mainnolink

mainrelease: $(release_obj) runner/main.o.release
	$(CC) $(FLAGS) -o ../build/Delta2 $^ $(LIBS) $(LIBNAMES) -O3

maindebug: $(obj) runner/main.o
	$(CC) $(FLAGS) -o ../build/Delta2.debug $^ $(LIBS) $(LIBNAMES) -O0

mainnolink: $(obj) runner/main.o

test: $(all_obj) runner/unit_test_runner.o 
	$(CC) $(FLAGS) -o ../build/test_Delta2 $^ $(LIBS) $(LIBNAMES)

ifeq (cleandebug, $(filter cleandebug,$(MAKECMDGOALS)))
 dont_include_deps = 1
endif

ifeq (cleanrelease, $(filter cleanrelease,$(MAKECMDGOALS)))
 dont_include_deps = 1
endif

ifeq (cleandep, $(filter cleandep,$(MAKECMDGOALS)))
 dont_include_deps = 1
endif

ifeq (clean, $(filter clean,$(MAKECMDGOALS)))
 dont_include_deps = 1
endif

ifndef dont_include_deps
 -include $(all_dep)   # include all dep files in the makefile
endif

# rule to generate a dep file by using the C preprocessor
# (see man cpp for details on the -MM and -MT options)
%.d: %.cpp
	$(CC) $(CFLAGS) $(FLAGS) $(DEFINES) $(INCLUDE) $< -MM -MT $(@:.d=.o) >$@

%.o.release: %.cpp
	$(CC) $(CFLAGS) -O3 -g $(FLAGS) $(DEFINES) -DNDEBUG $(INCLUDE) -o $@ -c $<

%.o: %.cpp
	$(CC) $(CFLAGS) -O0 -g $(FLAGS) $(DEFINES) $(INCLUDE) -o $@ -c $<

.PHONY: cleandebug
cleandebug:
	rm -f $(all_obj) ../build/Delta2.debug ./runner/main.o

.PHONY: cleanrelease
cleanrelease:
	rm -f $(release_obj) ../build/Delta2 ./runner/main.o.release

.PHONY: cleandep
cleandep:
	rm -f $(all_dep)

.PHONY: clean
clean: cleandebug cleanrelease cleandep

.PHONY: printsrc
printsrc:
	echo ${src}
	echo ${test_src}
