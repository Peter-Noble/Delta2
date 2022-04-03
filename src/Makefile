CC = dpcpp


# Excluded: $(wildcard ./timestepping/*.cpp)
strategy_src = $(wildcard ./strategies/broad_phase/*.cpp) $(wildcard ./strategies/contact_detection/*.cpp) $(wildcard ./strategies/contact_detection_continuous/*.cpp) $(wildcard ./strategies/contact_force/*.cpp) $(wildcard ./strategies/friction/*.cpp) $(wildcard ./strategies/PDE/*.cpp) $(wildcard ./strategies/time_step_size/*.cpp)
all_src = $(wildcard *.cpp) $(wildcard ./model/*.cpp) $(wildcard ./common/*.cpp) $(wildcard ./collision_detection/*.cpp) $(shell ls ./embree_common/**/*.cpp | grep -v "cmake" | grep -v "taskschedulerppl" | grep -v "taskschedulertbb") $(wildcard ./quickhull/*.cpp) $(wildcard ./quickhull/Structs/*.cpp) $(strategy_src)
test_src = $(wildcard test_*.cpp) $(wildcard ./model/test_*.cpp) $(wildcard ./common/test_*.cpp) $(wildcard ./timestepping/test_*.cpp) $(wildcard ./collision_detection/test_*.cpp)
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

CFLAGS = -std=c++17 -qopenmp
FLAGS = -march=native -g -qopenmp
DEFINES = -DTASKING_INTERNAL

# DEFINES = -DEXCLUDE_GOOGLE_BENCH -DLIKWID_PERFMON

#-I/ddn/home/wxmv95/catch2 -I/ddn/apps/Cluster-Apps/likwid/5.0.1/include/ 
# -I/home/peter/Documents/likwid-5.2.0/src/includes
# -L/home/peter/Documents/likwid-5.2.0
INCLUDE = -I/home/peter/Documents/eigen -I/opt/intel/oneapi/vtune/latest/include -I/home/peter/Documents/simde -I/home/peter/Documents/embree-3.13.2.x86_64.linux/include/embree3 -I/home/peter/Documents/embree-3.13.2.x86_64.linux/include -I/home/peter/Documents/Catch2/single_include -I/home/peter/Documents/libigl/include -I/home/peter/Documents/tetgen1.6.0 -I/home/peter/Documents/glad/include -I/usr/include
LIBS = -L/home/peter/Documents/embree-3.13.2.x86_64.linux/lib -lgmp
LIBNAMES = /opt/intel/oneapi/vtune/latest/sdk/lib64/libittnotify.a /home/peter/Documents/embree-3.13.2.x86_64.linux/lib/libembree3.so.3 /home/peter/Documents/tetgen1.6.0/libtet.a /home/peter/Documents/glad/src/libglad.a /usr/lib/x86_64-linux-gnu/libglfw.so.3.3

release: mainrelease
debug: maindebug
alldebug: maindebug test

mainrelease: $(release_obj) runner/main.o.release
	$(CC) $(FLAGS) $(LIBS) -o ../build/Delta2 $^ $(LIBNAMES)

maindebug: $(obj) runner/main.o
	$(CC) $(FLAGS) $(LIBS) -o ../build/Delta2.debug $^ $(LIBNAMES)

test: $(all_obj) runner/unit_test_runner.o 
	$(CC) $(FLAGS) $(LIBS) -o ../build/test_Delta2 $^ $(LIBNAMES)

 -include $(all_dep)   # include all dep files in the makefile

# rule to generate a dep file by using the C preprocessor
# (see man cpp for details on the -MM and -MT options)
%.d: %.cpp
	$(CC) $(CFLAGS) $(FLAGS) $(DEFINES) $(INCLUDE) $< -MM -MT $(@:.d=.o) >$@

%.o.release: %.cpp
	$(CC) $(CFLAGS) -O3 $(FLAGS) $(DEFINES) $(INCLUDE) -o $@ -c $<

%.o: %.cpp
	$(CC) $(CFLAGS) -O0 $(FLAGS) $(DEFINES) $(INCLUDE) -o $@ -c $<

.PHONY: clean
clean:
	cleandebug
	cleanrelease
	cleandep
.PHONY: cleandebug
cleandebug:
	rm -f $(all_obj) Delta2
.PHONY: cleanrelease
cleanrelease:
	rm -f $(release_obj) Delta2

.PHONY: cleandep
cleandep:
	rm -f $(all_dep)

.PHONY: printsrc
printsrc:
	echo ${src}
	echo ${test_src}
