# Include the configuration file
include config.mk

# Lists of files to be built
objs=butler_volmer.o
src=$(patsubst %.o,%.cc,$(objs) $(mg_objs))
execs=bv_test bv_snaps

all: $(execs)

# Include the file dependencies
-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) $(iflags) -MM $(src) >Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm $(objs) $(execs)

%.o: %.cc
	$(cxx) $(cflags) $(iflags) -c $<

bv_test: bv_test.cc $(objs)
	$(cxx) $(cflags) -o $@ $^

bv_snaps: bv_snaps.cc $(objs)
	$(cxx) $(cflags) -o $@ $^

.PHONY: clean depend
