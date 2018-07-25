# Load the common configuration file
include config.mk
objs=
src=$(patsubst %.o,%.cc,$(objs))
execs=bck_banded bck_thomas cn_thomas bck_thomas_adaptivex

all: $(execs)


bck_banded: bck_banded.cc
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lp_lflags)

bck_thomas: bck_thomas.cc
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lp_lflags)

cn_thomas: cn_thomas.cc
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lp_lflags)

bck_thomas_adaptivex: bck_thomas_adaptivex.cc
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lp_lflags)

clean:
	rm -f $(execs) $(objs) libf2d.a

.PHONY: clean depend
