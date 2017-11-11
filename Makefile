TARGETS= coloreddbg query validatemantis

ifdef D
	DEBUG=-g -DDEBUG
	OPT=
else
	DEBUG=
	OPT=-Ofast
endif

ifdef NH
	ARCH=
else
	ARCH=-msse4.2 -D__SSE4_2_
endif

ifdef P
	PROFILE=-pg -no-pie # for bug in gprof.
endif

CXX = g++ -std=c++11
CC = gcc -std=gnu11
LD= g++ -std=c++11

CXXFLAGS += -Wall $(DEBUG) $(PROFILE) $(OPT) $(ARCH) -m64 -I. -Isdsl/include -Iinclude \
-Wno-unused-result -Wno-strict-aliasing -Wno-unused-function -Wno-sign-compare

CFLAGS += -Wall $(DEBUG) $(PROFILE) $(OPT) $(ARCH) -m64 -I. \
-Wno-unused-result -Wno-strict-aliasing -Wno-unused-function -Wno-sign-compare \
-Wno-implicit-function-declaration

LDFLAGS += $(DEBUG) $(PROFILE) $(OPT) -lsdsl -lpthread -lboost_system \
-lboost_thread -lm -lz -lrt

#
# declaration of dependencies
#

all: $(TARGETS)

# dependencies between programs and .o files

coloreddbg:             coloreddbg.o cqf/gqf.o hashutil.o
query:             			query.o cqf/gqf.o hashutil.o
validatemantis:     		validatemantis.o cqf/gqf.o hashutil.o

# dependencies between .o files and .h files

coloreddbg.o: 		cqf/gqf.h hashutil.h util.h coloreddbg.h bitvector.h cqf.h
query.o: 					cqf/gqf.h hashutil.h util.h coloreddbg.h bitvector.h cqf.h kmer.h
validatemantis.o: cqf/gqf.h hashutil.h util.h coloreddbg.h bitvector.h cqf.h kmer.h
hashutil.o: 								hashutil.h

# dependencies between .o files and .cc (or .c) files

cqf/gqf.o: cqf/gqf.c cqf/gqf.h

#
# generic build rules
#

$(TARGETS):
	$(LD) $^ $(LDFLAGS) -o $@

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDE) $< -c -o $@

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE) $< -c -o $@

clean:
	rm -f *.o core cqf/gqf.o $(TARGETS)
