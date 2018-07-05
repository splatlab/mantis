TARGETS= mantis 

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

LOC_INCLUDE=include
LOC_SRC=src
OBJDIR=obj

CXXFLAGS += -Wall $(DEBUG) $(PROFILE) $(OPT) $(ARCH) -m64 -I. -I$(LOC_INCLUDE) \
-Wno-unused-result -Wno-strict-aliasing -Wno-unused-function -Wno-sign-compare

CFLAGS += -Wall $(DEBUG) $(PROFILE) $(OPT) $(ARCH) -m64 -I. -I$(LOC_INCLUDE)\
-Wno-unused-result -Wno-strict-aliasing -Wno-unused-function -Wno-sign-compare \
-Wno-implicit-function-declaration

LDFLAGS += $(DEBUG) $(PROFILE) $(OPT) -lsdsl -lpthread -lboost_system \
-lboost_thread -lm -lz -lrt

#
# declaration of dependencies
#

all: $(TARGETS)

# dependencies between programs and .o files
mantis:									$(OBJDIR)/kmer.o $(OBJDIR)/mantis.o $(OBJDIR)/validatemantis.o $(OBJDIR)/gqf.o $(OBJDIR)/hashutil.o $(OBJDIR)/query.o $(OBJDIR)/coloreddbg.o $(OBJDIR)/bitvector.o $(OBJDIR)/util.o  $(OBJDIR)/MantisFS.o

# dependencies between .o files and .h files
$(OBJDIR)/mantis.o:					$(LOC_SRC)/mantis.cc
$(OBJDIR)/MantisFs.o:       $(LOC_SRC)/MantisFS.cc $(LOC_INCLUDE)/MantisFS.h
$(OBJDIR)/util.o:           $(LOC_SRC)/util.cc $(LOC_INCLUDE)/util.h
$(OBJDIR)/bitvector.o:      $(LOC_SRC)/bitvector.cc $(LOC_INCLUDE)/bitvector.h
$(OBJDIR)/kmer.o:           $(LOC_SRC)/kmer.cc $(LOC_INCLUDE)/kmer.h
$(OBJDIR)/coloreddbg.o: 		$(LOC_INCLUDE)/cqf/gqf.h $(LOC_INCLUDE)/hashutil.h $(LOC_INCLUDE)/util.h $(LOC_INCLUDE)/coloreddbg.h $(LOC_INCLUDE)/bitvector.h $(LOC_INCLUDE)/cqf.h
$(OBJDIR)/query.o: 					$(LOC_INCLUDE)/cqf/gqf.h $(LOC_INCLUDE)/hashutil.h $(LOC_INCLUDE)/util.h $(LOC_INCLUDE)/coloreddbg.h $(LOC_INCLUDE)/bitvector.h $(LOC_INCLUDE)/cqf.h $(LOC_INCLUDE)/kmer.h
$(OBJDIR)/validatemantis.o: $(LOC_INCLUDE)/cqf/gqf.h $(LOC_INCLUDE)/hashutil.h $(LOC_INCLUDE)/util.h $(LOC_INCLUDE)/coloreddbg.h $(LOC_INCLUDE)/bitvector.h $(LOC_INCLUDE)/cqf.h $(LOC_INCLUDE)/kmer.h
$(OBJDIR)/hashutil.o: 			$(LOC_INCLUDE)/hashutil.h

# dependencies between .o files and .cc (or .c) files

$(OBJDIR)/gqf.o: $(LOC_SRC)/cqf/gqf.c $(LOC_INCLUDE)/cqf/gqf.h

#
# generic build rules
#

$(TARGETS):
	$(LD) $^ $(LDFLAGS) -o $@

$(OBJDIR)/%.o: $(LOC_SRC)/%.cc | $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $@ $<

$(OBJDIR)/%.o: $(LOC_SRC)/%.c | $(OBJDIR)
	$(CC) $(CFLAGS) $(INCLUDE) -c -o $@ $<

$(OBJDIR)/%.o: $(LOC_SRC)/cqf/%.c | $(OBJDIR)
	$(CC) $(CFLAGS) $(INCLUDE) -c -o $@ $<

$(OBJDIR):
	@mkdir -p $(OBJDIR)

clean:
	rm -f $(OBJDIR)/*.o core $(TARGETS)
