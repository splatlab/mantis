TARGETS= mantis 
MYLIBS=-L/home/paf2023/tools/lib \
-L/pbtech_mounts/softlib001/apps/EL6/spack/opt/spack/linux-rhel6-x86_64/gcc-6.3.0/binutils-2.28-phpi22i4vio3f55thjj5xxeupzmj463k/lib \
-L/pbtech_mounts/softlib001/apps/EL6/spack/opt/spack/linux-rhel6-x86_64/gcc-6.3.0/gettext-0.19.8.1-cnbmeuapkcyxq23tumt4esnuvijpl4w7/lib \
-L/pbtech_mounts/softlib001/apps/EL6/spack/opt/spack/linux-rhel6-x86_64/gcc-6.3.0/bzip2-1.0.6-x5hwxgowvbzmnowngidcjin5qicfg4fl/lib \
-L/pbtech_mounts/softlib001/apps/EL6/spack/opt/spack/linux-rhel6-x86_64/gcc-6.3.0/libxml2-2.9.4-ox4x2yk4rqlh3ag4wxxgzficxgwbcsy7/lib \
-L/pbtech_mounts/softlib001/apps/EL6/spack/opt/spack/linux-rhel6-x86_64/gcc-6.3.0/xz-5.2.3-npj3lafufjb7d57yhkthdo3qeeollgiq/lib \
-L/pbtech_mounts/softlib001/apps/EL6/spack/opt/spack/linux-rhel6-x86_64/gcc-6.3.0/python-2.7.14-4q4ykdbphrdjx2n4ztg5getsrjkl7dkz/lib \
-L/pbtech_mounts/softlib001/apps/EL6/spack/opt/spack/linux-rhel6-x86_64/gcc-6.3.0/ncurses-6.0-37ewjt6dvfa4vn3coyfnx5ejnfzg6d27/lib \
-L/pbtech_mounts/softlib001/apps/EL6/spack/opt/spack/linux-rhel6-x86_64/gcc-6.3.0/binutils-2.28-pkmmuwsdolh4vk2qrjvevqbrds55oxet/lib \
-L/pbtech_mounts/softlib001/apps/EL6/spack/opt/spack/linux-rhel6-x86_64/gcc-6.3.0/binutils-2.28-3yfqdmrhogpxaehsissamfd4bagokf3y/lib \
-L/pbtech_mounts/softlib001/apps/EL6/spack/opt/spack/linux-rhel6-x86_64/gcc-6.3.0/binutils-2.28-kamnyqrq5i6e6ghsccvnx4banzgdvn6r/lib \
-L/pbtech_mounts/softlib001/apps/EL6/spack/opt/spack/linux-rhel6-x86_64/gcc-4.4.6/binutils-2.28-t6r4t2u7mpl73qqhveutncj4gamajiip/lib \
-L/pbtech_mounts/softlib001/apps/EL6/spack/opt/spack/linux-rhel6-x86_64/gcc-4.4.6/binutils-2.28-a4bouodbguxjzgeoogr7qjqculofd4cx/lib \
-L/pbtech_mounts/softlib001/apps/EL6/spack/opt/spack/linux-rhel6-x86_64/gcc-6.3.0/zlib-1.2.11-7eta3q6r3ozz5vh7zaffgjajf4fwvxby/lib \
-L/pbtech_mounts/softlib001/apps/EL6/spack/opt/spack/linux-rhel6-x86_64/gcc-4.5.4/zlib-1.2.11-upeb67xkoclf5fgx7musxjw5v7gk7yw7/lib \
-L/pbtech_mounts/softlib001/apps/EL6/spack/opt/spack/linux-rhel6-x86_64/gcc-4.4.6/zlib-1.2.11-erzei6xutxne2tjssbmelgutsnhz6x6s/lib \
-L/pbtech_mounts/softlib001/apps/EL6/spack/opt/spack/linux-rhel6-x86_64/gcc-6.3.0/zlib-1.2.11-rjqs3cq27e2o4tbce7kzmictrih57avc/lib \
-L/pbtech_mounts/softlib001/apps/EL6/spack/opt/spack/linux-rhel6-x86_64/gcc-6.3.0/boost-1.65.1-oiwzkauf2kl22pu7j2vjupr3ffu3mzts/lib \
-L/pbtech_mounts/softlib001/apps/EL6/spack/opt/spack/linux-rhel6-x86_64/gcc-6.3.0/boost-1.65.1-5awddeormrqlmqq2frhjrdcu5szebebq/lib \
-L/pbtech_mounts/softlib001/apps/EL6/spack/opt/spack/linux-rhel6-x86_64/gcc-6.3.0/boost-1.65.1-sqviukvzkgl524rq2dtz6d5ytukaxlj3/lib \
-L/pbtech_mounts/softlib001/apps/EL6/spack/opt/spack/linux-rhel6-x86_64/gcc-6.3.0/boost-1.64.0-sx2zxhmitq565rmn3jutk67ffqrko6zx/lib \
-L/pbtech_mounts/softlib001/apps/EL6/spack/opt/spack/linux-rhel6-x86_64/gcc-4.4.6/boost-1.63.0-gpfsh3erjtf2xyksy66crnxzqib6qrva/lib \
-L/pbtech_mounts/softlib001/apps/EL6/spack/opt/spack/linux-rhel6-x86_64/gcc-6.3.0/boost-1.63.0-eu4aw5v6ewmqd5kuzieorttiljxcuscv/lib \
-L/pbtech_mounts/softlib001/apps/EL6/spack/opt/spack/linux-rhel6-x86_64/gcc-6.3.0/boost-1.54.0-7qaxqizxejlur5pr5mfecf2umfp52lmr/lib \
-L/pbtech_mounts/softlib001/apps/EL6/spack/opt/spack/linux-rhel6-x86_64/gcc-4.4.6/gcc-6.3.0-han2ovjpjdnuwkhtfhugdq66cjkte7q3/lib64 \
-L/pbtech_mounts/softlib001/apps/EL6/spack/opt/spack/linux-rhel6-x86_64/gcc-4.4.6/gcc-6.3.0-xnoykzftyh27xs6vqkwvtft2cce5yd4k/lib64 \
-L/pbtech_mounts/softlib001/apps/EL6/spack/opt/spack/linux-rhel6-x86_64/gcc-4.4.6/gcc-4.5.4-tmeuqohwhscsb2hegxgfjmtxy7ywelud/lib64 \

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

LDFLAGS += $(MYLIBS) $(DEBUG) $(PROFILE) $(OPT) -lsdsl -lpthread -lboost_system \
-lboost_thread -lm -lz -lrt

#
# declaration of dependencies
#

all: $(TARGETS)

# dependencies between programs and .o files
#mantis:									$(OBJDIR)/kmer.o $(OBJDIR)/mantis.o $(OBJDIR)/validatemantis.o $(OBJDIR)/gqf.o $(OBJDIR)/hashutil.o $(OBJDIR)/query.o $(OBJDIR)/queries.o $(OBJDIR)/coloreddbg.o $(OBJDIR)/bitvector.o $(OBJDIR)/util.o  $(OBJDIR)/MantisFS.o
mantis:									$(OBJDIR)/kmer.o $(OBJDIR)/mantis.o $(OBJDIR)/validatemantis.o $(OBJDIR)/gqf.o $(OBJDIR)/hashutil.o $(OBJDIR)/queries.o $(OBJDIR)/coloreddbg.o $(OBJDIR)/bitvector.o $(OBJDIR)/util.o  $(OBJDIR)/MantisFS.o

# dependencies between .o files and .h files
$(OBJDIR)/mantis.o:					$(LOC_SRC)/mantis.cc
$(OBJDIR)/MantisFs.o:       $(LOC_SRC)/MantisFS.cc $(LOC_INCLUDE)/MantisFS.h
$(OBJDIR)/util.o:           $(LOC_SRC)/util.cc $(LOC_INCLUDE)/util.h
$(OBJDIR)/bitvector.o:      $(LOC_SRC)/bitvector.cc $(LOC_INCLUDE)/bitvector.h
$(OBJDIR)/kmer.o:           $(LOC_SRC)/kmer.cc $(LOC_INCLUDE)/kmer.h
$(OBJDIR)/coloreddbg.o: 		$(LOC_INCLUDE)/cqf/gqf.h $(LOC_INCLUDE)/hashutil.h $(LOC_INCLUDE)/util.h $(LOC_INCLUDE)/coloreddbg.h $(LOC_INCLUDE)/bitvector.h $(LOC_INCLUDE)/cqf.h
#$(OBJDIR)/query.o: 					$(LOC_INCLUDE)/cqf/gqf.h $(LOC_INCLUDE)/hashutil.h $(LOC_INCLUDE)/util.h $(LOC_INCLUDE)/coloreddbg.h $(LOC_INCLUDE)/bitvector.h $(LOC_INCLUDE)/cqf.h $(LOC_INCLUDE)/kmer.h
$(OBJDIR)/queries.o: 					$(LOC_INCLUDE)/cqf/gqf.h $(LOC_INCLUDE)/hashutil.h $(LOC_INCLUDE)/util.h $(LOC_INCLUDE)/coloreddbg.h $(LOC_INCLUDE)/bitvector.h $(LOC_INCLUDE)/cqf.h $(LOC_INCLUDE)/kmer.h
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
