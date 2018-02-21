CVERSION = -std=c++0x
CC = g++ $(CVERSION)
SRCDIR = src
BUILDDIR = build
INCDIR = include
TARGET = bin

SRCEXT = cpp
SOURCES = $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS = $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
DEBUG = -g
CFLAGS =  -c $(DEBUG) -Wall
LFLAGS =  $(DEBUG)
LIB = -L lib
INC = -I include

COMMON_OBJS = $(BUILDDIR)/channel.o $(BUILDDIR)/function.o $(BUILDDIR)/sparsematrix.o $(BUILDDIR)/node.o  

OBJS1 = $(BUILDDIR)/rcm.o $(COMMON_OBJS)
OBJS2 = $(BUILDDIR)/hybrid.o $(BUILDDIR)/pcm.o $(COMMON_OBJS) 
OBJS3 = $(BUILDDIR)/hybrid_mac_inde.o $(BUILDDIR)/pcm.o $(COMMON_OBJS) 
OBJS4 = $(BUILDDIR)/rcm_mac.o $(COMMON_OBJS) 
OBJS5 = $(BUILDDIR)/hybrid_pam.o $(BUILDDIR)/pcm.o $(COMMON_OBJS)
OBJS6 = $(BUILDDIR)/hybrid_mac_synthetic.o $(BUILDDIR)/pcm.o $(COMMON_OBJS)

# For hybrid schemes on AWGN channels, including asymmetric sources and simplified methods.
rcm : $(OBJS1)
	$(CC) $(LFLAGS) $(OBJS1) -o $(TARGET)/$@
hybrid : $(OBJS2)
	$(CC) $(LFLAGS) $(OBJS2) -o $(TARGET)/$@

# For MAC channel independent decoding for the joint decoder.
hybrid_mac_inde : $(OBJS3)
	$(CC) $(LFLAGS) $(OBJS3) -o $(TARGET)/$@ 
$(TARGET)/rcm_mac : $(OBJS4)
	$(CC) $(LFLAGS) $(OBJS4) -o rcm_mac
$(TARGET)/hybrid_pam : $(OBJS5)
	$(CC) $(LFLAGS) $(OBJS5) -o hybrid_pam

# For MAC channel synthetic decoding.
hybrid_mac_syn : $(OBJS6)
	$(CC) $(LFLAGS) $(OBJS6) -o $(TARGET)/$@

$(BUILDDIR)/rcm.o : $(INCDIR)/node.h $(INCDIR)/sparsematrix.h $(INCDIR)/function.h point_to_point/rcm.cpp
	$(CC) $(CFLAGS) $(INC) point_to_point/rcm.cpp -o $@

$(BUILDDIR)/hybrid.o : $(INCDIR)/node.h $(INCDIR)/sparsematrix.h $(INCDIR)/function.h point_to_point/hybrid.cpp
	$(CC) $(CFLAGS) $(INC) point_to_point/hybrid.cpp -o $@

$(BUILDDIR)/hybrid_mac_inde.o : $(INCDIR)/node.h $(INCDIR)/sparsematrix.h $(INCDIR)/function.h hybrid_mac_independent/hybrid_mac_inde.cpp
	$(CC) $(CFLAGS) $(INC) hybrid_mac_independent/hybrid_mac_inde.cpp -o $@

$(BUILDDIR)/hybrid_mac_synthetic.o : $(INCDIR)/node.h $(INCDIR)/sparsematrix.h $(INCDIR)/function.h hybrid_mac_synthetic/hybrid_mac_synthetic.cpp
	$(CC) $(CFLAGS) $(INC) hybrid_mac_synthetic/hybrid_mac_synthetic.cpp -o $@

$(BUILDDIR)/hybrid_pam.o : node.h sparsematrix.h func.h hybrid_pam.cpp
	$(CC) $(CFLAGS) $(INC) hybrid_pam.cpp

$(BUILDDIR)/rcm_mac.o : node.h sparsematrix.h func.h rcm_mac.cpp
	$(CC) $(CFLAGS) $(INC) rcm_mac.cpp

$(BUILDDIR)/hybrid_mac_syn.o : $(INCDIR)/node.h $(INCDIR)/sparsematrix.h $(INCDIR)/function.h hybrid_mac_synthetic/hybrid_mac_syn.cpp

$(BUILDDIR)/channel.o : $(INCDIR)/function.h $(SRCDIR)/channel.cpp
	$(CC) $(CFLAGS) $(INC) $(SRCDIR)/channel.cpp -o $@ 

$(BUILDDIR)/function.o : $(INCDIR)/function.h $(SRCDIR)/function.cpp
	$(CC) $(CFLAGS) $(INC) $(SRCDIR)/function.cpp -o $@

$(BUILDDIR)/sparsematrix.o : $(INCDIR)/sparsematrix.h $(SRCDIR)/sparsematrix.cpp
	$(CC) $(CFLAGS) $(INC) $(SRCDIR)/sparsematrix.cpp -o $@

$(BUILDDIR)/node.o : $(INCDIR)/node.h $(SRCDIR)/node.cpp
	$(CC) $(CFLAGS) $(INC) $(SRCDIR)/node.cpp -o $@

$(BUILDDIR)/pcm.o : $(INCDIR)/function.h $(SRCDIR)/pcm.cpp
	$(CC) $(CFLAGS) $(INC) $(SRCDIR)/pcm.cpp -o $@

start:
	mkdir -p $(BUILDDIR) $(TARGET)	

clean:
	rm ./$(BUILDDIR)/*.o

