OBJS1 = hybrid_mac_inde.o channel.o function.o sparsematrix.o node.o 
OBJS2 = rcm_mac.o channel.o function.o sparsematrix.o node.o
OBJS3 = hybrid_pam.o channel.o function.o sparsematrix.o node.o
CVERSION = -std=c++0x
CC = g++ $(CVERSION)
SRCDIR = src
BUILDDIR = build
TARGET = bin

SRCEXT = cpp
SOURCES = $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS = $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
DEBUG = -g
CFLAGS =  -c $(DEBUG) -Wall
LFLAGS =  $(DEBUG)
LIB = -L lib
INC = -I include

hybrid_mac_inde : $(OBJS1)
	$(CC) $(LFLAGS) $(OBJS1) -o $(TARGET) hybrid_mac_inde
rcm_mac : $(OBJS2)
	$(CC) $(LFLAGS) $(OBJS2) -o rcm_mac
hybrid_pam : $(OBJS3)
	$(CC) $(LFLAGS) $(OBJS3) -o hybrid_pam

func.h : include/node.h
node.h : include/sparsematrix.h

hybrid_mac_inde.o : include/%.h hybrid_mac_inde/hybrid_mac_inde.cpp
	$(CC) $(CFLAGS) $(INC) hybrid_mac_inde/hybrid_mac_inde.cpp

hybrid_pam.o : node.h sparsematrix.h func.h hybrid_pam.cpp
	$(CC) $(CFLAGS) $(INC) hybrid_pam.cpp

rcm_mac.o : node.h sparsematrix.h func.h rcm_mac.cpp
	$(CC) $(CFLAGS) $(INC) rcm_mac.cpp

$(BUILDDIR)/%.o: include/%.h src/%.cpp
  $(CC) $(CFLAGS) $(INC) -c -o $@ $<

