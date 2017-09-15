OBJS1 = hybrid_mac.o channel.o function.o sparsematrix.o node.o
OBJS2 = rcm_mac_pure.o channel.o function.o sparsematrix.o node.o
OBJS3 = mac_synthetic2.o channel.o function.o sparsematrix.o node.o
OBJS4 = hybrid_mac_synthetic.o channel.o function.o sparsematrix.o node.o
OBJS5 = ldgm_mac.o channel.o function.o sparsematrix.o node.o
OBJS6 = mac_synthetic3.o channel.o function.o sparsematrix.o node.o
CC = g++
DEBUG = -g
CFLAGS =  -c $(DEBUG) -std=c++0x
LFLAGS =  $(DEBUG)

hybrid_mac : $(OBJS1)
	$(CC) $(LFLAGS) $(OBJS1) -o hybrid_mac
rcm_pure : $(OBJS2)
	$(CC) $(LFLAGS) $(OBJS2) -o rcm_pure
mac_synthetic2 : $(OBJS3)
	$(CC) $(LFLAGS) $(OBJS3) -o mac_synthetic2
mac_synthetic : $(OBJS4)
	$(CC) $(LFLAGS) $(OBJS4) -o mac_synthetic
ldgm_mac : $(OBJS5)
	$(CC) $(LFLAGS) $(OBJS5) -o ldgm_mac
mac_synthetic3 : $(OBJS6)
	$(CC) $(LFLAGS) $(OBJS6) -o mac_synthetic3

func.h : node.h
node.h : sparsematrix.h

hybrid_mac.o : node.h sparsematrix.h func.h hybrid_mac.cpp
	$(CC) $(CFLAGS) hybrid_mac.cpp

mac_synthetic2.o : node.h sparsematrix.h func.h mac_synthetic2.cpp
	$(CC) $(CFLAGS) mac_synthetic2.cpp

mac_synthetic3.o : node.h sparsematrix.h func.h mac_synthetic3.cpp
	$(CC) $(CFLAGS) mac_synthetic3.cpp

hybrid_mac_synthetic.o : node.h sparsematrix.h func.h hybrid_mac_synthetic.cpp
	$(CC) $(CFLAGS) hybrid_mac_synthetic.cpp

ldgm_mac.o : node.h sparsematrix.h func.h ldgm_mac.cpp
	$(CC) $(CFLAGS) ldgm_mac.cpp

rcm_mac_pure.o : node.h sparsematrix.h func.h rcm_mac_pure.cpp
	$(CC) $(CFLAGS) rcm_mac_pure.cpp

channel.o : func.h channel.cpp
	$(CC) $(CFLAGS) channel.cpp

function.o: func.h function.cpp
	$(CC) $(CFLAGS) function.cpp

sparsematrix.o: sparsematrix.h sparsematrix.cpp
	$(CC) $(CFLAGS) sparsematrix.cpp

node.o : node.h node.cpp
	$(CC) $(CFLAGS) node.cpp

clean:
	\rm *.o *~ hybrid 
