INCLUDE_DIRS = 
LIB_DIRS = 
#CC = icc
CC = gcc
CPP = g++

CDEFS=
CFLAGS= -g -Wall -fopenmp $(INCLUDE_DIRS) $(CDEFS)
LIBS=

PRODUCT= bisection newton_upd bruteroot regulafalsi_upd regulafalsi_upd2 bruteroot_upd train_bruteroot

HFILES= 
CFILES= newton_upd.c bruteroot regulafalsi_upd regulafalsi_upd2 bruteroot_upd train_bruteroot
CPPFILES= bisection.cpp

SRCS= ${HFILES} ${CFILES}
OBJS= ${CFILES:.c=.o}

all:	${PRODUCT}

clean:
	-rm -f *.o *.NEW *~
	-rm -f ${PRODUCT} ${DERIVED} ${GARBAGE}

bruteroot:	bruteroot.c
	$(CC) $(CFLAGS) -o $@ bruteroot.c -lm

bisection:	bisection.cpp
	$(CPP) $(CFLAGS) -o $@ bisection.cpp -lm

regulafalsi:	regulafalsi.c
	$(CC) $(CFLAGS) -o $@ regulafalsi.c -lm

newton_upd:	newton_upd.c
	$(CC) $(CFLAGS) -o $@ newton_upd.c -lm

regulafalsi_upd: regulafalsi_upd.c
	$(CC) $(CFLAGS) -o $@ regulafalsi_upd.c -lm

regulafalsi_upd2: regulafalsi_upd2.c
	$(CC) $(CFLAGS) -o $@ regulafalsi_upd2.c -lm

bruteroot_upd: bruteroot_upd.c
	$(CC) $(CFLAGS) -o $@ bruteroot_upd.c -lm
train_bruteroot: train_bruteroot.c
	$(CC) $(CFLAGS) -o $@ train_bruteroot.c -lm
