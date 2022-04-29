INCLUDE_DIRS = 
LIB_DIRS = 
#CC = icc
CC = gcc
CPP = g++

CDEFS=
CFLAGS= -g -Wall -fopenmp $(INCLUDE_DIRS) $(CDEFS)
LIBS=

PRODUCT= bisection newton_upd bruteroot regulafalsi_upd regulafalsi_upd2 bruteroot_upd  timeinterp_upd timeinterp_upd_lineup

HFILES= 
CFILES= newton_upd.c bruteroot regulafalsi_upd regulafalsi_upd2 bruteroot_upd  timeinterp_upd timeinterp_upd_lineup
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

timeinterp_upd: timeinterp_upd.c
	$(CC) $(CFLAGS) -o $@ timeinterp_upd.c -lm

timeinterp_upd_lineup: timeinterp_upd_lineup.c
	$(CC) $(CFLAGS) -o $@ timeinterp_upd_lineup.c -lm
