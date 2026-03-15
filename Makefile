CC=			gcc
CXX=		g++
CFLAGS=		-std=c99 -g -Wall -O3
CXXFLAGS=	$(CFLAGS)
CPPFLAGS=
INCLUDES=
LOBJS=		kommon.o knhx.o tree.o io.o msa.o model.o sfunc.o scfg2.o
AOBJS=
PROG=		phycfg
LIBS=		-lpthread -lz -lm

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address -ldl
endif

.SUFFIXES:.c .cpp .o
.PHONY:all clean depend

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

libphycfg.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

phycfg:libphycfg.a $(AOBJS) main.o
		$(CC) $(CFLAGS) $(AOBJS) main.o -o $@ -L. -lphycfg $(LIBS)

clean:
		rm -fr *.o a.out $(PROG) *~ *.a *.dSYM

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c *.cpp)

# DO NOT DELETE

io.o: phycfg.h kommon.h kseq.h
knhx.o: knhx.h
kommon.o: kommon.h
main.o: kommon.h phycfg.h ketopt.h
model.o: pcpriv.h phycfg.h kommon.h
msa.o: kommon.h phycfg.h
scfg2.o: pcpriv.h phycfg.h kommon.h
sfunc.o: pcpriv.h phycfg.h
tree.o: kommon.h knhx.h phycfg.h khashl.h
