OPT=-g -Wall

#TOT=rf_pendol lorenz_int rtbps_int cmani_rtbp

TOT=chou_comandes

tot : $(TOT)

# ==========
# Utilitats
# ==========

chou_comandes : chou_comandes.c rk78.o
	gcc -o chou_comandes $(OPT) chou_comandes.c rk78.o -lm

rf_pendol : rf_pendol.c pendol.o rk78.o
	gcc -o rf_pendol $(OPT) rf_pendol.c pendol.o rk78.o -lm

lorenz_int : lorenz_int.c lorenz.o flux.o rk78.o
	gcc -o lorenz_int $(OPT) lorenz_int.c lorenz.o flux.o rk78.o -lm

rtbps_int : rtbps_int.c rtbps.o flux.o rk78.o
	gcc -o rtbps_int $(OPT) rtbps_int.c rtbps.o flux.o rk78.o -lm

cmani_rtbp : cmani_rtbp.c cmani.o rtbps.o sislin.o flux.o rk78.o
	gcc -o cmani_rtbp $(OPT) cmani_rtbp.c cmani.o rtbps.o sislin.o flux.o rk78.o -lm

cmani_pendol : cmani_pendol.c cmani.o pendol.o sislin.o flux.o rk78.o
	gcc -o cmani_pendol $(OPT) cmani_pendol.c cmani.o pendol.o sislin.o flux.o rk78.o -lm

# ========
# Rutines
# ========

rk78.o : rk78.c
	gcc -c $(OPT) rk78.c

flux.o : flux.c
	gcc -c $(OPT) flux.c

sislin.o : sislin.c
	gcc -c $(OPT) sislin.c

cmani.o : cmani.c
	gcc -c $(OPT) cmani.c

pendol.o : pendol.c
	gcc -c $(OPT) pendol.c

lorenz.o : lorenz.c
	gcc -c $(OPT) lorenz.c

rtbps.o : rtbps.c
	gcc -c $(OPT) rtbps.c

#========
# Neteja
#========

clean : 
	rm -f *.o

realclean : clean
	rm -f $(TOT)
