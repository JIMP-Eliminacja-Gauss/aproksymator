CFLAGS=-ggdb -Wall -Wextra -pedantic -lm 

aprox: main.o splines.o points.o aproksymator_na_bazie.o gaus/libge.a
	$(CC) $(CFLAGS) -o aprox  main.o splines.o points.o aproksymator_na_bazie.o -L gaus -l ge

intrp: main.o splines.o points.o interpolator.o gaus/libge.a
	$(CC) $(CFLAGS) -o intrp  main.o splines.o points.o interpolator.o -L gaus -l ge

prosta: main.o splines.o points.o prosta.o
	$(CC) $(CFLAGS) -o prosta  main.o splines.o points.o prosta.o	

aproksymator_na_bazie.o: makespl.h points.h gaus/piv_ge_solver.h gaus/conj_grad_method.h
	$(CC) $(CFLAGS) -I gaus -c aproksymator_na_bazie.c

interpolator.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) $(CFLAGS) -I gaus -c interpolator.c

.PHONY: clean

clean:
	-rm *.o aprox intrp prosta
