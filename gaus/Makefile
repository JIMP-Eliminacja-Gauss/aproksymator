libge: matrix.o conj_grad_method.o
	ar rvs libge.a $^

matrix.o: matrix.c matrix.h
pivot.o: pivot.c matrix.h
piv_ge_solver.o: piv_ge_solver.c piv_ge_solver.h matrix.h
conj_grad_method.o: conj_grad_method.c conj_grad_method.h

.PHONY: clean

clean:
	-rm *.o libge.a
