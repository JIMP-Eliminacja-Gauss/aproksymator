#include "conj_grad_method.h"

int piv_ge_solver(matrix_t *eqs) { 
  if (eqs != NULL) {
      if (conj_grad_solver(eqs) == 0) 
          return 0;
      else 
          return 1;
  }
  else
    return 1;
}

